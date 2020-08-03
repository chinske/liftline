function liftlineui()
% LIFTLINEUI Lifting-Line Theory User Interface
% Copyright 2020 Christopher Chinske
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
% 
% LIFTLINEUI
% 

selected_aircraft = 'None';
selected_runcon = 'None';

run_flag_global = 1;

run_flag_local = 1;
while run_flag_local == 1 && run_flag_global == 1
    home
    disp(' ')
    disp('LIFTLINE')
    disp(' ')
    disp('------------------------------')
    disp(['Selected Aircraft: ',selected_aircraft])
    disp(['Selected Run Conditions: ',selected_runcon])
    disp('------------------------------')
    disp(' ')
    disp('Options:')
    disp('[1] Select aircraft script')
    disp('[2] Select run conditions script')
    disp('[c] Continue')
    disp('[e] Exit')
    disp(' ')
    
    choice = input('Select an option: ','s');
    if strcmpi(choice,'e')
        % exit
        run_flag_global = 0;
    elseif strcmpi(choice,'c')
        % continue
        run_flag_local = 0;
    elseif strcmpi(choice,'1')
        % select aircraft
        [aircraft_file,aircraft_path] = uigetfile();
        selected_aircraft = [aircraft_path,aircraft_file];
    elseif strcmpi(choice,'2')
        % select run conditions
        [runcon_file,runcon_path] = uigetfile();
        selected_runcon = [runcon_path,runcon_file];
    else
        % invalid input
    end
    
end

run_flag_local = 1;
while run_flag_local == 1 && run_flag_global == 1
    home
    disp('Options:')
    disp('[1] Grid Convergence Analysis')
    disp('[2] Find AoA for L = W')
    disp('[3] Run Case as Defined in Run Conditions Script')
    disp('[4] AoA Range')
    disp('[5] Spar Shear and Bending Moment')
    disp(' ')
    
    choice = input('Select an option: ','s');
    if strcmpi(choice,'1')
        % grid convergence analysis
        grid_convergence(selected_aircraft,selected_runcon)
        run_flag_global = 0;
    elseif strcmpi(choice,'2')
        % find AoA for L = W
        find_aoa(selected_aircraft,selected_runcon)
        run_flag_global = 0;
    elseif strcmpi(choice,'3')
        % run case as defined in run conditions script
        run_case(selected_aircraft,selected_runcon)
        run_flag_global = 0;
    elseif strcmpi(choice,'4')
        % AoA range
        run_aoa_range(selected_aircraft,selected_runcon)
        run_flag_global = 0;
    elseif strcmpi(choice,'5')
        % spar shear and bending moment
        run_flag_global = 0;
    else
        % invalid input
    end
    
end

end % end main function

% --------------------------------------------------
function grid_convergence(selected_aircraft,selected_runcon)

disp('Loading aircraft script...')
run(selected_aircraft)

disp('Loading run conditions script...')
run(selected_runcon)

disp('Recommended CL and CD significant figures: 5')
sigfig_CL = input('Desired CL significant figures: ');
sigfig_CD = input('Desired CD significant figures: ');

tol_CL = 0.5 * 10^(2-sigfig_CL);
tol_CD = 0.5 * 10^(2-sigfig_CD);

disp(['CL Tolerance: ',num2str(tol_CL)])
disp(['CD Tolerance: ',num2str(tol_CD)])

disp('Running grid convergence analysis...')
ncoef = 2;
perr_CL = 1;
perr_CD = 1;
CL_old = 10;
CD_old = 10;
while abs(perr_CL) >= tol_CL || abs(perr_CD) >= tol_CD
    
    ncoef = ncoef + 1;
    
    [theta,y,c,a0,alpha,alpha_zl,clmax_vec] = ...
        mpegrid(ybp,cbp,twista,clalp,alpzl,clmax,alpha_r,ncoef);
    An = mpesolve(theta,c,a0,alpha,alpha_zl,b);
    [Gamma,CL,CDv,cl] = mperesult(An,theta,c,U,b,S);
    
    CD = CDv;
    
    perr_CL = 100.*(CL - CL_old)./CL;
    perr_CD = 100.*(CD - CD_old)./CD;
    
    CL_old = CL;
    CD_old = CD;
    
    disp(['ncoef: ',num2str(ncoef), ...
        ' | perr_CL: ',num2str(abs(perr_CL)), ...
        ' | perr_CD: ',num2str(abs(perr_CD))])
    
end
disp(['Required ncoef: ',num2str(ncoef)])
end

% --------------------------------------------------
function find_aoa(selected_aircraft,selected_runcon)
% use the secant method to find AoA for L = W

disp('Loading aircraft script...')
run(selected_aircraft)

disp('Loading run conditions script...')
run(selected_runcon)

tol = 5E-6;
perr = 1;
alpha_r0 = 0;
alpha_r1 = 5;
it = 0;
while perr > tol
    
    it = it + 1;
    
    % eval for alpha_r0
    [theta,y,c,a0,alpha,alpha_zl,clmax_vec] = ...
        mpegrid(ybp,cbp,twista,clalp,alpzl,clmax,alpha_r0,ncoef);
    An = mpesolve(theta,c,a0,alpha,alpha_zl,b);
    [Gamma,CL,CDv,cl] = mperesult(An,theta,c,U,b,S);
    
    f0 = 0.5.*rho.*(U.^2).*S.*CL - W;
    
    % eval for alpha_r1
    [theta,y,c,a0,alpha,alpha_zl,clmax_vec] = ...
        mpegrid(ybp,cbp,twista,clalp,alpzl,clmax,alpha_r1,ncoef);
    An = mpesolve(theta,c,a0,alpha,alpha_zl,b);
    [Gamma,CL,CDv,cl] = mperesult(An,theta,c,U,b,S);
    
    f1 = 0.5.*rho.*(U.^2).*S.*CL - W;
    
    % apply the secant method
    alpha_r_next = alpha_r1 - ...
        (f1.*(alpha_r0-alpha_r1))./(f0-f1);
    
    alpha_r0 = alpha_r1;
    alpha_r1 = alpha_r_next;
    perr = 100.*abs((alpha_r1 - alpha_r0)./alpha_r1);
    
end

% eval at solution point
alpha_r = alpha_r1;
[theta,y,c,a0,alpha,alpha_zl,clmax_vec] = ...
    mpegrid(ybp,cbp,twista,clalp,alpzl,clmax,alpha_r,ncoef);
An = mpesolve(theta,c,a0,alpha,alpha_zl,b);
[Gamma,CL,CDv,cl] = mperesult(An,theta,c,U,b,S);

L = 0.5.*rho.*(U.^2).*S.*CL;

disp(['Required alpha_r: ',num2str(alpha_r)])
disp(['Input Weight: ',num2str(W)])
disp(['Computed Lift: ',num2str(L)])
disp(['Percent Error: ',num2str(100.*abs((L-W)./W))])
disp(['Iterations: ',num2str(it)])

end

% --------------------------------------------------
function run_case(selected_aircraft,selected_runcon)

disp('Loading aircraft script...')
run(selected_aircraft)

disp('Loading run conditions script...')
run(selected_runcon)

disp('Running case...')
run_all_mpe

% display results
disp('--------------------------------------------------')
disp(' ')

disp(['Wing lift coefficient: ',num2str(CL)])
disp(['Wing vortex-induced drag coefficient: ',num2str(CDv)])
disp(' ')

L = 0.5.*rho.*(U.^2).*S.*CL;
Dv = 0.5.*rho.*(U.^2).*S.*CDv;
disp(['Lift (N): ',num2str(L)])
disp(['Vortex-induced drag (N): ',num2str(Dv)])
disp(' ')

disp('--------------------------------------------------')

% plot results
s.y = y;
s.b = b;
s.Gamma = Gamma;
s.cl = cl;
s.clmax_vec = clmax_vec;

plotflag = [1 1 0 0];
mpeplot(s,plotflag)

end

% --------------------------------------------------
function run_aoa_range(selected_aircraft,selected_runcon)

disp('Loading aircraft script...')
run(selected_aircraft)

disp('Loading run conditions script...')
run(selected_runcon)

aoa_i = input('Input initial AoA (deg): ');
aoa_f = input('Input final AoA (deg):   ');
d_aoa = input('Input increment (deg):   ');

disp('Iterating over AoA range...')
% iterate over AoA range
i = 1;
for alpha_r = aoa_i:d_aoa:aoa_f
    
    run_all_mpe
    aoa_values(i) = alpha_r;
    CL_values(i) = CL;
    CDv_values(i) = CDv;
    
    % compute how much of the wing is stalled
    stall_ratio(i) = sum(cl > clmax_vec)./length(cl);
    
    i = i+1;
    
end

% --------------------------------------------------
% separate data based on stall ratio

disp('Preparing data to plot...')

% 50 % of the wing is stalled => full stall
sr_full_stall = 0.5;

% initialize plotdata indices
plotdata(1).ind = 1;
plotdata(2).ind = 1;
plotdata(3).ind = 1;

% check first value
stall_state = check_stall_state(stall_ratio(1),0,sr_full_stall);
plotdata(stall_state).CL_values = CL_values(1);
plotdata(stall_state).CDv_values = CDv_values(1);
plotdata(stall_state).aoa_values = aoa_values(1);
plotdata(stall_state).ind = plotdata(stall_state).ind + 1;

% check rest of values
for i = 2:length(stall_ratio)
    
    % check if stall_ratio is different than previous value
    % if different, carry previous plot forward by one point
    % this will prevent gaps in the plot
    if stall_ratio(i) ~= stall_ratio(i-1)
        stall_state = ...
            check_stall_state(stall_ratio(i-1),0,sr_full_stall);
        ind = plotdata(stall_state).ind;
        
        plotdata(stall_state).CL_values(ind) = CL_values(i);
        plotdata(stall_state).CDv_values(ind) = CDv_values(i);
        plotdata(stall_state).aoa_values(ind) = aoa_values(i);
        plotdata(stall_state).ind = plotdata(stall_state).ind + 1;
    end
    
    stall_state = ...
        check_stall_state(stall_ratio(i),0,sr_full_stall);
    ind = plotdata(stall_state).ind;
    
    plotdata(stall_state).CL_values(ind) = CL_values(i);
    plotdata(stall_state).CDv_values(ind) = CDv_values(i);
    plotdata(stall_state).aoa_values(ind) = aoa_values(i);
    plotdata(stall_state).ind = plotdata(stall_state).ind + 1;
    
end

% --------------------------------------------------
% plot results

figure
for i = 1:length(plotdata)
    
    if i == 1
        lspec = '-k';
        lwid = 1;
        dname = 'Normal';
    elseif i == 2
        lspec = '-b';
        lwid = 2;
        dname = 'Partially Stalled';
    elseif i == 3
        lspec = '-r';
        lwid = 2;
        dname = 'Fully Stalled';
    end
    
    plot(plotdata(i).aoa_values,plotdata(i).CL_values, ...
        lspec,'LineWidth',lwid,'DisplayName',dname)
    hold on
    xlabel('\alpha (deg)')
    ylabel('C_L')
    title('Wing Lift Coefficient')
    legend show
    legend('Location','southoutside')
    grid on
    grid minor
    
end

figure
for i = 1:length(plotdata)
    
    if i == 1
        lspec = '-k';
        lwid = 1;
        dname = 'Normal';
    elseif i == 2
        lspec = '-b';
        lwid = 2;
        dname = 'Partially Stalled';
    elseif i == 3
        lspec = '-r';
        lwid = 2;
        dname = 'Fully Stalled';
    end
    
    plot(plotdata(i).CDv_values,plotdata(i).CL_values, ...
        lspec,'LineWidth',lwid,'DisplayName',dname)
    hold on
    xlabel('C_{Dv}')
    ylabel('C_L')
    title('Wing Drag Polar')
    legend show
    legend('Location','southoutside')
    grid on
    grid minor
    
end

end

% --------------------------------------------------
function stall_state = ...
    check_stall_state(stall_ratio,val_no_stall,val_full_stall)

if stall_ratio == val_no_stall
    stall_state = 1;
elseif stall_ratio > val_no_stall && stall_ratio < val_full_stall
    stall_state = 2;
elseif stall_ratio >= val_full_stall
    stall_state = 3;
end

end