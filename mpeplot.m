function mpeplot(s,plotflag)
% MPEPLOT Plot MPE Results
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
% MPEPLOT(S,PLOTFLAG) plots MPE results.
% 
% S is a structure array that contains the fields required for the
% specified plots.
% 
% PLOTFLAG is a vector that specifies which plots to generate.  The
% elements of PLOTFLAG specify the following plots.
% 
%   1. Circulation
%   2. Section Lift Coefficient
%   3. Shear Diagram
%   4. Bending Moment Diagram
% 
% For example, PLOTFLAG = [1 0 1 1] would plot Circulation, Shear Diagram,
% and Bending Moment Diagram.
% 
% The required fields in S for each plot are listed below.
% 
% Circulation : y,b,Gamma
% Section Lift Coefficient : y,b,cl,clmax_vec
% Shear Diagram : y,b,V
% Bending Moment Diagram : y,b,V,M
% 

plabel = input('Enter plot label: ','s');

if plotflag(1) == 1
    y = s.y;
    b = s.b;
    Gamma = s.Gamma;
    
    figure(1)
    plot(y/(b/2),Gamma)
    xlabel('y/(b/2)')
    ylabel('Gamma (m^2/s)')
    title('Circulation')
    grid on
end

if plotflag(2) == 1
    y = s.y;
    b = s.b;
    cl = s.cl;
    clmax_vec = s.clmax_vec;
    
    figure(2)
    plot(y/(b/2),cl)
    hold on
    plot(y/(b/2),clmax_vec,'-r')
    hold off
    xlabel('y/(b/2)')
    ylabel('Section Lift Coefficient')
    title('Section Lift Coefficient')
    grid on
end

if plotflag(3) == 1
    y = s.y;
    b = s.b;
    V = s.V;
    
    figure(3)
    plot(y/(b/2),V,'DisplayName',plabel)
    xlabel('y/(b/2)')
    ylabel('Shear (N)')
    title('Shear Diagram')
    legend show
    legend('Location','best');
    grid on
    hold on
end

if plotflag(4) == 1
    y = s.y;
    b = s.b;
    M = s.M;
    
    figure(4)
    plot(y/(b/2),M,'DisplayName',plabel)
    xlabel('y/(b/2)')
    ylabel('Bending Moment (N*m)')
    title('Bending Moment Diagram')
    legend show
    legend('Location','best')
    grid on
    hold on
end