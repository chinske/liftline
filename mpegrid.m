function [theta,y,c,a0,alpha,alpha_zl,clmax_vec] = ...
    mpegrid(ybp,cbp,twista,clalp,alpzl,clmax,alpha_r,ncoef)
% MPEGRID Monoplane Equation Grid
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
% [theta,y,c,a0,alpha,alpha_zl,clmax_vec] = ...
% MPEGRID(ybp,cbp,twista,clalp,alpzl,clmax,alpha_r,ncoef) returns vectors
% of coordinates and section values at N spanwise locations from
% 
%   0 <= y <= b/2.
% 
% MPEGRID outputs are typically used as inputs to MPESOLVE.
% 

% convert twist angle (twista) to rad
twista = twista.*pi./180;

% convert root angle of attack (alpha_r) to rad
alpha_r = alpha_r.*pi./180;

% convert airfoil parameters to 1/rad (clalp) and rad (alpzl)
clalp = clalp.*180./pi;
alpzl = alpzl.*pi./180;

% use ybp to get span
b = 2*ybp(length(ybp));

% build theta and y vectors
theta = linspace(pi/2,pi,ncoef);
y = (b/2)*(-cos(theta));

% numerically, cos(pi/2) evaluates to a very small, non-zero number
% so, force y = 0 at root
y(1) = 0;

% --------------------------------------------------
% build vectors, start at panel 1
panel = 1;

% initialize twist at root and tip of panel
twist_r = 0;
twist_t = twista(panel);

for i = 1:length(theta)
    
    % check for breakpoint crossing (i.e., next panel)
    % if next panel, re-initialize twist at root and tip of panel
    if y(i) > ybp(panel+1)
        panel = panel + 1;
        twist_r = twista(panel-1);
        twist_t = twista(panel);
    end
    
    % ------------------------------
    % build chord vector
    
    % local (panel) root and tip values
    y_r = ybp(panel);
    y_t = ybp(panel+1);
    c_r = cbp(panel);
    c_t = cbp(panel+1);
    
    c(i) = c_r + ((c_t-c_r)./(y_t-y_r)).*(y(i)-y_r);
    
    % ------------------------------
    % build section lift-curve slope vector
    clalp_r = clalp(panel,1);
    clalp_t = clalp(panel,2);
    
    a0(i) = clalp_r + ((clalp_t-clalp_r)./(y_t-y_r)).*(y(i)-y_r);
    
    % ------------------------------
    % build section zero lift angle of attack vector
    alpzl_r = alpzl(panel,1);
    alpzl_t = alpzl(panel,2);
    
    alpha_zl(i) = alpzl_r + ((alpzl_t-alpzl_r)./(y_t-y_r)).*(y(i)-y_r);
    
    % ------------------------------
    % build clmax vector
    clmax_r = clmax(panel,1);
    clmax_t = clmax(panel,2);
    
    clmax_vec(i) = clmax_r + ((clmax_t-clmax_r)./(y_t-y_r)).*(y(i)-y_r);
    
    % ------------------------------
    % build angle of attack vector
    twist(i) = twist_r + ((twist_t-twist_r)./(y_t-y_r)).*(y(i)-y_r);
    alpha(i) = alpha_r + twist(i);
    
end