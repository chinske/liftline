function [V,M,M_root] = spareq(Gamma,y,rho,U,varargin)
% SPAREQ Spar Equilibrium
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
% [V,M,M_root] = SPAREQ(Gamma,y,rho,U)
% [V,M,M_root] = SPAREQ(Gamma,y,rho,U,ploads,uloads)
% 
% Returns shear and bending moment for the cantilever wing.  Optional input
% arguments ploads and uloads specify point loads and uniform distributed
% loads along the wing, respectively.
% 

if nargin == 4
    % ploads and uloads not provided
    ploads = [0,0];
    uloads = [0,y(length(y)),0];
elseif nargin == 6
    ploads = varargin{1};
    uloads = varargin{2};
else
    error('Invalid number of input arguments.')
end

% lift per unit span
lift = rho.*U.*Gamma;

% initialize V and M
V = zeros(size(y));
M = zeros(size(y));

% get semi-span from y
n = length(y);
b2 = y(n);

% ----------
% test case
% p0 = 1;
% lift = p0.*(1 - (y./b2).^3);
% ----------

% check ploads y-values strictly increasing
for i = 2:size(ploads,1)
    if ploads(i,1) <= ploads(i-1,1)
        error('ploads y-values must be strictly increasing.')
    end
end

% check uloads that y2 is always greater than y1
for i = 1:size(uloads,1)
    if uloads(i,2) <= uloads(i,1)
        error('Check uloads.  y2 must always be greater than y1.')
    end
end

% get number of ploads
n_ploads = size(ploads,1);
index_next_pload = n_ploads;

% get number of uloads
n_uloads = size(uloads,1);

% check for ploads at tip (e.g., tip tank)
% if last entry of ploads is at tip, apply at tip
if ploads(n_ploads,1) >= b2
    V(n) = -ploads(n_ploads,2);
    index_next_pload = index_next_pload - 1;
end

% integrate from tip to root for V
for i = n-1:-1:1
    
    % y-values
    a = y(i+1); % a = outboard / tip
    b = y(i);   % b = inboard / root
    
    % lift per unit span values
    lift_a = lift(i+1);
    lift_b = lift(i);
    
    % check for ploads within interval b <= y < a
    sum_local_ploads = 0;
    flag = 1;
    while (flag == 1) && (index_next_pload > 0)
        if ploads(index_next_pload,1) < a && ...
                ploads(index_next_pload,1) >= b
            sum_local_ploads = sum_local_ploads + ...
                ploads(index_next_pload,2);
            index_next_pload = index_next_pload - 1;
        else
            flag = 0;
        end
    end
    
    % check for uloads within interval b <= y <= a
    % if any overlap exists, apply uniform load to entire interval
    sum_local_uloads = 0;
    for ii = 1:n_uloads
        y1 = uloads(ii,1);
        y2 = uloads(ii,2);
        if ~(b>y2) && ~(a<y1)
            % overlap exists
            sum_local_uloads = sum_local_uloads + ...
                uloads(ii,3);
        end
    end
    
    fa = lift_a + sum_local_uloads;
    fb = lift_b + sum_local_uloads;
    I = (b-a).*(fa + fb)./2;
    V(i) = V(i+1) + I - sum_local_ploads;
end

% integrate from tip to root for M
for i = n-1:-1:1
    
    % y-values
    a = y(i+1); % a = outboard / tip
    b = y(i);   % b = inboard / root
    
    % shear values
    fa = V(i+1);
    fb = V(i);
    
    I = (b-a).*(fa + fb)./2;
    M(i) = M(i+1) + I;
end

M_root = M(1);