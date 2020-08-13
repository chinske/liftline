function [Gamma,CL,CDv,cl] = mperesult(An,theta,c,U,b,S)
% MPERESULT Monoplane Equation Result
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
% [Gamma,CL,CDv,cl] = MPERESULT(An,theta,c,U,b,S) returns aerodynamic
% results based on the output of MPESOLVE.
% 

% compute spanwise circulation distribution
Gamma = zeros(size(An));
for i = 1:length(An)
    for j = 1:length(An)
        n = 2*j-1;
        Gamma(i) = Gamma(i) + An(j).*sin(n.*theta(i));
    end
end
Gamma = 2.*b.*U.*Gamma;

% compute wing lift coefficient
CL = An(1).*pi.*(b.^2)./S;

% compute wing vortex-induced drag coefficient
CDv_sum = 0;
for i = 1:length(An)
    n = 2*i-1;
    CDv_sum = CDv_sum + n.*An(i).^2;
end
CDv = pi.*((b.^2)./S).*CDv_sum;

% compute section lift coefficient
cl = 2.*Gamma./(U.*c);