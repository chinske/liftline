function An = mpesolve(theta,c,a0,alpha,alpha_zl,b)
% MPESOLVE Monoplane Equation Solver
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
% An = MPESOLVE(theta,c,a0,alpha,alpha_zl,b) solves the monoplane equation
% for symmetrical loading.  Returns 1-by-n vector An of odd term Fourier
% sine series coefficients, An = [A1, A3, ... A(2n-1)].
% 
% theta    : 1-by-n vector, transformed spanwise coordinate (rad)
% c        : 1-by-n vector, chord values
% a0       : 1-by-n vector, section lift-curve slope values (per rad)
% alpha    : 1-by-n vector, section angle of attack values (rad)
% alpha_zl : 1-by-n vector, section zero lift angle of attack values (rad)
% b        : span

n = length(theta);

% remove tip values
% Gamma(theta=pi) = 0 => trivial row in matrix equation
theta = theta(1:n-1);
c = c(1:n-1);
a0 = a0(1:n-1);
alpha = alpha(1:n-1);
alpha_zl = alpha_zl(1:n-1);

% --------------------------------------------------
% solve matrix equation

mu = c.*a0./(4.*b);

bvec = mu.*(alpha - alpha_zl).*sin(theta);

for i = 1:length(theta)
    for j = 1:length(theta)
        A(i,j) = sin((2.*j-1).*theta(i)).* ...
                 ((2.*j-1).*mu(i) + sin(theta(i)));
    end
end

bvec = bvec';
An = A\bvec;
An = An';

% --------------------------------------------------
% include trivial An value
An(n) = 0;