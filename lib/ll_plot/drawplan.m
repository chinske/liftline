function drawplan(ybp,cbp,sweepa,twista)
% DRAWPLAN Draw Planform
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
% DRAWPLAN(ybp,cbp,sweepa,twista) plots the planform.
% 

% find LE x-values
xle = 0;
for i = 2:length(ybp)
    
    % panel taper ratio and aspect ratio
    tr = cbp(i)./cbp(i-1);
    ar = 4.*(ybp(i)-ybp(i-1))./(cbp(i-1)+cbp(i));
    
    % panel tangent of sweep angle, LE (TSALE)
    tsale = tand(sweepa(i-1)) + ((1-tr)./(ar.*(1+tr)));
    
    xle(i) = xle(i-1) - (ybp(i)-ybp(i-1)).*tsale;
end

% find TE x-values
xte = xle - cbp;

% plot wing planform
x = [xle,fliplr(xte)];
y = [ybp,fliplr(ybp)];
figure
plot(y,x,'-k')
axis equal
xlabel('y')
ylabel('x')

% plot twist vs. span
twist = [0,twista];
figure
plot(ybp/ybp(length(ybp)),twist)
xlabel('y/(b/2)')
ylabel('Twist (deg)')