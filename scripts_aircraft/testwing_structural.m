% --- EXAMPLE: ploads
% two 25 N stores, located at y = 1 m and y = 2 m
% one 25 N tip tank
% 
% ploads = [1,-25; ...
%           2,-25; ...
%           b/2,-25];
% uloads = [0,b/2,0];   % no uloads, but must be provided as an input

% --- EXAMPLE: uloads
% equipment from y = 0 m to 1 m / 10 Newtons per meter
% fuel tank from y = 1 m to 2 m / 20 Newtons per meter
% 
ploads = [0,0];         % no ploads, but must be provided as an input
uloads = [0,1,-10;1,2,-20];