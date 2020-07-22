% Test Wing

W = 700;

AR = 9;
b = 4.572;
S = (b^2)/AR;
lambda = 0.4;
cr = 0.726;

% --------------------------------------------------
% breakpoints (meters)
ybp = [0,b/2];

% chord at breakpoints (meters)
cbp = [cr,lambda*cr];

% sweep angles (deg)
sweepa = 0;

% c/4 twist angles (deg)
twista = 0;

% --------------------------------------------------
% airfoil parameters

% airfoil section lift-curve slope (per deg)
clalp = [2*pi,2*pi]*pi/180;

% airfoil section zero lift angle of attack (deg)
alpzl = [-1.2,-1.2];

% airfoil section maximum lift coefficient
clmax = [0.84,0.84];