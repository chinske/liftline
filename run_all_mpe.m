[theta,y,c,a0,alpha,alpha_zl,clmax_vec] = ...
mpegrid(ybp,cbp,twista,clalp,alpzl,clmax,alpha_r,ncoef);
An = mpesolve(theta,c,a0,alpha,alpha_zl,b);
[Gamma,CL,CDv,cl] = mperesult(An,theta,c,U,b,S);