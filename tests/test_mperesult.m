function fail_count = test_mperesult()

testname1 = 'mperesult';
fail_count = 0;

% --------------------------------------------------
% check against Nelson handout
testname2 = 'Nelson';
error_limit = 0.0001;

addpath scripts_aircraft
addpath scripts_cases
testwing
runcon

[theta,y,c,a0,alpha,alpha_zl,clmax_vec] = ...
    mpegrid(ybp,cbp,twista,clalp,alpzl,clmax,alpha_r,ncoef);
An = mpesolve(theta,c,a0,alpha,alpha_zl,b);
[Gamma,CL,CDv,cl] = mperesult(An,theta,c,U,b,S);

CL_ref = 0.4655;
CDv_ref = 0.0078;

if max(abs(CL - CL_ref)) < error_limit
    print_message(testname1,testname2,'CL Error','Pass')
else
    print_message(testname1,testname2,'CL Error','Fail')
    fail_count = fail_count + 1;
end

if max(abs(CDv - CDv_ref)) < error_limit
    print_message(testname1,testname2,'CDv Error','Pass')
else
    print_message(testname1,testname2,'CDv Error','Fail')
    fail_count = fail_count + 1;
end

end