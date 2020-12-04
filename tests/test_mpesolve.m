function fail_count = test_mpesolve()

testname1 = 'mpesolve';
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

An_ref(1) = 0.0165;
An_ref(2) = 0.0001;
An_ref(3) = 0.0009;
An_ref(4) = -0.0001;
An_ref(5) = 0;

if max(abs(An - An_ref)) < error_limit
    print_message(testname1,testname2,'An Error','Pass')
else
    print_message(testname1,testname2,'An Error','Fail')
    fail_count = fail_count + 1;
end

end