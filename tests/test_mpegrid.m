function fail_count = test_mpegrid()

testname1 = 'mpegrid';
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

theta_ref(1) = 1.5708;
theta_ref(2) = 1.9635;
theta_ref(3) = 2.3562;
theta_ref(4) = 2.7489;
theta_ref(5) = pi;

if max(abs(theta - theta_ref)) < error_limit
    print_message(testname1,testname2,'THETA Error','Pass')
else
    print_message(testname1,testname2,'THETA Error','Fail')
    fail_count = fail_count + 1;
end

c_ref(1) = 0.7260;
c_ref(2) = 0.5593;
c_ref(3) = 0.4180;
c_ref(4) = 0.3236;
c_ref(5) = lambda*cr;

if max(abs(c - c_ref)) < error_limit
    print_message(testname1,testname2,'Chord Error','Pass')
else
    print_message(testname1,testname2,'Chord Error','Fail')
    fail_count = fail_count + 1;
end

end