fail_count = 0;

fail_count = fail_count + test_mpegrid();
fail_count = fail_count + test_mpesolve();
fail_count = fail_count + test_mperesult();

disp(['Failed Tests: ',num2str(fail_count)])