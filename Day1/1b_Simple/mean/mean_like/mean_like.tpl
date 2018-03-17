DATA_SECTION
  init_int nobs;
  init_vector Y(1,nobs);
  init_int eof;

   !!    if(eof!=999){
   !!       cout << "Error reading data file." <<endl;
   !!      ad_exit(1);
   !!     }
          
  
PARAMETER_SECTION
  init_number mu;
  init_number log_sigma;
  number sigma;
  vector sqdev(1,nobs);
  objective_function_value f;
  
PROCEDURE_SECTION
  sqdev = square(mu-Y);
  // sigma = sqrt(sqdev / nobs);
  sigma = exp(log_sigma);
  f = (0.5 * nobs * log(2 * M_PI)) + nobs * log(sigma) + sum(sqdev)/(2 * square(sigma));

REPORT_SECTION
  report << "Observed" << endl;
  report << Y << endl;
  report << "Mean" << endl;
  report << mu << endl;
  report << "Sigma" << endl;
  report << sigma << endl;
  report << "NLL" << endl;
  report <<  f << endl;
