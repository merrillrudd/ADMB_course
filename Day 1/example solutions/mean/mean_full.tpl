DATA_SECTION
  init_int nobs;
  init_vector Y(1,nobs);
  init_int eof;
  
  LOCAL_CALCS
  if(eof!=999){
          cout << "Error reading data file." <<endl;
          ad_exit(1);
  }
          
  
PARAMETER_SECTION
  init_number mu;
  init_number logsigma;
  number rss;
  objective_function_value f;
  
PROCEDURE_SECTION
  rss = sum(square(Y-mu));
  f = (0.5 * nobs * log(2 * M_PI)) + (nobs *logsigma) + rss/(2*square(exp(logsigma)));
  
REPORT_SECTION
  report << "Observed" << endl;
  report << Y << endl;
  report << "Average" << endl;
  report << mu << endl;
  report << "Sigma" << endl;
  report << exp(logsigma) << endl;
  report << "NLL" << endl;
  report <<  f << endl;
