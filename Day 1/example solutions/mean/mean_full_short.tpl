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
  number rss;
  objective_function_value f;
  
PROCEDURE_SECTION
  rss = sum(square(Y-mu));
  f = (0.5 * nobs * log(2 * M_PI)) + (0.5 * nobs * log(rss)) - (0.5 * nobs * log(nobs)) + (nobs/2);
  
REPORT_SECTION
  report << "Observed" << endl;
  report << Y << endl;
  report << "Average" << endl;
  report << mu << endl;
  report << "NLL" << endl;
  report <<  f << endl;
