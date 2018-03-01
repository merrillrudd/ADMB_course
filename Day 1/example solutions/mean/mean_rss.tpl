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
  objective_function_value f;
  
PROCEDURE_SECTION
  f = sum(square(Y-mu));
  
REPORT_SECTION
  report << "Observed" << endl;
  report << Y << endl;
  report << "Average" << endl;
  report << mu << endl;
  report << "NLL" << endl;
  report <<  f << endl;
