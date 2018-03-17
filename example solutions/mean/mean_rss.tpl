DATA_SECTION
  init_int nobs;
  init_vector Y(1,nobs);
  
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
