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
//  vector sqdev(1,nobs);
  objective_function_value f;
  
PROCEDURE_SECTION
  // sqdev = square(mu-Y);
  // f = sum(sqdev);
  f = sum(square(Y - mu));

REPORT_SECTION
  report << "Observed" << endl;
  report << Y << endl;
  report << "Mean" << endl;
  report << mu << endl;
  report << "NLL" << endl;
  report <<  f << endl;
