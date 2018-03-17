DATA_SECTION
  init_int nobs;
  init_vector Y(1,nobs);
  init_vector x(1,nobs);
  init_int eof;
  
  LOCAL_CALCS
  if(eof!=999){
          cout << "Error reading data file." <<endl;
          ad_exit(1);
  }
          
  
PARAMETER_SECTION
  init_number a;
  init_number b;
  vector pred_Y(1,nobs);
  vector sqdev(1,nobs);
  sdreport_number aa;
  objective_function_value f;
  
PROCEDURE_SECTION
  aa = a;
  pred_Y = a*x+b;
  sqdev = square(pred_Y-Y);
  f = sum(sqdev);

REPORT_SECTION
  report << "Observed" << endl;
  report << Y << endl;
  report << "Predicted" << endl;
  report << pred_Y << endl;
  report << "Parameter a" << endl;
  report << a << endl;
  report << "Parameter b" << endl;
  report << b << endl; 
  report << "NLL" << endl;
  report <<  f << endl;

  ofstream f1("Param_est");
  f1 << "Parameter a" << endl;
  f1 << a << endl;
  f1 << "Parameter b" << endl;
  f1 << b << endl;
