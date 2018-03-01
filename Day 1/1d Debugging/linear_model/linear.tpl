DATA_SECTION
  init_int n;
  init_vector y(1,n);
  init_vector x(1,n);
  init_int eof;
  
  LOCAL_CALCS
  if(eof!=999){
          cout << "Error reading data file." <<endl;
          ad_exit(1);
  }
  if(eof==999) cout << "Data read correctly!" << endl;
          
  
PARAMETER_SECTION
  init_number a;
  init_number b;
  init_number logSigma;
  sdreport_number sigmasq;
  sdreport_vector ypred(2,n);
  objective_function_value nll;
  
PROCEDURE_SECTION
  sigmasq = exp(2*logSigma);
  ypred = a*x + b;
  nll = 0.5 * (n * log(2 * M_PI * sigmasq) + sum(square(y - ypred)) / sigmasq);

REPORT_SECTION
  report << "Observed" << endl;
  report << y << endl;
  report << "Predicted" << endl;
  report << ypred << endl;
  report << "NLL" << endl;
  report << nll << endl;
