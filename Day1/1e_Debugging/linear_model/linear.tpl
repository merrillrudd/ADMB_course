DATA_SECTION
  init_int DEBUG_FLAG;
  init_int n;
  init_vector y(1,n);
  init_vector x(1,n);
  
PARAMETER_SECTION
  init_number a;
  init_number b;
  init_number logSigma;
  number sigma;
  sdreport_number sigmasq;
  sdreport_vector ypred(2,n);
  objective_function_value nll;
  
PROCEDURE_SECTION
  sigma = exp(logSigma);
  if(DEBUG_FLAG != 0) cout << "sigma " << sigma << endl;
  
  sigmasq = square(sigma);
  if(DEBUG_FLAG !=0) cout << "sigmasq" << sigmasq << endl;

  ypred = a*x + b;
  if(DEBUG_FLAG != 0) cout << "ypred" << ypred << endl;
  
  nll = 0.5 * (n * log(2 * M_PI * sigmasq))  + (sum(square(y - ypred)) / sigmasq);
  
REPORT_SECTION
  report << "Observed" << endl;
  report << y << endl;
  report << "Predicted" << endl;
  report << ypred << endl;
  report << "NLL" << endl;
  report << nll << endl;
