DATA_SECTION
  init_int nobs;
  init_matrix dat(1,nobs,1,2);// read in data
  vector age(1,nobs);
  vector len(1,nobs);
  
PARAMETER_SECTION
  init_bounded_number linf(20,100);
  init_number t0;
  init_number k;
  init_number logSigma;
  sdreport_vector Lpred(1,nobs);
  number RSS;
  number sigma;
  // set up parameters from .xlsx
  objective_function_value nll;

PRELIMINARY_CALCS_SECTION
 // get age and length data
 age = column(dat,1);
 len = column(dat,2);
 // set initial value for Linf
 linf = max(len);  

PROCEDURE_SECTION
  Lpred = linf * (1 - mfexp(-k * (age - t0))) ;
  RSS = sum(square(len-Lpred));
  sigma = exp(logSigma);
  nll = 0.5 * nobs * log(2.0 * M_PI) + (nobs * logSigma) + (RSS / (2.0*square(sigma)));
  
REPORT_SECTION
  report << "Observed" << endl;
  report << len <<  endl;
  report << "Predicted" << endl;
  report << Lpred << endl;
  report << "log sigma" << endl;
  report << logSigma << endl;
  report << "RSS" << endl;
  report << RSS << endl;
  report << "NLL" << endl;
  report << nll << endl;
