DATA_SECTION
  init_int nobs;
  init_vector Y(1,nobs);         
  
PARAMETER_SECTION
  init_number mu;
  init_number logsigma;
  sdreport_number sigma;
  number rss;
  objective_function_value f;
  
PROCEDURE_SECTION
  sigma = exp(logsigma);
  rss = sum(square(Y-mu));
  f = (0.5 * nobs * log(2 * M_PI)) +
    (nobs *logsigma) + rss/(2*square(sigma));
  
REPORT_SECTION
  report << "Observed" << endl;
  report << Y << endl;
  report << "Average" << endl;
  report << mu << endl;
  report << "Sigma" << endl;
  report << sigma << endl;
  report << "NLL" << endl;
  report <<  f << endl;

