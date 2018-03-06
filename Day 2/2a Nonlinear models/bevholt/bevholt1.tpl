DATA_SECTION
  init_int nobs;
  init_matrix data(1,nobs,1,2);
  vector S(1,nobs);
  !! S = column(data,1);
  vector Robs(1,nobs)
  !! Robs = column(data,2);

PARAMETER_SECTION
  init_number Rmax;
  init_number Shalf;
  init_number logSigma;
  sdreport_number sigma;
  sdreport_vector Rpred(1,nobs);
  number RSS;
  objective_function_value nll;

PROCEDURE_SECTION
  Rpred = Rmax * elem_div(S,S+Shalf);
  RSS = sum(square(log(Robs)-log(Rpred)));
  sigma = exp(logSigma);
  nll = 0.5 * nobs *log(2.0 * M_PI) + (nobs * logSigma) + (RSS/(2.0 * square(sigma)));

REPORT_SECTION
  report << "Observed_rec" << endl;
  report << Robs << endl;
  report << "Predicted_rec" << endl;
  report << Rpred << endl;
  report << "Observed_spawn" << endl;
  report << S << endl;
  
