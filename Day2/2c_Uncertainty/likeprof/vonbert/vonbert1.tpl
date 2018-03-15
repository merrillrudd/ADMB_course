DATA_SECTION
  init_int nobs;
  init_matrix data(1,nobs,1,2);
  vector ages(1,nobs);
  vector Lobs(1,nobs);

PARAMETER_SECTION
  init_number Linf;
  init_number K;
  init_number t0;
  init_number logSigma;
  sdreport_number sigma;
  sdreport_vector Lpred(1,nobs);
  likeprof_number Linf_prof;
  likeprof_number K_prof;
  number RSS;
  objective_function_value nll;

PRELIMINARY_CALCS_SECTION
  ages = column(data,1);
  Lobs = column(data,2);
  Linf = 1.1 * max(Lobs);  // sensible starting value

PROCEDURE_SECTION
  Linf_prof = Linf;
  K_prof = K;
  Lpred = Linf * (1.0 - exp(-K * (ages - t0)));
  RSS = sum(square(Lobs-Lpred));
  sigma = exp(logSigma);
  nll = 0.5 * nobs * log(2.0 * M_PI) + (nobs * logSigma) + (RSS / (2.0*square(sigma)));
  
