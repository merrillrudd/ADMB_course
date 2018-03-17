DATA_SECTION
  init_int nobs;
  init_matrix data(1,nobs,1,2);
  vector SB_obs(1,nobs);
  vector Rec_obs(1,nobs);

  !! SB_obs = column(data,1);
  !! Rec_obs = column(data,2);

PARAMETER_SECTION
  init_number alpha;
  init_number beta;
  init_number logSigma
  sdreport_vector Rec_pred(1,nobs);
  number sigma;
  number RSS;
  objective_function_value nll;

// PRELIMINARY_CALCS_SECTION
//   Rmax = 1.1 * max(Rec_obs);
//   S50 = 0.5 * max(SB_obs);
//   logSigma = log(1);

PROCEDURE_SECTION
  Rec_pred = elem_div(alpha * SB_obs, 1 + beta*SB_obs);
  RSS = sum(square(Rec_obs - Rec_pred));
  sigma = exp(logSigma);
  nll = 0.5 * nobs * log(2.0 * M_PI) + (nobs * logSigma) + (RSS / (2.0*square(sigma)));

REPORT_SECTION
  report << "nll" << endl;
  report << nll << endl;
  report << "SB_obs" << endl;
  report << SB_obs << endl;
  report << "Rec_obs" << endl;
  report << Rec_obs << endl;
  report << "Rec_pred" << endl;
  report << Rec_pred << endl;
  
  
