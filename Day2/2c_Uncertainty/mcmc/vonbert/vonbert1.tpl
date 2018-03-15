DATA_SECTION
  init_int nobs;
  init_matrix data(1,nobs,1,2);
  vector ages(1,nobs);
  vector Lobs(1,nobs);

  //used for mcmc
  int iter;
  !! iter = 0;

PARAMETER_SECTION
  init_number Linf;
  init_number K;
  init_number t0;
  init_number logSigma;
  sdreport_number sigma;
  sdreport_vector Lpred(1,nobs);
  number RSS;
  objective_function_value nll;

PRELIMINARY_CALCS_SECTION
  ages = column(data,1);
  Lobs = column(data,2);
  Linf = 1.1 * max(Lobs);  // sensible starting value

PROCEDURE_SECTION
  Lpred = Linf * (1.0 - exp(-K * (ages - t0)));
  RSS = sum(square(Lobs-Lpred));
  sigma = exp(logSigma);
  nll = 0.5 * nobs * log(2.0 * M_PI) + (nobs * logSigma) + (RSS / (2.0*square(sigma)));

  if(mceval_phase()){
    mcmc_output();
  }

FUNCTION mcmc_output
	if(iter==0)
	{
		ofstream myfile("refpar.mcmc");
		myfile<<"Linf\t K"<<endl;
	}
	iter++;
	ofstream myfile("refpar.mcmc",ios::app);
	myfile<< Linf << "\t" << K << endl;
