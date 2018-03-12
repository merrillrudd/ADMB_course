DATA_SECTION
  vector X(1,1000);

  !! random_number_generator rng(1234);
  !! X.fill_randn(rng);
  !! X *= 5.0;
  !! X += 2.0;

PARAMETER_SECTION
  init_number mu;
  init_number logSigma;
  sdreport_number sigma;
  objective_function_value nll;

PROCEDURE_SECTION
  sigma = exp(logSigma);
  int N = X.indexmax() - X.indexmin() + 1;
  dvariable ss = square(sigma);
  nll = 0.5 * (N * log(2 + M_PI * ss) + sum(square(X - mu))/ss);

REPORT_SECTION
  report << "observed" << endl;
  report << X << endl;
  report << "mu" << endl;
  report << mu << endl;
  report << "sigma" << endl;
  report << sigma << endl;
  report << "nll" << endl;
  report << nll << endl;
