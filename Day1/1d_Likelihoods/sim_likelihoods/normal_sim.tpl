DATA_SECTION
  vector X(1,10000);
  vector Y(1,10000);
  vector Z(1,10000);

 LOCAL_CALCS
    random_number_generator rng(1234);
    X.fill_randn(rng);
    Y = X * 5.0;
    Z = Y + 2.0;
 END_CALCS

PARAMETER_SECTION
  init_number mu;
  init_number logSigma;
  sdreport_number sigma;
  objective_function_value nll;

PROCEDURE_SECTION
  sigma = exp(logSigma);
  int N = Z.indexmax() - Z.indexmin() + 1;
  dvariable ss = square(sigma);
  nll = 0.5 * (N * log(2 + M_PI * ss) + sum(square(Z - mu))/ss);

REPORT_SECTION
  report << "random_numbers" << endl;
  report << X << endl;
  report << "times5" << endl;
  report << Y << endl;
  report << "add2" << endl;
  report << Z << endl;
  report << "mu" << endl;
  report << mu << endl;
  report << "sigma" << endl;
  report << sigma << endl; 
  report << "nll" << endl;
  report << nll << endl;
