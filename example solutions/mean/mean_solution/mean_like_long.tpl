// Likelihood (long)
DATA_SECTION
  init_int n
  init_vector y(1,n)

PARAMETER_SECTION
  init_number mu
  init_number logSigma
  objective_function_value f

PROCEDURE_SECTION
  f = 0.5*n*log(2.0*M_PI) + n*logSigma + sum(square(y-mu))/(2.0*square(exp(logSigma)));
