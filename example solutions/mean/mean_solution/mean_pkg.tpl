// Likelihood (package)
GLOBALS_SECTION
  #include "dnorm.cpp"

DATA_SECTION
  init_int n
  init_vector y(1,n)

PARAMETER_SECTION
  init_number mu
  init_number logSigma
  objective_function_value f

PROCEDURE_SECTION
  f = dnorm(y, mu, exp(logSigma));
