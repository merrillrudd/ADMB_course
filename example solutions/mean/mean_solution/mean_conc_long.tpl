// Concentrated (full)
DATA_SECTION
  init_int n
  init_vector y(1,n)

PARAMETER_SECTION
  init_number mu
  objective_function_value f

PROCEDURE_SECTION
  f = 0.5*n*log(2.0*M_PI) + 0.5*n*log(sum(square(y-mu))) - 0.5*n*log(n) + n/2;
