// Concentrated (short)
DATA_SECTION
  init_int n
  init_vector y(1,n)

PARAMETER_SECTION
  init_number mu
  objective_function_value f

PROCEDURE_SECTION
  f = 0.5 * n * log(sum(square(y-mu)));
