// Residual sum of squares
DATA_SECTION
  init_int n
  init_vector y(1,n)

PARAMETER_SECTION
  init_number mu
  objective_function_value f

PROCEDURE_SECTION
  f = sum(square(y-mu));
