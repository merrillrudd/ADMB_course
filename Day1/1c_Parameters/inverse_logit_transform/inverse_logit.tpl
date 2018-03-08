DATA_SECTION
  init_int dim;
  init_vector X(1,dim);
  
PARAMETER_SECTION
  init_vector a(1,dim-1);
  sdreport_vector p(1,dim);
  objective_function_value nll;
  
PROCEDURE_SECTION
  p = a2p(a);
  nll = -gammln(sum(X)+1.0) + sum(gammln(X+1.0))
        - sum(elem_prod(X, log(p)));
  
FUNCTION dvar_vector a2p(const dvar_vector& a)
  dvar_vector p(1,dim);
  dvar_vector expa = exp(a);
  p(1,dim-1) = expa / (1 + sum(expa));
  p(dim) = 1 - sum(p(1,dim-1));
  return p;

