DATA_SECTION
  // read in data

PARAMETER_SECTION
  // set up parameters from .xlsx
  objective_function_value nll;

PRELIMINARY_CALCS_SECTION
 // get age and length data
 // set initial value for Linf

PROCEDURE_SECTION
  Lpred = 
  RSS = 
  sigma = 
  nll = 0.5 * nobs * log(2.0 * M_PI) + (nobs * logSigma) + (RSS / (2.0*square(sigma)));
  
