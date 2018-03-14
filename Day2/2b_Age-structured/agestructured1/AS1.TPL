DATA_SECTION
 init_int Nyear;                                 // Number of years
 init_int Nage;                                  // Last age-class
 init_number Mval;                               // Natural mortality
 init_vector Weight(0,Nage);                     // Weight-at-age
 init_number SigCatch;                           // CV for catch data
 init_number SigCPUE;                            // Sigma for log(CPUE)
 init_number Omega;                              // Weight for proportion data
 init_matrix CatchCPUE(1,Nyear,0,2)              // Catch and cpue data
 init_matrix Propn(1,Nyear,-1,Nage)              // Proportion data
 vector Catch(1,Nyear);                          // Catch data
 vector CPUE(1,Nyear);                           // CPUE data
 
 !! Catch = column(CatchCPUE,1);                 // Extract the catch data
 !! CPUE = column(CatchCPUE,2);                  // Extract the CPUE data

PARAMETER_SECTION
 init_vector logN1(1,Nage)                // Log of abundance at age year 1
 init_vector logR(1,Nyear)                // log of recruits in all years
 init_bounded_number Sel50(0,Nage,1)             // Age-at-50%-selectivity
 init_bounded_number Sel95(0,Nage,1)             // Age-at-95%-selectivity
 init_vector logF(1,Nyear,1)                     // Log of fishing mortality
 init_number logq(1)                             // Log of catchability

 matrix N(1,Nyear,0,Nage)                      // Numbers-at-age
 vector S(0,Nage)                                // Selectivity-at-age
 matrix F(1,Nyear,0,Nage)                        // F-matrix
 matrix Z(1,Nyear,0,Nage)                        // Z-matrix
 
 vector Catch_pred(1,Nyear);                      // Predicted catch (weight);
 vector CPUE_pred(1,Nyear);                       // Predicted CPUE (weight);
 matrix Propn_pred(1,Nyear,0,Nage);               // Predicted catch-at-age
 vector Bio(1,Nyear);                             // vulnerable biomass
 
 number NLL1;                                   // Catch data
 number NLL2;                                   // CPUE data
 number NLL3;                                   // Proportion data
  objective_function_value objn;

// // ==========================================================================

PROCEDURE_SECTION
// int Year;

 // Set up the selectivity pattern
 Select();

 // Project the model forward and compute various outputs
 Numbers();

 // observation model - predict catch, cpue, and catch-at-age
 Predict_Observations();
 
 // Compute the likelihood
 Likelihood();

 objn = NLL1 + NLL2 + NLL3;
  cout << objn << endl;



// ---------------------------------------------------------------------------- 

FUNCTION Select
 int Age;
 
 for (Age=0;Age<=Nage;Age++){
   S(Age) = 1.0 / (1 + exp(-log(19) * (Age-Sel50) / (Sel95-Sel50)));
  }

//alternative:
//for(int Age=0;Age<=Nage;Age++){}

// in R
//for(Age in 0:Nage){}



// ---------------------------------------------------------------------------- 

FUNCTION Numbers
 int Age,Year;
 
 // Clear the N matrix
 N.initialize();

 // Compute the F matrix
 for (Year=1;Year<=Nyear;Year++){
  for (Age=0;Age<=Nage;Age++){
   F(Year,Age) = mfexp(logF(Year)) * S(Age);
  }
 }
 Z = F + Mval;  

 // Insert the abundance from ages 1-Nage in the first year
 for (Age=1;Age<=Nage;Age++){
  N(1,Age) = mfexp(logN1(Age));
 }
 // Insert the recruits age=0 for all years
 for (Year=1;Year<=Nyear;Year++){
  N(Year,0) = mfexp(logR(Year)); 
 }
 
 // Project the whole N matrix
 for (Year=1;Year<Nyear;Year++){
  for (Age=0;Age<Nage;Age++){
   N(Year+1,Age+1) = N(Year,Age) * mfexp(-Z(Year,Age));
  }
 }

FUNCTION Predict_Observations
 int Year, Age;
 
 // Compute the predicted exploitable biomass, catch-at-age and catch
 Propn_pred.initialize();
 for (Year=1;Year<=Nyear;Year++){
   Bio(Year) = 0;
   Catch_pred(Year) = 0;
   
   for (Age=0;Age<=Nage;Age++){
     Propn_pred(Year,Age) = F(Year,Age) / Z(Year,Age) * N(Year,Age) *
                                 (1.0-mfexp(-Z(Year,Age)));
     Catch_pred(Year) += Weight(Age) * Propn_pred(Year,Age);
     Bio(Year) += Weight(Age) * S(Age) * N(Year,Age);
    }
    
   CPUE_pred(Year) = mfexp(logq) * Bio(Year); 
   Propn_pred(Year) /= sum(Propn_pred(Year));
   //  Propn_pred(Year) = Propn_pred(Year) /  sum(Propn_pred(Year));
  }

// ---------------------------------------------------------------------------- 

FUNCTION Likelihood
 int Year,Age;
 
 // Catch data
 // normal likelihood
 NLL1 = 0;
 for (Year=1;Year<=Nyear;Year++){
  NLL1 += square( (Catch(Year)-Catch_pred(Year))/Catch_pred(Year));
 }
 NLL1 = NLL1 / (2.0*square(SigCatch));
 // NLL1 /= (2.0*square(SigCatch));
 
 // CPUE data
 // lognormal likelihood
 NLL2 = 0;
 for (Year=1;Year<=Nyear;Year++){
  NLL2 += square( log(CPUE(Year)) - log(CPUE_pred(Year)) );
 }
 NLL2 = NLL2 / (2.0*square(SigCPUE)); 
 
 // Catch-at-age data
 // multinomial likelihood
 NLL3 = 0;
 for (Year=1;Year<=Nyear;Year++){
  for (Age=0;Age<=Nage;Age++){
   if (Propn(Year,Age) >0)
    NLL3 += Propn(Year,Age) * log(Propn_pred(Year,Age) / Propn(Year,Age));
  }
 }
 NLL3 = -1*Omega*NLL3;  
 
// ==========================================================================

REPORT_SECTION
 report << "jnll"  << " " << "nll_catch"  << " " << "nll_cpue"  << " " << "nll_caa"  << endl;
 report << "jnll nll_catch nll_cpue nll_caa" << endl; 
 report << objn << " " << NLL1 << " " << NLL2 << " " << NLL3 << endl;
 report << "Selex" << endl;
 report << S << endl;
 report << "F" << endl;
 report << mfexp(logF) << endl;

