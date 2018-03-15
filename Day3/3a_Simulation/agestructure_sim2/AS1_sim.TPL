DATA_SECTION
 int sim;
 int rseed;
 LOCAL_CALCS
   sim = 0;
   rseed = 0;
   int on,opt;
   if((on=option_match(ad_comm::argc,ad_comm::argv,"-sim",opt))>-1){
     sim=1;
     rseed=atoi(ad_comm::argv[on+1]);
   }
 END_CALCS

 init_adstring datafile;
 init_adstring controlfile;
 
 !!ad_comm::change_datafile_name(datafile);
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

 !!ad_comm::change_datafile_name(controlfile);
 init_vector sim_logN(1,Nage);
 init_vector sim_logR(1,Nyear);
 init_number sim_Sel50;
 init_number sim_Sel95;
 init_vector sim_logF(1,Nyear);
 init_number sim_logq;
 init_number sim_sigmaProc;
 init_number sim_sigmaObs;


PARAMETER_SECTION
 init_vector logN1(1,Nage)                // Log of abundance at age year 1
 init_vector logR(1,Nyear)                // log of recruits in all years
 init_bounded_dev_vector recdevs(1,Nyear,-10.0,10.0,2); // recruitment deviations
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

 vector Catch_obs(1,Nyear);                      // Observed catch (weight);
 vector CPUE_obs(1,Nyear);                       // Observed CPUE (weight);
 matrix Propn_obs(1,Nyear,0,Nage);               // Observed catch-at-age
 
 vector Bio(1,Nyear);                             // vulnerable biomass
 
 number NLL1;                                   // Catch data
 number NLL2;                                   // CPUE data
 number NLL3;                                   // Proportion data
  objective_function_value objn;

// // ==========================================================================

PRELIMINARY_CALCS_SECTION
 // set starting values without simulation
 if(sim==0){
  logN1 = 1;
  logR = 1;
  Sel50 = 3;
  Sel95 = 4;
  logF = 1;
  logq = log(0.05);

  for(int Year=1;Year<=Nyear;Year++){
    Catch_obs(Year) = Catch(Year);
    CPUE_obs(Year) = CPUE(Year);
    for(int Age=0;Age<=Nage;Age++){
      Propn_obs(Year,Age) = Propn(Year,Age);
    }
  }
 }

 // run simulation mode
 if(sim){
   run_simulation();
 }

 
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



// ---------------------------------------------------------------------------- 
FUNCTION run_simulation
  random_number_generator rng(rseed);

  dvector recdevs(1,Nyear);          // recruitment deviations
  dvector obsdevs(1,Nyear);          // observation deviations
  recdevs.fill_randn(rng);           // fill devs with standard random normal(0,1)
  obsdevs.fill_randn(rng);
  recdevs *= sim_sigmaProc;
  obsdevs *= sim_sigmaObs;

  Sel50 = sim_Sel50;               // rewrite with sim values
  Sel95 = sim_Sel95;
  Select();

  // Compute the F matrix
  for (int Year=1;Year<=Nyear;Year++){
   for (int Age=0;Age<=Nage;Age++){
    F(Year,Age) = mfexp(logF(Year) + obsdevs(Year)) * S(Age);
   }
  }
  Z = F + Mval;
  cout << obsdevs << endl;

 // Insert the abundance from ages 1-Nage in the first year
  for (int Age=1;Age<=Nage;Age++){
   N(1,Age) = mfexp(logN1(Age));
  }
  // Insert the recruits age=0 for all years
  for (int Year=1;Year<=Nyear;Year++){
   N(Year,0) = mfexp(logR(Year) + recdevs(Year)); 
  }
 
  // Project the whole N matrix
  for (int Year=1;Year<Nyear;Year++){
   for (int Age=0;Age<Nage;Age++){
    if(Age<(Nage-1)) N(Year+1,Age+1) = N(Year,Age) * mfexp(-Z(Year,Age));
    if(Age==(Nage-1)) N(Year+1,Age+1) = N(Year,Age) * mfexp(-Z(Year,Age)) +  N(Year,Age+1) * mfexp(-Z(Year,Age+1));
   }
  }

  
  // Compute the predicted exploitable biomass, catch-at-age and catch
  for (int Year=1;Year<=Nyear;Year++){
    Bio(Year) = 0;
    Catch_obs(Year) = 0;
   
    for (int Age=0;Age<=Nage;Age++){
      Propn_obs(Year,Age) = F(Year,Age) / Z(Year,Age) * N(Year,Age) *
                                 (1.0-mfexp(-Z(Year,Age)));
      Catch_obs(Year) += Weight(Age) * Propn_obs(Year,Age);
      Bio(Year) += Weight(Age) * S(Age) * N(Year,Age);
     }
    
    CPUE_obs(Year) = mfexp(sim_logq) * Bio(Year); 
    Propn_obs(Year) /= sum(Propn_obs(Year));
   }
  
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
   if(Age<(Nage-1)) N(Year+1,Age+1) = N(Year,Age) * mfexp(-Z(Year,Age));
   if(Age==(Nage-1)) N(Year+1,Age+1) = N(Year,Age) * mfexp(-Z(Year,Age)) +  N(Year,Age+1) * mfexp(-Z(Year,Age+1));
  }
//  for (Age=0;Age<Nage;Age++){
//   N(Year+1,Age+1) = N(Year,Age) * mfexp(-Z(Year,Age));
//  }
//  N(Year+1,Nage) += N(Year,Nage) * mfexp(-Z(Year,Nage));
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
  NLL1 += square( (Catch_obs(Year)-Catch_pred(Year))/Catch_pred(Year));
 }
 NLL1 = NLL1 / (2.0*square(SigCatch));
 // NLL1 /= (2.0*square(SigCatch));
 
 // CPUE data
 // lognormal likelihood
 NLL2 = 0;
 for (Year=1;Year<=Nyear;Year++){
  NLL2 += square( log(CPUE_obs(Year)) - log(CPUE_pred(Year)) );
 }
 NLL2 = NLL2 / (2.0*square(SigCPUE)); 
 
 // Catch-at-age data
 // multinomial likelihood
 NLL3 = 0;
 for (Year=1;Year<=Nyear;Year++){
  for (Age=0;Age<=Nage;Age++){
   if (Propn_obs(Year,Age) >0)
    NLL3 += Propn_obs(Year,Age) * log(Propn_pred(Year,Age) / Propn_obs(Year,Age));
  }
 }
 NLL3 = -1*Omega*NLL3;  
 
// ==========================================================================

REPORT_SECTION
 report << "jnll"  << " " << "nll_catch"  << " " << "nll_cpue"  << " " << "nll_caa"  << endl;
 report << objn << " " << NLL1 << " " << NLL2 << " " << NLL3 << endl;
 report << "Selex" << endl;
 report << S << endl;
 report << "F" << endl;
 report << mfexp(logF) << endl;
 report << "Recruits" << endl;
 report << mfexp(logR) << endl;
 report << "VulBio" << endl;
 report << Bio << endl;

