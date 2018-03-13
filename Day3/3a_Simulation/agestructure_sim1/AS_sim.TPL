DATA_SECTION
// create a command line option
 //-sim rseed to simulate data
 
 int sim;
 int rseed;
 LOCAL_CALCS
   int on = 0;
   rseed = 0;
   sim = 0;
   if(ad_comm::argc > 1){
     int on = 0;
     if( (on=option_match(ad_comm::argc,ad_comm::argv,"-sim")) > -1){
       sim = 1;
       rseed = atoi(ad_comm::argv[on + 1]);
     }
   }
 END_CALCS

 init_adstring datafile
 init_adstring controlfile

 !! ad_comm::change_datafile_name(datafile);
// read in data file
 init_int Nyear;                                 // Number of years
 init_int Nage;                                  // Last age-class
 init_number Mval;                               // Natural mortality
 init_vector Weight(1,Nage);                     // Weight-at-age
 init_number logR0;
 init_number SigR;
 init_number SigF;
 init_number SigCatch;                           // CV for catch data
 init_number SigCPUE;                            // Sigma for log(CPUE)
 init_number Omega;                              // Weight for proportion data
 init_matrix CatchCPUE(1,Nyear,0,2);              // Catch and cpue data
 init_matrix Propn(1,Nyear,1,Nage);              // Proportion data
 init_int eof;
 vector Catch(1,Nyear);
 vector CPUE(1,Nyear);
 
       !! Catch = column(CatchCPUE,1);                 // Extract the catch data
       !! CPUE = column(CatchCPUE,2);                  // Extract the CPUE data

 LOCAL_CALCS
   if(eof != 999){
      cout<< "Data read incorrectly" << endl;
      ad_exit(1);
   }
 END_CALCS

 !! ad_comm::change_datafile_name(controlfile);
// read in control file
 init_int Rphz;                                 // Number of years
 init_int Sel50phz;                                  // Last age-class
 init_int Sel95phz;
 init_int Fphz;
 init_int qphz;
 init_int eofc;
 
 LOCAL_CALCS
   if(eofc != 999){
      cout<< "Control file read incorrectly" << endl;
      ad_exit(1);
   }
 END_CALCS


PARAMETER_SECTION
 init_vector logR(1,Nyear,Rphz);                                // Log of recruits
 init_bounded_number Sel50(0,Nage,Sel50phz);             // Age-at-50%-selectivity
 init_bounded_number Sel95(0,Nage,Sel95phz);             // Age-at-95%-selectivity
 init_vector logF(1,Nyear,Fphz);                     // Log of fishing mortality
 init_number logq(qphz);                               // Log of catchability

 vector Catch_obs(1,Nyear);
 vector CPUE_obs(1,Nyear);
 matrix Propn_obs(1,Nyear,1,Nage);

 matrix N(1,Nyear,1,Nage);                      // Numbers-at-age
 vector Nt(1,Nyear);
 vector S(1,Nage);                                // Selectivity-at-age
 matrix F(1,Nyear,1,Nage);                        // F-matrix
 matrix Z(1,Nyear,1,Nage);                        // Z-matrix
 vector R(1,Nyear);                               //recruitment
// vector SB(1,Nyear);                              // spawning biomass
// number SB0;                                      // equilibrium spawning biomass
 vector log_rec_devs(1,Nyear);                      // log recruitment deviations
 vector log_obs_devs(1,Nyear);
 vector Ftrue(1,Nyear);
 vector Rtrue(1,Nyear);
 
 vector Catch_pred(1,Nyear);                      // Predicted catch (weight);
 vector CPUE_pred(1,Nyear);                       // Predicted CPUE (weight);
 matrix Propn_pred(1,Nyear,1,Nage);               // Predicted catch-at-age
 vector VulBio(1,Nyear);                             // vulnerable biomass
 
 number Like1;                                   // Catch data
 number Like2;                                   // CPUE data
 number Like3;                                   // Proportion data
 objective_function_value objn;

 sdreport_number Bratio;                         // Relative depletion between current and initial year
 
// ==========================================================================
PRELIMINARY_CALCS_SECTION
  // rewrite obs
  Simulate(sim);

PROCEDURE_SECTION
 int Age, Year;
  
 // cout << sim << endl;
  //cout << Catch_obs << endl;
 // Set up the selectivity pattern
 Select();

 // bring sim flag back to zero before running EM
 sim = 0;
 // Project the model forward and compute various outputs
 Numbers();
 
 // observation model - predict catch, cpue, and catch-at-age
 Predict_Observations();
 
 // Compute the likelihood
 Likelihood();

 objn = Like1 + Like2 + Like3;

 // Statistic to summarize
 Bratio = VulBio(Nyear)/VulBio(1);

FUNCTION void Simulate(const int& sim)
  int Year, Age;

  // generate random normal deviates for recruitment and fishing  processes
  random_number_generator rng(rseed);
  dvector epsilon_obs_devs(1,Nyear);
  dvector epsilon_rec_devs(1,Nyear);

  epsilon_rec_devs.fill_randn(rng);
  epsilon_obs_devs.fill_randn(rng);

  if(sim == 0){
    SigR = 0;
    SigF = 0;
  }

  log_rec_devs = dvar_vector(epsilon_rec_devs * SigR - 0.5 * square(SigR));
  log_obs_devs = dvar_vector(epsilon_obs_devs * SigF - 0.5 * square(SigF));

  Select();
  Numbers();
  Predict_Observations();

  Ftrue = exp(logF + log_obs_devs);
  Rtrue = exp(logR + log_rec_devs);
  
  // rewrite observed catch with observation error
  // rewrite CPUE and Propn with process error
  if(sim != 0){
  for(Year=1; Year<=Nyear; Year++){
    Catch_obs(Year) = Catch_pred(Year);
    CPUE_obs(Year) = CPUE_pred(Year);

    for(Age=1; Age<=Nage; Age++){
      Propn_obs(Year,Age) = Propn_pred(Year,Age);
    }
  }
  }
  if(sim == 0){
  for(Year=1; Year<=Nyear; Year++){
    Catch_obs(Year) = Catch(Year);
    CPUE_obs(Year) = CPUE(Year);

    for(Age=1; Age<=Nage; Age++){
      Propn_obs(Year,Age) = Propn(Year,Age);
    }
  }
  }
  

  ofstream simfile("AS_sim.sim");
  simfile << "F" << endl;
  simfile << Ftrue  << endl;
  simfile << "R" << endl;
  simfile << Rtrue << endl;
  simfile << "Nt" << endl;
  simfile << Nt << endl;
  simfile << "VulBio" << endl;
  simfile << VulBio << endl;
  simfile << "Catch_obs" << endl;
  simfile << Catch_obs << endl;
  simfile << "CPUE_obs" << endl;
  simfile << CPUE_obs << endl;
  simfile << "Propn_obs" << endl;
  simfile << Propn_obs << endl;
  simfile << "Selex" << endl;
  simfile << S << endl;
  simfile << "log_R_devs" << endl;
  simfile << log_rec_devs << endl;
  simfile << "log_F_devs" << endl;
  simfile << log_obs_devs << endl;
 
// ---------------------------------------------------------------------------- 

FUNCTION Select
 int Age;
 
 for (Age=1;Age<=Nage;Age++){
   S(Age) = 1.0 / (1 + exp(-log(19) * ((Age-1)-Sel50) / (Sel95-Sel50)));
  }

// ---------------------------------------------------------------------------- 

FUNCTION Numbers
 int Age,Year;


 // Compute the F matrix
 for (Year=1;Year<=Nyear;Year++){
  for (Age=1;Age<=Nage;Age++){
   if(sim!=0) F(Year,Age) = mfexp(logF(Year) + log_obs_devs(Year)) * S(Age);
   if(sim==0) F(Year,Age) = mfexp(logF(Year)) * S(Age);
  }
 }
 Z = F + Mval;  

 // Clear the N matrix
 N.initialize();
 Nt.initialize();
 
 // Insert the recruitments
 if(sim!=0) R(1) = exp(logR(1) + log_rec_devs(1));
 if(sim==0) R(1) = exp(logR(1));
 N(1,1) = R(1);
 for (Age=2;Age<=Nage;Age++){
  N(1,Age) = N(1,Age-1) * exp(-Z(1,Age-1));
  Nt(1) += N(1,Age);
 }
 //SB0 = 0;
// SB.initialize();
 //for(Age=1; Age<=Nage; Age++){
  // SB0 += exp(logR0) * Maturity(Age) * Weight(Age);
   //SB(1) += N(1,Age) * Maturity(Age) * Weight(Age);
 //}
 
 // Project the whole N matrix
 for (Year=2;Year<=Nyear;Year++){
 // R(Year) = ((4 * steep * exp(logR0) * SB(Year-1)) / (SB0 * (1-steep) + SB(Year-1) * (5*steep - 1))) * exp(log_rec_devs(Year));
  if(sim!=0) R(Year) = exp(logR(Year) + log_rec_devs(Year));
  if(sim==0) R(Year) = exp(logR(Year));
  for (Age=1;Age<=Nage;Age++){
    if(Age == 1) N(Year,Age) = R(Year);
    if(Age > 1 & Age < Nage) N(Year,Age) = N(Year-1,Age-1) * mfexp(-Z(Year-1,Age-1));
    if(Age == Nage) N(Year,Age) = N(Year-1,Age-1) * exp(-Z(Year-1,Age-1)) + N(Year-1,Age) * exp(-Z(Year-1,Age));
    if(Age > 1) Nt(Year) += N(Year,Age);
  }
 }

FUNCTION Predict_Observations
 int Year, Age;
 
 // Compute the predicted exploitable biomass, catch-at-age and catch
 Propn_pred.initialize();
 for (Year=1;Year<=Nyear;Year++){
   VulBio(Year) = 0;
   Catch_pred(Year) = 0;
   for (Age=1;Age<=Nage;Age++){
     Propn_pred(Year,Age) = F(Year,Age) / Z(Year,Age) * N(Year,Age) * (1.0-mfexp(-Z(Year,Age)));
     Catch_pred(Year) += Weight(Age) * Propn_pred(Year,Age);
     VulBio(Year) += Weight(Age) * S(Age) * N(Year,Age); // * mfexp(-Z(Year,Age)/2.0);
    } 
   CPUE_pred(Year) = exp(logq) * VulBio(Year); 
   Propn_pred(Year) /= sum(Propn_pred(Year));
  }

// ---------------------------------------------------------------------------- 

FUNCTION Likelihood
 int Year,Age;
 
 // Catch data
 // normal likelihood
 Like1 = 0;
 for (Year=1;Year<=Nyear;Year++){
  Like1 += square( (Catch_obs(Year)-Catch_pred(Year))/Catch_pred(Year));
 }
 Like1 = Like1 / (2.0*square(SigCatch));
 
 // CPUE data
 // lognormal likelihood
 Like2 = 0;
 for (Year=1;Year<=Nyear;Year++){
  Like2 += square( log(CPUE_obs(Year)) - log(CPUE_pred(Year)) );
 }
 Like2 = Like2 / (2.0*square(SigCPUE)); 
 
 // Catch-at-age data
 // multinomial likelihood
 Like3 = 0;
 for (Year=1;Year<=Nyear;Year++){
  for (Age=1;Age<=Nage;Age++){
   if (Propn_obs(Year,Age) >0)
    Like3 += Propn_obs(Year,Age) * log(Propn_pred(Year,Age) / Propn_obs(Year,Age));
  }
 }
 Like3 = -1*Omega*Like3;  
 
// ==========================================================================

REPORT_SECTION
 report << "jnll"  << " " << "nll_catch"  << " " << "nll_cpue"  << " " << "nll_caa"  << endl;
 report << objn << " " << Like1 << " " << Like2 << " " << Like3 << endl;
 report << "Selex" << endl;
 report << S << endl;
 report << "F" << endl;
 report << exp(logF) << endl;
 report << "Ftrue" << endl;
 report << Ftrue << endl;
 report << "Recruits" << endl;
 report << R << endl;
 report << "Rtrue" << endl;
 report << Rtrue << endl;
 report << "VulBio" << endl;
 report << VulBio << endl;
 report << "Nt" << endl;
 report << Nt << endl;
 report << "Catch_obs" << endl;
 report << Catch_obs << endl;
 report << "Catch_pred" << endl;
 report << Catch_pred << endl;
 report << "CPUE_obs" << endl;
 report << CPUE_obs << endl;
 report << "CPUE_pred" << endl;
 report << CPUE_pred << endl;
 report << "Propn_obs" << endl;
 report << Propn_obs << endl;
 report << "Propn_pred" << endl;
 report << Propn_pred << endl;
 report << "Weight" << endl;
 report << Weight << endl;
 

