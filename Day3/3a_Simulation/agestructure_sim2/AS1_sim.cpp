#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <AS1_sim.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
   sim = 0;
   rseed = 0;
   int on,opt;
   if((on=option_match(ad_comm::argc,ad_comm::argv,"-sim",opt))>-1){
     sim=1;
     rseed=atoi(ad_comm::argv[on+1]);
   }
  datafile.allocate("datafile");
  controlfile.allocate("controlfile");
ad_comm::change_datafile_name(datafile);
  Nyear.allocate("Nyear");
  Nage.allocate("Nage");
  Mval.allocate("Mval");
  Weight.allocate(0,Nage,"Weight");
  SigCatch.allocate("SigCatch");
  SigCPUE.allocate("SigCPUE");
  Omega.allocate("Omega");
  CatchCPUE.allocate(1,Nyear,0,2,"CatchCPUE");
  Propn.allocate(1,Nyear,-1,Nage,"Propn");
  Catch.allocate(1,Nyear);
  CPUE.allocate(1,Nyear);
 Catch = column(CatchCPUE,1);                 // Extract the catch data
 CPUE = column(CatchCPUE,2);                  // Extract the CPUE data
ad_comm::change_datafile_name(controlfile);
  sim_logN.allocate(1,Nage,"sim_logN");
  sim_logR.allocate(1,Nyear,"sim_logR");
  sim_Sel50.allocate("sim_Sel50");
  sim_Sel95.allocate("sim_Sel95");
  sim_logF.allocate(1,Nyear,"sim_logF");
  sim_logq.allocate("sim_logq");
  sim_sigmaProc.allocate("sim_sigmaProc");
  sim_sigmaObs.allocate("sim_sigmaObs");
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  logN1.allocate(1,Nage,"logN1");
  logR.allocate(1,Nyear,"logR");
  recdevs.allocate(1,Nyear,-10.0,10.0,2,"recdevs");
  Sel50.allocate(0,Nage,1,"Sel50");
  Sel95.allocate(0,Nage,1,"Sel95");
  logF.allocate(1,Nyear,1,"logF");
  logq.allocate(1,"logq");
  N.allocate(1,Nyear,0,Nage,"N");
  #ifndef NO_AD_INITIALIZE
    N.initialize();
  #endif
  S.allocate(0,Nage,"S");
  #ifndef NO_AD_INITIALIZE
    S.initialize();
  #endif
  F.allocate(1,Nyear,0,Nage,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  Z.allocate(1,Nyear,0,Nage,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  Catch_pred.allocate(1,Nyear,"Catch_pred");
  #ifndef NO_AD_INITIALIZE
    Catch_pred.initialize();
  #endif
  CPUE_pred.allocate(1,Nyear,"CPUE_pred");
  #ifndef NO_AD_INITIALIZE
    CPUE_pred.initialize();
  #endif
  Propn_pred.allocate(1,Nyear,0,Nage,"Propn_pred");
  #ifndef NO_AD_INITIALIZE
    Propn_pred.initialize();
  #endif
  Catch_obs.allocate(1,Nyear,"Catch_obs");
  #ifndef NO_AD_INITIALIZE
    Catch_obs.initialize();
  #endif
  CPUE_obs.allocate(1,Nyear,"CPUE_obs");
  #ifndef NO_AD_INITIALIZE
    CPUE_obs.initialize();
  #endif
  Propn_obs.allocate(1,Nyear,0,Nage,"Propn_obs");
  #ifndef NO_AD_INITIALIZE
    Propn_obs.initialize();
  #endif
  Bio.allocate(1,Nyear,"Bio");
  #ifndef NO_AD_INITIALIZE
    Bio.initialize();
  #endif
  NLL1.allocate("NLL1");
  #ifndef NO_AD_INITIALIZE
  NLL1.initialize();
  #endif
  NLL2.allocate("NLL2");
  #ifndef NO_AD_INITIALIZE
  NLL2.initialize();
  #endif
  NLL3.allocate("NLL3");
  #ifndef NO_AD_INITIALIZE
  NLL3.initialize();
  #endif
  objn.allocate("objn");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
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
 
}

void model_parameters::userfunction(void)
{
  objn =0.0;
 // Set up the selectivity pattern
 Select();
 // Project the model forward and compute various outputs
 Numbers();
 // observation model - predict catch, cpue, and catch-at-age
 Predict_Observations();
 // Compute the likelihood
 Likelihood();
 objn = NLL1 + NLL2 + NLL3;
}

void model_parameters::run_simulation(void)
{
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
  dvar_vector Ftrue(1,Nyear);
  dvar_vector Rtrue(1,Nyear);
  for(int Year=1;Year<=Nyear;Year++){
    Ftrue(Year) = mfexp(logF(Year)) * mfexp(obsdevs(Year));
    Rtrue(Year) = mfexp(logR(Year)) * mfexp(recdevs(Year));
  }
  // Compute the F matrix
  for (int Year=1;Year<=Nyear;Year++){
   for (int Age=0;Age<=Nage;Age++){
    F(Year,Age) = Ftrue(Year) * S(Age);
   }
  }
  Z = F + Mval;
 // Insert the abundance from ages 1-Nage in the first year
  for (int Age=1;Age<=Nage;Age++){
   N(1,Age) = mfexp(logN1(Age));
  }
  // Insert the recruits age=0 for all years
  for (int Year=1;Year<=Nyear;Year++){
   N(Year,0) = Rtrue(Year); 
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
   ofstream sim("AS1.sim");
   sim << "Ftrue" << endl;
   sim << Ftrue << endl;
   sim << "Rtrue" << endl;
   sim << Rtrue << endl; 
}

void model_parameters::Select(void)
{
 int Age;
 for (Age=0;Age<=Nage;Age++){
   S(Age) = 1.0 / (1 + exp(-log(19) * (Age-Sel50) / (Sel95-Sel50)));
  }
}

void model_parameters::Numbers(void)
{
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
 }
}

void model_parameters::Predict_Observations(void)
{
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
}

void model_parameters::Likelihood(void)
{
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
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
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
 report << "Catch_obs" << endl;
 report << Catch_obs << endl;
 report << "Catch_pred" << endl;
 report << Catch_pred << endl;
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
