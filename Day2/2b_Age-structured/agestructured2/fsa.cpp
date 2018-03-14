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
#include <fsa.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  minAge.allocate("minAge");
  maxAge.allocate("maxAge");
  minYear.allocate("minYear");
  maxYear.allocate("maxYear");
  catchNo.allocate(minAge,maxAge,minYear,maxYear,"catchNo");
  stockMeanWeight.allocate(minAge,maxAge,minYear,maxYear,"stockMeanWeight");
  propMature.allocate(minAge,maxAge,minYear,maxYear,"propMature");
  M.allocate(minAge,maxAge,minYear,maxYear,"M");
  minAgeS.allocate("minAgeS");
  maxAgeS.allocate("maxAgeS");
  minYearS.allocate("minYearS");
  maxYearS.allocate("maxYearS");
  surveyTime.allocate("surveyTime");
  survey.allocate(minAgeS,maxAgeS,minYearS,maxYearS,"survey");
  logC.allocate(minAge,maxAge,minYear,maxYear);
 logC=log(catchNo);
  logSurvey.allocate(minAgeS,maxAgeS,minYearS,maxYearS);
 logSurvey=log(survey);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  logN1Y.allocate(minAge,maxAge,"logN1Y");
  logN1A.allocate(minYear+1,maxYear,"logN1A");
  logFY.allocate(minYear,maxYear,"logFY");
  logFA.allocate(minAge,maxAge-3,"logFA");
  logVarLogCatch.allocate("logVarLogCatch");
  logQ.allocate(minAgeS,maxAgeS,"logQ");
  logVarLogSurvey.allocate("logVarLogSurvey");
  nll.allocate("nll");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  F.allocate(minAge,maxAge,minYear,maxYear,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  logN.allocate(minAge,maxAge,minYear,maxYear,"logN");
  #ifndef NO_AD_INITIALIZE
    logN.initialize();
  #endif
  predLogC.allocate(minAge,maxAge,minYear,maxYear,"predLogC");
  #ifndef NO_AD_INITIALIZE
    predLogC.initialize();
  #endif
  predLogSurvey.allocate(minAgeS,maxAgeS,minYearS,maxYearS,"predLogSurvey");
  #ifndef NO_AD_INITIALIZE
    predLogSurvey.initialize();
  #endif
  ss.allocate("ss");
  #ifndef NO_AD_INITIALIZE
  ss.initialize();
  #endif
  logFAextend.allocate(minAge,maxAge,"logFAextend");
  #ifndef NO_AD_INITIALIZE
    logFAextend.initialize();
  #endif
  ssb.allocate(minYear,maxYear,"ssb");
  fbar24.allocate(minYear,maxYear,"fbar24");
}

void model_parameters::userfunction(void)
{
  nll =0.0;
  logFAextend.initialize();
  logFAextend(minAge,maxAge-3)=logFA;  
  F=outer_prod(exp(logFAextend),exp(logFY));
  cout << F << endl;
  nll = 1;
  // logN.colfill(minYear,logN1Y);
  // for(int y=minYear+1; y<=maxYear; ++y){
  //   logN(minAge,y)=logN1A(y);
  //   for(int a=minAge+1; a<=maxAge; ++a){
  //     logN(a,y)=logN(a-1,y-1)-F(a-1,y-1)-M(a-1,y-1);
  //   }
  // }
  // predLogC=log(F)-log(F+M)+log(1.0-exp(-F-M))+logN;
  // ss=exp(logVarLogCatch);
  // int N=(maxYear-minYear+1)*(maxAge-minAge+1);
  // nll=0.5*(N*log(2.0*M_PI*ss)+sum(square(logC-predLogC))/ss);
  // ss=exp(logVarLogSurvey);
  // for(int y=minYearS; y<=maxYearS; ++y){
  //   for(int a=minAgeS; a<=maxAgeS; ++a){
  //     predLogSurvey(a,y)=logQ(a)-(F(a,y)+M(a,y))*surveyTime+logN(a,y);
  //     nll+=0.5*(log(2.0*M_PI*ss)+square(logSurvey(a,y)-predLogSurvey(a,y))/ss);
  //   }
  // }
  // ssb=colsum(elem_prod(elem_prod(exp(logN),propMature),stockMeanWeight));
  // fbar24=colsum(F.sub(2,4))/3.0;
}

void model_parameters::preliminary_calculations(void){
#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::report(const dvector& gradients){}

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
