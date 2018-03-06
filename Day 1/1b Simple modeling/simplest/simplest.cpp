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
#include <simplest.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  nobs.allocate("nobs");
  Y.allocate(1,nobs,"Y");
  x.allocate(1,nobs,"x");
  eof.allocate("eof");
  if(eof!=999){
          cout << "Error reading data file." <<endl;
          ad_exit(1);
  }
          
  
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  a.allocate("a");
  b.allocate("b");
  pred_Y.allocate(1,nobs,"pred_Y");
  #ifndef NO_AD_INITIALIZE
    pred_Y.initialize();
  #endif
  sqdev.allocate(1,nobs,"sqdev");
  #ifndef NO_AD_INITIALIZE
    sqdev.initialize();
  #endif
  aa.allocate("aa");
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::userfunction(void)
{
  f =0.0;
  aa = a;
  pred_Y = a*x+b;
  sqdev = square(pred_Y-Y);
  f = sum(sqdev);
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
  report << "Observed" << endl;
  report << Y << endl;
  report << "Predicted" << endl;
  report << pred_Y << endl;
  report << "Parameter a" << endl;
  report << a << endl;
  report << "Parameter b" << endl;
  report << b << endl; 
  report << "NLL" << endl;
  report <<  f << endl;
  ofstream f1("Param_est");
  f1 << "Parameter a" << endl;
  f1 << a << endl;
  f1 << "Parameter b" << endl;
  f1 << b << endl;
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
