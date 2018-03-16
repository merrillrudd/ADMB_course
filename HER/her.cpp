#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
  #include <admodel.h>
  #include <string.h>
  #include <time.h>
  #undef COUT
  #define COUT(object) cout << #object "\n" << setw(6) \
  << setprecision(3) << setfixed() << object << endl;
    
    
    
  
    
    
    
    
    
  
    
  
  
  
  
  //system("cp ham.rep run1.rep");
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <her.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
    int on = 0;
    retro_yrs = 0;
    rseed  = 0;
    b_simulation_flag = 0;
    if (ad_comm::argc > 1)
    {
      int on = 0;
      if ( (on=option_match(ad_comm::argc,ad_comm::argv,"-sim")) > -1 )
      {
        b_simulation_flag = 1;
        rseed = atoi(ad_comm::argv[on+1]);
      }
      if ( (on=option_match(ad_comm::argc,ad_comm::argv,"-retro")) > -1 )
      {
        retro_yrs = atoi(ad_comm::argv[on+1]);
        cout<<"|—————————————————————————————————————————————————|\n";
        cout<<"| Implementing Retrospective analysis             |\n";
        cout<<"|—————————————————————————————————————————————————|\n";
        cout<<"| Number of retrospective years = "<<retro_yrs<<endl;
      }
    }
  DEBUG_FLAG.allocate("DEBUG_FLAG");
  DataFile.allocate("DataFile");
  ControlFile.allocate("ControlFile");
 ad_comm::change_datafile_name(DataFile);
  dat_syr.allocate("dat_syr");
  dat_nyr.allocate("dat_nyr");
  mod_syr.allocate("mod_syr");
  mod_nyr.allocate("mod_nyr");
  sage.allocate("sage");
  nage.allocate("nage");
  rec_syr = mod_syr + sage;
  age.allocate(sage,nage);
 age.fill_seqadd(sage,1);
  nFecBlocks.allocate("nFecBlocks");
  nFecBlockYears.allocate(1,nFecBlocks,"nFecBlockYears");
  fec_slope.allocate(1,nFecBlocks,"fec_slope");
  fec_inter.allocate(1,nFecBlocks,"fec_inter");
  data_ct_raw.allocate(dat_syr,dat_nyr,1,3,"data_ct_raw");
  data_sp_waa.allocate(dat_syr,dat_nyr,sage-1,nage,"data_sp_waa");
  data_cm_waa.allocate(dat_syr,dat_nyr,sage-1,nage,"data_cm_waa");
  data_cm_comp.allocate(dat_syr,dat_nyr,sage-1,nage,"data_cm_comp");
  data_sp_comp.allocate(dat_syr,dat_nyr,sage-1,nage,"data_sp_comp");
  data_egg_dep.allocate(dat_syr,dat_nyr,1,3,"data_egg_dep");
  data_mileday.allocate(dat_syr,dat_nyr,1,3,"data_mileday");
  avg_sp_waa.allocate(sage,nage);
    int n = data_sp_waa.rowmax() - data_sp_waa.rowmin() + 1;
    avg_sp_waa = colsum(data_sp_waa)(sage,nage) / n;
  Eij.allocate(mod_syr,mod_nyr,sage,nage);
    int iyr = mod_syr;
    
    // for each fecundity block,
    // do the regression to calculate eggs
    // on the condition the year is less than the year in this block
    // iteratively go through and do calculation until the number of iyrs are in the fecundity block
    for(int h = 1; h <= nFecBlocks; h++){
      do{
        Eij(iyr) = 1.e-6 *
                (data_sp_waa(iyr)(sage,nage) * fec_slope(h) - fec_inter(h));
        iyr ++;
      }while(iyr <= nFecBlockYears(h));
    }
  dat_eof.allocate("dat_eof");
 if(dat_eof != 999){cout<<"Error reading data file, aborting."<<endl; exit(1);}
 if(dat_eof == 999){cout << "Data read correctly!" << endl;}
 ad_comm::change_datafile_name(ControlFile);
 n_theta = 6;
  theta_DM.allocate(1,n_theta,1,7,"theta_DM");
  theta_ival.allocate(1,n_theta);
  theta_lb.allocate(1,n_theta);
  theta_ub.allocate(1,n_theta);
  theta_phz.allocate(1,n_theta);
  theta_iprior.allocate(1,n_theta);
  theta_p1.allocate(1,n_theta);
  theta_p2.allocate(1,n_theta);
 theta_ival = column(theta_DM,1);
 theta_lb  = column(theta_DM,2);
 theta_ub  = column(theta_DM,3);
 theta_phz = ivector(column(theta_DM,4));
 theta_iprior = ivector(column(theta_DM,5));
 theta_p1 = column(theta_DM,6);
 theta_p2 = column(theta_DM,7);
  nMatBlocks.allocate("nMatBlocks");
  maturity_cont.allocate(1,nMatBlocks,1,4,"maturity_cont");
  mat_a50.allocate();
  mat_a95.allocate();
  mat_phz.allocate(1,nMatBlocks);
  nMatBlockYear.allocate(1,nMatBlocks);
    mat_a50 = column(maturity_cont,1);
    mat_a95 = column(maturity_cont,2);
    mat_phz = ivector(column(maturity_cont,3));
    nMatBlockYear = ivector(column(maturity_cont,4));
  mort_type.allocate("mort_type");
  mort_dev_phz.allocate("mort_dev_phz");
  nMortBlocks.allocate("nMortBlocks");
  nMortBlockYear.allocate(1,nMortBlocks,"nMortBlockYear");
 nSlxCols = 9;
  nSlxBlks.allocate("nSlxBlks");
  selex_cont.allocate(1,nSlxBlks,1,nSlxCols,"selex_cont");
  nSelType.allocate(1,nSlxBlks);
  nslx_phz.allocate(1,nSlxBlks);
  nslx_rows.allocate(1,nSlxBlks);
  nslx_cols.allocate(1,nSlxBlks);
  nslx_syr.allocate(1,nSlxBlks);
  nslx_nyr.allocate(1,nSlxBlks);
    // selectivity model (see control file for options)
    nSelType = ivector(column(selex_cont,2));
    // estimation phase for selectivity
    nslx_phz = ivector(column(selex_cont,7));
    // start of each selectivity time block
    nslx_syr = ivector(column(selex_cont,8));
    // end of each selectivity time block
    nslx_nyr = ivector(column(selex_cont,9));
    // determine dimensions for log_slx_pars ragged object.
    for(int h = 1; h <= nSlxBlks; h++){
      nslx_rows(h) = 1;
      switch(nSelType(h)){
        case 1: // logistic 2-parameters
          nslx_cols = int(2);
        break;
      }
    }
  nMiscCont.allocate("nMiscCont");
  dMiscCont.allocate(1,nMiscCont,"dMiscCont");
  data_catch.allocate(dat_syr,dat_nyr,1,3);
    // include catch scalar dMiscCont(1)
    data_catch = data_ct_raw;
    for( int i = dat_syr; i <= dat_nyr; i++ ) {
      data_catch(i,2) = dMiscCont(1) * data_ct_raw(i,2);
    }
    if(dMiscCont(2)) cout<<"Condition model on Ft"<<endl;
  ctl_eof.allocate("ctl_eof");
    if(ctl_eof != 999){
      cout<<"Error reading control file, aborting."<<ctl_eof<<endl; 
      exit(1);
    }
    if(ctl_eof == 999){
      cout << "Control file read correctly!" << endl;
    }
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  nll_total.allocate("nll_total");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::userfunction(void)
{
  nll_total =0.0;
  nll_total = 1;
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
