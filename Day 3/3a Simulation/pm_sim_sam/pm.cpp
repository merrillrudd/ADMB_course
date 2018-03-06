#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	#define REPORT(object) \
	report << #object "\n" \
	<< object << endl;
	#include <admodel.h>
	#include <time.h>
	#include <stats.cxx>
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;
	
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <pm.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
		sim=0; //default 0 read in data
		rseed=0; //default rseed value
		int on,opt;
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-sim",opt))>-1) //option match is a true false
		{
			sim=1;
			rseed=atoi(ad_comm::argv[on+1]);
		}
  datafile.allocate("datafile");
  ctlfile.allocate("ctlfile");
 ad_comm::change_datafile_name(datafile);
  syr.allocate("syr");
  nyr.allocate("nyr");
  nage.allocate("nage");
  linf.allocate("linf");
  vonbk.allocate("vonbk");
  a.allocate("a");
  b.allocate("b");
  iyr.allocate(syr,nyr,"iyr");
  yt.allocate(syr,nyr,"yt");
  ct.allocate(syr,nyr,"ct");
  A.allocate(syr,nyr,1,nage,"A");
  eof.allocate("eof");
		if(eof != 999)
		{
			cout<<"Error reading data.\n Fix it."<<endl;
			exit(1);
		}
 ad_comm::change_datafile_name(ctlfile);
  npar.allocate("npar");
  ipar.allocate(1,npar,"ipar");
  eofc.allocate("eofc");
		if(eofc != 999)
		{
			cout<<"Error reading control file.\n Fix it."<<endl;
			exit(1);
		}
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  log_k.allocate("log_k");
  r.allocate(0.001,1.000,"r");
  log_q.allocate("log_q");
  log_isig.allocate(3,"log_isig");
  log_itau.allocate(3,"log_itau");
  xt.allocate(syr,nyr,-15.,15.,2,"xt");
		log_k = log(ipar(1));
		r = ipar(2);
		log_q = log(ipar(3));
		log_isig=log(1./ipar(4));
		log_itau=log(1./ipar(5));
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  k.allocate("k");
  #ifndef NO_AD_INITIALIZE
  k.initialize();
  #endif
  q.allocate("q");
  #ifndef NO_AD_INITIALIZE
  q.initialize();
  #endif
  fpen.allocate("fpen");
  #ifndef NO_AD_INITIALIZE
  fpen.initialize();
  #endif
  sig.allocate("sig");
  #ifndef NO_AD_INITIALIZE
  sig.initialize();
  #endif
  tau.allocate("tau");
  #ifndef NO_AD_INITIALIZE
  tau.initialize();
  #endif
  pt.allocate(syr,nyr+1,"pt");
  #ifndef NO_AD_INITIALIZE
    pt.initialize();
  #endif
  it.allocate(syr,nyr,"it");
  #ifndef NO_AD_INITIALIZE
    it.initialize();
  #endif
  nu.allocate(syr,nyr,"nu");
  #ifndef NO_AD_INITIALIZE
    nu.initialize();
  #endif
  depletion.allocate("depletion");
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
	if(sim)
	{
		run_simulation_model();
	}
	
}

void model_parameters::userfunction(void)
{
  f =0.0;
	//*** Main function calls  ***//
	// 1) population_dynamics
		population_dynamics();
	// 2) observation_model
		observation_model();
	// 3) calc_objective_function
		calc_objective_function();
	//****************************//
}

void model_parameters::population_dynamics(void)
{
	{
	//transorm parameters
	k = mfexp(log_k);
	//initialize population at unfished
	pt(syr) = 1.0;
	fpen = 0.0;
	int i;
	for(i=syr; i<=nyr; i++)
	{
		pt(i+1)=posfun((pt(i)+r*pt(i)*(1.-pt(i))-ct(i)/k)*mfexp(xt(i)),0.001,fpen);
	}
	depletion = pt(nyr);
	//cout<<pt<<endl;
	//exit(1);
	//TODO come back and add process errors
	}
}

void model_parameters::observation_model(void)
{
	{
		q = mfexp(log_q);
		it = q*k*pt(iyr); //NB iyr is an ivector
		nu = log(yt) - log(it);
	}
}

void model_parameters::calc_objective_function(void)
{
	{
		dvar_vector likevec(1,8);
		//likelihoods
		//obs residuals
		tau = sqrt(1./mfexp(log_itau));
		likevec(1) = dnorm(nu,tau);
		// process residuals
		sig=sqrt(1./mfexp(log_isig));
		likevec(2)=dnorm(xt,sig);
		//penlty for 0 biomass
		likevec(3)=1.e5*fpen;
		//priors
		//lognormal priors for r,k
		likevec(4)=dlnorm(k,8,0.25);
		likevec(5)=dlnorm(r,-1.38,0.51);
		//uniform prior on q
		likevec(6)=-log(q);
		likevec(7)=dgamma(1.0/(tau*tau),1.71,0.0086);
		likevec(8)=dgamma(1.0/(sig*sig),3.79,0.0102);;
		// total objective function
		f=sum(likevec);
		//add penalties
		if(fpen>0)cout<<"Fpen = "<<fpen<<endl;
	}
}

void model_parameters::run_simulation_model(void)
{
	random_number_generator rng(rseed);
	dvector tmp_nu(syr,nyr);
	dvector tmp_xt(syr,nyr);
	tmp_nu.fill_randn(rng);
	tmp_xt.fill_randn(rng);
	tau = sqrt(1./mfexp(log_itau));
	sig = sqrt(1./mfexp(log_isig));
	xt=tmp_xt*sig;
	population_dynamics();
	observation_model();
	nu=tmp_nu*tau;
	yt=value(elem_prod(it,exp(nu)));
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
	REPORT(fpen);
	REPORT(pt);
	//extracting values from dvar_vectors
	//placing into a dvector.
	dvector bt = value(pt(iyr)*k);
	REPORT(bt);
	REPORT(iyr);
	REPORT(nu);
	REPORT(sig);
	REPORT(tau);
	REPORT(r);
	REPORT(k);
	REPORT(q);
	REPORT(ct);
	REPORT(yt);
	REPORT(it);
	double fmsy=value(r/2.);
	double msy=value(r*k/4.);
	double tac=value(fmsy*pt(nyr+1)*k);
	REPORT(fmsy);
	REPORT(msy);
	REPORT(tac);
}

void model_parameters::final_calcs()
{
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;
	cout<<endl<<endl<<"*******************************************"<<endl;
	cout<<"--Start time: "<<ctime(&start)<<endl;
	cout<<"--Finish time: "<<ctime(&finish)<<endl;
	cout<<"--Runtime: ";
	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
	cout<<"*******************************************"<<endl;
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

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
