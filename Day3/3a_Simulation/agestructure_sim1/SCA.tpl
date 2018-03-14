//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
//Programer: Robert Ahrens													 
//Purpose:SCA simulator-assessment model.											 
//Notes: 				 
//							 
//																 
//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

DATA_SECTION
	int sim;
	int rseed;
	
	LOCAL_CALCS
		sim=0;
		rseed=0;
		int on,opt;
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-sim",opt))>-1)
		{
			sim=1;
			rseed=atoi(ad_comm::argv[on+1]);
		}
	END_CALCS

	init_adstring datafile;
	init_adstring ctlfile;
	//change to the new data file
	!!ad_comm::change_datafile_name(datafile);
	//for mcmc stuff
	init_int syr;//first year
	init_int eyr;//last year
	init_int nage;//number of ages
	init_number fwalpha;//ford-walfor alpha
	init_number fwrho;//ford-walford rho
	init_number wla;//W-L a
	init_number wlb;//W-L b
	init_vector yt(syr,eyr);//cpue
	init_vector ct(syr,eyr);//catch
	init_matrix pat(syr,eyr,1,nage);//proportions at age
	//!!cout<<pat<<endl;
	//!!ad_exit(1);
	init_number iro;// initial guess for for equilibrium unfished recruitment
	init_number icr;// initial guess for compensation ratio note steepness=cr/(cr+4)
	init_number irbar;// inital guess for average recuitment ofver the timeseries
	init_number ifbar;// inital guess for average fishing mortality rate
	init_number iahat;// inital guess for age at 50% vulberability to gear
	init_number ighat;//inital guess at how vunerability changes with age
	init_int eof;
	int iter;
	!!iter=0;

	LOCAL_CALCS
		if(eof!=999)
		{
			cout<<"Error reading data.\n Fix it."<<endl;
			ad_exit(1);
		}
	END_CALCS
	vector age(1,nage);//ages
	vector la(1,nage);//lengths
	vector wa(1,nage);//weights
	vector fa(1,nage);//fecundity
	number m;//natural norality
	number vbk;//
	number linf;//
	LOCAL_CALCS
		vbk=-log(fwrho);
		linf=fwalpha/(1-fwrho);
		m=1.2*vbk;
		age.fill_seqadd(1,1);
		la[1]=fwalpha;
		for (int i=2;i<=nage;i++)la[i]=fwalpha+fwrho*la[i-1];//ford walford growth
		wa=wla*pow(la,wlb);//legth to weight cm-kg
		//weight at age * maturity at age
		//fa=1./(1.+exp(-(age-log(3.)/vbk)/(0.1*log(3.)/vbk)));
		//fa=elem_prod(wa,fa);
		fa=elem_prod(wa,plogis(age,log(3.)/vbk,0.1*log(3.)/vbk));
	END_CALCS
	//load in data for simulations
	!!ad_comm::change_datafile_name(ctlfile);
	init_number simsig;
	init_number simtau;
	init_number simqo;
	init_number simao;
	init_number simro;
	init_number simcr;
	init_number simrbar;
	init_number simahat;
	init_number simghat;
	init_vector simF(syr,eyr);
	init_int eofc;
	LOCAL_CALCS
		if(eofc!=999)
		{
			cout<<"Error reading control file.\n Fix it."<<endl;
			ad_exit(1);
		}
	END_CALCS

PARAMETER_SECTION
	init_number log_ro(4);//Unfished recruitment
	init_number log_cr(5);//compensation ratio
	init_number log_rbar;//average recruitment over the time series
	init_number log_fbar;//average fishing mortality rate over the time series
	init_bounded_number ahat(0,nage);//Age at 50% vulnerability
	init_bounded_number ghat(0,5);//Shape of logistic vunerability
	init_number ptdev;//precision of total variance
	init_bounded_number rho(0,1);//Proportion of total variance that is process error
	init_bounded_dev_vector wt(syr-nage,eyr,-10.,10.,2);//Recruitment deviation from r_bar for all inital ages and all recruitments
	init_bounded_dev_vector log_ft_dev(syr,eyr,-10.,10.,3);//Fishing mortality deviations from F_bar
	objective_function_value nll;//Total negative log likelihood
	!!log_ro=log(iro);//set some inital values or use a pin file
	!!log_cr=log(icr);
	!!log_rbar=log(irbar);
	!!log_fbar=log(ifbar);
	!!rho=0.5;
	!!ptdev=log(1./0.08);
	!!ahat=iahat;
	!!ghat=ighat;
	number ro;//declare variable for back transformation
	number cr;
	number rbar;
	number fbar;
	number so;//additional variables for S-R fit to Beverton-Holt so=a beta=b
	number beta;
	number q;//y=qVB from exp(Z_bar)
	number fmsy;//MRP to be calculated at the end
	number msy;
	number bmsy;
	number sig;//process error
	number tau;//observation error
	sdreport_number tdev;

	vector va(1,nage);//vulnerability or selectivity
	vector ft(syr,eyr);//fishing mortality rates
	vector bt(syr,eyr+1);//Vulnerable  Biomass
	vector ct_hat(syr,eyr);// Estimates total catch
	vector ct_resid(syr,eyr);// Catch residual for likelihood
	vector yt_resid(syr,eyr);// Cpue residuals for likelihood
	vector rt(syr+1,eyr); // Predicted recruits using Beverton-Holt
	vector rt_resid(syr+1,eyr);//Deviations from Bevertin-Holt and predicted Recruitment each year

	matrix Nt(syr,eyr+1,1,nage);//Numbers at age
	matrix Zt(syr,eyr,1,nage);//Total mortality at age M+F*vul
	matrix pahat(syr,eyr,1,nage);//Predicted proprtions at age

PRELIMINARY_CALCS_SECTION
  if(sim)
  {
  	run_data_simulation();
  }

PROCEDURE_SECTION
	initialization();
	statedynamics();
	observation_model();
	stock_recruit_model();
	objective_function();

	if(mceval_phase())
	{ 
		forecast();
		get_CF(value(fmsy),value(msy),value(bmsy));
		mcmc_output();
	}
	if(last_phase())
	{
		forecast();
		get_CF(value(fmsy),value(msy),value(bmsy));
	} 

FUNCTION initialization
	va=plogis(age,ahat,ghat);//vulnerability
	//va=1./(1.+mfexp(-(age-ahat)/ghat));
	ft=mfexp(log_fbar+log_ft_dev);//Total fishing mortality
	Zt=m+outer_prod(ft,va);//Total mortality by age

FUNCTION statedynamics
	dvar_vector lxo=pow(mfexp(-m),age-1.);//unfished surviorship
	lxo(nage)/=(1.-mfexp(-m));//add the plus group
	Nt(syr,1)=mfexp(log_rbar+wt(syr-1));//initalize numbers at age 1 in the first year as rbar and deviation
	for(int j=2;j<=nage;j++) Nt(syr,j)=mfexp(log_rbar+wt(syr-j))*lxo(j);//initalize numbers at age in the first year as rbar and deviation and survivorship unfished
	for(int i=syr;i<=eyr;i++) //for (1 in sry:eyr)R equivalent
	{
		Nt(i+1,1)=mfexp(log_rbar+wt(i));//Recruitment each year as rbar and deviation
		Nt(i+1)(2,nage)=++elem_prod(Nt(i)(1,nage-1),mfexp(-Zt(i)(1,nage-1)));//N(t+1,a+1)=N(t,a)exp(-Z(t,a)) note ++ shifts 1-(nage-1) to 2-nage could also use the .shift function
		Nt(i+1,nage)+=Nt(i,nage)*mfexp(-Zt(i,nage));//Plus group calcuation += just adds 
	}
	bt=Nt*elem_prod(wa,va);//Vulnerable biomass using matrix multiplication matrix(t,a)*vector(a)=vector(t) 

FUNCTION observation_model
	dvar_matrix C(syr,eyr,1,nage);//catch at age
	dvar_matrix F(syr,eyr,1,nage);//F at age
	F=outer_prod(ft,va);//F at age as F(t)*vulnerability at age
	C=elem_prod(elem_div(F,Zt),elem_prod(1.-mfexp(-Zt),Nt));//Baranov catch equation F/Z*N(1-exp(-Z))
	ct_hat=C*wa;//Total annual catch as matrix(t,a)*vector(a)=vector(t)
	//catch residuals
	ct_resid=log(ct)-log(ct_hat);//catch residuals 
	//cpue residuals ala waters and ludwig 1994
	yt_resid=log(yt)-log(bt(syr,eyr));//Z values by year
	q=mfexp(mean(yt_resid));//q as exp(Z_bar)
	yt_resid-=mean(yt_resid);//Z residuals as Z values by year - Zbar
	//predicted proportions in the catch
	for(int i=syr;i<=eyr;i++)pahat(i)=C(i)/sum(C(i));// proportions at age

FUNCTION stock_recruit_model
	ro=mfexp(log_ro);//UNfished recruitment
	cr=mfexp(log_cr)+1.;//Compensation ratio
	dvar_vector lxo=pow(mfexp(-m),age-1.);//Survivorship
	lxo(nage)/=(1.-mfexp(-m));//Plus group
	dvariable phieo=lxo*fa;//unfished eggs per recruit 
	so=cr/phieo;//alpha as we know it relative max egg to recruit survival scaled as fecundity is scaled
	beta=(cr-1.)/(ro*phieo);//scale parameter of B-H
	dvar_vector sbt=Nt*fa; //NUmbers*fecundity matrix(t,a)*vector(a)=vector(t)
	dvar_vector nmr=so*sbt(syr,eyr-1); //for each year so*eggs
	dvar_vector den=(1.+beta*sbt(syr,eyr-1));//for each year 1+b*eggs
	dvar_vector tmp_rt=++elem_div(nmr,den);//for each year so*eggs/(1+b*eggs) with index shift so tmp_rt has index of (syr+1,eyr) to match rt 
	rt=column(Nt.sub(syr+1,eyr),1);//estimated recruitment from r_bar and deviations .sub extracts sub matrix of Nt and extract the first column 
	rt_resid=log(tmp_rt)-log(rt);//residuals as predicted -observed(estimated)

FUNCTION objective_function 
	dvar_vector nll_vec(1,4);//vector for likelihood components
	tdev=sqrt(1./mfexp(ptdev));//transform total precision to total variance
	sig=sqrt(rho*1./mfexp(ptdev));//standard deviation (process error) as a proprtion of the total variance
	tau=sqrt((1.-rho)*1./mfexp(ptdev));//standard deviation (observation error) as a proprtion of the total variance
	nll_vec.initialize();
	nll_vec(1)=dnorm(ct_resid,0.05);//catch residuals sd must be specified can not be estimated
	nll_vec(2)=dnorm(yt_resid,tau);//observation residuals from cpue time series
	nll_vec(3)=dnorm(rt_resid,sig);//process residuals from stock recruitment fit
	double tau2;
	nll_vec(4)=dmvlogistic(pat,pahat,tau2);//using multivariate logistic for proportions at age tau2 is returned from the function as the MLE estimate of the variance
	dvar_vector p_vec(1,5);// hey there are priors
	p_vec.initialize();// create p_vec
	dvariable h=cr/(4.+cr);//convert compensation ratio to steepness
	p_vec(1)=dbeta((h - 0.2)/0.8,2.0,2.0);// put a weak beta prior on steepness
	
	if(last_phase())
	{
		p_vec(3)=dnorm(wt,2);//week priors on the deviaitons
		p_vec(4)=dnorm(log_ft_dev,2);
	}
	else
	{
		p_vec(3)=100.*norm2(wt); //stong priors preventing big deviations
		p_vec(4)=100.*norm2(log_ft_dev);

	}
	p_vec(5)=dbeta(rho,50,50);// strong beta prior on the proportion of total variance as process errror
	nll=sum(nll_vec)+sum(p_vec);//total posterior as likelihoods and priors

FUNCTION mcmc_output
	if(iter==0)
	{
		ofstream ofs("refpar.mcmc");
		ofs<<"fmsy\t bmsy\t msy\t b/bmsy\t f/fmsy"<<endl;
		//ofs<<"r\t k\t q\t sig\t"<<endl;
	}
	iter++;
	double fratio=value(ft(eyr)/fmsy);
	double bratio=value(Nt(eyr)*wa/bmsy);
	ofstream ofs("refpar.mcmc",ios::app);
	ofs<<fmsy<<"\t"<<bmsy<<"\t"<<msy<<"\t"<<bratio<<"\t"<<fratio<<endl;

FUNCTION forecast

FUNCTION run_data_simulation
	random_number_generator rng(rseed);//create the random number generator object
	dmatrix C(syr,eyr,1,nage); //Catch at age
	dvector tmp(syr-nage,eyr);//anomalies
	dvector eps(syr,eyr);//anomalies
	dvector simqt(syr,eyr);//catchability over time
	tmp.fill_randn(rng);//fill with standard random normal(0,1) 
	eps.fill_randn(rng);
	wt=tmp*simsig;//recruitment anomalies 
	eps*=simtau;//observation anomlaies
	log_ro=log(simro);//simulation assumes B-H recruitment
	log_cr=log(simcr);
	log_rbar=log(simrbar);
	ahat=simahat;//selectivity
	ghat=simghat;
	ro=mfexp(log_ro);//back transform
	cr=mfexp(log_cr);
	rbar=mfexp(log_rbar);
	dvector lxo=pow(mfexp(-m),age-1.);//unfished survivoarship so simulation is starting at unfished 
	lxo(nage)/=(1.-mfexp(-m));//plus group
	double phieo=lxo*fa;//eggs per recruit
	double phibo=lxo*wa;//biomass per recruit
	double Bo=value(ro)*phibo;//Initial biomass
	so=cr/phieo;//B-H alpha
	beta=(cr-1.)/(ro*phieo);//B-H Beta
	//selectivity
	va=plogis(age,ahat,ghat);//logistic selectivity
	//Make some fish
	//initial numbers
	Nt(syr,1)=mfexp(log_rbar+wt(syr-1));//inital age 1 based on r_bar
	for(int j=2;j<=nage;j++)Nt(syr,j)=mfexp(log_rbar+wt(syr-j))*lxo(j);//inital number at age based on r_bar
	ft=simF;//F time series from control
	for(int i=syr;i<=eyr;i++)
	{
		dvector ba=value(elem_prod(Nt(i),wa));//biomass at age 
		Zt(i)=m+ft(i)*va;//Z at age
		//update numbers
		double sbt=value(Nt(i)*fa);//spawning biomass value is used because Nt is a dvar_matrix
		simqt(i)=simqo/(simao+(1.-simao)*sum(ba)/Bo); // q can be density dependent control from control file
		Nt(i+1,1)=so*sbt/(1.+beta*sbt)*mfexp(wt(i));// assume B-H recruitment 
		Nt(i+1)(2,nage)=++elem_prod(Nt(i)(1,nage-1),mfexp(-Zt(i)(1,nage-1)));//NUmbers at age update
		Nt(i+1,nage)+=Nt(i,nage)*mfexp(-Zt(i,nage));//PLus group
		dvector zttmp=value(Zt(i));//cfreate a temp vector for Z because Z is dvar_matrix and operations below are on dvectors
		C(i)=elem_prod(elem_div(value(ft(i)*va),zttmp),elem_prod(1.-mfexp(-zttmp),value(Nt(i))));//Baranov catch equations
		pat(i)=rmvlogistic(C(i)(1,nage),0.3,rseed+i);//samples proprtions at age from random multivariate logistic

	}
	ct=C*wa;// create catch ate age data
	dvar_vector wbar(syr,eyr);//create a mean weight vector just in case
	for (int i=syr;i<=eyr;i++) wbar(i)=ct(i)/sum(C(i)); // mean weight as biomass/numbers
	bt=Nt*elem_prod(wa,va); //biomass at age
	yt=elem_prod(simqt,value(elem_prod(bt(syr,eyr),mfexp(eps))));//cpue is measuring vulnerable biomass
	wt=0;
	//use this code to create data files for more basic methods.
	/*
	ofstream ofs("DD.dat");
	ofs<<"#syr"<<endl;
	ofs<<syr<<endl;
	ofs<<"#eyr"<<endl;
	ofs<<eyr<<endl;
	ofs<<"#surv"<<endl;
	ofs<<mfexp(-m)<<endl;
	ofs<<"#rho"<<endl;
	ofs<<0.73<<endl;
	ofs<<"#alpha"<<endl;
	ofs<<11<<endl;
	ofs<<"#agek"<<endl;
	ofs<<3<<endl;
	ofs<<"#nage"<<endl;	
	ofs<<nage<<endl;
	ofs<<"#wa"<<endl;
	ofs<<wa<<endl;
	ofs<<"#yt"<<endl;
	ofs<<yt<<endl;
	ofs<<"#ct"<<endl;
	ofs<<ct<<endl;
	ofs<<"#wbar"<<endl;
	ofs<<wbar<<endl;
	ofs<<"#eof"<<endl;
	ofs<<999<<endl;
	*/
	//cout<<yt<<endl;
	//cout<<""<<endl;
	//cout<<ct<<endl;
	//cout<<""<<endl;
	//cout<<pat<<endl;
	//ad_exit(1);



FUNCTION void calc_partials(const double& fe, double& phie, double& phif, double& phiq, double& dphif_df, double& dphiq_df, double& dRe_df)
	//Use this function to calculate the partial derivatives (Table 2 in Martell et al. 2008)
	//Arguments: fe=fishing rate
	int i;
	dvector lx=(pow(exp(-m),age-1.));
	lx(nage)/=(1.-exp(-(m)));
	dvector lz(1,nage);
	dvector za=value(m+fe*va);
	dvector sa=1.-exp(-za);
	dvector qa=elem_prod(elem_div(value(va),za),sa);
	double dlz_df=0;

	lz[1]=1.0; 
	dphiq_df=0; dphif_df=0;
	phie=(sum(elem_prod(lx,fa)));
	for(i=1; i<=nage; i++)
	{
		if(i>1) lz[i]=lz[i-1]*exp(-za[i-1]);
		if(i>1) dlz_df=dlz_df*exp(-za[i-1]) - lz[i-1]*value(va[i-1])*exp(-za[i-1]);
		if(i==nage){ //6/11/2007 added plus group.
			lz[i]/=(1.-mfexp(-za[i]));
			//dlz_df=dlz_df*mfexp(-za[i-1]) - lz[i-1]*va[i-1]*mfexp(-za[i-1])/(1.-mfexp(-za[i]))
			dlz_df=value(dlz_df/(1.-mfexp(-za[i]))
					-lz[i-1]*mfexp(-za[i-1])*va[i]*mfexp(-za[i])
			/((1.-mfexp(-za[i]))*(1.-mfexp(-za[i]))));
		}	
		dphif_df=dphif_df+(fa[i])*dlz_df;
		dphiq_df=dphiq_df+(wa[i]*qa[i]*dlz_df+(lz[i]*wa[i]*value(va[i]*va[i]))/za[i]*(exp(-za[i])-sa[i]/za[i]));
	}
	phif=sum(elem_prod(lz,(fa)));
	phiq=sum(elem_prod(elem_prod(lz,(wa)),qa));
	dRe_df=value(ro/(cr-1.))*phie/square(phif)*dphif_df;
	//dphif_df=sum(elem_prod(fa,dlz_df));
	//dvector t2=elem_div(elem_prod(elem_prod(lz,value(va)),wa),za);
	//dvector t3=exp(-za)-elem_div(sa,za);
	//dphiq_df=sum(elem_prod(elem_prod(wa,qa),dlz_df)+elem_prod(t2,t3));



FUNCTION void get_CF(double& fe, double& msy,double& bmsy)
	//This function uses Newton-Raphson method to iteratively solve for F*
	//Then calculates C* given F* (See eq 1.3 in Martell 2008 CJFAS)
	int iter;
	double dy,ddy,re;
	double phie,phif,phiq,dphif_df,dphiq_df,dRe_df;
	fe=(m);  //initial guess for Fmsy
	for(iter= 1; iter<=50; iter++)
	{
		calc_partials(fe,phie,phif,phiq,dphif_df,dphiq_df,dRe_df);
		re=value(ro*(cr-phie/phif)/(cr-1.));
		dy=re*phiq+fe*phiq*dRe_df+fe*re*dphiq_df;
		ddy=phiq*dRe_df+re*dphiq_df;
		//Newton update
		fe=fe-dy/ddy;
		if(sfabs(dy)<1.e-10)break;
		//cout<<"Fe dy\t"<<fe<<" "<<dy<<" "<<fe-dy/ddy<<endl;
	}
	msy=fe*re*phiq;
	bmsy=re*phif;
	//cout<<"Fe "<<fe<<endl;



REPORT_SECTION
	
	//double fmsy;
	//double msy;
	//double bmsy;
	//if(last_phase())
	//{
	//	get_CF(fmsy,msy,bmsy);
	//}
	REPORT(ro);
	REPORT(cr);
	REPORT(mfexp(log_rbar));
	REPORT(mfexp(log_rbar+wt));
	REPORT(mfexp(log_fbar));
	REPORT(mfexp(log_fbar+log_ft_dev));
	REPORT(ahat);
	REPORT(ghat);
	REPORT(ct);
	REPORT(ct_hat);
	REPORT(fmsy);
	REPORT(msy);
	REPORT(bmsy);

TOP_OF_MAIN_SECTION
	
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);


GLOBALS_SECTION
	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << endl;

	#include <admodel.h>
	#include <time.h>
	#include <contrib.h>//IF you have ADMB-11
        #include <stats.cxx>
        #include<statsLib.h>
        time_t start,finish;
	long hour,minute,second;
	double elapsed_time;

FINAL_SECTION
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;
	cout<<"*******************************************"<<endl;
	cout<<"--Start time: "<<ctime(&start)<<endl;
	cout<<"--Finish time: "<<ctime(&finish)<<endl;
	cout<<"--Runtime: ";
	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
	cout<<"*******************************************"<<endl;


