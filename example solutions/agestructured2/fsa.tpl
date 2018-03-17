DATA_SECTION
  init_int minAge
  init_int maxAge
  init_int minYear
  init_int maxYear
  init_matrix catchNo(minAge,maxAge,minYear,maxYear)
  init_matrix stockMeanWeight(minAge,maxAge,minYear,maxYear)
  init_matrix propMature(minAge,maxAge,minYear,maxYear)
  init_matrix M(minAge,maxAge,minYear,maxYear)

  init_int minAgeS
  init_int maxAgeS
  init_int minYearS
  init_int maxYearS
  init_number surveyTime
  init_matrix survey(minAgeS,maxAgeS,minYearS,maxYearS)

  init_int rec_model;

  matrix logC(minAge,maxAge,minYear,maxYear)
  !! logC=log(catchNo);
  matrix logSurvey(minAgeS,maxAgeS,minYearS,maxYearS)
  !! logSurvey=log(survey);

PARAMETER_SECTION
  init_vector logN1Y(minAge,maxAge)
  init_vector logN1A(minYear+1,maxYear)
  init_vector logFY(minYear,maxYear)
  init_vector logFA(minAge,maxAge-3)
  init_vector logVarLogCatch(1,2)
  init_vector logQ(minAgeS,maxAgeS)
  init_vector logQnew(minAgeS,maxAgeS)
  init_number logVarLogSurvey
  init_bounded_number Steep(0.2,0.99);

  objective_function_value nll

  matrix F(minAge,maxAge,minYear,maxYear)
  matrix logN(minAge,maxAge,minYear,maxYear)
  matrix predLogC(minAge,maxAge,minYear,maxYear)
  matrix predLogSurvey(minAgeS,maxAgeS,minYearS,maxYearS)
  number ss
  vector logFAextend(minAge,maxAge)
  sdreport_vector ssb(minYear,maxYear)

  number alpha
  number beta

PROCEDURE_SECTION
  logFAextend.initialize();
  logFAextend(minAge,maxAge-3)=logFA;  

  F=outer_prod(exp(logFAextend),exp(logFY));

  logN.colfill(minYear,logN1Y);

  ssb.initialize();
  for(int a=minAge;a<=maxAge;a++){
    ssb(minYear) += exp(logN(a,minYear)) * propMature(a,minYear) * stockMeanWeight(a,minYear);
  }
  alpha = (ssb(minYear)*(1-Steep))/(4*Steep*exp(logN(minAge,minYear)));
  beta = (5*Steep-1)/(4*Steep*exp(logN(minAge,minYear)));
  for(int y=minYear+1; y<=maxYear; ++y){
    if(rec_model==1) logN(minAge,y) = ssb(y-1)/(alpha + beta*ssb(y-1));  //logN1A(y);
    if(rec_model==2) logN(minAge,y) = alpha * ssb(y-1) * mfexp(-beta * ssb(y-1));
    for(int a=minAge+1; a<=maxAge; ++a){
      logN(a,y)=logN(a-1,y-1)-F(a-1,y-1)-M(a-1,y-1);
      if(a==maxAge){
        logN(a,y)=log(exp(logN(a,y))+exp(logN(a,y-1)-F(a,y-1)-M(a,y-1)));
      }
    }
    for(int a=minAge;a<=maxAge;a++){
      ssb(y) += exp(logN(a,y)) * propMature(a,y) * stockMeanWeight(a,y);
    }
  }

  predLogC=log(F)-log(F+M)+log(1.0-exp(-F-M))+logN;

  nll = 0;
  for(int y=minYear; y<=maxYear; y++){
    for(int a=minAge; a<=maxAge; a++){
      if(a==minAge) ss = exp(logVarLogCatch(1));
      if(a > minAge) ss = exp(logVarLogCatch(2));
      nll += 0.5 * (log(2.0 * M_PI * ss) + square(logC(a,y) - predLogC(a,y))/ss);
    }
  }
  // ss=exp(logVarLogCatch);
  // int N=(maxYear-minYear+1)*(maxAge-minAge+1);
  // nll=0.5*(N*log(2.0*M_PI*ss)+sum(square(logC-predLogC))/ss);

  ss=exp(logVarLogSurvey);
  for(int y=minYearS; y<=maxYearS; ++y){
    for(int a=minAgeS; a<=maxAgeS; ++a){
      if(y<2000) predLogSurvey(a,y)=logQ(a)-(F(a,y)+M(a,y))*surveyTime+logN(a,y);
      if(y>=2000) predLogSurvey(a,y)=logQnew(a) - (F(a,y)+M(a,y)) * surveyTime + logN(a,y);
      nll+=0.5*(log(2.0*M_PI*ss)+square(logSurvey(a,y)-predLogSurvey(a,y))/ss);
    }
  }
  //ssb=colsum(elem_prod(elem_prod(exp(logN),propMature),stockMeanWeight));

REPORT_SECTION
  report << nll << endl; 
