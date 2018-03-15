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

  matrix logC(minAge,maxAge,minYear,maxYear)
  !! logC=log(catchNo);
  matrix logSurvey(minAgeS,maxAgeS,minYearS,maxYearS)
  !! logSurvey=log(survey);

PARAMETER_SECTION
  init_vector logN1Y(minAge,maxAge)       // number of individuals year 1
  init_vector logN1A(minYear+1,maxYear)   // number of individauls first age
  init_vector logFY(minYear,maxYear)      // fishing mortality each year
  init_vector logFA(minAge,maxAge-3)      // fishing mortality for most ages
  init_number logVarLogCatch
  init_vector logQ(minAgeS,maxAgeS)
  init_number logVarLogSurvey

  objective_function_value nll

  matrix F(minAge,maxAge,minYear,maxYear)
  matrix logN(minAge,maxAge,minYear,maxYear)
  matrix predLogC(minAge,maxAge,minYear,maxYear)
  matrix predLogSurvey(minAgeS,maxAgeS,minYearS,maxYearS)
  number ss
  vector logFAextend(minAge,maxAge)
  sdreport_vector ssb(minYear,maxYear)

PROCEDURE_SECTION
  logFAextend.initialize();
  logFAextend(minAge,maxAge-3)=logFA;  

  F=outer_prod(exp(logFAextend),exp(logFY));

  // initialize N matrix for first year
  logN.colfill(minYear,logN1Y);
  // for(int a=minAge;a<=maxAge;a++){
  //   logN(a,minYear) = logN1Y(a);
  // }


  for(int y=minYear+1; y<=maxYear; ++y){
    logN(minAge,y)=logN1A(y);
    for(int a=minAge+1; a<=maxAge; ++a){
      logN(a,y)=logN(a-1,y-1)-F(a-1,y-1)-M(a-1,y-1);
    }
  }

  predLogC=log(F)-log(F+M)+log(1.0-exp(-F-M))+logN;

  ss=exp(logVarLogCatch);
  int N=(maxYear-minYear+1)*(maxAge-minAge+1);
  nll=0.5*(N*log(2.0*M_PI*ss)+sum(square(logC-predLogC))/ss);

  ss=exp(logVarLogSurvey);
  for(int y=minYearS; y<=maxYearS; ++y){
    for(int a=minAgeS; a<=maxAgeS; ++a){
      predLogSurvey(a,y)=logQ(a)-(F(a,y)+M(a,y))*surveyTime+logN(a,y);
      nll+=0.5*(log(2.0*M_PI*ss)+square(logSurvey(a,y)-predLogSurvey(a,y))/ss);
    }
  }

   // for(int y=minYear;y<=maxYear;y++){
   //    for(int a=minAge;a<=maxAge;a++)
   //         ssb(y) += exp(logN(a,y)) * propMature(a,y) * weight(a,y);
   //  }
   // }
  ssb=colsum(elem_prod(elem_prod(exp(logN),propMature),stockMeanWeight));

