DATA_SECTION
  // catch and life history info
  init_int minAge;
  init_int maxAge;
  init_int minYear;
  init_int maxYear;
  init_matrix catchNo(minAge,maxAge,minYear,maxYear);
  init_matrix stockMeanWeight(minAge,maxAge,minYear,maxYear);
  init_matrix propMature(minAge,maxAge,minYear,maxYear);
  init_matrix M(minAge,maxAge,minYear,maxYear);

  // survey info
  init_int minAgeS;
  init_int maxAgeS;
  init_int minYearS;
  init_int maxYearS;
  init_number surveyTime;
  init_matrix survey(minAgeS,maxAgeS,minYearS,maxYearS);

  // log transform catch and survey data
  matrix logC(minAge,maxAge,minYear,maxYear)
  !! logC=log(catchNo);
  matrix logSurvey(minAgeS,maxAgeS,minYearS,maxYearS)
  !! logSurvey=log(survey);

PARAMETER_SECTION
  init_vector logN1Y(minAge,maxAge);                               // numbers year 1
  init_vector logN1A(minYear+1,maxYear);                           // recruits age-1
  init_vector logFY(minYear,maxYear);                              // fishing mortality by year
  init_vector logFA(minAge,maxAge-3);                              // fishing mortality by age - estimate only first few ages
  init_number logVarLogCatch;                                      // log variance of the log catch
  init_vector logQ(minAgeS,maxAgeS);                               // log catchability coefficient for survey
  init_number logVarLogSurvey;                                     // log variance of log survey

  objective_function_value nll;

  matrix F(minAge,maxAge,minYear,maxYear);                        // Fishing mortality by age over time
  matrix logN(minAge,maxAge,minYear,maxYear);                     // log Numbers at age over time
  matrix predLogC(minAge,maxAge,minYear,maxYear);                 // predicted log catch at age over time
  matrix predLogSurvey(minAgeS,maxAgeS,minYearS,maxYearS);        // predicted log survey at age over time
  number ss;                                                      // variance for observation error                                          
  vector logFAextend(minAge,maxAge);                              // log fishing mortality at age for all ages
  sdreport_vector ssb(minYear,maxYear);                           // spawning stock biomass

PROCEDURE_SECTION
  logFAextend.initialize();
  logFAextend(minAge,maxAge-3) = logFA;  

  F = outer_prod(exp(logFAextend), exp(logFY));  

  logN.colfill(minYear,logN1Y);
  for(int y=minYear+1; y<=maxYear; ++y){
    logN(minAge,y)=logN1A(y);
    for(int a=minAge+1; a<=maxAge; ++a){
      logN(a,y)=logN(a-1,y-1)-F(a-1,y-1)-M(a-1,y-1);
    }
  }

  predLogC = log(F)-log(F+M)+log(1.0-exp(-F-M))+logN;

  ss = exp(logVarLogCatch);
  int N = (maxYear-minYear+1)*(maxAge-minAge+1);
  nll = 0.5 * (N * log(2.0 * M_PI * ss) + sum(square(logC-predLogC))/ss);

  ss = exp(logVarLogSurvey);
  for(int y=minYearS; y<=maxYearS; ++y){
    for(int a=minAgeS; a<=maxAgeS; ++a){
      predLogSurvey(a,y) = logQ(a)-(F(a,y)+M(a,y)) * surveyTime + logN(a,y);
      nll += 0.5 * (log(2.0 * M_PI * ss) + square(logSurvey(a,y)-predLogSurvey(a,y))/ss);
    }
  }
  ssb = colsum(elem_prod(elem_prod(exp(logN),propMature),stockMeanWeight));

