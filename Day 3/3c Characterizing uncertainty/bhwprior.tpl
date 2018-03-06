GLOBALS_SECTION
  #include <fstream>
  using namespace std;
  ofstream clogf("program.log");
  #define TRACE(obj) clogf<<"line "<<__LINE__<<", file "<<__FILE__<<", "<<#obj" =\n " \
                          <<obj<<endl<<endl; 

DATA_SECTION
  init_int nR
  init_int nC
  init_matrix obs(1,nR,1,nC)
  vector ssb(1,nR)
  !! ssb=column(obs,1);
  vector logR(1,nR)
  !! logR=column(obs,2);
 
  !! ad_comm::change_datafile_name("prior.cor");
  
  init_int idxA
  init_adstring nameA
  init_number estA
  init_number sdA
  init_number corrAA

  init_int idxB
  init_adstring nameB
  init_number estB
  init_number sdB
  init_number corrBA
  init_number corrBB

  vector meanAB(1,2)
  matrix covAB(1,2,1,2)
 LOC_CALCS
  meanAB(1)=estA;
  meanAB(2)=estB;
  covAB(1,1)=square(sdA);
  covAB(1,2)=sdA*sdB*corrBA;
  covAB(2,1)=covAB(1,2);
  covAB(2,2)=square(sdB);
  TRACE(meanAB);
  TRACE(covAB);
 END_CALCS

PARAMETER_SECTION
  init_number loga;
  init_number logb;
  init_number logSigma;
  number sigmaSq;
  vector pred(1,nR);
  vector vecAB(1,2);

  objective_function_value nll;
 
PROCEDURE_SECTION
  sigmaSq=exp(2.0*logSigma);
  pred=loga+log(ssb)-log(1+exp(logb)*ssb);
  nll=0.5*(nR*log(2*M_PI*sigmaSq)+sum(square(logR-pred))/sigmaSq);

  vecAB(1)=loga;
  vecAB(2)=logb;
  dvar_vector diff=vecAB-meanAB;
  nll+=0.5*(log(2.0*M_PI)*2.0+log(det(covAB))+diff*inv(covAB)*diff);

  