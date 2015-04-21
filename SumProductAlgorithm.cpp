// (c) daniel ford, daniel.jb.ford@gmail.com

// sum-product algorithm for MQL beast problem

#include<cmath>
#include<stdlib.h>
#include "C:/Users/ghost/Desktop/robotics/code/dfrobot/matrix/libdfr-matrix.h"
#include "C:/Users/ghost/Desktop/robotics/code/dfrobot/random/libdfr-rv.h"

// definitions and enums
#define BINS        6
#define ANG_UNIFORM 1./BINS
#define SEX         2
#define BODY        2
enum {male, female};

// create and initialize variable and factor storage 
Matrix S(SEX); Matrix B(BODY); Matrix Bf(BODY);   // 2x1
Matrix RIAA(BINS); Matrix RIAAf(BINS);            // 6x1
Matrix LOAA(BINS*2); Matrix LIAA(BINS); Matrix ROAA(BINS*2);
Matrix LOAAf(BINS*2); Matrix LIAAf(BINS); Matrix ROAAf(BINS*2);
Matrix RILA(BINS); Matrix RILAf(BINS);
Matrix ROLA(BINS*2); Matrix LILA(BINS); Matrix LOLA(BINS*2);
Matrix ROLAf(BINS*2); Matrix LOLAf(BINS*2); Matrix LILAf(BINS);

Matrix Fa(BODY); Matrix Fb(BINS*2); Matrix Fc(BINS);
Matrix Fd(BINS*2); Matrix Fe(BINS*2); Matrix Ff(BINS*2);
Matrix Fg(BINS); Matrix Fh(BINS*2); Matrix Fi(BINS);
Matrix Fj(BINS*2); Matrix Fk(BINS*2); Matrix Fl(BINS);
Matrix Fm(BINS*2); Matrix Fn(BINS*2); Matrix Fo(BINS);
Matrix Fp(BODY); Matrix Fq(BINS); Matrix Fr(BINS);

Matrix bodyTrue(BODY, BODY); 
Matrix bodySex(BODY, BODY);
Matrix ILsexM(BINS); Matrix ILsexF(BINS); Matrix ILsex(BINS,2);
Matrix IAsexM(BINS); Matrix IAsexF(BINS); Matrix IAsex(BINS,2);
Matrix obsTrue(BINS,BINS); 
Matrix inner(BINS, BINS); 
Matrix outIn(BINS, BINS*2); 
Matrix outer(BINS*2, BINS*2); 

int trueSex = 0, calcSex = 0;

// prototypes
void initLeaves();
void sampleBeast();
void leavesToRoot();
void rootToLeaves();
int invCDF(const Matrix& pdf);
double binToAngle(const double& low, const double& high, const int& bin);
double binToAngle2Q(const double& low, const double& high, const int& bin);
Matrix percents(const double& p_m, const int& actual, const int& end=0);

int main(int argc, char** argv){
  
  if( argv[1] == NULL ){
    cout << "usage is p6 [number of runs]" << endl;
    exit(-1);
  }
  int rightCount = 0; int runs = atoi(argv[1]);
  rvSeed();

  Matrix p7plot(10,2);
  
  // load distributions from files
  bodyTrue.fromFile("bodyTrue.dat");  
  bodySex.I();
  ILsexM.fromFile("ILsexM.dat"); ILsexF.fromFile("ILsexF.dat");
  ILsexM = ILsexM.norm(); ILsexF = ILsexF.norm(); ILsex = ILsexM.conc(ILsexF,COLS);
  IAsexM.fromFile("IAsexM.dat"); IAsexF.fromFile("IAsexF.dat");
  IAsexM = IAsexM.norm(); IAsexF = IAsexF.norm(); IAsex = IAsexM.conc(IAsexF,COLS);
  obsTrue.fromFile("angPDF.dat");
  inner.fromFile("inner.dat");
  outIn.fromFile("outerIn.dat");
  outer.fromFile("outer.dat");
  /*
  bodyTrue.print(); bodySex.print(); ILsex.print(); IAsex.print();
  obsTrue.print(); inner.print(); outIn.print(); outer.print();
  */
  for(int i=0; i<runs; i++){
    // initialize leaves (actually just zeroing observation arrays in this case)
    initLeaves();
    
    // get sampled creature
    sampleBeast();
    
    // compute messages from leaves to root
    leavesToRoot();
    
    // compute messages from root to leaves
    // don't need to do this since we only need p(S) and S is the root
    //rootToLeaves();  
    
    // for marginal of S, pick highest probability  

    if( trueSex == male)
      cout << "true sex: male" << endl;
    else
      cout << "true sex: female" << endl;
      
    if(S[0][0] > S[1][0]){
      cout << "calc sex: male" << endl;
      calcSex = male;
      if( trueSex == calcSex )
        rightCount++;
    }
    else if(S[0][0] < S[1][0]){
      cout << "calc sex: female" << endl;
      calcSex = female;
      if( trueSex == calcSex )
        rightCount++;
    }
    else if(S[0][0] == S[1][0])
      cout << "Houston, we have a problem" << endl;  

    cout << "posterior (M F):" << endl;
    (S.T()).print();      
    
    p7plot = percents(S[0][0], trueSex);
    
  }
  
  cout << endl << "*********************" << endl;
  cout << "Accuracy over " << runs << " runs: " << double(rightCount)*100/runs << endl;

  p7plot = percents(2, 1, 1);  
  p7plot.print();
  p7plot.printFile("p7.dat");
 
  return 0;
  
}

// compute and pass messages from leaves to root
void leavesToRoot(){

  // compute first set of factor message (Fa - Fi)
  Fa = bodyTrue*Bf;
  Fb = obsTrue*RIAAf;
  Fc = obsTrue*LIAAf;
  Fd = outer*LOAAf;
  Fe = outer*ROAAf;
  Ff = outer*ROLAf;
  Fg = obsTrue*LILAf;
  Fh = outer*LOLAf;
  Fi = obsTrue*RILAf;
  /*
  cout << "1st factors" << endl;
  Fa.print(); Fb.print(); Fc.print(); Fd.print(); Fe.print();
  Ff.print(); Fg.print(); Fh.print(); Fi.print();
  */
  // compute second set of node messages (LOAA, ROAA, ROLA, LOLA)
  LOAA = Fd; ROAA = Fe; ROLA = Ff; LOLA = Fh;
  /*
  cout << "LOAA, ROAA, ROLA, LOLA" << endl;
  LOAA.print(); ROAA.print(); ROLA.print(); LOLA.print();
  */
  // compute second set of factor messages (Fj - Fk)
  Fj = outIn*LOAA;
  Fk = outIn*LOLA;

  //Fj = Fj.norm(); Fk = Fk.norm();
  /*
  cout << "2nd factors" << endl;
  Fj.print(); Fk.print();
  */
  // compute third set of node messages (LIAA, LILA)
  LIAA = Fc.emult(Fj);
  LILA = Fg.emult(Fk);

  //LIAA = LIAA.norm(); LILA = LILA.norm();
  /*
  cout << "LIAA, LILA" << endl;
  LIAA.print(); LILA.print();
  */
  // compute third set of factor messages (Fl - Fo)
  // note: properly set up, dot products perform the multiplication and marginalization needed for s-p
  Fl = inner*LIAA;
  Fm = outIn*ROAA;
  Fn = outIn*ROLA;
  Fo = inner*LILA;

  //Fl = Fl.norm(); Fm = Fm.norm(); Fn = Fn.norm(); Fo = Fo.norm();
  /*
  cout << "3rd factors" << endl;
  Fl.print(); Fm.print(); Fn.print(); Fo.print();
  */
  // compute fourth set of node messages (B, RIAA, RILA)
  // note: need to do element-wise multiplication to combine factor messages incoming to a node
  B = Fa; 
  RIAA = (Fb.emult(Fl)).emult(Fm); 
  RILA = (Fn.emult(Fo)).emult(Fi);

  //RIAA = RIAA.norm(); RILA = RILA.norm();
  /*
  cout << "RIAA, RILA" << endl;
  RIAA.print(); RILA.print();
  */
  // compute fourth set of factor messages (Fp - Fr)
  Fp = bodySex*B;
  Fq = IAsex.T()*RIAA;
  Fr = ILsex.T()*RILA;

  /*
  cout << "4th factors" << endl;
  Fp.print(); Fq.print(); Fr.print();
  */
  
  // compute root node (S) [actually messages from factors to S]
  S = ((Fp.emult(Fq)).emult(Fr)).norm();

}

// compute and pass messages from leaves to root
// not needed at the moment, S is root so we have the full marginal
void rootToLeaves(){

}

void initLeaves(){

  /*
  cout << "values to params" << endl;
  (Bf.T()).print(); (RIAAf.T()).print(); (LIAAf.T()).print(); (LOAAf.T()).print();
  (ROAAf.T()).print(); (ROLAf.T()).print(); (LILAf.T()).print(); (LOLAf.T()).print();
  (RILAf.T()).print();
  */

  // zero out observation arrays b/c they're global, need to be reset each iteration
  Bf.zeros(); RIAAf.zeros(); LIAAf.zeros(); LOAAf.zeros(); ROAAf.zeros();
  ROLAf.zeros(); LILAf.zeros(); LOLAf.zeros(); RILAf.zeros();
  
}

void sampleBeast(){

  int params = 9;

  double p_male = 2./3;
  double p_circToHex = 0.35;
  double p_hexToCirc = 0.5;

  enum {circbod, hexbod};
  enum {shape, liaa, loaa, lila, lola, riaa, roaa, rila, rola};
  
  Matrix beast(params);
  Matrix pdf(5,6); pdf.fromFile("PDFtable.dat");
  Matrix pdfAngNoise(6,6); pdfAngNoise.fromFile("angPDF.dat");
  
  // extract, normalize, transpose angle distributions
  Matrix pdfRIAA_f(pdf.slice(0,0,0,(pdf.cols()-1)));
  pdfRIAA_f = (pdfRIAA_f.norm()).T();
  
  Matrix pdfRIAA_m(pdf.slice(1,1,0,(pdf.cols()-1)));
  pdfRIAA_m = (pdfRIAA_m.norm()).T();
  
  Matrix pdfRILA_f(pdf.slice(2,2,0,(pdf.cols()-1)));
  pdfRILA_f = (pdfRILA_f.norm()).T();  
  
  Matrix pdfRILA_m(pdf.slice(3,3,0,(pdf.cols()-1)));
  pdfRILA_m = (pdfRILA_m.norm()).T();
  
  Matrix pdfOuter(pdf.slice(4,4,0,(pdf.cols()-1)));
  pdfOuter = (pdfOuter.norm()).T();
  
  // generate true shape
  if( rvStdUniform() <= p_male )
    { beast[shape][0] = circbod; trueSex = 0;}
  else
    { beast[shape][0] = hexbod; trueSex = 1;}

  // sample true RIAA and RILA given sex
  if( beast[shape][0] == circbod ){
    beast[riaa][0] = invCDF(pdfRIAA_m);
    beast[rila][0] = invCDF(pdfRILA_m); 
  }
  else if( beast[shape][0] == hexbod ){
    beast[riaa][0] = invCDF(pdfRIAA_f);
    beast[rila][0] = invCDF(pdfRILA_f);
  }
  
  // calculate true LIAA and LILA bins
  Matrix liaaPDF = inner.slice(beast[riaa][0], beast[riaa][0], 0, (inner.cols()-1)); liaaPDF = liaaPDF.T();
  beast[liaa][0] = invCDF(liaaPDF);
  Matrix lilaPDF = inner.slice(beast[rila][0], beast[rila][0], 0, (inner.cols()-1)); lilaPDF = lilaPDF.T();
  beast[lila][0] = invCDF(lilaPDF);
  
  // sample true absolute outer limb bins  
  beast[loaa][0] = invCDF(pdfOuter)+beast[liaa][0];
  beast[lola][0] = invCDF(pdfOuter)+beast[lila][0];
  beast[roaa][0] = invCDF(pdfOuter)+beast[riaa][0];
  beast[rola][0] = invCDF(pdfOuter)+beast[rila][0];
  /*
  cout << "True sample..." << endl;
  beast.print();
  */
  
  // observe RIAA, LIAA, RILA, LILA
  Matrix pdfA = pdfAngNoise.slice(beast[riaa][0],beast[riaa][0],0,(pdfAngNoise.cols()-1)); pdfA = pdfA.T();
  beast[riaa][0] = invCDF(pdfA);
  pdfA = pdfAngNoise.slice(beast[liaa][0],beast[liaa][0],0,(pdfAngNoise.cols()-1)); pdfA = pdfA.T();
  beast[liaa][0] = invCDF(pdfA);
  pdfA = pdfAngNoise.slice(beast[rila][0],beast[rila][0],0,(pdfAngNoise.cols()-1)); pdfA = pdfA.T();
  beast[rila][0] = invCDF(pdfA);
  pdfA = pdfAngNoise.slice(beast[lila][0],beast[lila][0],0,(pdfAngNoise.cols()-1)); pdfA = pdfA.T();
  beast[lila][0] = invCDF(pdfA);

  // sample relative outer limb angles
  // note: needed to add small regularization terms (0.01) to the outer observation matrix to avoid NaN problems
  Matrix outPDF(BINS*2, BINS*2); outPDF.fromFile("outer.dat");
  outPDF.slice(beast[loaa][0],beast[loaa][0],0,(outPDF.cols()-1));
  outPDF = outPDF.T();
  beast[loaa][0] = invCDF(outPDF);

  outPDF.slice(beast[lola][0],beast[lola][0],0,(outPDF.cols()-1));
  outPDF = outPDF.T();
  beast[lola][0] = invCDF(outPDF);
  
  outPDF.slice(beast[roaa][0],beast[roaa][0],0,(outPDF.cols()-1));
  outPDF = outPDF.T();
  beast[roaa][0] = invCDF(outPDF);
  
  outPDF.slice(beast[rola][0],beast[rola][0],0,(outPDF.cols()-1));
  outPDF = outPDF.T();
  beast[rola][0] = invCDF(outPDF);

  // observe shape and change accordingly
  if( beast[shape][0] == circbod ){
    if( rvStdUniform() <= p_circToHex )
      beast[shape][0] = hexbod;
  }
  if( beast[shape][0] == hexbod ){
    if( rvStdUniform() <= p_hexToCirc )
      beast[shape][0] = circbod;
  }
  //*/
  /*
  cout << "Observed sample..." << endl;
  beast.print();
  */
  // assign sampled values to observed params

  Bf[beast[shape][0]][0] = 1;
  RIAAf[beast[riaa][0]][0] = 1;
  LIAAf[beast[liaa][0]][0] = 1;
  LOAAf[beast[loaa][0]][0] = 1;
  ROAAf[beast[roaa][0]][0] = 1;
  ROLAf[beast[rola][0]][0] = 1;
  LILAf[beast[lila][0]][0] = 1;
  LOLAf[beast[lola][0]][0] = 1;
  RILAf[beast[rila][0]][0] = 1;
  /*
  cout << "values to params" << endl;
  (Bf.T()).print(); (RIAAf.T()).print(); (LIAAf.T()).print(); (LOAAf.T()).print();
  (ROAAf.T()).print(); (ROLAf.T()).print(); (LILAf.T()).print(); (LOLAf.T()).print();
  (RILAf.T()).print();
  */
}

// sample from given pdf using inverse CDF routine
// returns index of original pdf, not the value at that index
int invCDF(const Matrix& pdf)
{

	double tempSamp = 0.;
	int j = 0;

	// integrate to generate CDF
	Matrix CDF(pdf.rows(),1);
	CDF[0][0] = pdf[0][0];
	for(int i=0;i<CDF.rows()-1;i++)
	{
    CDF[i+1][0] = CDF[i][0] + pdf[i+1][0];	
  }
  
	// choose sample from uniform generator according to CDF 
  j = 0;
	tempSamp = rvStdUniform(0, 0.999);		
	while( CDF[j][0] < tempSamp )
    j++;
	
	return j;
	
}

// converts bin number to random angle within that bin range
// assumes 6 bins of 15 degrees each
double binToAngle(const double& low, const double& high, const int& bin){

  int binSize = 15; int binNum = 6;
  double angle = 0;
  
  angle = fmod(rvStdUniform(low, high), binNum) + bin*binSize;
  
  return angle;

}

// special exception for 180-degree outer angles
// converts bin number to random angle within that bin range
// assumes 12 bins of 15 degrees each
double binToAngle2Q(const double& low, const double& high, const int& bin){

  const int binSize = 15; const int binNum = 12;
  double angle = 0;
  
  angle = fmod(rvStdUniform(low, high), binNum) + bin*binSize;
  
  return angle;

}

// track percentage correctness vs. posterior probabilities
Matrix percents(const double& p_m, const int& actual, const int& end){

  static Matrix results(10,2);
  static int mCount = 0;
  
  if( actual == male ){
    
    mCount++;
    
    if( (0 <= p_m) && (p_m < 0.1) )
        results[0][1]++;

    if( (0.1 <= p_m) && (p_m < 0.2) )
        results[1][1]++;              
      
    if( (0.2 <= p_m) && (p_m < 0.3) )
        results[2][1]++;                  
    
    if( (0.3 <= p_m) && (p_m < 0.4) )
        results[3][1]++;  

    if( (0.4 <= p_m) && (p_m < 0.5) )
        results[4][1]++;

    if( (0.5 <= p_m) && (p_m < 0.6) )
        results[5][1]++;        
    
    if( (0.6 <= p_m) && (p_m < 0.7) )
        results[6][1]++;  

    if( (0.7 <= p_m) && (p_m < 0.8) )
        results[7][1]++;
    
    if( (0.8 <= p_m) && (p_m < 0.9) )
        results[8][1]++;

    if( (0.9 <= p_m) && (p_m < 1.0) )
        results[9][1]++;
  }
  
  if( mCount != 0 && end )
    results = results*100/mCount;  
    
  results[0][0] = 0.05;
  results[1][0] = 0.15;
  results[2][0] = 0.25;
  results[3][0] = 0.35;
  results[4][0] = 0.45;
  results[5][0] = 0.55;
  results[6][0] = 0.65;
  results[7][0] = 0.75;
  results[8][0] = 0.85;
  results[9][0] = 0.95;
  
  return results;
  
}