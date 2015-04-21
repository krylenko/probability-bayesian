// (c) daniel ford, daniel.jb.ford@gmail.com

// simple Metropolis-Hastings sampler

#include<stdlib.h>
#include<cmath>
#include<sys/time.h>

#include "C:/Users/ghost/Desktop/robotics/code/dfrobot/matrix/libdfr-matrix.h"
#include "C:/Users/ghost/Desktop/robotics/code/dfrobot/random/libdfr-rv.h"

#define DIM   2
#define GRID  50
#define X     1
#define Y     0

// fn prototypes
double secret(const Matrix& in);
//void totUp(const Matrix& Z);

int main(int argc, char** argv){

  // check cmd line and save iterations
  if( argc != 2 ){
    cout << "usage is p3 [iterations]" << endl;
    exit(-1);
  }
  int iter = atoi(argv[1]);

  // seed PRNG w/ current microsecond
  struct timeval tv;
  gettimeofday(&tv,NULL);
  rvSeed(tv.tv_usec);
  
  // init count array
  Matrix counts(GRID+1,GRID+1);
  
  // init Q function
  Matrix mu(DIM); 
  Matrix sigma(DIM,DIM); sigma.I(); //sigma = sigma*3.;
  Gaussian Q(mu,sigma);
  
  // init random location
  Matrix Z(DIM);
  Z[X][0] = (int)rvStdUniform(-25,25);
  Z[Y][0] = (int)rvStdUniform(-25,25);
  cout << "started at " << Z[X][0] << ", " << Z[Y][0] << endl;
  
  // run sampler
  for(int i=0;i<iter;++i){
  
    // sample Q(Z|Z_prev)
    Q.mu = Z; Matrix newZ(DIM); 
    newZ[X][0] = rvGaussian(Z[X][0],1); newZ[Y][0] = rvGaussian(Z[Y][0],1); 
    
    // compute p(accept)
    double A = 0;
    double pRatio = secret(newZ)/secret(Z);
    if( pRatio < 1 )
      A = pRatio;
    else
      A = 1;
    
    // choose emission and increment proper count
    int Xnew = (int)newZ[X][0]; int Ynew = (int)newZ[Y][0]; 
    if( (newZ[X][0] >= -25 && newZ[X][0] <= 25 ) && 
        (newZ[Y][0] >= -25 && newZ[Y][0] <= 25 ) ){ 
      double pick = rvStdUniform();
      if( (A == 1) || (pick <= A) ){ 
        counts[Xnew+(GRID/2)][(GRID/2)-Ynew] += 1;
        Z = newZ;
      }
      else{
        counts[(int)Z[X][0]+(GRID/2)][(GRID/2)-(int)Z[Y][0]] += 1;
      }
    }
  
  }

  // save counts to file
  char file[100];
  sprintf(file, "p3-%d.dat", iter);
  counts.printFile(file);
  
  return 0;

}

// distribution to be explored
double secret(const Matrix& in){

  int i = (int)in[X][0], j = (int)in[Y][0];
  double p = 0;

  if ((i >= -5) && (i <= -3) && (j >= -4) && (j <= 4)) 
    p = 100;
  else if ((i >= -3) && (i <= 7) && (j == -2)) 
    p = 100;
  else if ((i >= -3) && (i <= 7) && (j == 2)) 
    p = 100;
  else
    p = 2000.0 / pow((10 + i*i+j*j),2);
    
  return p;  
}
