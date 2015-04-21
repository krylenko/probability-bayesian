// (c) daniel ford, daniel.jb.ford@gmail.com

// simple Markov chain Monte Carlo sampler

#include<stdlib.h>
#include<cmath>
#include<sys/time.h>

#include "C:/Users/ghost/Desktop/robotics/code/dfrobot/matrix/libdfr-matrix.h"
#include "C:/Users/ghost/Desktop/robotics/code/dfrobot/random/libdfr-rv.h"

//#define N         20
#define OFFSET    10
#define power     100
#define THRESHOLD 1e-10

int main(int argc, char** argv){

  if(argc != 2)
    exit(-1);
  const int N = atoi(argv[1]);

  Matrix T(N,N); Matrix newT(N,N);
  Matrix uniform(N); uniform.fill(1./N);
  Matrix stationary(N);
  int time = 0;
  
  // generate transition matrix
  for(int j=0;j<N;++j){
    for(int i=0;i<N;++i){
      T[i][j] = rvGaussian(fabs(i-j),5) + OFFSET; // add offset to avoid negative "probabilities"
    }
  }
  T = T.norm(COLS);
  //T.printFile("p4a.dat");   // save matrix for part a)
  T.printFile("p4d.dat");   // save matrix for part d)

  // find stationary probability vector
  
  //cout << " **** initial" << endl;
  //uniform.print();
  
  cout << endl << " **** step 1" << endl;
  stationary = T*uniform; time += 1;
  //stationary.print();
  
  double R = (uniform-stationary).mag()/(uniform.mag()+stationary.mag());
  cout << "R " << R << endl;
  
  for(int k=2;k<=power;++k){
    if(R < THRESHOLD)
      break;
    cout << " **** step " << k << endl;
    Matrix old = stationary;
    stationary = T*stationary;
    R = (old-stationary).mag()/(old.mag()+stationary.mag());
    cout << "R " << R << endl;
    //stationary.print();
    time += 1;
  }  
  
  cout << endl << "N " << N << " time " << time << endl;
  
  return 0;

}