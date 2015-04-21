// (c) daniel ford, daniel.jb.ford@gmail.com

// EM for Gaussian mixtures
// full covariance matrix learning

#include<stdlib.h>
#include<sys/time.h>
#include<cmath>
#include<vector>

#include "C:/Users/ghost/Desktop/robotics/code/dfrobot/matrix/libdfr-matrix.h"
#include "C:/Users/ghost/Desktop/robotics/code/dfrobot/random/libdfr-rv.h"

enum{x,y}; enum{OFF,ON};
#define HOLDFLG 1     
#define HOLDOUT 10    // how often to hold points out
#define N       2000
#define DIM     2

// function prototypes
Matrix E_step(const int& K, const Matrix& data, vector<Gaussian>& clusters,
              Matrix& resp, Matrix& weights, Matrix& LL, Matrix& LL_HO);
void M_steps(const int& K, const Matrix& data, vector<Gaussian>& clusters,
             const Matrix& resp, Matrix& weights, const Matrix& respTot);
void partition(const int& K, const Matrix& data, const Matrix& resp);

int main(int argc, char** argv){

  // check input for proper args
  if(argc != 4){
    cout << "usage is p2 [K] [data filename] [num runs]" << endl;
    exit(-1);
  }

  // init variables and arrays
  const int K = atoi(argv[1]);
  const int RUNS = atoi(argv[3]);
  Matrix data(N,DIM);
  data.fromFile(argv[2]);
  data.printFile("p2-init.dat");
 
  // seed PRNG w/ current microsecond
  struct timeval tv;
  gettimeofday(&tv,NULL);
  rvSeed(tv.tv_usec);
 
  // init clusters
  vector<Gaussian> clusters;
  clusters.reserve(K);
  for(int k=0;k<K;++k){
      Matrix mu(DIM); Matrix sigma(DIM,DIM);
      int startidx = (int)rvStdUniform(0,N);
      mu[x][0] = data[startidx][x];
      mu[y][0] = data[startidx][y];
      sigma.I();
      clusters.push_back(Gaussian(mu,sigma));
  }
      
  // init responsibilities and weights
  Matrix resp(N,K); resp.fill(1./K);
  Matrix weights(K); weights.fill(1./K);
  
  // init log likelihoods 
  Matrix logLike(RUNS,2);
  Matrix logLikeHeldOut(RUNS,2);
  
  // actual EM iterations
  for(int i=0;i<RUNS;i++){
    
    // calculate log likelihood and expectation
    Matrix respTot = E_step(K, data, clusters, resp, weights, logLike, logLikeHeldOut);  
    
    // maximize parameters
    M_steps(K, data, clusters, resp, weights, respTot);
 
  }

  // save log likelihoods for graphing/checking
  logLike.printFile("LL-p2.dat");
  logLikeHeldOut.printFile("LL-HO-p2.dat");
  
  // spit out separate files for colored plot
  partition(K, data, resp);
 
  return 0;
  
}

// calculate responsibilities
Matrix E_step(const int& K, const Matrix& data, vector<Gaussian>& clusters,
              Matrix& resp, Matrix& weights, Matrix& LL, Matrix& LL_HO){

  static int runIdx = 0;
  Matrix respTot(K);
  LL[runIdx][1] = 0; LL_HO[runIdx][1] = 0;
  
  // loop over data
  for(int i=0;i<N;++i){
  
    #if HOLDFLG   
    // hold out every tenth data point and calc log likelihood
    double totalKho = 0;
    if( i%HOLDOUT == 0 ){
      Matrix heldOut(DIM);
      heldOut[x][0] = data[i][x]; heldOut[y][0] = data[i][y];      
      for(int k=0;k<K;++k){
        totalKho += weights[k][0]*clusters[k].pX(heldOut);
      }
      ++i;
    }
    #endif
    
    // extract transposed data point
    Matrix point(DIM);
    point[x][0] = data[i][x]; point[y][0] = data[i][y]; 

    // calc cluster normalization factor
    double totalK = 0;
    for(int k=0;k<K;++k)
      totalK += weights[k][0]*clusters[k].pX(point);

    // calc current cluster's responsibilities 
    for(int k=0;k<K;++k){
      double thisK = weights[k][0]*clusters[k].pX(point);
      if(totalK != 0)
        resp[i][k] = thisK/totalK;
      else resp[i][k] = 0;
    }
    // save log likelihoods
    LL[runIdx][0] = LL_HO[runIdx][0] = runIdx;
    LL[runIdx][1] += (1./N)*log(totalK);

    #if HOLDFLG
    // have to check mod of prev i since we've already incremented
    if( (i-1)%HOLDOUT == 0 )
      LL_HO[runIdx][1] += ((double)HOLDOUT/(double)N)*log(totalKho);
    #endif
    
  }
  
  // calculate responsibility sum per cluster
  for(int k=0;k<K;++k)
    respTot[k][0] = resp.sum(COLS,k);
  
  runIdx++;

  // print out some stuff to watch while it's running
  //resp.print();
  //respTot.print();  
  for(int k=0;k<K;++k){
    cout << "*** " << k+1 << endl;
    clusters[k].mu.print(); clusters[k].sigma.print();
  }
  
  return respTot;
  
}

// update parameters based on responsibilities
void M_steps(const int& K, const Matrix& data, vector<Gaussian>& clusters, const Matrix& resp, Matrix& weights, const Matrix& respTot){

  double Xvar = 0, Yvar = 0;
  
  // loop over clusters
  for(int k=0;k<K;++k){
  
    clusters[k].mu.zeros(); 
    clusters[k].sigma.zeros(); 
  
    // loop over data and calculate relevant sums for this cluster
    for(int i=0;i<N;++i){

      // update parameters based on non-held-out points
      #if HOLDFLG
      if(i%HOLDOUT != 0){
      #endif
        Matrix point(DIM);
        point[x][0] = data[i][x]; point[y][0] = data[i][y]; 
    
        // mean
        clusters[k].mu[x][0] += data[i][x]*resp[i][k]/respTot[k][0];
        clusters[k].mu[y][0] += data[i][y]*resp[i][k]/respTot[k][0];       
        
        // variance
        clusters[k].sigma = clusters[k].sigma + ((point-clusters[k].mu)*((point-clusters[k].mu).T()))*resp[i][k];
      #if HOLDFLG
      }
      #endif
      
    }
    // normalize and update parameters
    clusters[k].sigma = clusters[k].sigma/respTot[k][0];    
    weights[k][0] = respTot[k][0]/N;
  }  

}

// separate data into K cluster files using the responsibilities
void partition(const int& K, const Matrix& data, const Matrix& resp){

  // set up variables and arrays
  int flag = 0;
  vector<Matrix> groups;
  groups.reserve(K);
  for(int k=0;k<K;++k){
      groups.push_back(Matrix(1,DIM)); }  
  
  // loop over data to extract each point to its maximally responsible cluster
  for(int i=0;i<N;++i){
  
    #if HOLDFLG
    if(i%HOLDOUT!=0){
    #endif
      // simple argmax to select cluster
      double max = resp[i][0];
      int argmax = 0;
      for(int k=0;k<K;++k){
        if( resp[i][k] > max )
          argmax = k;
      }
    
      // add point to that cluster's file
      Matrix point(1,DIM);
      point[0][x] = data[i][x];
      point[0][y] = data[i][y];
      groups[argmax] = groups[argmax].conc(point,ROWS);
    #if HOLDFLG
    }
    #endif
  }
  
  // print partition files
  for(int k=0;k<K;++k){
    // slice off (0,0) point from matrix init
    groups[k] = groups[k].slice(1,(groups[k].rows()-1),0,(groups[k].cols()-1));  
    char file[23];
    sprintf(file,"p2-%i-clusters-%i.dat", K, k+1);
    groups[k].printFile(file);
  }
  
}