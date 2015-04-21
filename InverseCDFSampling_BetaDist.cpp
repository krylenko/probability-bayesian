// (c) daniel ford, daniel.jb.ford@gmail.com

// sample from B(2,5) w/ 10, 100, 1000, 10000 samples

//absolute paths in DOS
#include "C:/Users/ghost/Desktop/robotics/code/dfrobot/matrix/libdfr-matrix.h"
#include "C:/Users/ghost/Desktop/robotics/code/dfrobot/random/libdfr-rv.h"

// relative paths in DOS
//#include "../matrix/libdfr-matrix.h"
//#include "../random/libdfr-rv.h"

#include <cmath>
#include <stdlib.h>

#define PI	3.141592653589793
#define Gc	sqrt(2*PI)

int factorial(int n);
double betaCalc(const int& alpha, const int& beta, const double& x);

int main(int argc, char** argv)
{

	int alpha = 2, beta = 5;
	
	double step = 0.01;		
	int length = 1;					

	int samples = atoi(argv[1]);
	int binNum = 50;

	rvSeed();						// init random number gens

	// generate actual Beta distribution for comparison	
	Matrix x(0, double(length), step);
	Matrix betas(x);	
	betas.zeros();
	betas = (x.conc(betas,1)).T();
	
	for(int i=0;i<betas.rows();++i)
		betas[i][1] = betaCalc(alpha, beta, betas[i][0]);

	// generate samples and create histograms
	Matrix betaSamp = invCDFsample(samples, step, betas);
	Matrix histSamp = rvBin(betaSamp.T(), 0,
							double(length), binNum);  						

	betas.printFile("p13_real.dat", NOBLANKS);
	histSamp.printFile("p13_hist.dat", NOBLANKS);
	
	return 0;
	
}

int factorial(int n)
{
	if( n<0 ) return -1;
	if( n == 0 ) return 1;
	return n*factorial(n-1);
}

double betaCalc(const int& alpha, const int& beta, const double& x)
{
	return (factorial((alpha+beta)-1)/
			(factorial(alpha-1)*factorial(beta-1))) *
		   pow(x, (alpha-1)) * pow((1-x),(beta-1)) ;
}
	
	