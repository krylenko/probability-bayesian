// (c) daniel ford, daniel.jb.ford@gmail.com

// sample from N(0,1) w/ 10, 100, 1000, 10000 samples

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

int main(int argc, char** argv)
{

	double gaussSigma = 1., gaussMean = 0.;
	double step = 0.1;		
	int length = 5;					// actual length is L*2
									// assumes symmetry around 0,0

	int samples = atoi(argv[1]);
	int binNum = 50;

	rvSeed();						// init random number gens

	// generate actual Gaussian for comparison	
	Matrix x(-double(length), double(length), step);
	Matrix gauss(x);	
	gauss.zeros();
	gauss = (x.conc(gauss,1)).T();
	
	for(int i=0;i<gauss.rows();++i)
		gauss[i][1] = (1/(sqrt(gaussSigma)*Gc)) *
						exp(-(pow((gauss[i][0]-gaussMean),2) /
						(2*gaussSigma)));

	// generate samples and create histograms
	Matrix gaussSamp = invCDFsample(samples, step, gauss);
	Matrix histSamp = rvBin(gaussSamp.T(), -double(length),
							double(length), binNum);  					

	gauss.printFile("p12_real.dat", NOBLANKS);
	histSamp.printFile("p12_hist.dat", NOBLANKS);
	
	return 0;
	
}
	
	