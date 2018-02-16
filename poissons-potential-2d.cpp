#include <iostream>
#include <vector> /* Using vectors instead of arrays */
#include <sstream>
#include <fstream>
typedef std::vector<double> double_vec; /* Simplification of code */

int main(int argc, char* argv[])
{
	/*	Variable Initialization	*/
	std::stringstream sN, sTop, sRight, sBottom, sLeft;
	double Top, Right, Bottom, Left;
	int N, Num;
	
	/*	Argument Handling (Needs Improvemnt)	*/
	Num = 100;
	sN << argv[2];
	sTop << argv[3];
	sRight << argv[4];
	sBottom << argv[5];
	sLeft << argv[6];
	if(argc == 8){
		std::stringstream sNum;
		sNum << argv[7];
		sNum >> Num;
	}
	sN >> N;
	sTop >> Top;
	sRight >> Right;
	sLeft >> Left;
	sBottom >> Bottom;

	/* Phi vector initialization	*/
	std::vector<double_vec> phi(N+1,double_vec(N+1, 0));	
	for(int i = 0; i < N+1; i++){
		phi[0][i] = Top;
		phi[N][i] = Bottom;
		phi[i][N] = Right;
		phi[i][0] = Left;
	}

	/* Steping Action	*/
	const double pi = 3.14159265359;
	const double omega = 2 / (1 + pi / N); /* Used as a lazy value */
	std::vector<double_vec> phi_new(N+1,double_vec(N+1)); /* */
	for(int count = 0; count < Num; count++){
		for(int i = 1; i < N; i++){
			for(int j = 1; j < N; j++){
				phi_new[i][j] = phi[i][j] + (omega / 4) * (phi[i+1][j] + phi_new[i-1][j] + phi[i][j+1] + phi_new[i][j-1] - 4 * phi[i][j]);
			}
		}
		std::swap(phi, phi_new); /* Using swap which is faster than reassigning values again */
	}

	/*	Data file output	*/
	std::ofstream datafile;
	datafile.open(argv[1]);
	datafile << "# This is a data file of potential over a plane\n";
	datafile << "# Variables: N=" << N << ", Top=" << Top << ", Right=" << Right << ", Bottom=" << Bottom << ", Left=" << Left << ", cycles=" << Num << "\n\n";
	for(int i = 0; i < N+1; i++){
		for(int j = 0; j < N+1; j++){
			datafile << i << "	" << j << "	" << phi[i][j] << "\n";
		}
		datafile << "\n";
	}
	datafile.close();
	return 0;
}
