#include <iostream>
#include <vector> /* Using vectors instead of arrays */
#include <sstream>
#include <fstream>

typedef std::vector<double> double_vec; /* Simplification of code */
typedef std::vector<double_vec> double_vec_vec;

int main(int argc, char* argv[])
{
	/*	Variable Initialization	*/
	std::stringstream sN, sTop, sRight, sBottom, sLeft, sFront, sBack;
	double Top, Right, Bottom, Left, Front, Back;
	int N, Num;
	
	/*	Argument Handling (Needs Improvemnt)	*/
	Num = 100;
	sN << argv[2];
	sTop << argv[3];
	sRight << argv[4];
	sBottom << argv[5];
	sLeft << argv[6];
	sFront << argv[7];
	sBack << argv[8];
	if(argc == 10){
		std::stringstream sNum;
		sNum << argv[9];
		sNum >> Num;
	}
	
	sN >> N;
	sTop >> Top;
	sRight >> Right;
	sLeft >> Left;
	sBottom >> Bottom;
	sFront >> Front;
	sBack >> Back;

	if( N < 4){std::cout << "N must be greater than three!\n"; return(1);}

	/* Phi vector initialization	*/
	std::vector<double_vec_vec> phi(N+1,double_vec_vec(N+1, double_vec(N+1,0)));
	for(int k = 1; k < N; k++){
		for(int i = 0; i < N+1; i++){
			phi[i][k][N] = Top;
			phi[i][k][0] = Bottom;
			phi[N][i][k] = Right;
			phi[0][i][k] = Left;
			phi[i][0][k] = Front;
			phi[i][N][k] = Back;
	}}
	std::cout << phi[1][N][N] << "\n";

	/* Steping Action	*/
	const double pi = 3.14159265359;
	const double omega = 2 / (1 + pi / N*N); /* Used as a lazy value */
	std::vector<double_vec_vec> phi_new(N+1,double_vec_vec(N+1, double_vec(N+1))); /* */
	for(int count = 0; count < Num; count++){
		for(int i = 1; i < N; i++){
			for(int j = 1; j < N; j++){
				for(int k = 1; k < N; k++){
					phi_new[i][j][k] = phi[i][j][k] + (omega / 6) * (phi[i+1][j][k] + phi_new[i-1][j][k] + phi[i][j+1][k] + phi_new[i][j-1][k] + phi[i][j][k+1] + phi_new[i][j][k-1] - 6 * phi[i][j][k]);

					// phi_new[i][j][k] = (1 / 6) * (phi[i+1][j][k] + phi_new[i-1][j][k] + phi[i][j+1][k] + phi_new[i][j-1][k] + phi[i][j][k+1] + phi_new[i][j][k-1]);
				}
			}
		}
		std::swap(phi, phi_new); /* Using swap which is faster than reassigning values again */
	}

	/*	Data file output	*/
	std::ofstream datafile;
	datafile.open(argv[1]);
	datafile << "# This is a data file of potential over a plane\n";
	datafile << "# Variables: N=" << N << ", Top=" << Top << ", Right=" << Right << ", Bottom=" << Bottom << ", Left=" << Left << ", Front=" << Front << ", Back=" << Back << ", cycles=" << Num << "\n\n";
	for(int i = 0; i < N+1; i++){
		for(int j = 0; j < N+1; j++){
			for(int k = 0; k < N+1; k++){
				datafile << i << "	" << j << "	" << k << "	" << phi[i][j][k] << "\n";
			}
		}
		datafile << "\n";
	}
	datafile.close();
	return 0;
}
