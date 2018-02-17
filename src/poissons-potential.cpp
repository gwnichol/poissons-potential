/*
 * */

#include <iostream>
#include <vector> /* Using vectors instead of arrays */
#include <sstream>
#include <fstream>
#include <cstring>
typedef std::vector<double> double_vec; /* Simplification of code */
typedef std::vector<double_vec> double_vec_vec;

void print_help_text(char* arg){
	std::cout << "Usage:\n" << arg << " [options]\n\nOptions:\n"
		<< "  -h,  --help         Shows the help prompt\n"
		<< "  -3D, --3demensions  Activates 3 spacial demensions\n"
		<< "  -Ef, --field        Outputs electric field approximation instead of voltage\n"
		<< "  -s,  --size         Number of points on one edge of the surface (>2). Default: 25\n"
		<< "  -l,  --length       Length of one edge of the surface in meters. Default: 1.0\n"
		<< "                      Has no effect without field approximation\n"
		<< "  -c,  --count        The number of iterations that will be performed. Default: 100\n"
		<< "  -b,  --boundaries   Sets the boundary conditions (V) Takes 4 or 6 values\n"
		<< "                      2D: Top Right Bottom Left. Default: 0 0 0 0\n"
		<< "                      3D: Top Right Bottom Left Front Back. Default: 0 0 0 0 0 0\n"
		<< "  -f,  --filename     Sets the file name. Default: \"data.dat\"\n";
}

int main(int argc, char* argv[])
{
	/*	Variable Initialization	*/

	double Top, Right, Bottom, Left, Front, Back, Length, Delta;
	int N, Num;
	std::string filename;
	const double pi = 3.14159265359;
	bool CUBE = 0, EFIELD = 0;

	/* Default Values */
	Top = 1;
	Right = 1;
	Bottom = 1;
	Left = 1;
	Front = 1;
	Back = 1;
	N = 25;
	Num = 100;
	filename = "data.dat";
	Length = 1.0;

	/*	Argument Handling (Needs Improvemnt)	*/
	for (int i = 0; i < argc; i++){
		if((std::strcmp(argv[i], "-3D") == 0) || (std::strcmp(argv[i], "--3demensions") == 0)){
			CUBE = 1;
		}
		if((std::strcmp(argv[i], "-Ef") == 0) || (std::strcmp(argv[i], "--field") == 0)){
			EFIELD = 1;
		}
	}
	for (int i = 0; i < argc; i++){
		if((std::strcmp(argv[i], "-h") == 0) || (std::strcmp(argv[i], "--help") == 0)){
			std::cout << "This program uses Poisson's equation to map potential in 2-space when given boundary conditions\nIt can then output an electric field aproximation if desired\n";
			print_help_text(argv[0]);
			return 1;
		} else if ((std::strcmp(argv[i], "-s") == 0) || (std::strcmp(argv[i], "--size") == 0)){
			if(i + 1 < argc){
				std::stringstream sN;
				sN << argv[i + 1];
				sN >> N;
				}
			if(not (i + 1 < argc) || (N < 3)){
				std::cout << "Argument, " << argv[i] << ", needs an integer size value greater than 2.\n";
				print_help_text(argv[0]);
				return 1;
			}
		} else if ((std::strcmp(argv[i], "-c") == 0) || (std::strcmp(argv[i], "--count") == 0)){
			if(i + 1 < argc){
				std::stringstream sNum;
				sNum << argv[i + 1];
				sNum >> Num;
			}
			if( not (i + 1 < argc) || (Num < 1)){
				std::cout << "Argument,  " << argv[i] << ", needs an integer iteration count.\n";
				print_help_text(argv[0]);
				return 1;
			}
		} else if ((std::strcmp(argv[i], "-b") == 0) || (std::strcmp(argv[i], "--boundaries") == 0)){
			if( i + 4 < argc){
				std::stringstream sTop, sRight, sBottom, sLeft;
				sTop << argv[i + 1];
				sRight << argv[i + 2];
				sBottom << argv[i + 3];
				sLeft << argv[i + 4];
				sTop >> Top;
				sBottom >> Bottom;
				sLeft >> Left;
				sRight >> Right;
				if((i + 6 < argc) & (CUBE)){
					std::stringstream sFront, sBack;
					sFront << argv[i + 5];
					sBack << argv[i + 6];
					sFront >> Front;
					sBack >> Back;
				}}
			if((not (i + 6 < argc)) && (CUBE)){
				std::cout << "Argument, " <<  argv[i] << ", need six integer boundary values directly following it when 3D is active.\n";
				print_help_text(argv[0]);
				return 1;
			}
			if((not (i + 4 < argc)) && (not (CUBE))){
				std::cout << "Argument, " << argv[i] << ", needs four integer boundary values directly following it.\n";
				print_help_text(argv[0]);
				return 1;
			}
		} else if ((std::strcmp(argv[i], "-f") == 0) || (std::strcmp(argv[i], "--datafile") == 0)){
			if(i + 1 < argc){
				filename = argv[i + 1];
			} else {
				std::cout << "Argument, " << argv[i] << ", needs a filename following it.\n";
				print_help_text(argv[0]);
				return 1;
			}
		} else if((std::strcmp(argv[i], "-l") == 0) || (std::strcmp(argv[i], "--length") == 0)){
			if(i + 1 < argc){
				std::stringstream sLength;
				sLength << argv[i + 1];
				sLength >> Length;
			}
			if( (not (i + 1 < argc)) || (Length == 0) ){
				std::cout << "Argument, " << argv[i] << ", needs a number greater than zero following it.\n";
				print_help_text(argv[0]);
				return 1;
			}
		}
	}
	
	std::cout << "Initializing Vector\n";
	/* Phi vector initialization	*/
	if( not CUBE ){
	/* --------------- 2D Creation ------------ */
	std::vector<double_vec> phi(N+1,double_vec(N+1, 0));
	for(int i = 0; i < N+1; i++){
		phi[i][N] = Top;
		phi[i][0] = Bottom;
		phi[N][i] = Right;
		phi[0][i] = Left;
	}

	std::cout << "Steping\n";
	/* Steping Action	*/

	const double omega = 2 / (1 + pi / N); /* Used as a relaxation constant */
	std::vector<double_vec> phi_new(N+1,double_vec(N+1)); /* */
	for(int count = 0; count < Num; count++){
		for(int i = 1; i < N; i++){
			for(int j = 1; j < N; j++){
				phi_new[i][j] = phi[i][j] + (omega / 4) * (phi[i+1][j] + phi_new[i-1][j] + phi[i][j+1] + phi_new[i][j-1] - 4 * phi[i][j]);
			}
		}
		std::swap(phi, phi_new); /* Using swap which is faster than reassigning values again */
	}

	std::cout << "Writing to " << filename << ".\n";
	/*	Data file output	*/
	std::ofstream datafile;
	datafile.open(filename);
	if(not(EFIELD)){
		datafile << "# This is a data file of potential over a plane\n";
		datafile << "# Variables: N=" << N << ", Top=" << Top << ", Right=" << Right << ", Bottom=" << Bottom << ", Left=" << Left << ", cycles=" << Num << "\n\n";
		datafile << "# X	Y	V\n";
		for(int i = 0; i < N+1; i++){
			for(int j = 0; j < N+1; j++){
				datafile << i << "	" << j << "	" << phi[i][j] << "\n";
			}
			datafile << "\n";
		}
	} else {
		datafile << "# This is a data file of electric field over a plane\n"
			<< "# Variables: N=" << N << ", Top=" << Top << ", Right=" << Right << ", Bottom=" << Bottom << ", Left=" << Left << ", cycles=" << Num << ", Length=" << Length << "\n\n"
			<< "# X	Y	Ex	Ey\n";
		Delta = Length / N;
		for(int i =1; i < N; i++){
			for(int j = 1; j < N; j++){
				datafile << i * Delta << "	" << j * Delta << "	" << (phi[i][j] - phi[i+1][j]) / Delta << "	" << (phi[i][j] - phi[i][j+1]) / Delta << "\n";
			}
			datafile << "\n";
		}
	}
	datafile.close();}
	else{
	/* --------------------------------------------- */
	/* ---------------- 3D Creation ---------------- */
	
	std::vector<double_vec_vec> phi(N+1,double_vec_vec(N+1, double_vec(N+1,0)));
	for(int k = 1; k < N; k++){ /* Sets boundary conditions */
		for(int i = 0; i < N + 1; i++){
			phi[i][k][N] = Top;
			phi[i][k][0] = Bottom;
			phi[N][i][k] = Right;
			phi[0][i][k] = Left;
			phi[i][0][k] = Front;
			phi[i][N][k] = Back;
	}}

	/* Steping Action	*/
	const double omega = 2 / (1 + pi / N*N);	/* Used as a relaxation constant: N^2 instead of N : Found it works better */
	std::vector<double_vec_vec> phi_new(N+1,double_vec_vec(N+1, double_vec(N+1))); /* */
	for(int count = 0; count < Num; count++){
		for(int i = 1; i < N; i++){
			for(int j = 1; j < N; j++){
				for(int k = 1; k < N; k++){
					phi_new[i][j][k] = phi[i][j][k] + (omega / 6) * (phi[i+1][j][k] + phi_new[i-1][j][k] + phi[i][j+1][k] + phi_new[i][j-1][k] + phi[i][j][k+1] + phi_new[i][j][k-1] - 6 * phi[i][j][k]); /* Incorporates the third spacial demension */
				}
			}
		}
		std::swap(phi, phi_new); /* Using swap which is faster than reassigning values again */
	}

	/*	Data file output	*/
	std::ofstream datafile;
	datafile.open(filename);
	if(not(EFIELD)){
		datafile << "# This is a data file of potential over a plane\n";
		datafile << "# Variables: N=" << N << ", Top=" << Top << ", Right=" << Right << ", Bottom=" << Bottom << ", Left=" << Left << ", Front=" << Front << ", Back=" << Back << ", cycles=" << Num << "\n\n";
		datafile << "# X	Y	Z	V\n";
		for(int i = 0; i < N+1; i++){
			for(int j = 0; j < N+1; j++){
				for(int k = 0; k < N+1; k++){
					datafile << i << "	" << j << "	" << k << "	" << phi[i][j][k] << "\n";
				}
			}
			datafile << "\n";
		}
	} else {
		datafile << "# This is a data file of approximate electric field vectors\n"
			<< "# Variable: N=" << N << ", Top=" << Top << ", Right=" << Right << ", Bottom=" << Bottom << ", Left=" << Left << ", Front=" << Front << ", Back=" << Back << ", cycles=" << Num << ", Length=" << Length << "\n\n"
			<< "# X	Y	Z	Ex	Ey	Ez\n";
		Delta = Length / N;
		for(int i = 1; i < N; i++){
			for(int j = 1; j < N; j++){
				for(int k = 1; k < N; k++){
					datafile << i * Delta << "	" << j * Delta << "	" << k * Delta << "	" << (phi[i][j][k] - phi[i+1][j][k]) / Delta << "	" << (phi[i][j][k] - phi[i][j+1][k]) / Delta << "	" << (phi[i][j][k] - phi[i][j][k+1]) / Delta << "\n";
				}
			}
			datafile << "\n";
		}
	}
	datafile.close();
	}
	/* -------------------------------------------------- */
	return 0;
}
