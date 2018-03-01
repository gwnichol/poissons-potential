/*
 * */

#include <iostream> /* System IO */
#include <vector> /* Using vectors instead of arrays */
#include <string> /* Argument Handling */
#include <fstream>

typedef std::vector<double> double_vec; /* Simplification of code */
int main(){
	/* Default Values */
	int x, y;
	double potentialVoltage = 4;
	int N = 500;
	int Num = 1000;
	std::string filename = "data.dat";
	int numLines = 4;
	double lines[numLines][4] = {
		{0.2, 0.2, 0.2, 0.8},
		{0.2, 0.8, 0.8, 0.8},
		{0.2, 0.2, 0.8, 0.2},
		{0.2, 0.8, 0.2, 0.8},
	};

	/* --------------- 2D Creation ------------ */
	std::vector<double_vec> phi(N+1,double_vec(N+1, 0));
	std::vector<std::vector<bool>> iter(N+1, std::vector<bool>(N+1, 1));
	// Create initialization lines
	for (int i = 0; i < numLines; i++){
		double len = (N+1)*((((lines[i][2] - lines[i][0])*(lines[i][2] - lines[i][0]))+((lines[i][3] - lines[i][1])*(lines[i][3] - lines[i][1]))),0.5);
		int dist = 2*len; 
		std::cout << len << " " << dist << "\n";
		for( int t = 1; t < dist; t++){
			x = int ((N+1)*(t*(lines[i][2] - lines[i][0])/dist + lines[i][0]));
			y = int ((N+1)*(t*(lines[i][3] - lines[i][1])/dist + lines[i][1]));
			iter[x][y] = 0;
			phi[x][y] = potentialVoltage;
		}
	}
	std::cout << "Steping\n";
	/* Steping Action	*/

	const double pi = 3.14159265359;
	const double omega = 2 / (1 + pi / N); /* Used as a lazy value */
	std::vector<double_vec> phi_new(N+1,double_vec(N+1)); /* */
	for(int count = 0; count < Num; count++){
		for(int i = 1; i < N; i++){
			for(int j = 1; j < N; j++){
				if(iter[i][j]){
				phi_new[i][j] = phi[i][j] + (omega / 4) * (phi[i+1][j] + phi_new[i-1][j] + phi[i][j+1] + phi_new[i][j-1] - 4 * phi[i][j]);
			}}
		}
		std::swap(phi, phi_new); /* Using swap which is faster than reassigning values again */
	}


	std::cout << "Writing to " << filename << ".\n";
	/*	Data file output	*/
	std::ofstream datafile;
	datafile.open(filename);
	datafile << "# This is a data file of potential over a plane\n";
	datafile << "# X	Y	V\n";
	for(int i = 0; i < N+1; i++){
		for(int j = 0; j < N+1; j++){
			datafile << i << "	" << j << "	" << phi[i][j] << "\n";
		}
		datafile << "\n";
	}
	return 0;
}
