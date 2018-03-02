/*
 * */

#include <iostream> /* System IO */
#include <vector> /* Using vectors instead of arrays */
#include <fstream> /* Data file output */
#include <string> /* Argument Handling */
#include <cmath>

typedef std::vector<double> double_vec; /* Simplification of code */
typedef std::vector<double_vec> double_vec_vec;

void print_help_text(char* arg){
	std::cout << "Usage:\n" << arg << " [options]\n\nOptions:\n"
		<< "  -h,  --help         Shows the help prompt\n"
		<< "  -3D, --3demensions  Activates 3 spacial demensions\n"
		<< "  -stl --STLOBJECT    Activates usage of STL files in 3D\n"
		<< "                      Needs 3D activated to work\n"
		<< "  -Ef, --field        Outputs electric field approximation instead of voltage\n"
		<< "  -n,  --size         Number of points on one edge of the surface (>2). Default: 25\n"
		<< "  -l,  --length       Length of one edge of the surface in meters. Default: 1.0\n"
		<< "                      Has no effect without field approximation\n"
		<< "  -c,  --count        The number of iterations that will be performed. Default: 100\n"
		<< "  -b,  --boundaries   Sets the boundary conditions (V) Takes 4 or 6 values\n"
		<< "                      2D: Top Right Bottom Left. Default: 0 0 0 0\n"
		<< "                      3D: Top Right Bottom Left Front Back. Default: 0 0 0 0 0 0\n"
		<< "  -f,  --filename     Sets the file name. Default: \"data.dat\"\n"
		<< "\n"
		<< "  STL OBJECT INSERTION\n"
		<< "    * Starts with \"[\" and ends with \"]\"\n"
		<< "    * Takes arguments within braces\n"
		<< "  -t,                 Sets the filename of STL File\n"
		<< "  -s,                 Sets the scale of the file : 0 < s < 1\n"
		<< "  -x, -y, -z          Sets the x, y, and z shift of the object : 0 < x,y,z < 1\n";
}

struct STL {
	std::string name;
	double scale;
	double x, y, z;
	double volt;
    STL(std::string p_name, double p_scale, double p_x, double p_y, double p_z, double p_volt)
        : name(std::move(p_name)), scale(p_scale), x(p_x), y(p_y), z(p_z), volt(p_volt){}
};

int main(int argc, char* argv[])
{
	/*	Variable Initialization	*/

	double Top, Right, Bottom, Left, Front, Back, Length, Delta;
	int N, Num;
	std::string filename, stlFile;
	const double pi = 3.14159265359;
	bool CUBE = 0, EFIELD = 0, STLOBJECT = 0;

	/* Default Values */
	Top = 0;
	Right = 0;
	Bottom = 0;
	Left = 0;
	Front = 0;
	Back = 0;
	N = 25;
	Num = 100;
	filename = "data.dat";
	stlFile = "shape.stl";
	Length = 1.0;
	std::vector<STL> objects;

	/*	Argument Handling (Needs Improvemnt)	*/
	for (int i = 0; i < argc; i++){
		if((std::string(argv[i]) == "-3D") || (std::string(argv[i]) == "--3demensions")){
			CUBE = 1;
		}
		if((std::string(argv[i]) == "-Ef") || (std::string(argv[i]) == "--field")){
			EFIELD = 1;
		}
		if((std::string(argv[i]) == "-stl") || (std::string(argv[i]) == "--stl-object")){
			STLOBJECT = 1;
		}
	}
	for (int i = 0; i < argc; i++){
		if((std::string(argv[i]) == "-h") || (std::string(argv[i]) == "--help")){
			std::cout << "This program uses Poisson's equation to map potential in 2-space when given boundary conditions\nIt can then output an electric field aproximation if desired\n";
			print_help_text(argv[0]);
			return 1;
		} else if ((std::string(argv[i]) == "-n") || (std::string(argv[i]) == "--number")){
			if (i + 1 < argc)
			{
				try 
				{
					N = std::stoi( argv[i + 1], nullptr );
				}
				catch ( const std::invalid_argument& ia)
				{
					std::cout << "Unable to convert, " << argv[i + 1] << ", into an integer!\n";
					print_help_text(argv[0]);
					return 1;
				}
			}
			if (not (i + 1 < argc) || (N < 3))
			{
				std::cout << "Argument, " << argv[i] << ", needs an integer size value greater than 2.\n";
				print_help_text( argv[0] );
				return 1;
			}
		} 
		else if ((std::string(argv[i]) == "-c") || (std::string(argv[i]) == "--count"))
		{
			if(i + 1 < argc){
				try 
				{
					Num = std::stoi( argv[i + 1], nullptr );
				}
				catch ( const std::invalid_argument& ia)
				{
					std::cout << "Unable to convert, " << argv[i + 1] << ", into an integer!\n";
					print_help_text( argv[0] );
					return 1;
				}
			}
			if( not (i + 1 < argc) || (Num < 1)){
				std::cout << "Argument,  " << argv[i] << ", needs an integer iteration count.\n";
				print_help_text(argv[0]);
				return 1;
			}
		} 
		else if ( (std::string( argv[i] ) == "-b") || (std::string(argv[i]) == "--boundaries") )
		{
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
			if( i + 4 < argc)
			{
				try 
				{
					Top = std::stod( argv[i + 1], nullptr );
				}
				catch ( const std::invalid_argument& ia)
				{
					std::cout << "Error: \"" << argv[i + 1] << "\" is not a number!\n";
					print_help_text( argv[0] );
					return 1;
				}
				try 
				{
					Right = std::stod( argv[i + 2], nullptr );
				}
				catch ( const std::invalid_argument& ia)
				{
					std::cout << "Error: \"" << argv[i + 2] << "\" is not a number!\n";
					print_help_text( argv[0] );
					return 1;
				}
				try 
				{
					Bottom = std::stod( argv[i + 3], nullptr );
				}
				catch ( const std::invalid_argument& ia)
				{
					std::cout << "Error: \"" << argv[i + 3] << "\" is not a number!\n";
					print_help_text( argv[0] );
					return 1;
				}
				try 
				{
					Left = std::stod( argv[i + 4], nullptr );
				}
				catch ( const std::invalid_argument& ia)
				{
					std::cout << "Error: \"" << argv[i + 4] << "\" is not a number!\n";
					print_help_text( argv[0] );
					return 1;
				}
				if((i + 6 < argc) & (CUBE))
				{
					try
					{
						Front = std::stod( argv[i + 5], nullptr);
					}
					catch ( const std::invalid_argument& ia)
					{
						std::cout << "Error: " << argv[i + 5] << " is not a number!\n";
						print_help_text( argv[0] );
						return 1;
					}
					try
					{
						Back = std::stod( argv[i + 6], nullptr);
					}
					catch ( const std::invalid_argument& ia)
					{
						std::cout << "Error: " << argv[i + 5] << " is not a number!\n";
						print_help_text( argv[0] );
						return 1;
					}
					i = i + 2;
				}
				i = i + 4;	
			}
		} else if ((std::string(argv[i]) == "-f") || (std::string(argv[i]) == "--datafile")){
			if(i + 1 < argc){
				filename = argv[i + 1];
			} else {
				std::cout << "Argument, " << argv[i] << ", needs a filename following it.\n";
				print_help_text(argv[0]);
				return 1;
			}
		} else if((std::string(argv[i]) == "-l") || (std::string(argv[i]) == "--length")){
			if(i + 1 < argc){
				try
				{
					Length = std::stod( argv[i + 1], nullptr);
				}
				catch ( const std::invalid_argument& ia)
				{
					std::cout << "Error: \"" << argv[i + 1] << "\" is not a number!\n";
					print_help_text(argv[0]);
					return 1;
				}
			}
			if( (not (i + 1 < argc)) || (Length == 0) ){
				std::cout << "Argument, " << argv[i] << ", needs a number greater than zero following it.\n";
				print_help_text(argv[0]);
				return 1;
			}
		} else if((std::string(argv[i]) == "-stl") || (std::string(argv[i]) == "--stl-object")){
			if(i + 2 < argc){
				if(std::string(argv[i+1]) == "["){
					double scale = 0.5, x = 0.5, y = 0.5, z = 0.5, volt = 1;
					std::string name = "object.stl";
					for(int i_stl = i + 1; i_stl < argc; i_stl++){
						std::string arg_stl = argv[i_stl];
						if((arg_stl == "-s") | (arg_stl == "--scale")){
							if(i_stl + 1 < argc){
								try
								{
									scale = std::stod( argv[i_stl + 1], nullptr );
								}
								catch ( const std::invalid_argument& ia)
								{
								std::cout << "Unable to convert, " << argv[i_stl + 1] << ", into an number!\n";
								print_help_text(argv[0]);
								return 1;
								}
								if((scale <= 0) | (scale >= 1)){
									std::cout << "Error: " << argv[i_stl] << " needs to fit 0 < " << argv[i_stl] << " < 1!\n";
									return 1;
								}
							} else {
								std::cout << "Argument, " << argv[i_stl] << ", needs a number greater than zero and less than one following it.\n";
								print_help_text(argv[0]);
								return 1;
							}
						} else if (arg_stl == "-x"){
							if(i_stl + 1 < argc){
								try
								{
									x = std::stod( argv[i_stl + 1], nullptr );
								}
								catch ( const std::invalid_argument& ia)
								{
								std::cout << "Unable to convert, " << argv[i_stl + 1] << ", into an number!\n";
								print_help_text(argv[0]);
								return 1;
								}
								if((x <= 0) | (x >= 1)){
									std::cout << "Error: " << argv[i_stl] << " needs to fit 0 < " << argv[i_stl] << " < 1!\n";
									print_help_text(argv[0]);
									return 1;
								}
							} else {
								std::cout << "Argument, " << argv[i_stl] << ", needs a number greater than zero and less than one following it.\n";
								print_help_text(argv[0]);
								return 1;
							}
						} else if (arg_stl == "-y"){
							if(i_stl + 1 < argc){
								try
								{
									y = std::stod( argv[i_stl + 1], nullptr );
								}
								catch ( const std::invalid_argument& ia)
								{
								std::cout << "Unable to convert, " << argv[i_stl + 1] << ", into an number!\n";
								print_help_text(argv[0]);
								return 1;
								}
								if((y <= 0) | (y >= 1)){
									std::cout << "Error: " << argv[i_stl] << " needs to fit 0 < " << argv[i_stl] << " < 1!\n";
									print_help_text(argv[0]);
									return 1;
								}
							} else {
								std::cout << "Argument, " << argv[i_stl] << ", needs a number greater than zero and less than one following it.\n";
								print_help_text(argv[0]);
								return 1;
							}

						} else if (arg_stl == "-z"){
							if(i_stl + 1 < argc){
								try
								{
									z = std::stod( argv[i_stl + 1], nullptr );
								}
								catch ( const std::invalid_argument& ia)
								{
								std::cout << "Unable to convert, " << argv[i_stl + 1] << ", into an number!\n";
								print_help_text(argv[0]);
								return 1;
								}
								if((z <= 0) | (z >= 1)){
									std::cout << "Error: " << argv[i_stl] << " needs to fit 0 < " << argv[i_stl] << " < 1!\n";
									print_help_text(argv[0]);
									return 1;
								}
							} else {
								std::cout << "Argument, " << argv[i_stl] << ", needs a number greater than zero and less than one following it.\n";
								print_help_text(argv[0]);
								return 1;
							}
						} else if (arg_stl == "-v"){
							if(i_stl + 1 < argc){
								try
								{
									volt = std::stod( argv[i_stl + 1], nullptr );
								}
								catch ( const std::invalid_argument& is)
								{
								std::cout << "Unable to convert, " << argv[i_stl + 1] << ", into an number!\n";
								print_help_text(argv[0]);
								return 1;
								}
							} else {
								std::cout << "Argument, " << argv[i_stl] << ", needs a number following it.\n";
								print_help_text(argv[0]);
								return 1;
							}
						} else if (arg_stl == "-t"){
							if(i_stl + 1 < argc){
								name = argv[i_stl + 1];
							} else {
								std::cout << "Argument, " << argv[i_stl] << ", needs a filename following it.\n";
								print_help_text(argv[0]);
								return 1;
							}
						} else if (arg_stl == "]"){
							std::cout << "Creating Object: " << name << "\n";
							objects.emplace_back(name, scale, x, y, z, volt);
							break;
						}
					}
				} else {
				}
			}
		}



	}
	
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
	else if(not(STLOBJECT)){
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
	} else if(STLOBJECT){ /* STL File Input ---------------------------------------*/
		std::cout << "You chose to use an stl file\n";
		std::vector<std::vector<std::vector<bool>>> iter(N+1, std::vector<std::vector<bool>>(N+1, std::vector<bool>(N+1, 1)));
		std::vector<double_vec_vec> phi(N+1,double_vec_vec(N+1, double_vec(N+1,0)));
		char header[80];
		uint32_t numT;
		float normal[3];
		float vertex[3];
		uint16_t attrib;
		std::ofstream datafile;
		datafile.open(filename);
		datafile << "# This is a data file of potential over a plane\n";
		datafile << "# Variables: N=" << N << ", Top=" << Top << ", Right=" << Right << ", Bottom=" << Bottom << ", Left=" << Left << ", Front=" << Front << ", Back=" << Back << ", cycles=" << Num << "\n";
		datafile.close();

		for(unsigned int num_object = 0; num_object < objects.size(); num_object++){
		std::cout << "Started reading STL\n";
		std::ifstream file (objects[num_object].name, std::ios::in|std::ios::binary);
		if(file.is_open())
        {
                file.read(header,80);
                file.read(reinterpret_cast<char *>(&numT), sizeof(numT));
				std::vector<std::vector<std::vector<double>>> triangles(numT, std::vector<std::vector<double>>(3, std::vector<double>(3)));
				for(unsigned int i = 0; i < numT; i++){
					file.read(reinterpret_cast<char *>(&normal), sizeof(normal));
					for(int j = 0; j < 3; j++){
					file.read(reinterpret_cast<char *>(&vertex), sizeof(normal));
					triangles[i][j][0] = vertex[0];
					triangles[i][j][1] = vertex[1];
					triangles[i][j][2] = vertex[2];
					}
					file.read(reinterpret_cast<char *>(&attrib), sizeof(attrib));
				}
        file.close();
		std::cout << "Finished reading STL\n";

		double x_max, x_min, y_max, y_min, z_max, z_min = 0;
		for(unsigned int i = 0; i < numT; i++){
			for(int j = 0; j < 3; j++){
				if(triangles[i][j][0] <= x_min){x_min = triangles[i][j][0];}
				if(triangles[i][j][0] >= x_max){x_max = triangles[i][j][0];}
				if(triangles[i][j][1] <= y_min){y_min = triangles[i][j][1];}
				if(triangles[i][j][1] >= y_max){y_max = triangles[i][j][1];}
				if(triangles[i][j][2] <= z_min){z_min = triangles[i][j][2];}
				if(triangles[i][j][2] >= z_max){z_max = triangles[i][j][2];}
			}
		}
		double dx, dy, dz;
		double rescale = 1;
		dx = x_max - x_min;
		dy = y_max - y_min;
		dz = z_max - z_min;
	
		if((dx >= dy) & (dx >= dz)){rescale = 1/dx;}
		else if((dy >= dx) & (dy >= dz)){rescale = 1/dy;}
		else if((dz >= dy) & (dz >= dx)){rescale = 1/dz;}

		dx = N * dx * rescale * objects[num_object].scale;
		dy = N * dy * rescale * objects[num_object].scale;
		dz = N * dz * rescale * objects[num_object].scale;

		for(unsigned int i = 0; i < numT; i++){
			for(int j = 0; j < 3; j++){
				if(x_min < 0){triangles[i][j][0] = triangles[i][j][0] - x_min;}
				triangles[i][j][0] = N * rescale * objects[num_object].scale * triangles[i][j][0] + N * objects[num_object].x - dx/2;
				if(y_min < 0){triangles[i][j][1] = triangles[i][j][1] - y_min;}
				triangles[i][j][1] = N * rescale * objects[num_object].scale * triangles[i][j][1] + N * objects[num_object].y - dy/2;
				if(z_min < 0){triangles[i][j][2] = triangles[i][j][2] - z_min;}
				triangles[i][j][2] = N * rescale * objects[num_object].scale * triangles[i][j][2] + N * objects[num_object].z - dz/2;
			}
		}

		double side_dx, side_dy, side_dz, side_len, side_x, side_y, side_z, hypo_x, hypo_y, hypo_z, hypo_len, hypo_dx, hypo_dy, hypo_dz;
		for(unsigned int i = 0; i < numT; i++){
			side_dx = (triangles[i][1][0] - triangles[i][0][0]);
			side_dy = (triangles[i][1][1] - triangles[i][0][1]);
			side_dz = (triangles[i][1][2] - triangles[i][0][2]);
			side_len = std::pow((side_dx*side_dx + side_dy*side_dy + side_dz*side_dz),0.5);
			for(int side_t = 0; side_t < side_len; side_t++){
				side_x = side_dx*(side_t/side_len) + triangles[i][0][0];
				side_y = side_dy*(side_t/side_len) + triangles[i][0][1];
				side_z = side_dy*(side_t/side_len) + triangles[i][0][2];
				hypo_dx = (triangles[i][2][0] - side_x);
				hypo_dy = (triangles[i][2][1] - side_y);
				hypo_dz = (triangles[i][2][2] - side_z);

				hypo_len = std::pow((hypo_dx*hypo_dx + hypo_dy*hypo_dy + hypo_dz*hypo_dz),0.5);
				for(int hypo_t = 0; hypo_t < hypo_len; hypo_t++)
				{
					hypo_x = hypo_dx*(hypo_t/hypo_len) + side_x;
					hypo_y = hypo_dy*(hypo_t/hypo_len) + side_y;
					hypo_z = hypo_dz*(hypo_t/hypo_len) + side_z;
					try{
						iter.at(hypo_x).at(hypo_y).at(hypo_z) = 0;
						phi.at(hypo_x).at(hypo_y).at(hypo_z) = objects[num_object].volt;
					}
					catch(std::out_of_range){
					}
				}
			}	
		}
		std::ofstream datafile;
		datafile.open(filename, std::ofstream::app);
		datafile << "# Using: STLFile=" << objects[num_object].name << ", x_shift=" << N * objects[num_object].x << ", y_shift=" << N * objects[num_object].y << ", z_shift=" << N * objects[num_object].z << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << "\n\n";
		datafile.close();
		}}
	std::cout << "Started Stepping\n";

	/* Steping Action	*/
	const double omega = 2 / (1 + pi / N*N);	/* Used as a relaxation constant: N^2 instead of N : Found it works better */
	std::vector<double_vec_vec> phi_new(N+1,double_vec_vec(N+1, double_vec(N+1))); /* */
	for(int count = 0; count < Num; count++){
		for(int i = 1; i < N; i++){
			for(int j = 1; j < N; j++){
				for(int k = 1; k < N; k++){
					if(iter[i][j][k]){
					phi_new[i][j][k] = phi[i][j][k] + (omega / 6) * (phi[i+1][j][k] + phi_new[i-1][j][k] + phi[i][j+1][k] + phi_new[i][j-1][k] + phi[i][j][k+1] + phi_new[i][j][k-1] - 6 * phi[i][j][k]); /* Incorporates the third spacial demension */
				}}
			}
		}
		std::swap(phi, phi_new); /* Using swap which is faster than reassigning values again */
	}

	std::cout << "Writing to file\n";
	/*	Data file output	*/
	datafile.open(filename, std::ofstream::app);
	datafile << "# X	Y	Z	V\n";
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			for(int k = 0; k < N; k++){
				datafile << i << "	" << j << "	" << k << "	" << phi[i][j][k] << "\n";
			}
		}
		datafile << "\n";
	}
	} /* End of STL */
	/* -------------------------------------------------- */
	return 0;
}

