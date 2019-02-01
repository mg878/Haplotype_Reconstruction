#include <algorithm>
#include <iostream>
#include <list>
#include <string>
#include <array>
#include <fstream>
#include <math.h> 
#include <vector>
#include <sstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <numeric>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_pow_int.h>
#include <time.h> 
#include <stdio.h>
#include <unistd.h>
using namespace std;


int main(int argc, char* argv[])
{

	vector<string> table;
	ifstream bottleneck;
	//Change the directory below
	bottleneck.open(argv[1]);
	while (!bottleneck.eof())
	{
		string line;
		getline(bottleneck, line);
		table.push_back(line);
	}
	bottleneck.close();
	
	double row = table.size() - 1;
	chdir(argv[2]);
	ofstream myfile;
	myfile.open(argv[3]);
	cout << row << endl;
	myfile << row << endl;
	myfile.close();


	return 0;
} 