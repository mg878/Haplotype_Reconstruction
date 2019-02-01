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
	
	vector < vector<string> > rows;
	for (unsigned int i = 0; i < table.size(); i++)
	{
		vector<string> row;
		istringstream iss(table[i]);
		string term;
		while (iss >> term)
		{
			row.push_back(term);
		}
		if (!row.empty())
		{
			rows.push_back(row);
		}
	}

	vector<double> f_after;
	vector<string> temp;
	for (unsigned int i = 0; i < rows.size(); i++)
	{
		temp = (rows[i]);
		unsigned int sz = temp.size();
		string temp = rows[i][sz-2];
		double num_1 = atof(temp.c_str());
		candidates.push_back(temp);
	}
	
	
	double row = table.size() - 1;
	chdir(argv[2]);
	ofstream myfile;
	myfile.open(argv[3]);
	cout << row << endl;
	myfile << row << endl;
	myfile.close();


	return 0;
} 