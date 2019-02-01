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

struct container_is_empty_t {
	template<class C>
	bool operator()(C const& c)const {
		return c.empty();
	}
	template<class T, size_t N>
	bool operator()(T(&)[N])const {
		return N > 0;
	}
};
static container_is_empty_t const container_is_empty;

int main(int argc, char* argv[])
{
	
	vector<string> table_2;
	ifstream test_2;
	//Change the directory below
	test_2.open(argv[1]);
	while (!test_2.eof())
	{
		string line_2;
		getline(test_2, line_2);
		table_2.push_back(line_2);
	}
	test_2.close();
	vector < vector<string> > rows_2;
	for (unsigned int i = 0; i < table_2.size(); i++)
	{
		vector<string> row_2;
		istringstream iss(table_2[i]);
		string term_2;
		while (iss >> term_2)
		{
			row_2.push_back(term_2);
		}
		if (!row_2.empty())
		{
			rows_2.push_back(row_2);
		}
	}

	vector<string> candidates;
	vector<string> temp_2;
	for (unsigned int i = 1; i < 7; i++)
	{
		temp_2 = (rows_2[i]);
		unsigned int sz = temp_2.size();
		string temp_2 = rows_2[i][0];
		candidates.push_back(temp_2);
	}
	
	vector<double> f_before;
	vector<double> f_after;
	vector<string> temp;
	for (unsigned int i = 1; i < 7; i++)
	{
		temp = (rows_2[i]);
		unsigned int sz = temp.size();
		string temp_1 = rows_2[i][sz-3];
		double num_1 = atof(temp_1.c_str());
		f_before.push_back(num_1);
		string temp_2 = rows_2[i][sz-1];
		double num_2 = atof(temp_2.c_str());
		f_after.push_back(num_2);
	}
	

	chdir(argv[2]);
	ofstream myfile;
	myfile.open(argv[3]);
	for (unsigned int i = 0; i < f_after.size(); i++)
	{
		cout << i+1 << "\t" << candidates[i] << "\t" << f_before[i] << "\t" << f_after[i] << endl;
		myfile << i+1 << "\t" << candidates[i] << "\t" << f_before[i] << "\t" << f_after[i] << endl;
	}
	myfile.close();

	return 0;
}
