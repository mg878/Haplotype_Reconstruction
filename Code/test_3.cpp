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

//Normalizing the frequencies after changing them randomly
void NormaliseFreqs(vector<double> &init_freqs)
{
	double tot = 0;
	for (unsigned int i = 0; i<init_freqs.size(); i++) {
		tot = tot + init_freqs[i];
	}
	for (unsigned int i = 0; i<init_freqs.size(); i++) {
		init_freqs[i] = init_freqs[i] / tot;
	}
}


int main(int argc, char* argv[])
{
	
	string name1 = argv[1];
	if (FILE *file1 = fopen(name1.c_str(), "r")) 
	{
			fclose(file1);
		vector<string> table;
		ifstream test;
		//Change the directory below
		test.open(argv[1]);
		while (!test.eof())
		{
			string line;
			getline(test, line);
			table.push_back(line);
		}
		test.close();
		
		vector < vector<string> > rows_2;
		for (unsigned int i = 0; i < table.size(); i++)
		{
			vector<string> row_2;
			istringstream iss(table[i]);
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

		vector<double> f_before;
		vector<double> f_after;
		vector<string> temp;
		vector<string> candidates;
		for (unsigned int i = 0; i < rows_2.size(); i++)
		{
			temp = (rows_2[i]);
			unsigned int sz = temp.size();
			string temp_1 = rows_2[i][sz-2];
			double num_1 = atof(temp_1.c_str());
			if (num_1 < 1e-10)
			{
				num_1 = 0;
			}
			string temp_2 = rows_2[i][sz-1];
			double num_2 = atof(temp_2.c_str());
			if (num_2 < 1e-10)
			{
				num_2 = 0;
			}
			if (num_1 < 1e-7 && num_2 > 1e-5)
			{
			cout << "mutation: " << num_1 << "\t" << num_2 << endl;
			}
			else if (num_1 != num_2 && num_1 != 1)
			{
				f_before.push_back(num_1);
				f_after.push_back(num_2);
				candidates.push_back(rows_2[i][sz-3]);
			}
		}
		
		NormaliseFreqs(f_after);
		NormaliseFreqs(f_before);
		
		for (unsigned int i = 0; i < f_after.size(); i++)
		{
			cout << i+1 << "\t" << candidates[i] << "\t" << f_before[i] << "\t" << f_after[i] << endl;
		}
		
		chdir(argv[2]);
		if (f_after.size() > 0)
		{
			ofstream myfile;
			myfile.open("outcome_1.txt");
			for (unsigned int i = 0; i < f_after.size(); i++)
			{
				myfile << i+1 << "\t" << candidates[i] << "\t" << f_before[i] << "\t" << f_after[i] << endl;
			}
			myfile.close();
		}
	
		return 0;
	}
	
	else
	{
		return false;
	}
}
