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

double Average(vector<double> init_freqs)
{
	double tot = 0;
	for (unsigned int i = 0; i<init_freqs.size(); i++) {
		tot = tot + init_freqs[i];
	}
	
	return tot/(init_freqs.size());
}

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
	
	vector<double> real_freq_b;
	vector<double> real_freq_a;
	vector<string> temp;
	for (unsigned int i = 0; i < rows_2.size(); i++)
	{
		temp = (rows_2[i]);
		unsigned int sz = temp.size();
		string tot_before = rows_2[i][sz-2];
		double tb = atof(tot_before.c_str());
		string tot_after = rows_2[i][sz-1];
		double ta = atof(tot_after.c_str());
		string temp_1 = rows_2[i][sz-6];
		double num_1 = atof(temp_1.c_str());
		string temp_2 = rows_2[i][sz-5];
		double num_2 = atof(temp_2.c_str());
		string temp_3 = rows_2[i][sz-4];
		double num_3 = atof(temp_3.c_str());
		string temp_4 = rows_2[i][sz-3];
		double num_4 = atof(temp_4.c_str());
		real_freq_b.push_back(num_1/tb);
		real_freq_a.push_back(num_2/ta);
		real_freq_b.push_back(num_3/tb);
		real_freq_a.push_back(num_4/ta);
	}
	
		
	/*double average = accumulate( bottleneck_list.begin(), bottleneck_list.end(), 0.0)/bottleneck_list.size(); 
	chdir(argv[2]);
	ofstream myfile;
	myfile.open(argv[3]);
	myfile << 1/average << endl;
	cout << 1/average << endl;
	myfile.close();*/
	
	/*for (unsigned int i = 0; i < real_freq_b.size(); i++)
	{
		myfile << i + 1 << "\t" << candidates[i] << "\t" << real_freq_b[i] << "\t" << real_freq_a[i] << endl;
	}*/
	chdir(argv[2]);
	ofstream myfile;
	myfile.open(argv[3]);
	for (unsigned int i = 0; i < real_freq_a.size(); i++)
	{
		myfile << i << "\t" << real_freq_b[i] << "\t" << real_freq_a[i] << endl;
	}
	myfile.close();

	return 0;
} 