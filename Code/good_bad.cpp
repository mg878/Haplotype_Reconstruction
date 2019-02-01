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

	chdir(argv[1]);
	vector<string> table_1;
	ifstream test_1;
	test_1.open(argv[2]);
	while (!test_1.eof())
	{
		string line_1;
		getline(test_1, line_1);
		table_1.push_back(line_1);
	}
	test_1.close();
	vector < vector<string> > rows_1;
	for (unsigned int i = 0; i < table_1.size(); i++)
	{
		vector<string> row_1;
		istringstream iss(table_1[i]);
		string term_1;
		while (iss >> term_1)
		{
			row_1.push_back(term_1);
		}
		if (!row_1.empty())
		{
			rows_1.push_back(row_1);
		}
	}
	
	vector<string> hap_1;
	for (unsigned int i = 0; i < rows_1.size(); i++)
	{
		hap_1.push_back(rows_1[i][1]);		
	}
	
	vector<string> table_2;
	ifstream test_2;
	test_2.open(argv[3]);
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
	
	vector<string> hap_2;
	for (unsigned int i = 0; i < rows_2.size(); i++)
	{
		hap_2.push_back(rows_2[i][1]);		
	}
	
	int check = 0;
	for (unsigned int i = 0; i < hap_1.size(); i++)
	{
		//cout << "i" << "\t"  << i << "\t" << hap_1[i] << endl;
		for (unsigned int j = 0; j < hap_2.size(); j++)
		{
			//cout  << "j" << "\t" << j << "\t" << hap_2[j] << endl;
			if (hap_1[i] == hap_2[j])
			{
				check++;
				//cout << "************************************************" << endl;
				//cout << hap_1[i] << "\t" << hap_2[j] << endl;
			}
		}
	}
	ofstream myfile;
	myfile.open("number.txt");
	//cout << check << endl;
	//myfile << check << endl;
	cout << hap_1.size() << endl;
	myfile << hap_1.size() << endl;
	myfile.close();


	return 0;
} 