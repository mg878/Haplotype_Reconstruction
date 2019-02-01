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

//Calculating the Dirichlet Multinomial given C=noise parameter, n=reads, x=frequencies of haplotypes
double Likelihood(int c, vector<double> n, vector<double> x)
{
	double sum_1 = 0;
	for (unsigned int i = 0; i < n.size(); i++)
	{
		sum_1 += n[i];
	}
	double sum_2 = 0;
	for (unsigned int i = 0; i < n.size(); i++)
	{
		sum_2 += gsl_sf_lngamma(n[i] + 1);
	}
	double sum_3 = 0;
	for (unsigned int i = 0; i < x.size(); i++)
	{
		sum_3 += (c*x[i]);
	}
	double sum_4 = 0;
	for (unsigned int i = 0; i < x.size(); i++)
	{
		sum_4 += (n[i] + c * x[i]);
	}
	double sum_5 = 0;
	for (unsigned int i = 0; i < n.size(); i++)
	{
		sum_5 += gsl_sf_lngamma(n[i] + c * x[i]);
	}
	double sum_6 = 0;
	for (unsigned int i = 0; i < n.size(); i++)
	{
		sum_6 += gsl_sf_lngamma(c*x[i]);
	}
	return gsl_sf_lngamma(sum_1 + 1) - sum_2 + gsl_sf_lngamma(sum_3) - gsl_sf_lngamma(sum_4) + sum_5 - sum_6;
}

//Given a timepoint (i.e 0 for before and 1 for after transmission), this function finds the likelihood of a set of haplotypes with their contribution
//in matching with the partial reads and their corresponding frequencies.
double Dirichlet_m(double timepoint, vector<string> haplotypes, vector<double> freq, vector<vector<vector<int>>> &CONTRIBS, vector<vector<vector<string>>> &partitions)
{
	if (timepoint == 0)
	{
		timepoint = 3;
	}
	else if (timepoint == 1)
	{
		timepoint = 1;
	}
	vector<string> candidates = haplotypes;
	vector<double> q = freq;
	/*for (int i = 0; i < q.size(); i++)
	{
	cout << endl << endl << "BEGINNIG: " << q[i] << endl;
	}*/
	double ULTIMA_THULE = 0;
	int SIZE_1 = candidates.size();
	for (unsigned int i = 0; i<CONTRIBS.size(); i++)
	{
		vector<double> inf;
		vector<double> nn;
		vector<double> temp = q;
		temp.erase(remove(temp.begin(), temp.end() - 1, q[SIZE_1]));
		unsigned int check = 0;
		vector<double> majorCheckvec;
		double NUM1 = 0;
		double sum_nn = 0;
		unsigned int majorCheck1 = 0;
		vector<double> vaccum_holder;
		for (unsigned int j = 0; j<(CONTRIBS[i]).size(); j++)
		{
			double sum = 0;
			int majorCheck = 0;
			for (unsigned int k = 0; k<(CONTRIBS[i][j]).size(); k++)
			{
				if (!(CONTRIBS[i][j]).empty())
				{
					sum += q[(CONTRIBS[i][j][k])];
					majorCheck++;
					majorCheck1++;
					majorCheckvec.push_back(q[(CONTRIBS[i][j][k])]);
				}
			}
			if ((CONTRIBS[i][j]).empty())
			{
				check++;
				double SIZE1_1 = (partitions[i][j]).size();
				NUM1 = atoi((partitions[i][j]).at(SIZE1_1 - timepoint).c_str());
				sum_nn += NUM1;
				vaccum_holder.push_back(NUM1);
			}
			if (majorCheck != 0)
			{
				inf.push_back(sum);
				double SIZE = (partitions[i][j]).size();
				double NUM = atoi((partitions[i][j]).at(SIZE - timepoint).c_str());
				nn.push_back(NUM);
			}
		}
		if (check == 0 && majorCheck1 == temp.size())
		{
			inf.push_back(q[SIZE_1]);
			nn.push_back(0.0);
		}
		unsigned int check2 = 0;
		if (check == 0 && majorCheck1 != temp.size())
		{
			for (unsigned int m = 0; m < majorCheckvec.size(); m++)
			{
				temp.erase(remove(temp.begin(), temp.end() - 1, majorCheckvec[m]));
			}
			double temp_sum = 0;
			for (unsigned int g = 0; g < temp.size(); g++)
			{
				temp_sum += temp[g];
			}
			inf.push_back(temp_sum);
			nn.push_back(0.0);
			inf.push_back(q[SIZE_1]);
			nn.push_back(0.0);
			check2++;
		}
		unsigned int check3 = 0;
		if (check != 0 && majorCheck1 != temp.size() && majorCheck1 != 0 && check2 == 0)
		{
			for (unsigned int m = 0; m < majorCheckvec.size(); m++)
			{
				temp.erase(remove(temp.begin(), temp.end() - 1, majorCheckvec[m]));
			}
			double temp_sum = 0;
			for (unsigned int g = 0; g < temp.size(); g++)
			{
				temp_sum += temp[g];
			}
			double total = accumulate(vaccum_holder.begin(), vaccum_holder.end(), 0);
			inf.push_back(temp_sum);
			nn.push_back(0.0);
			inf.push_back(q[SIZE_1]);
			nn.push_back(total);
			check3++;
		}
		if (check != 0 && majorCheck1 != temp.size() && majorCheck1 == 0 && check2 == 0 && check3 == 0)
		{
			double temp_sum = 0;
			for (unsigned int g = 0; g < temp.size(); g++)
			{
				temp_sum += temp[g];
			}
			double total = accumulate(vaccum_holder.begin(), vaccum_holder.end(), 0);
			inf.push_back(q[SIZE_1]);
			nn.push_back(total);
			inf.push_back(temp_sum);
			nn.push_back(0.0);
		}

		if (check != 0 && majorCheck1 == temp.size() && check2 == 0 && check3 == 0)
		{
			inf.push_back(q[SIZE_1]);
			nn.push_back(sum_nn);
		}
		ULTIMA_THULE += Likelihood(200, nn, inf);
		temp = q;
	}
	return ULTIMA_THULE;
}

//For a set of haplotypes (i.e. vector candidates), this function produces their optimum frequencies given short-read data.
double OptimumFreqs(double num, vector<string> &candidates, vector<double> &init_freqs, vector<vector<vector<int>>> &CONTRIBS, vector<vector<vector<string>>> &partitions)
{
	gsl_rng_env_setup();
	gsl_rng *rgen = gsl_rng_alloc(gsl_rng_taus);
	int seed = (int)time(NULL);
	gsl_rng_set(rgen, seed);
	double init_size = init_freqs.size();
	vector<double> init_freqs_store = init_freqs;
	double L = -1e9;
	unsigned int check = 1;
	unsigned int max_it = 1e4;
	double changex = 1e-2;
	double L_store1 = -1e9;

	for (unsigned int it = 0; it < max_it; it++)
	{
		double r = 2 * (gsl_rng_uniform(rgen) - 0.5);
		int j = floor(init_freqs.size()*gsl_rng_uniform(rgen));
		init_freqs[j] = init_freqs[j] + (r*changex);

		for (unsigned int i = 0; i < init_size; i++)
		{
			if (init_freqs[i] < 1e-11)
			{
				init_freqs[i] = 1e-11;
			}
		}

		if (init_freqs[init_size - 1] > 1e-2)
		{
			init_freqs[init_size - 1] = 1e-2;
		}
		NormaliseFreqs(init_freqs);
		L = Dirichlet_m(num, candidates, init_freqs, CONTRIBS, partitions);

		if (L > L_store1)
		{
			L_store1 = L;
			init_freqs_store = init_freqs;
			check = 1;
		}
		if (L <= L_store1)
		{
			init_freqs = init_freqs_store;
			check++;
		}
		if (check >= 50)
		{
			break;
		}
	}
	return L_store1;
}



int main(int argc, char* argv[])
{
	string name1 = argv[1];
	if (FILE *file1 = fopen(name1.c_str(), "r")) 
	{
        fclose(file1);
		vector<string> table_2;
		ifstream test_2;
		//Change the directory below
		test_2.open(argv[2]);
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
		for (unsigned int i = 0; i < rows_2.size(); i++)
		{
			temp_2 = (rows_2[i]);
			unsigned int sz = temp_2.size();
			string temp_2 = rows_2[i][sz-3];
			candidates.push_back(temp_2);
		}
		
		vector<double> f_before;
		vector<double> f_after;
		vector<string> temp;
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
			f_before.push_back(num_1);
			string temp_2 = rows_2[i][sz-1];
			double num_2 = atof(temp_2.c_str());
			if (num_2 < 1e-10)
			{
				num_2 = 0;
			}
			f_after.push_back(num_2);
		}
		
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
		vector < vector<string> > rows;
		for (unsigned int i = 0; i < table.size(); i++)
		{
			vector<string> row;
			istringstream iss(table[i]);
			string term;
			unsigned int check_naser = 0;
			while (iss >> term)
			{
				row.push_back(term);
				if (term == "X")
				{
					check_naser++;
					row.clear();
					break;
				}
			}
			if (!row.empty())
			{
				rows.push_back(row);
			}
		}

		vector<string> loci;
		vector<string> tempak;
		for (unsigned int i = 0; i < rows.size(); i++)
		{
			tempak = (rows[i]);
			unsigned int sz = tempak.size();
			for (unsigned int j = 1; j < sz - 6; j++)
			{
				loci.push_back(rows[i][j]);
			}
		}
		
		int TOTAL_N = 0;
		for (unsigned int i=0; i<rows.size(); i++)
		{
			unsigned int sz = rows[i].size();
			string temp_1 = rows[i][sz-1];
			string temp_2 = rows[i][sz-3];
			int num_1 = atoi(temp_1.c_str());
			int num_2 = atoi(temp_2.c_str());
			TOTAL_N += num_1 + num_2;
		}
		
		vector<int> positions;
		for (unsigned int i = 0; i < loci.size(); i++)
		{
			int num = atoi(loci.at(i).c_str());
			positions.push_back(num);
		}


		sort(positions.begin(), positions.end());
		positions.erase(unique(positions.begin(), positions.end()), positions.end());


		vector<string> Fpositions;
		for (int i : positions)
		{
			Fpositions.push_back(to_string(i));
		}

		vector<vector<string>> locus = rows;
		vector<string> temp1;
		for (unsigned int i = 0; i < rows.size(); i++)
		{
			temp1 = (rows[i]);
			unsigned int sz1 = temp1.size();
			for (unsigned int j = 1; j < sz1 - 6; j++)
			{
				int count = -1;
				for (string k : Fpositions)
				{
					count++;
					if (rows[i][j] == k)
					{
						locus[i][j] = " ";
						locus[i][j] = to_string(count);
						break;
					}
				}
			}
		}

		vector<vector<vector<string>>> partitions;
		for (unsigned int k = 1; k <= Fpositions.size(); k++)
		{
			for (unsigned int i = 0; i < Fpositions.size(); i++)
			{
				vector<string> pick2;
				vector<vector<string>> partialPar2;
				for (unsigned int j = 0; j < locus.size(); j++)
				{
					pick2 = locus[j];
					string two = to_string(k);
					unsigned int count_checker = 0;
					if (pick2[0] == two)
					{
						for (unsigned int m = 0; m < k; m++)
						{
							if (pick2[m + 1] == to_string(i + m))
							{
								count_checker++;
							}
						}
						if (count_checker == k)
						{
							partialPar2.push_back(pick2);
						}
					}
				}
				if (!partialPar2.empty())
				{
					partitions.push_back(partialPar2);
				}
			}
		}

		//You can uncomment the following lines if you want to check whether the read-data looks OK and the partitioning happened correctly.
		int read_count = 0;
		for (unsigned int i = 0; i<partitions.size(); i++)
		{
			for (unsigned int j = 0; j<(partitions[i]).size(); j++)
			{
				for (unsigned int k = 0; k<(partitions[i][j]).size(); k++)
				{
					//cout << partitions[i][j][k] << "  ";
					read_count++;
				}
				//cout << endl << "--------------------------------------" << endl;
			}
			//cout << "************************************" << endl;
		}

		vector<string> tampa;
		for (unsigned int i = 0; i < (partitions).size(); i++)
		{
			double SIZE0 = (partitions[i]).size();
			for (unsigned int j = 0; j < SIZE0 - 1; j++)
			{
				double SIZE1 = (partitions[i][j]).size();
				double NUM1 = atoi((partitions[i][j]).at(SIZE1 - 1).c_str());
				double SIZE2 = (partitions[i][j + 1]).size();
				double NUM2 = atoi((partitions[i][j + 1]).at(SIZE2 - 1).c_str());
				if (NUM2 > NUM1)
				{
					tampa = partitions[i][j];
					partitions[i][j] = partitions[i][j + 1];
					partitions[i][j + 1] = tampa;
				}
			}
		}

		
		vector<double> freqs_before;
		vector<double> freqs_after;
		//You can change the cut-off limit in the log-likelihood improvement that you demand every time the computer adds a new haplotype 
		//In this case, we demand the improvement to be, at least, 5 log-likelihood units.
		//The following loop determines how many times the computer makes an attempt to find the optimum haplotypes and their corresponding frequencies.
		//In this case, the computer makes 20 attempts to optimize candidates.size() haplotypes.
		vector<double> init_freqs1;
		for (unsigned int i = 0; i < candidates.size() + 1; i++)
		{
			if (i == 0)
			{
				init_freqs1.push_back(1 - ((candidates.size())*1e-2));
			}
			else if (i != 0)
			{
				init_freqs1.push_back(1e-2);
			}
		}
		vector<double> init_freqs2;
		for (unsigned int i = 0; i < candidates.size() + 1; i++)
		{
			if (i == 0)
			{
				init_freqs2.push_back(1 - ((candidates.size())*1e-2));
			}
			else if (i != 0)
			{
				init_freqs2.push_back(1e-2);
			}
		}


		vector<vector<vector<int>>> CONTRIBS;
		for (unsigned int i = 0; i < partitions.size(); i++)
		{
			vector<vector<int>> contribs;
			for (unsigned int j = 0; j < (partitions[i]).size(); j++)
			{
				vector<int> Contribs;
				vector<string> tempa;
				tempa = partitions[i][j];
				unsigned int sz = tempa.size();
				int num = sz - 6;
				for (unsigned int s = 0; s < candidates.size(); s++)
				{
					int counter = 0;
					for (int k = 1; k < num; k++)
					{
						int tem = atoi((tempa[k]).c_str());
						if ((candidates[s])[tem] == (tempa[num])[k - 1])
						{
							counter++;
						}
					}
					int dum = atoi((tempa[0]).c_str());
					if (counter == dum)
					{
						Contribs.push_back(s);
					}
				}
				contribs.push_back(Contribs);
			}
			CONTRIBS.push_back(contribs);
		}

		gsl_rng_env_setup();
		gsl_rng *rgen = gsl_rng_alloc(gsl_rng_taus);
		int seed = (int)time(NULL);
		gsl_rng_set(rgen, seed);

		OptimumFreqs(0, candidates, init_freqs1, CONTRIBS, partitions);
		OptimumFreqs(1, candidates, init_freqs2, CONTRIBS, partitions);
		//L_store = L_store1 + L_store2;
		//cout << endl << L_store1 << " + " << L_store2 << " = " << L_store << endl;
		//cout << L_store << endl;
		
		/*for (unsigned int i = 0; i < f_after.size(); i++)
		{
			cout << init_freqs2[i]  << "\t" << f_after[i] << endl << init_freqs1[i] << "\t" << f_before[i] << endl;
		}
		
		cout << endl << "****************************************" << endl;*/
		//SetCurrentDirectory(argv[3]);
		chdir(argv[3]);
		ofstream myfile;
		myfile.open(argv[4]);
		for (unsigned int i = 0; i < f_after.size(); i++)
		{
			myfile << i << "\t" << init_freqs1[i] << "\t" << init_freqs2[i] << endl;
		}
		myfile.close();

		return 0;
	}
	
	else
	{
		return false;
	}
}
