#include <algorithm>
#include <iostream>
#include <list>
#include <string>
#include <array>
#include <fstream>
#include <math.h> 
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
	double L = -1e7;
	unsigned int check = 1;
	unsigned int max_it = 8e3;
	double changex = 1e-2;
	double L_store1 = -1e7;

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

		NormaliseFreqs(init_freqs);
		if (init_freqs[init_size - 1] > 1e-2)
		{
			init_freqs[init_size - 1] = 1e-2;
			double TOTAL = 0;
			for (unsigned int i = 0; i < init_size-1; i++)
			{
				TOTAL += init_freqs[i];
			}
			for (unsigned int i = 0; i < init_size-1; i++)
			{
				init_freqs[i] = init_freqs[i]*(0.99/TOTAL);
			}
		}
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
	vector<string> temp;
	for (unsigned int i = 0; i < rows.size(); i++)
	{
		temp = (rows[i]);
		unsigned int sz = temp.size();
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

	vector<vector<string>> shuffle;
	for (unsigned int i = 0; i < partitions.size(); i++)
	{
		for (unsigned int j = 0; j < (partitions[i]).size(); j++)
		{
			vector<string> shuffle_holder;
			for (unsigned int k = 0; k < (partitions[i][j]).size(); k++)
			{
				shuffle_holder.push_back(partitions[i][j][k]);
			}
			shuffle.push_back(shuffle_holder);
		}
	}

	string unknown;
	for (unsigned int i = 0; i<positions.size(); i++)
	{
		unknown.append("X");
	}
	vector<string> candidates;
	double L_new = -1e7;
	double L_store = -1e7;
	double L_preserve_1 = -1e7 - 300;
	double L_preserve_2 = -1e7 - 1;
	double loop = 0;
	vector<double> all_likelihoods;
	vector<double> haplotypes;
	vector<double> freqs_before;
	vector<double> freqs_after;
	double BIC_check = 2;
	//You can change the cut-off limit in the log-likelihood improvement that you demand every time the computer adds a new haplotype 
	//In this case, we demand the improvement to be, at least, 5 log-likelihood units.
	while (BIC_check > 0)
	{
		loop++;
		candidates.push_back(unknown);
		unsigned int hap_size = candidates.size();
		L_preserve_1 = L_preserve_2;
		vector<double> temp_likelihoods;
		vector<vector<double>> freqs_1;
		vector<vector<double>> freqs_2;
		vector<vector<string>> temp_candidates;
		int candidate_size1 = candidates.size();
		//The following loop determines how many times the computer makes an attempt to find the optimum haplotypes and their corresponding frequencies.
		//In this case, the computer makes 20 attempts to optimize candidates.size() haplotypes.
		double sacred_check = 0;
		for (unsigned int sample = 0; sample < 10; sample++)
		{
			//?could make the condition below better?
			L_new = L_preserve_2 - 100 / (2*candidate_size1);
			L_store = -1e5;
			candidates[candidate_size1 - 1] = unknown;
			double candidate_size = candidates.size();
			double BAD_COUNT = 0;
			//This defines the cut-off on when to stop making further attempts to optimize the haplotypes
			//This depends on loop=number of haplotypes at any time and read_count=total number of partial-reads available.
			double THRESHOLD = 2 * loop * read_count;
			double attempts = 1;
			vector<double> init_freqs1;
			gsl_rng_env_setup();
			gsl_rng *rgen = gsl_rng_alloc(gsl_rng_taus);
			int seed = (int)time(NULL);
			gsl_rng_set(rgen, seed);
			for (unsigned int i = 0; i < candidates.size() + 1; i++)
			{
			
				init_freqs1.push_back(gsl_rng_uniform(rgen));			
			}
			
			NormaliseFreqs(init_freqs1);
			vector<double> init_freqs1_store = init_freqs1;
			
			vector<double> init_freqs2;
			for (unsigned int i = 0; i < candidates.size() + 1; i++)
			{
			
				init_freqs2.push_back(gsl_rng_uniform(rgen));			
			}
			
			NormaliseFreqs(init_freqs2);
			vector<double> init_freqs2_store = init_freqs2;
			
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


			for (unsigned int jjj = 0; jjj < 150; jjj++)
			{
				double randomX = floor(shuffle.size()*gsl_rng_uniform(rgen));
				string replaceX;
				vector<double> posX;
				double loci_numberX = atoi((shuffle[randomX][0]).c_str());
				for (unsigned int i = 0; i < loci_numberX; i++)
				{
					double shuf_sizeX = shuffle[randomX].size();
					replaceX.push_back(shuffle[randomX][shuf_sizeX - 6][i]);
					posX.push_back(atoi((shuffle[randomX][i + 1]).c_str()));
				}
				string sub_hapX = candidates[candidate_size - 1].substr(posX[0], posX.size());
				candidates[candidate_size - 1].replace(posX[0], posX.size(), replaceX);
			}
			double L_store1;
			double L_store2;
			while (true)
			{
				double j = floor(candidates.size()*gsl_rng_uniform(rgen));
				double random = floor(shuffle.size()*gsl_rng_uniform(rgen));
				string REPLACE;
				vector<double> pos;
				double loci_number = atoi((shuffle[random][0]).c_str());
				for (unsigned int i = 0; i < loci_number; i++)
				{
					double shuf_size = shuffle[random].size();
					REPLACE.push_back(shuffle[random][shuf_size - 6][i]);
					pos.push_back(atoi((shuffle[random][i + 1]).c_str()));
				}
				string sub_hap = candidates[j].substr(pos[0], pos.size());
				candidates[j].replace(pos[0], pos.size(), REPLACE);
				attempts++;

				vector<vector<vector<int>>> CONTRIBS;
				for (unsigned int i = 0; i < partitions.size(); i++)
				{
					vector<vector<int>> contribs;
					for (unsigned int j = 0; j < (partitions[i]).size(); j++)
					{
						vector<int> Contribs;
						vector<string> tempa;
						tempa = partitions[i][j];
						int sz = tempa.size();
						int num = sz - 6;
						for (unsigned int s = 0; s < candidates.size(); s++)
						{
							unsigned int counter = 0;
							for (int k = 1; k < num; k++)
							{
								int tem = atoi((tempa[k]).c_str());
								if ((candidates[s])[tem] == (tempa[num])[k - 1])
								{
									counter++;
								}
							}
							unsigned int dum = atoi((tempa[0]).c_str());
							if (counter == dum)
							{
								Contribs.push_back(s);
							}
						}
						contribs.push_back(Contribs);
					}
					CONTRIBS.push_back(contribs);
				}

				L_store1 = OptimumFreqs(0, candidates, init_freqs1, CONTRIBS, partitions);
				L_store2 = OptimumFreqs(1, candidates, init_freqs2, CONTRIBS, partitions);
				L_store = L_store1 + L_store2;

				if (L_store > L_new)
				{
					//cout << endl << L_store1 << " + " << L_store2 << " = " << L_store << "  " << L_new << endl;
					L_new = L_store;
					//cout << endl << "---->IMPROVE!" << endl;
					BAD_COUNT = 1;
					for (unsigned int i = 0; i < init_freqs1.size(); i++)
					{
						init_freqs1_store[i] = init_freqs1[i];
					}
					for (unsigned int i = 0; i < init_freqs2.size(); i++)
					{
						init_freqs2_store[i] = init_freqs2[i];
					}
				}
				else if (L_store <= L_new)
				{
					//cout << endl << L_store1 << " + " << L_store2 << " = " << L_store << "  " << L_new << endl;
					candidates[j].replace(pos[0], pos.size(), sub_hap);
					//cout << endl << "---->FAIL!" << endl;
					BAD_COUNT++;
				}
				if (BAD_COUNT > THRESHOLD || shuffle.size() == 0)
				{
					break;
				}
			}
			//Uncomment the following if you would like to see the output of every attempt the computer makes to find the optimum set

			/*cout << "sample " << sample << " is: " << L_new << endl;
			for (int i = 0; i < init_freqs1.size(); i++)
			{
			cout << init_freqs1[i] << "  ";
			}
			cout << endl;
			for (int i = 0; i < init_freqs1.size(); i++)
			{
			cout << init_freqs2[i] << "  ";
			}
			cout << endl << endl;
			for (int i = 0; i < candidates.size(); i++)
			{
			cout << candidates[i] << "   ";
			}
			cout << endl;*/
			if (loop > 1)
			{
				if (candidates[candidate_size-1] != candidates[candidate_size-2])
				{
					temp_likelihoods.push_back(L_new);
					freqs_1.push_back(init_freqs1_store);
					freqs_2.push_back(init_freqs2_store);
					temp_candidates.push_back(candidates);
				}
				else if (candidates[candidate_size-1] == candidates[candidate_size-2] && sacred_check < 5)
				{
					sample = 4;
					sacred_check ++;
				}
			}
			else if (loop == 1)
			{
				temp_likelihoods.push_back(L_new);
				freqs_1.push_back(init_freqs1_store);
				freqs_2.push_back(init_freqs2_store);
				temp_candidates.push_back(candidates);
			}
		}
		double max = temp_likelihoods[0];
		double max_it = 0;
		for (unsigned int mmm = 1; mmm < temp_likelihoods.size(); mmm++)
		{
			if (temp_likelihoods[mmm] > max)
			{
				max = temp_likelihoods[mmm];
				max_it = mmm;
			}
		}
		cout << endl << "number of haplotypes = " << hap_size;
		cout << endl;
		cout << temp_likelihoods[max_it];
		cout << endl;
		for (unsigned int i = 0; i < temp_candidates[max_it].size(); i++)
		{
			cout << temp_candidates[max_it][i] << '\t' << freqs_1[max_it][i] << '\t' << freqs_2[max_it][i] << endl;
		}
		all_likelihoods.push_back(temp_likelihoods[max_it]);
		L_preserve_2 = temp_likelihoods[max_it];
		BIC_check = -((-2*L_preserve_2)+hap_size*log (TOTAL_N)) + ((-2*L_preserve_1)+(hap_size-1)*log (TOTAL_N));
		if (BIC_check>0)
		{
			ofstream myfile;
			myfile.open("outcome_1.txt");
			for (unsigned int i = 0; i < temp_candidates[max_it].size(); i++)
			{
				myfile << i + 1 << "\t" << temp_candidates[max_it][i] << "\t" << freqs_1[max_it][i] << "\t" << freqs_2[max_it][i] << endl;
			}
			myfile.close();

		}
	}


	system("PAUSE");
	return 0;
}
