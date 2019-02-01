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

//Normalise the frequency of haplotypes to 1
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

//Calculate the Dirichlet Multinomial given c=noise parameter, n=reads, x=frequencies of haplotypes
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

/*Given a timepoint (0 for the donor and 1 for the recipient), 
this function finds the likelihood of a set of haplotypes with 
their frequencies (freq) and their contribution to the reads of 
a given partial haplotype set (CONTRIBS) with full partitioning 
of the partial sets (partitions) and inferred noise parameter c*/
double Dirichlet_m(double timepoint, vector<string> haplotypes, vector<double> freq, vector<vector<vector<int>>> &CONTRIBS, vector<vector<vector<string>>> &partitions, double c)
{
	//In Multi_locus_trajectories.out, go two columns back from the last column on the right --> first timepoint 
	if (timepoint == 0)
	{
		timepoint = 3;
	}
	//In Multi_locus_trajectories.out, go zero columns back from the last column on the right --> last timepoint 
	else if (timepoint == 1)
	{
		timepoint = 1;
	}
	vector<string> candidates = haplotypes;
	vector<double> q = freq;
	double ULTIMA_THULE = 0; //For a candidate haplotype set, this parameter calculates the sum of the likelihoods for each of the i partial haplotype sets
	int candidates_size = candidates.size();
	for (unsigned int i = 0; i<CONTRIBS.size(); i++)
	{
		vector<double> inf; //contains the sum of the frequencies of all contributing candidate haplotypes for each of the j members of the ith partial haplotype set
		vector<double> nn; //contains the corresponding total number of reads for each member of inf vector
		vector<double> temp = q;
		temp.erase(remove(temp.begin(), temp.end() - 1, q[candidates_size]));//erase the unkknown haplotype qx out of the frequency list 
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
				//add the frequency of all the haplotypes that share the same partial haplotype of interest
				if (!(CONTRIBS[i][j]).empty())
				{
					sum += q[(CONTRIBS[i][j][k])];//add the frequency of all haplotypes that contribute to the ith partial haplotype set (CONTRIBS[i][j][k] = ith partial haplotype set, jth element, kth contributing haplotype)
					majorCheck++;
					majorCheck1++;
					majorCheckvec.push_back(q[(CONTRIBS[i][j][k])]);
				}
			}
			//if none of the candidate haplotypes match with the partial haplotype of interest (i.e. do not 'contribute' to that partial haplotype set), the corresponding reads for that partial haplotype goes to the unknown haplotype qx)
			if ((CONTRIBS[i][j]).empty())
			{
				check++;
				double SIZE1_1 = (partitions[i][j]).size();
				NUM1 = atoi((partitions[i][j]).at(SIZE1_1 - timepoint).c_str());
				sum_nn += NUM1;
				vaccum_holder.push_back(NUM1);
			}
			//if, at least, there were one candidate haplotype that contributed to the ith partial haplotype set, uptade nn and inf
			if (majorCheck != 0)
			{
				inf.push_back(sum);
				double SIZE = (partitions[i][j]).size();
				double NUM = atoi((partitions[i][j]).at(SIZE - timepoint).c_str());
				nn.push_back(NUM);
			}
		}
		//if the candidate haplotypes covered all the partial haplotypes in the ith set, then there are zero reads left to be attributed to qx
		if (check == 0 && majorCheck1 == temp.size())
		{
			inf.push_back(q[candidates_size]);
			nn.push_back(0.0);
		}
		unsigned int check2 = 0;
		
		//if some (but not all) of the candidate haplotypes covered the entire partial haplotype set i, the remaining candidates would correspond to zero reads (so as the qx haplotype)
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
			inf.push_back(q[candidates_size]);
			nn.push_back(0.0);
			check2++;
		}
		
		//if some (but not all) of the candidate haplotypes did not contribute to the ith set, then they correspond to zero reads and the remaining unidentified partial haplotypes go to qx
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
			inf.push_back(q[candidates_size]);
			nn.push_back(total);
			check3++;
		}
		//if non of the candidate haplotypes contributed to ith set, then they all correspond to zero reads and everything else would be classified as qx
		if (check != 0 && majorCheck1 != temp.size() && majorCheck1 == 0 && check2 == 0 && check3 == 0)
		{
			double temp_sum = 0;
			for (unsigned int g = 0; g < temp.size(); g++)
			{
				temp_sum += temp[g];
			}
			double total = accumulate(vaccum_holder.begin(), vaccum_holder.end(), 0);
			inf.push_back(q[candidates_size]);
			nn.push_back(total);
			inf.push_back(temp_sum);
			nn.push_back(0.0);
		}

		//if all the candidate haplotypes covered some of the partial haplotypes, then all the remaining reads correspond to qx.
		if (check != 0 && majorCheck1 == temp.size() && check2 == 0 && check3 == 0)
		{
			inf.push_back(q[candidates_size]);
			nn.push_back(sum_nn);
		}
		ULTIMA_THULE += Likelihood(c, nn, inf);
		temp = q;
	}
	return ULTIMA_THULE;
}

//For a candidate haplotype set, this function finds an optimal frequency given the short-read data collected from Multi_locus_trajectories.out file
double OptimumFreqs(double timepoint, vector<string> &candidates, vector<double> &init_freqs, vector<vector<vector<int>>> &CONTRIBS, vector<vector<vector<string>>> &partitions, double c)
{
	gsl_rng_env_setup();
	gsl_rng *rgen = gsl_rng_alloc(gsl_rng_taus);
	int seed = (int)time(NULL);
	gsl_rng_set(rgen, seed);
	double init_size = init_freqs.size();
	vector<double> init_freqs_store = init_freqs;
	double L = -1e7; //initial (extremely low) value for the likelihood
	unsigned int check = 1;
	unsigned int max_it = 8e3; //maximum number of attempts for optimising the frequencies
	double changex = 1e-2; //magnitude of the incremental (random) change in frequency
	double L_store1 = -1e7; //improved likelihood after frequency optimisation

	for (unsigned int it = 0; it < max_it; it++)
	{
		double r = 2 * (gsl_rng_uniform(rgen) - 0.5);//random direction between -1 to +1
		int j = floor(init_freqs.size()*gsl_rng_uniform(rgen));
		init_freqs[j] = init_freqs[j] + (r*changex);

		for (unsigned int i = 0; i < init_size; i++)
		{
			if (init_freqs[i] < 1e-11)
			{
				//changed this from 1e-11 to 0, so maybe change it back if something went wrong
				init_freqs[i] = 0;
			}
		}

		NormaliseFreqs(init_freqs); //frequencies should add up to one after being randomly changed
		if (init_freqs[init_size - 1] > 1e-2)
		{
			init_freqs[init_size - 1] = 1e-2; //the frequency of the 'X' haplotype, qx, could only go up to 1%
			double TOTAL = 0;
			for (unsigned int i = 0; i < init_size-1; i++)
			{
				TOTAL += init_freqs[i];
			}
			for (unsigned int i = 0; i < init_size-1; i++)
			{
				init_freqs[i] = init_freqs[i]*(0.99/TOTAL);//in this case, the frequency of the rest of the haplotypes should add up to 99%
			}
		}
		L = Dirichlet_m(timepoint, candidates, init_freqs, CONTRIBS, partitions, c);//calculate the Dirichlet likelihood after frequency optimisation step 

		if (L > L_store1)//if the likelihood improved, store it and go from there in the next step
		{
			L_store1 = L;
			init_freqs_store = init_freqs;
			check = 1;
		}
		if (L <= L_store1)//if there was no improvement for over 50 consecutive attempts, then we have reached a maximum likelihood peak 
		{
			init_freqs = init_freqs_store;
			check++;
		}
		if (check >= 50)//50 is found heuristically from testing multiple simulations
		{
			break;
		}
	}
	return L_store1;
}

int main(int argc, char* argv[])
{
	string multi_locus_file = argv[1]; //first input is the Multi_locus_trajectories.out file
	double c = atoi(argv[2]); //second input is the inferred noise parameter C
	
	if (FILE *file1 = fopen(multi_locus_file.c_str(), "r")) //check whether the Multi_locus_trajectories.out file exists
	{
		fclose(file1);
		vector<string> table;
		ifstream reads;
		reads.open(multi_locus_file);
		while (!reads.eof())//read the content of Multi_locus_trajectories.out and store it in table
		{
			string line;
			getline(reads, line);
			table.push_back(line);
		}
		reads.close();
		
		vector < vector<string> > rows; //save each row of Multi_locus_trajectories.out
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
		
		vector<string> loci; //collect the loci numbers from Multi_locus_trajectories.out
		vector<string> temp;
		for (unsigned int i = 0; i < rows.size(); i++)
		{
			temp = (rows[i]);
			unsigned int sz = temp.size();
			for (unsigned int j = 1; j < sz - 6; j++) //according to the current format of SAMFIRE, from the second column all the way to the seventh last column (from the right) contains loci numbers
			{
				loci.push_back(rows[i][j]);
			}
		}
		
		int TOTAL_N = 0; //calculating the total number of reads for all the partial haplotypes -- this quantity is required for BIC calculations 
		for (unsigned int i=0; i<rows.size(); i++)
		{
			unsigned int sz = rows[i].size();
			string temp_1 = rows[i][sz-1];
			string temp_2 = rows[i][sz-3];
			int num_1 = atoi(temp_1.c_str());
			int num_2 = atoi(temp_2.c_str());
			TOTAL_N += num_1 + num_2;
		}
		
		vector<int> positions;//The following ~30 lines of code sorts the order of loci from 0 to the last variant locus and store them in Fpositions 
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

		
		vector<vector<vector<string>>> partitions; //partitions the partial haplotypes from Multi_locus_trajectories.out into groups with matching loci numbers
		for (unsigned int k = 1; k <= Fpositions.size(); k++)
		{
			for (unsigned int i = 0; i < Fpositions.size(); i++)
			{
				vector<string> pick;
				vector<vector<string>> partialPar;
				for (unsigned int j = 0; j < locus.size(); j++)
				{
					pick = locus[j];
					string pos_num = to_string(k);
					unsigned int check_match = 0;
					if (pick[0] == pos_num)
					{
						for (unsigned int m = 0; m < k; m++)
						{
							if (pick[m + 1] == to_string(i + m))
							{
								check_match++;
							}
						}
						if (check_match == k)
						{
							partialPar.push_back(pick);
						}
					}
				}
				if (!partialPar.empty())
				{
					partitions.push_back(partialPar);
				}
			}
		}
		
		vector<string> tampa;//Here we sort each element of a partition set from largest to smallest with respect to the total number of reads
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

		//Uncomment the following lines if you want to see how the jth member of the ith partial haplotype set looks like
		for (unsigned int i = 0; i<partitions.size(); i++)
		{
			for (unsigned int j = 0; j<(partitions[i]).size(); j++)
			{
				for (unsigned int k = 0; k<(partitions[i][j]).size(); k++)
				{
					//cout << partitions[i][j][k] << "  ";
				}
				//cout << endl << "--------------------------------------" << endl;
			}
			//cout << "**************************************" << endl;
		}
		
		vector<vector<string>> shuffle_partitions; //this stores all the information of every partial haplotype in one list which is later going to be used as a shuffling set for optimising a set of candidate haplotypes
		for (unsigned int i = 0; i < partitions.size(); i++)
		{
			for (unsigned int j = 0; j < (partitions[i]).size(); j++)
			{
				vector<string> shuffle_partitions_holder;
				for (unsigned int k = 0; k < (partitions[i][j]).size(); k++)
				{
					shuffle_partitions_holder.push_back(partitions[i][j][k]);
				}
				shuffle_partitions.push_back(shuffle_partitions_holder);
			}
		}
		
		//our 'default' candidate haplotype has no information in it, i.e. it starts with unknown nucleotide X at each position
		string unknown;
		for (unsigned int i = 0; i<positions.size(); i++)
		{
			unknown.append("X");
		}
			
		int row_numbers = rows.size();//number of rows in Multi_locus_trajectories.out
		vector<string> candidates;//a list of candidate haplotypes; starting from 1 haplotype and continue optimising/appending new ones until BIC becomes negative
		double L_new = -1e7;
		double L_store = -1e7;
		double L_preserve_1 = -1e7 - 300;
		double L_preserve_2 = -1e7 - 1;
		double loop = 0;
		vector<double> all_likelihoods;
		vector<double> haplotypes;
		vector<double> freqs_before;
		vector<double> freqs_after;
		double BIC_check=2;
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
			for (unsigned int sample = 0; sample < 20; sample++)
			{
				//?could make the condition below better?
				L_new = L_preserve_2 - 100 / (2*candidate_size1);
				L_store = -1e5;
				candidates[candidate_size1 - 1] = unknown;
				double candidate_size = candidates.size();
				double BAD_COUNT = 0;
				//This defines the cut-off on when to stop making further attempts to optimize the haplotypes
				//This depends on loop=number of haplotypes at any time and read_count=total number of partial-reads available.
				double THRESHOLD = 2 * loop * row_numbers;
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
					double randomX = floor(shuffle_partitions.size()*gsl_rng_uniform(rgen));
					string replaceX;
					vector<double> posX;
					double loci_numberX = atoi((shuffle_partitions[randomX][0]).c_str());
					for (unsigned int i = 0; i < loci_numberX; i++)
					{
						double shuf_sizeX = shuffle_partitions[randomX].size();
						replaceX.push_back(shuffle_partitions[randomX][shuf_sizeX - 6][i]);
						posX.push_back(atoi((shuffle_partitions[randomX][i + 1]).c_str()));
					}
					string sub_hapX = candidates[candidate_size - 1].substr(posX[0], posX.size());
					candidates[candidate_size - 1].replace(posX[0], posX.size(), replaceX);
				}
				double L_store1;
				double L_store2;
				while (true)
				{
					double j = floor(candidates.size()*gsl_rng_uniform(rgen));
					double random = floor(shuffle_partitions.size()*gsl_rng_uniform(rgen));
					string REPLACE;
					vector<double> pos;
					double loci_number = atoi((shuffle_partitions[random][0]).c_str());
					for (unsigned int i = 0; i < loci_number; i++)
					{
						double shuf_size = shuffle_partitions[random].size();
						REPLACE.push_back(shuffle_partitions[random][shuf_size - 6][i]);
						pos.push_back(atoi((shuffle_partitions[random][i + 1]).c_str()));
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

					L_store1 = OptimumFreqs(0, candidates, init_freqs1, CONTRIBS, partitions, c);
					L_store2 = OptimumFreqs(1, candidates, init_freqs2, CONTRIBS, partitions, c);
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
					if (BAD_COUNT > THRESHOLD || shuffle_partitions.size() == 0)
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
					int candidates_check = 0;
					for (unsigned int kappa = 0; kappa < candidate_size - 1; kappa++)
					{	
						if (candidates[candidate_size-1] != candidates[kappa])
						{
							candidates_check++;
						}
					}
					
					int really_bad = 0;
					for (unsigned int kappa = 0; kappa < candidate_size; kappa++)
					{	
						if (init_freqs1_store[kappa] < 1e-9 && init_freqs2_store[kappa] < 1e-9)
						{
							really_bad++;
						}
					}
					
					if (candidates_check == candidate_size - 1 && really_bad == 0)
					{
						temp_likelihoods.push_back(L_new);
						freqs_1.push_back(init_freqs1_store);
						freqs_2.push_back(init_freqs2_store);
						temp_candidates.push_back(candidates);
					}
					else if (candidates_check != candidate_size - 1 && really_bad != 0)
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
		
		return 0;
	}
	
	else
	{
		cout << "No Multi_locus_trajectories.out file found" << endl;
		return false;
	}
}