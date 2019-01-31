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
frequency (freq) and their contribution to the total reads of 
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
		ULTIMA_THULE += Likelihood(c, nn, inf);
		temp = q;
	}
	return ULTIMA_THULE;
}
