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

//Normalise the frequency of a vector to 1
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

