/*
	STOCHASTIC LOTKA-VOLTERRA MODEL WITH MULTIDIMENSIONAL NICHE SPACES
	by Tommaso Biancalani <tommasob@mit.edu>

	Supporting Material to the paper:
	"Framework for analyzing ecological trait-based models in multidimensional niche spaces" - Phys. Rev. E 91, 052107

	ABSTRACT
	We develop a theoretical framework for analyzing ecological models with a multidimensional niche space. Our approach relies on the fact that ecological niches are described by sequences of symbols, which allows us to include multiple phenotypic traits. Ecological drivers, such as competitive exclusion, are modeled by introducing the Hamming distance between two sequences. We show that a suitable transform diagonalizes the community interaction matrix of these models, making it possible to predict the conditions for niche differentiation and, close to the instability onset, the asymptotically long time population distributions of niches. We exemplify our method using the Lotka-Volterra equations with an exponential competition kernel.

	INSTALL
 	The code requires the GSL (GNU Scientific Library). On linux, I simply compile it with:
 	$ gcc lotka-volterra.c -lgsl -lgslcblas -lm -Wall -o "lotka-volterra"
 */

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/time.h>
#include <time.h>
#include <sys/select.h>
#include <math.h>

// Headers for special functions
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_exp.h>

// Headers for Mersenne-Twister random generator
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// Global variables
#define R_CONST		1.			// Competition kernel length scale
#define	SIGMA		2			// Competition kernel exponent
#define	V_CONST		100.		// System size of single compartment
#define	T_SAMP		1. 			// Sampling time (outputs every T_SAMP)
#define	TOT_TIME	100. 		// Evolve the system until t/V = TOT_TIME

// Sequence space definitions (sequences are called 'genomes')
typedef char GENOME_S;
typedef unsigned int GENOME_N;
#define	C_CONST	5.52486			// Competition kernel normalization constant
#define L 	3					// Phenotype length
#define NUM_GENOMES	54			// Total number of genome in the systems
const int DELTA[L] = {3,9,2};  	// It specifies the sequence space
GENOME_S *genome_0 = "000";
GENOME_S genome_db[NUM_GENOMES][L+1];
#define A(k1,k2) (k1)+2*(k2)  	// reaction k1 in genome k2
FILE *out;
unsigned int *org;				// org[I]: number of organisms with genome I

// GSL rn definitions (Mersenne-Twister)
const gsl_rng_type * rgT;
gsl_rng * rg;
double rn;
unsigned long int seed;
/**********************************************************/


GENOME_S *genome_num2str(GENOME_N I)
{
/*
	Helper for converting genomes from numbers (GENOME_N) to string (GENOME_S).
	Return converted genome in string format.
*/
	int l;
	GENOME_S *tbl = "0123456789abcdefghijklmnopqrstuvwxyz";
	GENOME_S *out;
	out = (GENOME_S *) malloc(sizeof(GENOME_S)*L);
	memcpy(out, genome_0, sizeof(GENOME_S)*L);
	int len = L-1;
	for(l=0; 1; l++)
	{
		 out[len--] = tbl[I % DELTA[l]];
		 if((I /= DELTA[l]) == 0)
		 	break;
	}
	return out;
}


unsigned int Hamming_dist(GENOME_S *genome_I, GENOME_S *genome_J)
{
/*
	Return Hamming distance between genomes I and J
*/
	int l;
	unsigned int ham_dist=L;
	for(l=0; l<L; l++)
		if(genome_I[l] - genome_J[l] == 0)
			ham_dist--;
	return ham_dist;
}


double G_IJ(GENOME_S *genome_I, GENOME_S *genome_J) {
/*
	Return competition rate G on genome I due to genome J.
*/
	double arg, ham_dist;
	ham_dist = ((double)Hamming_dist(genome_I, genome_J)) / R_CONST;
	arg = gsl_pow_int(ham_dist, SIGMA);
	return gsl_sf_exp(-arg);
}


int gillespie_evolve()
{
/*
	Print model dynamic on standard output. Dynamic is computed using Gillespie's algorithm. Return zero when time reaches TOT_TIME.
*/
	GENOME_N I, J;		// Genome indexes
	int k;				// Reaction index
	long int step; 		// Time iteration index
	int sample_step;	// Sampled time iteration index
	double tau, t; 		// Stochastic times
	double *T_rates; 	// Transition rates
	double T_sum, Tr2, T_tmp;
	double rn; 			// rand number

	// Allocate memory for transition rates
	T_rates = (double *)malloc(sizeof(double) * NUM_GENOMES * 2);

	// Main Loop
	t=0.; step=0; sample_step=0;
	while(t < TOT_TIME)
	{
		// Output as time reaches next sampled time
		if(t > (double) (T_SAMP * sample_step))
		{
			for(I=0; I<NUM_GENOMES; I++)
				fprintf(out, "%.3f ", org[I]/V_CONST);
			fprintf(out, "\n");
			sample_step++;
		}

		// Compute transiton rates
		for (I=0; I<NUM_GENOMES; I++)
		{
			T_rates[A(0, I)] = org[I]; 	// reaction 0 = birth
			T_rates[A(1, I)] = 0.; 		// reaction 1 = death
			for(J=0; J<NUM_GENOMES; J++)
				T_rates[A(1,I)] += G_IJ(genome_db[I], genome_db[J]) * org[J];
			T_rates[A(1,I)] *= org[I] / (V_CONST * C_CONST);
		}

		// Compute total transition rate
		T_sum = 0.;
		for(I=0; I<NUM_GENOMES; I++)
			T_sum += T_rates[A(0, I)] + T_rates[A(1, I)];

		// Generate next reaction time (tau)
		rn = gsl_rng_uniform_pos(rg);
		tau = (1./T_sum) * log(1./rn);

		// Update time
		t += tau;

		// Generate what reaction occurs and on which species
		Tr2 = gsl_rng_uniform_pos(rg) * T_sum;
		T_tmp = 0.;
		for (I=0; I<NUM_GENOMES; I++)
		{
			for(k=0; k<2; k++)
			{
				T_tmp += T_rates[A(k,I)];
				if(T_tmp >= Tr2)  goto gotten; // this is a fine use of goto :-)
			}
		}
		gotten:

		// Execute reaction i on species k
		switch(k) // Wich reaction?
		{
			// BIRTH
			case 0:
				org[I]++;
			break;

			// DEATH
			case 1:
				org[I]--;
			break;
		}
		// update step
		step = step+1;
	} // main loop
	return 0;
} // gillespe_evolve()


int main(int argc, char * argv[])
{
/*
	Main function. Arguments are not processed. The model parameters are set by the constants and globals defined at the beginning of the code.
*/
	GENOME_N I; // Genome index

	// Initialize random generator
	srand(time(NULL));
	seed = (unsigned long int)rand();
	rgT = gsl_rng_default; // default = Mersenne Twister
	rg = gsl_rng_alloc (rgT);
	gsl_rng_set (rg, seed);

	// Initialize organisms, e.g. org[23] = number of organisms with genome 23
	org = (unsigned int *)malloc(sizeof(unsigned int) * NUM_GENOMES);

	// Initialize system around homogeneous state
	for (I=0; I<NUM_GENOMES; I++)
		org[I] = V_CONST;

	// Set file descriptor for output
	out = stdout;

	// Initialize genome space
	for(I=0; I<NUM_GENOMES; I++)
		strcpy(genome_db[I], genome_num2str(I));

	// Compute model dynamic
	gillespie_evolve();

	// Free descriptors and quit
	fclose(out);
	gsl_rng_free(rg);
	return 0;
} // main()
