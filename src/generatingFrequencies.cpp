#include "generatingFrequencies.h"
#include <queue>
#include <cstdio>
#include <assert.h>

using namespace std;
using namespace boost::math;


int generate_inputs(int seed, int num_trees, int num_nodes, int num_samples, int num_mutations, double dirichlet_alpha, const char *input_file_location, const char *gamma_matrix_file_location, int deterministic_tree_index)
{
	
	int true_tree_index;	

	if(deterministic_tree_index != -1)
	{
		true_tree_index = deterministic_tree_index;
		if ( (true_tree_index < 0) || (true_tree_index >= num_trees) )
			fprintf(stderr, "Deterministic tree index out of range");
		assert( (true_tree_index >= 0) && (true_tree_index < num_trees));
	}
	else
	{
		boost::mt19937 gen_tree_index(seed);
		boost::random::uniform_int_distribution<> dist_tree_index(0, num_trees-1);
		true_tree_index = dist_tree_index(gen_tree_index);
	}

	ofstream input_file(input_file_location);
	if (input_file == NULL)
		fprintf(stderr, "Can not create file %s  \n", input_file_location);
	assert(input_file != NULL);

	input_file << "True_tree_index: " << true_tree_index << endl;
	input_file << "Num_nodes: " << num_nodes << endl;
	input_file << "Num_samples: " << num_samples << endl;
	input_file << "Seed: " << seed << endl << endl << endl;



	/*******************************************************************
	 Alpha frequencies are generated and then printed to file.
	********************************************************************/

	double alpha[num_nodes][num_samples];
	boost::mt19937 rng(seed);
	boost::gamma_distribution<> pdf(dirichlet_alpha);
	boost::variate_generator<boost::mt19937&, boost::gamma_distribution<> > generator(rng, pdf);
	
	for(int s = 0; s < num_samples; s++)
	{
		double sample_sum = 0.0;
		
		for(int i = 0; i < num_nodes; i++)
		{
			double sample = generator();
			alpha[i][s] = sample;
			sample_sum += sample;
		}
		
		for(int i = 0; i < num_nodes; i++)
		{
			alpha[i][s] = alpha[i][s] / sample_sum;
		}
	}
	
	input_file << "Alpha_freqencies:" << endl << endl << endl;
	for(int s=0; s<num_samples; s++)
	{
		for(int i=0; i<num_nodes; i++)
		{
			input_file << setw(10) << alpha[i][s] << " ";
		}
		
		input_file << endl;
	}
	
	input_file << endl << endl;

	
		
	/*************************************************************
	 Here we read gamma matrix.
	*************************************************************/

	int gamma[num_nodes][num_nodes];
	for (int i = 0; i < num_nodes; i++)
		for (int j = 0; j < num_nodes; j++)
			gamma[i][j]=0;

	ifstream gamma_matrix_file(gamma_matrix_file_location);
	if(gamma_matrix_file == NULL)
		fprintf(stderr, "Can not open file %s \n", gamma_matrix_file_location);
	assert(gamma_matrix_file != NULL);

	int read_int;
        gamma_matrix_file >> read_int >> read_int;


        for(int skip = 0; skip < true_tree_index; skip++)
                for(int i = 0; i < num_nodes; i++)
                        for(int j = 0; j < num_nodes; j++)
                                gamma_matrix_file >> read_int;


        for(int i=0; i<num_nodes;i++)
                for(int j=0; j<num_nodes; j++)
                        gamma_matrix_file >> gamma[i][j];

        gamma_matrix_file.close();
	

	/*****************************************************************
	 Here comes calculation and writing of q frequencies
	*****************************************************************/

	double q[num_nodes][num_samples];
	for(int s = 0; s < num_samples; s++)
	{
	
		for(int i = 0; i < num_nodes; i++)
		{
			q[i][s] = 0;

			for (int k = 0; k < num_nodes; k++)
			{
				q[i][s] += gamma[i][k] * alpha[k][s];
			}
		}

	}		

	input_file << "Q_freqencies:" << endl << endl << endl;
	for(int s=0; s<num_samples; s++)
	{
		for(int i=0; i<num_nodes; i++)
		{
			input_file << setw(10) << q[i][s] << " ";
		}
		
		input_file << endl;
	}
	
	input_file << endl << endl;
	

	/****************************************************************
	   Next we generate and  print locations of all mutations.
 	****************************************************************/
	 
	input_file << "Num_mutations: " << num_mutations << endl << endl << endl;
	
	typedef boost::mt19937 RNGType;
	RNGType rng_random(seed);
	boost::uniform_int<> mut_location(0 , num_nodes-1);
	boost::variate_generator<RNGType, boost::uniform_int<> > uni_integer(rng_random, mut_location);
	
	for(int i = 0; i < num_mutations; i++)
	{
		input_file << uni_integer() << " ";
	}
	input_file << endl;
	
	input_file.close();

	return true_tree_index;
}




void generate_frequencies(double error_rate , const char *input_file_location, const char *frequencies_file_location){
	
	int true_tree_index;
	int num_nodes;
	int num_samples;
	int seed;
	string skip_string;

	ifstream input_file(input_file_location);
	if(input_file == NULL)
		fprintf(stderr, "Can not open file %s \n", input_file_location);
	assert(input_file != NULL);

	ofstream frequencies_file(frequencies_file_location);
	if(frequencies_file == NULL)
		fprintf(stderr, "Can not open file %s \n ", frequencies_file_location); 	
	assert(frequencies_file != NULL);

	input_file >> skip_string >> true_tree_index >> skip_string >> num_nodes >> skip_string >> num_samples >> skip_string >> seed;
	
	double alpha[num_nodes][num_samples];
	input_file >> skip_string;
	for(int s = 0; s < num_samples; s++)
		for(int i = 0; i < num_nodes; i++)
			input_file >> alpha[i][s];
	
	double q[num_nodes][num_samples];
	input_file >> skip_string;
	for(int s = 0; s < num_samples; s++)
		for(int i=0; i < num_nodes; i++)
			input_file >> q[i][s];
	
	int num_mutations;
	input_file >> skip_string >> num_mutations;

	int locations_mutations[num_mutations];
	for(int i=0; i<num_mutations; i++)
		input_file >> locations_mutations[i];
	
	frequencies_file  << "Num_mutations: " << num_mutations << endl;
	frequencies_file  << "Num_samples: " << num_samples << endl;
	frequencies_file  << "Error_rate: " << error_rate << endl << endl << endl;
	
	input_file.close();

	typedef boost::mt19937 RNGType;
	RNGType rng(seed);
	boost::normal_distribution<> nd(0.0, error_rate);
	boost::variate_generator<RNGType, boost::normal_distribution<> > uni_gaussian(rng, nd);
	
	
	int counter = 0;
	double er;
	

	for (int i = 0; i < num_mutations; i++)
	{
		for(int s=0; s<num_samples; s++)
		{
			double frequency;
			
				if(error_rate != 0)
				{
					er = uni_gaussian();
					frequency = q[ locations_mutations[i] ][s] + er;
				}
				else
				{
					frequency = q[ locations_mutations[i] ][s];
				}

				if (frequency < 0)
				{
					frequency = 0;
				}
			
				if (frequency > 1)
				{
					frequency = 1;
				}
			
			frequencies_file << setw(10) << frequency << " ";
		}
		
		frequencies_file << endl;
	}
	
	input_file.close();
	frequencies_file.close();

}
