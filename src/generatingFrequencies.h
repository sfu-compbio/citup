#ifndef GENERATINGFREQUENCIES_H_INCLUDED
#define GENERATINGFREQUENCIES_H_INCLUDED

#include <iostream>
#include <fstream>


#include <boost/random.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/math/distributions/normal.hpp>

typedef boost::minstd_rand base_generator_type;
using namespace std;

int  generate_inputs(int seed, int num_trees, int num_nodes, int num_samples, int num_mutations, double dirichlet_alpha, const char *input_file_location, const char *gamma_matrix_file_location, int deterministic_tree_index);
void generate_frequencies(double error_rate, const char *input_file_location, const char *frequencies_file_location);

#endif // GENERATINGFREQUENCIES_H_INCLUDED
