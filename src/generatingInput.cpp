/*
 *
 * File:   main.cpp
 * Author: smalikic
 *
 * Created on June 9, 2013, 8:40 PM
 *
 */

#include <cstdlib>
#include <sstream>
#include <iostream>
#include <queue>
#include <cstring>
#include <string>
#include <algorithm>
#include <assert.h>
#include <ctime>
#include "generatingFrequencies.h"
#include "tclap/CmdLine.h"


using namespace std;


int main(int argc, char** argv) {
	
	
	TCLAP::CmdLine cmd("Command...", ' ', "1");
	
	TCLAP::ValueArg<int> seedArg("z","seed","Seed for random number generator",false,time(0),"integer",cmd);
	TCLAP::ValueArg<int> numNodes("n","numberNodes","Defines the number of nodes of the tree",false,7,"integer",cmd);
	TCLAP::ValueArg<int> deterministicTreeIndexArg("p","deterministicTreeIndex","Use if you want deterministic tree index",false,-1,"integer",cmd);
	TCLAP::ValueArg<int> numSamples("s","numberSamples","Defines the number of samples we have",false,1,"integer",cmd);
	TCLAP::ValueArg<int> numMutations("m","numberMutations","Defines the number of mutations that each sample contains",false,100,"integer",cmd);
	TCLAP::ValueArg<double> errorRate("e","errorRate","Defines the maximum error frequency, should be between 0 and 1",false,0.10,"real",cmd);
	TCLAP::ValueArg<double> dirichletAlphaArg("d","dirichletAlpha","Dirichlet parameter for leaf frequencies",false,5.0,"real",cmd);
	TCLAP::ValueArg<string> frequencies_file_arg("f","frequencies_file","Defines the name of frequencies file.",false,"./Frequencies.txt","string",cmd);
	TCLAP::ValueArg<string> input_file_arg("i","input_file","Defines the name of  input file.",false,"./Input.txt","string",cmd);
	TCLAP::ValueArg<string> gamma_adj_folder_arg("g","gamma_adjacency_folder_location","Defines location of folder that contains gamma and adj matrices files",false,"./GammaAdjMatrices/GammaAdjMatrices_Full","string",cmd);
	
	cmd.parse(argc,argv);
	
	
	int     seed                       =    seedArg.getValue();
	int     num_nodes                  =    numNodes.getValue();
	int 	deterministic_tree_index   =	deterministicTreeIndexArg.getValue();
	int     num_samples                =    numSamples.getValue();
	int     num_mutations              =    numMutations.getValue();
	double  error_rate                 =    errorRate.getValue();
	double  dirichletAlpha             =    dirichletAlphaArg.getValue();
	string  input_file_location        =    input_file_arg.getValue();
	string  frequencies_file_location  =    frequencies_file_arg.getValue();
	string  gamma_adj_folder_location  =	gamma_adj_folder_arg.getValue();

	
	ostringstream gamma_file_string_location;
	gamma_file_string_location << gamma_adj_folder_location << "/GammaMatrix" << num_nodes << ".txt";
	
	string gamma_matrix_file_location = gamma_file_string_location.str();	
 
	ifstream gamma_matrix_file(gamma_matrix_file_location.c_str());
	if(gamma_matrix_file == NULL)
	{
		cerr << "Error: Unable to open file " << gamma_matrix_file_location << endl;
		exit(1);
	}

	int skip_num_nodes, num_trees;
	gamma_matrix_file >> skip_num_nodes >> num_trees;
	gamma_matrix_file.close();

	generate_inputs(seed, num_trees , num_nodes, num_samples, num_mutations, dirichletAlpha, input_file_location.c_str(), gamma_matrix_file_location.c_str(), deterministic_tree_index);
	
	generate_frequencies(error_rate, input_file_location.c_str(), frequencies_file_location.c_str());
	
	
    return 0;
	
}

