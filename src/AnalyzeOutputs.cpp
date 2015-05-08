#include <iostream>
#include <vector>
#include <fstream>

#include <math.h>
#include <assert.h>
#include <sstream>
#include <cstdio>
#include <cstring>
#include <tclap/CmdLine.h>

using namespace std;


#define INF 100000000    //just infinity
#define MAX_MUT 5000
#define MAX_SAMPLES 10
#define MAX_NODES 1000
#define N MAX_NODES      //max number of vertices in one part

int cost[N][N];          //cost matrix
int n, max_match;        //n workers and n jobs
int lx[N], ly[N];        //labels of X and Y parts
int xy[N];               //xy[x] - vertex that is matched with x,
int yx[N];               //yx[y] - vertex that is matched with y
bool S[N], T[N];         //sets S and T in algorithm
int slack[N];            //as in the algorithm description
int slackx[N];           //slackx[y] such a vertex, that
// l(slackx[y]) + l(y) - w(slackx[y],y) = slack[y]
int previous[N];         //array for memorizing alternating paths

int 	input_mut_locations[MAX_MUT];
int 	output_mut_locations[MAX_MUT];
int     input_gamma_matrix[MAX_NODES][MAX_NODES];
int     output_gamma_matrix[MAX_NODES][MAX_NODES];
int 	input_adjacency_matrix[MAX_NODES][MAX_NODES];
int	output_adjacency_matrix[MAX_NODES][MAX_NODES];
int	input_num_nodes;
int 	output_num_nodes;
int     input_tree_index;
int	output_tree_index;
int 	num_mutations;
int 	read_int, skip_int;
int	tot_same_node_input;
int	tot_same_node_preserved;
int	tot_direct_neigbours_input;
int	tot_direct_neighbors_preserved;
int	tot_non_direct_ancestors_input;
int	tot_non_direct_ancestors_preserved;
int	tot_non_ancestors_input;
int	tot_non_ancestors_preserved;
double 	read_double;
double	skip_double;
string  skip_string;
double  input_alpha[MAX_NODES][MAX_SAMPLES];
double  output_alpha[MAX_NODES][MAX_SAMPLES];
double	input_q[MAX_NODES][MAX_SAMPLES];
double	output_q[MAX_NODES][MAX_SAMPLES];
int	current;
int  	num_samples;
string	gamma_adj_folder_location;
string  input_file_location;
string  output_file_location;


void 	ReadInput(const char *input_file_location, const char *output_file_location);
void    trap_ReadInput(const char *input_file_location, const char *output_file_location);
void    ReadGammaMatrix(int gamma[MAX_NODES][MAX_NODES], int num_leaves, int tree_index);
void    ReadAdjacencyMatrix(int adjacency[MAX_NODES][MAX_NODES], int num_nodes, int tree_index);
int	FirstMeasure();
double  SecondMeasure();
int 	ThirdMeasure();
void	FourthMeasure();
void    ShowTableComparison();
void    InitializeToZeroes();
int 	max(int a, int b);


void	init_labels();
void 	augment();
void 	update_labels();
void 	add_to_tree(int x, int prevx);
void	update_labels();
void 	add_to_tree(int x, int prevx);
int  	hungarian();



int main(int argc, char** argv)
{	

	TCLAP::CmdLine cmd("Command...", ' ', "1");
	
	TCLAP::ValueArg<string> gamma_adj_folder_arg("g","gamma_adjacency_folder_location","Defines location of folder that contains files with gamma and adjacency matrices",false,"./GammaAdjMatrices/GammaAdjMatrices_Full","string",cmd);
	TCLAP::ValueArg<string> input_file_arg("i","input_file","Defines exact name of input file",false,"./Input.txt","string",cmd);
	TCLAP::ValueArg<string> output_file_arg("o","output_file","Defines exact name of output file",false,"./Results.txt","string",cmd);
	TCLAP::ValueArg<string> stats_filename_arg("s","stats_file","Output statistics filename",false,"./Stats.txt","string",cmd);
	TCLAP::SwitchArg trap_results_arg("t","trap_results", "Set this if your output is obtained via trap", cmd);
	

	cmd.parse(argc,argv);
	
	
	input_file_location 		  =   input_file_arg.getValue();
	output_file_location	 	  =   output_file_arg.getValue();
	gamma_adj_folder_location         =   gamma_adj_folder_arg.getValue();
	string stats_filename   	  =   stats_filename_arg.getValue();	
	bool trap_results		  =   trap_results_arg.getValue();

	InitializeToZeroes();

	if(trap_results)
	{
		trap_ReadInput(input_file_location.c_str(), output_file_location.c_str());
		
		ReadGammaMatrix(input_gamma_matrix, input_num_nodes, input_tree_index);
		ReadAdjacencyMatrix(input_adjacency_matrix, input_num_nodes, input_tree_index);
	}
	else
	{
		ReadInput(input_file_location.c_str(), output_file_location.c_str());

		ReadGammaMatrix(input_gamma_matrix, input_num_nodes, input_tree_index);
		ReadAdjacencyMatrix(input_adjacency_matrix, input_num_nodes, input_tree_index);
		ReadGammaMatrix(output_gamma_matrix, output_num_nodes, output_tree_index);	
		ReadAdjacencyMatrix(output_adjacency_matrix, output_num_nodes, output_tree_index);
	}

	hungarian();
	
	
	ShowTableComparison();


	int num_misplaced_mutations 	  =   FirstMeasure();
	double average_freq_error 	  =   SecondMeasure();
	int num_nonpreserved_edges        =   ThirdMeasure();
	FourthMeasure();	
	cout << "1. SCORE (NUMBER OF MISPLACED MUTATIONS) IS : " << num_misplaced_mutations << endl;
	cout << "2. SCORE (AVERAGE FREQ. ERROR) IS : " << average_freq_error << endl;
	cout << "3. SCORE (NUMBER OF NON PRESERVED ADJACENCIES) IS : " << num_nonpreserved_edges << endl;
	cout << "MUTATIONS AT SAME NODE (INPUT vs PRESERVED): " << tot_same_node_input << "  " << tot_same_node_preserved << endl;
	cout << "MUTATIONS THAT ARE DIRECT NEIGHHBORS(INPUT vs PRESERVED " << tot_direct_neigbours_input << "  " << tot_direct_neighbors_preserved << endl;
	cout << "MUTATIONS WITH ANCESTRY RELATIONSHIP, NOT DIRECT, (INPUT vs PRESERVED) " << tot_non_direct_ancestors_input << "  " << tot_non_direct_ancestors_preserved << endl;
	cout << "MUTATIONS WITHOUT ANCESTRY RELATIONSHIP (INPUT vs OUTPUT) " << tot_non_ancestors_input << "  " << tot_non_ancestors_preserved << endl;
	
	int correct_tree = 0;
	if (input_num_nodes == output_num_nodes && num_nonpreserved_edges == 0)
	{
		correct_tree = 1;
	}

	ofstream stats_file(stats_filename.c_str());
	stats_file << "num_misplaced_mutations\t" << num_misplaced_mutations << endl;
	stats_file << "average_freq_error\t" << average_freq_error << endl;
	stats_file << "num_nonpreserved_edges\t" << num_nonpreserved_edges << endl;
	stats_file << "correct_tree\t" << correct_tree << endl;
	stats_file << "num_nodes_predicted\t" << output_num_nodes << endl;
	stats_file << "same_node_simulated\t" << tot_same_node_input << endl;
	stats_file << "same_node_predicted\t" << tot_same_node_preserved << endl;
	stats_file << "neighbours_simulated\t" << tot_direct_neigbours_input << endl;
	stats_file << "neighbours_predicted\t" << tot_direct_neighbors_preserved << endl;
	stats_file << "ancestors_simulated\t" << tot_non_direct_ancestors_input << endl;
	stats_file << "ancestors_predicted\t" << tot_non_direct_ancestors_preserved << endl;
	stats_file << "non_ancestors_simulated\t" << tot_non_ancestors_input << endl;
	stats_file << "non_ancestors_predicted\t" << tot_non_ancestors_preserved << endl;
	stats_file.close();

}



/********************************************************************************************************************************************************/


void ReadInput(const char *input_file_location_arg, const char *output_file_location_arg)
{
	
	ifstream input_file(input_file_location_arg, ifstream::in);
	if(input_file == NULL)
		fprintf(stderr, "Can not open input file %s \n", input_file_location_arg);
	assert(input_file != NULL);

	ifstream output_file(output_file_location_arg);
	if(output_file == NULL)
		fprintf(stderr, "Can not open output file %s \n", output_file_location_arg);
	assert(output_file != NULL);



	input_file >> skip_string >> input_tree_index;
	input_file >> skip_string >> input_num_nodes;
	if(input_num_nodes >= MAX_NODES)
		fprintf(stderr, "Number of nodes in %s is equal to %d that is >= MAX_NODES. Please change MAX_NODES defined in AnalyzeOutputs.cpp file. \n", input_file_location_arg, input_num_nodes);
	assert(input_num_nodes < MAX_NODES);

	input_file >> skip_string >> num_samples;
	if(num_samples >= MAX_SAMPLES)
		fprintf(stderr, "Number of samples in %s is equal to %d that is >= MAX_SAMPLES. Please change MAX_SAMPLES defined in AnalyzeOutputs.cpp file. \n", input_file_location_arg, num_samples);
	assert(num_samples < MAX_SAMPLES);


	int seed;
	input_file >> skip_string >> seed;
	
	double output_objective_value;
	output_file >> skip_string >> output_objective_value;

	int output_num_clusters; 
	output_file >> skip_string >> output_num_clusters;

	output_file >> skip_string >> output_tree_index;

	output_file >> skip_string >> output_num_nodes;
	if(output_num_nodes >= MAX_NODES)
		fprintf(stderr, "Number of nodes in %s is equal to %d that is >= MAX_NODES. Please change MAX_NODES defined in AnalyzeOutputs.cpp file. \n", output_file_location_arg, output_num_nodes);
	assert(output_num_nodes < MAX_NODES);


	n = max (input_num_nodes, output_num_nodes);

	int output_num_samples;
	output_file >> skip_string >> output_num_samples;
	if(output_num_samples >= MAX_SAMPLES)
		fprintf(stderr, "Number of samples in %s is equal to %d that is >= MAX_SAMPLES. Please change MAX_SAMPLES defined in AnalyzeOutputs.cpp file. \n", output_file_location_arg, output_num_samples);
	assert(output_num_samples < MAX_SAMPLES);


	int output_num_mutations; 
	output_file >> skip_string >> output_num_mutations;
	if(output_num_mutations > MAX_MUT)
		fprintf(stderr, "Number of mutations in %s is equal to %d that is >= MAX_MUT. Please change MAX_MUT defined in AnalyzeOutputs.cpp file. \n", output_file_location_arg, output_num_mutations );
	assert(output_num_mutations < MAX_MUT);


	string output_status_value; 
	output_file >> skip_string >> output_status_value;

	double output_cplex_time; 
	output_file >> skip_string >> output_cplex_time;

	int output_cplex_threads;
	output_file >> skip_string >> output_cplex_threads;

	int output_cplex_memory;
	output_file >> skip_string >> output_cplex_memory;
	
	//quadratic error reading
	output_file >> skip_string >> skip_string;
	

	// Input alpha frequencies:
	input_file >> skip_string;
	for(int s = 0; s < num_samples; s++)
		for(int i = 0; i < input_num_nodes; i++)
			input_file >> input_alpha[i][s];


	// Input q frequencies
	input_file >> skip_string;
	for(int s = 0; s < num_samples; s++)
		for(int i = 0; i < input_num_nodes; i++)
			input_file >> input_q[i][s];

	// Output alpha frequencies
	output_file >> skip_string;
	for(int s = 0; s < num_samples; s++)
		for(int i = 0; i < output_num_nodes; i++)
			output_file >> output_alpha[i][s];
		
	// Output q frequencies
	output_file >> skip_string;
	for(int s = 0; s < num_samples; s++)
		for(int i = 0; i < output_num_nodes; i++)
			output_file >> output_q[i][s];
	
		
	input_file >> skip_string >> num_mutations;
	if(output_num_mutations != num_mutations)
		fprintf(stderr, "ERROR! Number of mutations in %s is different from number of mutations in %s.! that is >= MAX_MUT.\n", output_file_location_arg, input_file_location_arg);
	assert(output_num_mutations == num_mutations);


	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			cost[i][j] = 0;
	
	
	vector<vector<int> > input_values;
	vector<vector<int> > result_values;

	input_values.resize(n);
	result_values.resize(n);

	for(int i = 0; i < num_mutations; i++)
	{
		input_file >> read_int;
		input_values[read_int].push_back(i);
		input_mut_locations[i] = read_int;
	}
	
	for(int i = 0; i < num_mutations; i++)
	{
		output_file >> read_int;
		result_values[read_int].push_back(i);
		output_mut_locations[i] = read_int;
	}
	
	
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			
			for(int k = 0; k < input_values[i].size(); k++)
			{
				current = input_values[i].at(k);
				
				for(int t = 0; t < result_values[j].size(); t++)
				{
					if(current == result_values[j].at(t))
					{
						cost[i][j]++;
						break;
					}
				}
			}
			
			cost[i][j] =-(input_values[i].size() + result_values[j].size() - 2*cost[i][j]);
			
		}
	}
	
	input_file.close();
	output_file.close();
	
}


/************************************************************************************************************************************/


void trap_ReadInput(const char *input_file_location_arg, const char *output_file_location_arg)
{
	
	ifstream input_file(input_file_location_arg, ifstream::in);
	if(input_file == NULL)
		fprintf(stderr, "Can not open input file %s \n", input_file_location_arg);
	assert(input_file != NULL);

	ifstream output_file(output_file_location_arg);
	if(output_file == NULL)
		fprintf(stderr, "Can not open output file %s \n", output_file_location_arg);
	assert(output_file != NULL);



	input_file >> skip_string >> input_tree_index;

	input_file >> skip_string >> input_num_nodes;
	if(input_num_nodes >= MAX_NODES)
		fprintf(stderr, "Number of nodes in %s is equal to %d that is >= MAX_NODES. Please change MAX_NODES defined in AnalyzeOutputs.cpp file. \n", input_file_location_arg, input_num_nodes);
	assert(input_num_nodes < MAX_NODES);

	input_file >> skip_string >> num_samples;
	if(num_samples >= MAX_SAMPLES)
		fprintf(stderr, "Number of samples in %s is equal to %d that is >= MAX_SAMPLES. Please change MAX_SAMPLES defined in AnalyzeOutputs.cpp file. \n", input_file_location_arg, num_samples);
	assert(num_samples < MAX_SAMPLES);

	int seed;
	input_file >> skip_string >> seed;
	

	/***************************************************************
		Reading gamma and adjacency matrices and num_nodes
		for trap outputs files
	***************************************************************/
	// This can be worked on in future (isomoprhism can be checked,
	// even if labellings are different

	output_tree_index = -1;

	double output_objective_value;
	output_file >> skip_string >> output_objective_value;
	
	int output_num_edges;
	output_file >> skip_string >> output_num_edges;
	assert(output_num_nodes == output_num_edges + 1);

	for(int i=0; i<output_num_edges; i++)
	{	
		int vertex_1, vertex_2;
		output_file >> vertex_1 >> vertex_2;
		output_adjacency_matrix[vertex_1][vertex_2]=1;
		output_adjacency_matrix[vertex_2][vertex_1]=1;
	}
	

	output_file >> skip_string >> output_num_nodes;
	if(output_num_nodes >= MAX_NODES)
		fprintf(stderr, "Number of nodes in %s is equal to %d that is >= MAX_NODES. Please change MAX_NODES defined in AnalyzeOutputs.cpp file. \n", output_file_location_arg, output_num_nodes);
	assert(output_num_nodes < MAX_NODES);

	n = max (input_num_nodes, output_num_nodes);

	for (int i=0; i<output_num_nodes; i++)
		for(int j=0; j<output_num_nodes; j++)
			output_file >> output_gamma_matrix[i][j];
	
	/***************************************************************/
	
	// Input alpha frequencies:
	input_file >> skip_string;
	for(int s = 0; s < num_samples; s++)
		for(int i = 0; i < input_num_nodes; i++)
			input_file >> input_alpha[i][s];


	// Input q frequencies
	input_file >> skip_string;
	for(int s = 0; s < num_samples; s++)
		for(int i = 0; i < input_num_nodes; i++)
			input_file >> input_q[i][s];

	int output_num_samples;
	output_file >> skip_string >> output_num_samples;
	if(output_num_samples >= MAX_SAMPLES)
		fprintf(stderr, "Number of samples in %s is equal to %d that is >= MAX_SAMPLES. Please change MAX_SAMLES defined in AnalyzeOutputs.cpp file. \n", output_file_location_arg, output_num_samples);
	assert(output_num_samples < MAX_SAMPLES);

	if(output_num_samples != num_samples)
		fprintf(stderr, "Something is wrong. Number of samples is different in %s and %s", input_file_location_arg, output_file_location_arg);
	assert(num_samples == output_num_samples);



	// Output alpha frequencies
	output_file >> skip_string;
	for(int s = 0; s < num_samples; s++)
		for(int i = 0; i < output_num_nodes; i++)
			output_file >> output_alpha[i][s];
		
	// Output q frequencies
	output_file >> skip_string;
	for(int s = 0; s < num_samples; s++)
		for(int i = 0; i < output_num_nodes; i++)
			output_file >> output_q[i][s];
	
		
	input_file >> skip_string >> num_mutations;
	if(num_mutations >= MAX_MUT)
		fprintf(stderr, "Number of mutations in %s is equal to %d that is >= MAX_MUT. Please change MAX_MUT defined in AnalyzeOutputs.cpp file. \n", input_file_location_arg, num_mutations);
	assert(num_mutations < MAX_MUT);


	
	int output_num_mutations;
	output_file >> skip_string >> output_num_mutations;
	if(output_num_mutations != num_mutations)
		fprintf(stderr, "Something is wrong. Num_mutaions is different in %s and %s", input_file_location_arg, output_file_location_arg);
	assert(num_mutations == output_num_mutations);


	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			cost[i][j] = 0;
	
	
	vector<vector<int> > input_values;
	vector<vector<int> > result_values;

	input_values.resize(n);
	result_values.resize(n);

	for(int i = 0; i < num_mutations; i++)
	{
		input_file >> read_int;
		input_values[read_int].push_back(i);
		input_mut_locations[i] = read_int;
	}
	
	for(int i = 0; i < num_mutations; i++)
	{
		output_file >> read_int;
		result_values[read_int].push_back(i);
		output_mut_locations[i] = read_int;
	}
	
	
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			
			for(int k = 0; k < input_values[i].size(); k++)
			{
				current = input_values[i].at(k);
				
				for(int t = 0; t < result_values[j].size(); t++)
				{
					if(current == result_values[j].at(t))
					{
						cost[i][j]++;
						break;
					}
				}
			}
			
			cost[i][j] = -(input_values[i].size() + result_values[j].size() - 2 * cost[i][j]);
			
		}
	}
	
	input_file.close();
	output_file.close();
	
}


/*******************************************************************************************************************************************/



void  ReadGammaMatrix(int gamma_matrix[MAX_NODES][MAX_NODES], int num_nodes, int tree_index)
{
	
	ostringstream gamma_location;
	gamma_location << gamma_adj_folder_location << "/GammaMatrix" << num_nodes << ".txt";	
	string gamma_matrix_location = gamma_location.str();
	
	ifstream gamma_matrix_file(gamma_matrix_location.c_str());
	if(gamma_matrix_file == NULL)
		fprintf(stderr, "Can not open gamma file %s \n", gamma_matrix_location.c_str());
	assert(gamma_matrix_file != NULL);
	
	int num_trees;
	gamma_matrix_file >> skip_int >> num_trees;
	
	if( (tree_index < 0) || (tree_index >= num_trees) )
	{
		fprintf(stderr, "ERROR! In reading gamma matrix. Tree index out of range!\n");
		assert( (tree_index >=0) && (tree_index < num_trees) );
	}
	
	for(int skip_tree = 0; skip_tree < tree_index; skip_tree++)
		for(int i = 0; i < num_nodes; i++)
			for(int j = 0; j < num_nodes; j++)
				gamma_matrix_file >> skip_int;
	
	
	for(int i=0; i<num_nodes;i++)
		for(int j=0; j<num_nodes; j++)
			gamma_matrix_file >> gamma_matrix[i][j];
				
	gamma_matrix_file.close();
}






void  ReadAdjacencyMatrix(int adjacency_matrix[MAX_NODES][MAX_NODES], int num_nodes, int tree_index)
{

	ostringstream adjacency_location;
	adjacency_location << gamma_adj_folder_location << "/AdjacencyMatrix" << num_nodes << ".txt";	
	string adjacency_matrix_location = adjacency_location.str();

	ifstream adjacency_matrix_file(adjacency_matrix_location.c_str());
	if(adjacency_matrix_file == NULL)
		fprintf(stderr, "Can not open adjacency file %s \n", adjacency_matrix_location.c_str());
	assert(adjacency_matrix_file != NULL);
	
	
	int num_trees;
	
	adjacency_matrix_file >> skip_int >> num_trees;	

	if( (tree_index < 0) || (tree_index >= num_trees) )
	{
		fprintf(stderr, "ERROR! In reading adjacency matrix. Tree index out of range!\n");
		assert( (tree_index >=0) && (tree_index < num_trees) );
	}
		



	int num_edges;	
	for(int skip_tree = 0; skip_tree < tree_index; skip_tree++)
	{	
		adjacency_matrix_file >> num_edges;
		for(int i = 0; i < num_edges; i++)
			adjacency_matrix_file >> skip_int >> skip_int;
	}



	int read_vertex_1, read_vertex_2;
	adjacency_matrix_file >> num_edges;	
	for(int i=0; i<num_edges; i++)
	{
		adjacency_matrix_file >> read_vertex_1 >> read_vertex_2;
		
		adjacency_matrix[read_vertex_1][read_vertex_2] = 1;
		adjacency_matrix[read_vertex_2][read_vertex_1] = 1;
	}
	

	adjacency_matrix_file.close();
}



/*******************************************************************************
	
	This function returns number of misplaced mutations.
	If mutation i is on node j in original tree but not 
	at node xy[i] then it is misplaced.

********************************************************************************/

int FirstMeasure()
{

	int num_misplaced_mutations = 0;
	for(int i = 0; i < num_mutations; i++)
	{
		if(output_mut_locations[i] != xy[input_mut_locations[i]])
		{
			num_misplaced_mutations++;
		}
	}

	return num_misplaced_mutations;	
}



/******************************************************************************

	This function returns average error among frequencies assigned to nodes.
	For example, if node i has frequency 0.15 and its corresponding node in
	output tree has frequency 0.12 then they introduced 0.03 error. 
	This number is normalized by num_nodes * num_samples, where num_nodes 
	is number of nodes in smaller tree.

******************************************************************************/

double SecondMeasure()
{
	
	double score = 0;
	
	if(input_num_nodes < output_num_nodes)
	{
		for(int s=0; s<num_samples; s++)
			for(int i=0; i< input_num_nodes; i++)
				score += fabs(output_q[xy[i]][s] - input_q[i][s]);
		
		score = score/(num_samples*input_num_nodes);
	}
	else
	{
		for(int s=0; s<num_samples; s++)
			for(int i=0; i< output_num_nodes; i++)
				score += fabs(output_q[i][s] - input_q[yx[i]][s]);
		
		score = score/(num_samples*output_num_nodes);
	}
	
	
	return score;
}



/*******************************************************************************

	In this measure we are counting number of non preserved adjacencies.
	Namely, if (i,j) is an edge in original tree then we are testing 
	whether (xy[i], xy[j]) is and edge in output tree.
	Here we also compare only smaller of input and output trees against
	the other one.

*******************************************************************************/

int ThirdMeasure()
{

	int score = 0;

	if(input_num_nodes > output_num_nodes)
	{
		for(int i = 0; i < output_num_nodes; i++)
			for(int j = 0; j < output_num_nodes; j++)
			{
				if(output_adjacency_matrix[i][j] && output_gamma_matrix[i][j] )
					if( !(input_adjacency_matrix[yx[i]][yx[j]] && input_gamma_matrix[yx[i]][yx[j]]) )
						score += 1;
			}
	}
	
	else
	{
		for(int i = 0; i < input_num_nodes; i++)
			for(int j = 0; j < input_num_nodes; j++)
			{
				if(input_adjacency_matrix[i][j] && input_gamma_matrix[i][j])
					if( !(output_adjacency_matrix[xy[i]][xy[j]] && output_gamma_matrix[xy[i]][xy[j]]) )
						score += 1;
			}	
	}

	return score;
}


void FourthMeasure()
{
	for(int mut1=0; mut1 < num_mutations; mut1++)
		for(int mut2 = mut1+1; mut2 < num_mutations; mut2++)
		{
			int loc1_input  = input_mut_locations[mut1];
			int loc2_input  = input_mut_locations[mut2];
			int loc1_output = output_mut_locations[mut1];
			int loc2_output = output_mut_locations[mut2];
			
			if (loc1_input == loc2_input)
			{
				tot_same_node_input += 1;
				if (loc1_output == loc2_output)
					tot_same_node_preserved += 1;
			} 
				
			else if ( input_gamma_matrix[loc1_input][loc2_input] && input_adjacency_matrix[loc1_input][loc2_input])
			{
				tot_direct_neigbours_input += 1;
				if (output_gamma_matrix[loc1_output][loc2_output] && output_adjacency_matrix[loc1_output][loc2_output])
	tot_direct_neighbors_preserved += 1;
			}

			else if ( input_gamma_matrix[loc2_input][loc1_input] && input_adjacency_matrix[loc2_input][loc1_input])
			{
				tot_direct_neigbours_input += 1;
				if ( output_gamma_matrix[loc2_output][loc1_output] && output_adjacency_matrix[loc2_output][loc1_output])
	tot_direct_neighbors_preserved += 1;
			}

			else if (input_gamma_matrix[loc1_input][loc2_input])
			{
				tot_non_direct_ancestors_input += 1;
				if (output_gamma_matrix[loc1_output][loc2_output])
					tot_non_direct_ancestors_preserved += 1;
			}
			
			else if (input_gamma_matrix[loc2_input][loc1_input])
			{
				tot_non_direct_ancestors_input += 1;
				if (output_gamma_matrix[loc2_output][loc1_output] == 1)
					tot_non_direct_ancestors_preserved += 1;
			}
			
			else if( (input_gamma_matrix[loc1_input][loc2_input] == 0) && (input_gamma_matrix[loc2_input][loc1_input] == 0) )
			{
				tot_non_ancestors_input += 1;
				if( (output_gamma_matrix[loc1_output][loc2_output] == 0) && (output_gamma_matrix[loc2_output][loc1_output] == 0) )
	tot_non_ancestors_preserved += 1;
			}

		}
}




void ShowTableComparison()
{


	cout << "ORIGINAL TREE - RUNNING TREE MAPPINGS" << endl << endl;
	for (int i = 0; i < n; i++)
		cout << setw(5) << i << setw(5) << xy[i] << endl << endl;
	

	cout << endl << "INPUT NODES ALPHA " << endl << endl;
	for(int s=0; s<num_samples; s++)
	{
		for(int i=0; i<n; i++)
				cout << setw(10) << input_alpha[i][s] << " ";
		
		cout << endl;
	}
	
	cout << endl << "OUTPUT NODES ALPHA " << endl << endl;
	for(int s=0; s<num_samples; s++)
	{
		for(int i=0; i<n; i++)
				cout << setw(10)  << output_alpha[xy[i]][s] << " ";
		
		cout << endl;
	}
	


	cout << endl << endl << endl;


	
	cout << "INPUT NODES q" << endl << endl;
	for(int s=0; s<num_samples; s++)
	{
		for(int i=0; i<n; i++)
			cout << setw(10) << input_q[i][s] << " ";
		
		cout << endl;
	}
	
	cout << endl << "OUTPUT NODES q " << endl << endl;
	for(int s=0; s<num_samples; s++)
	{
		for(int i=0; i<n; i++)
			cout << setw(10) << output_q[xy[i]][s] << " ";
		
		cout << endl;
	}
	
	cout << endl << endl << endl;
}






void InitializeToZeroes()
{
	for(int s=0; s<MAX_SAMPLES; s++)
		for(int i=0; i<MAX_NODES; i++)
		{
			output_q[i][s]=0;
			input_q[i][s]=0;
		}
	
	for(int s=0; s<MAX_SAMPLES; s++)
		for(int i=0; i<MAX_NODES; i++)
		{
			output_alpha[i][s]=0;
			input_alpha[i][s]=0;
		}
	

	for(int i=0; i<MAX_NODES; i++)
		for(int j=0; j<MAX_NODES; j++)
		{
			input_adjacency_matrix[i][j]=0;
			output_adjacency_matrix[i][j]=0;
		}

	tot_same_node_input			=	0;
	tot_same_node_preserved			=	0;
	tot_direct_neigbours_input		=	0;		
	tot_direct_neighbors_preserved		=	0;
	tot_non_direct_ancestors_input		=	0;
	tot_non_direct_ancestors_preserved	=	0;
	tot_non_ancestors_input			=	0;
	tot_non_ancestors_preserved		=	0;

}




int max(int a, int b)
{
	if(a>b) return a;
	return b;
}








void init_labels()
{
    memset(lx, 0, sizeof(lx));
    memset(ly, 0, sizeof(ly));
    for (int x = 0; x < n; x++)
        for (int y = 0; y < n; y++)
            lx[x] = max(lx[x], cost[x][y]);
}


void augment()                         //main function of the algorithm
{
    if (max_match == n) return;        //check wether matching is already perfect
    int x, y, root;                    //just counters and root vertex
    int q[N], wr = 0, rd = 0;          //q - queue for bfs, wr,rd - write and read
	//pos in queue
    memset(S, false, sizeof(S));       //init set S
    memset(T, false, sizeof(T));       //init set T
    memset(previous, -1, sizeof(previous));    //init set previous - for the alternating tree
    for (x = 0; x < n; x++)            //finding root of the tree
        if (xy[x] == -1)
        {
            q[wr++] = root = x;
            previous[x] = -2;
            S[x] = true;
            break;
        }
	
    for (y = 0; y < n; y++)            //initializing slack array
    {
        slack[y] = lx[root] + ly[y] - cost[root][y];
        slackx[y] = root;
    }
	
	
	//second part of augment() function
    while (true)                                                        //main cycle
    {
        while (rd < wr)                                                 //building tree with bfs cycle
        {
            x = q[rd++];                                                //current vertex from X part
            for (y = 0; y < n; y++)                                     //iterate through all edges in equality graph
                if (cost[x][y] == lx[x] + ly[y] &&  !T[y])
                {
                    if (yx[y] == -1) break;                             //an exposed vertex in Y found, so
					//augmenting path exists!
                    T[y] = true;                                        //else just add y to T,
                    q[wr++] = yx[y];                                    //add vertex yx[y], which is matched
					//with y, to the queue
                    add_to_tree(yx[y], x);                              //add edges (x,y) and (y,yx[y]) to the tree
                }
            if (y < n) break;                                           //augmenting path found!
        }
        if (y < n) break;                                               //augmenting path found!
		
        update_labels();                                                //augmenting path not found, so improve labeling
        wr = rd = 0;
        for (y = 0; y < n; y++)
			//in this cycle we add edges that were added to the equality graph as a
			//result of improving the labeling, we add edge (slackx[y], y) to the tree if
			//and only if !T[y] &&  slack[y] == 0, also with this edge we add another one
			//(y, yx[y]) or augment the matching, if y was exposed
            if (!T[y] &&  slack[y] == 0)
            {
                if (yx[y] == -1)                                        //exposed vertex in Y found - augmenting path exists!
                {
                    x = slackx[y];
                    break;
                }
                else
                {
                    T[y] = true;                                        //else just add y to T,
                    if (!S[yx[y]])
                    {
                        q[wr++] = yx[y];                                //add vertex yx[y], which is matched with
						//y, to the queue
                        add_to_tree(yx[y], slackx[y]);                  //and add edges (x,y) and (y,
						//yx[y]) to the tree
                    }
                }
            }
        if (y < n) break;                                               //augmenting path found!
    }
	
    if (y < n)                                                          //we found augmenting path!
    {
        max_match++;                                                    //increment matching
        //in this cycle we inverse edges along augmenting path
        for (int cx = x, cy = y, ty; cx != -2; cx = previous[cx], cy = ty)
        {
            ty = xy[cx];
            yx[cy] = cx;
            xy[cx] = cy;
        }
        augment();                                                      //recall function, go to step 1 of the algorithm
    }
}//end of augment() function



void update_labels()
{
    int x, y, delta = INF;             //init delta as infinity
    for (y = 0; y < n; y++)            //calculate delta using slack
        if (!T[y])
            delta = min(delta, slack[y]);
    for (x = 0; x < n; x++)            //update X labels
        if (S[x]) lx[x] -= delta;
    for (y = 0; y < n; y++)            //update Y labels
        if (T[y]) ly[y] += delta;
    for (y = 0; y < n; y++)            //update slack array
        if (!T[y])
            slack[y] -= delta;
}




void add_to_tree(int x, int prevx)
//x - current vertex,prevx - vertex from X before x in the alternating path,
//so we add edges (prevx, xy[x]), (xy[x], x)
{
    S[x] = true;                    //add x to S
    previous[x] = prevx;                //we need this when augmenting
    for (int y = 0; y < n; y++)    //update slacks, because we add new vertex to S
        if (lx[x] + ly[y] - cost[x][y] < slack[y])
        {
            slack[y] = lx[x] + ly[y] - cost[x][y];
            slackx[y] = x;
        }
}




int hungarian()
{
    int ret = 0;                      //weight of the optimal matching
    max_match = 0;                    //number of vertices in current matching
    memset(xy, -1, sizeof(xy));
    memset(yx, -1, sizeof(yx));
    init_labels();                    //step 0
    augment();                        //steps 1-3
    for (int x = 0; x < n; x++)       //forming answer there
        ret += cost[x][xy[x]];
    return ret;
}



