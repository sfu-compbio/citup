#include <stdio.h>
#include <string.h>
#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <ilconcert/iloenv.h>
#include <ilcplex/cplex.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <vector>
#include <assert.h>
#include "tclap/CmdLine.h"

#define MAX_KMEANS_ITERATIONS 5000
#define MAX_CLUSTERS 500
#define MAX_MUTATIONS 2000
#define MAX_SAMPLES 10
#define MAX_TREES 50
#define MAX_NODES 25
#define CPLEX_TIMELIM 1
#define CPLEX_TRELIM 4
#define CPLEX_THREADS 18
 
using namespace std;

string frequencies_file_location;
string output_file_location;
string clustered_mut_file_location;
string gamma_adj_folder_location;
string original_freq_file_location;
string objective_values_file_location;
double f[MAX_CLUSTERS][MAX_SAMPLES];
double original_f[MAX_MUTATIONS][MAX_SAMPLES];
double skip_double;
string skip_string;
int    skip_int;
int    num_clusters;
int    *mutations_clusters;
int    *num_points_in_cluster;
int    gamma_matrix[MAX_NODES][MAX_NODES];
int    num_nodes;
int    num_mutations;
int    tree_index;
int    num_trees;
int    num_samples;
double cplex_time_limit;
int    cplex_trelimit;
int    cplex_threads;
bool   node_has_mutation;
bool   quadratic_error;

void Read_clustered_mut_file(const char *clustered_mut_file_location);
int  SolveInstance(int tree_index, const char *output_file_location);




int main(int argc, char** argv){
	
	
	TCLAP::CmdLine cmd("Command...", ' ', "1");
	
	TCLAP::ValueArg<int> number_nodes("n", "numberNodes", "Defines the number of nodes our tree has", false, 5, "integer", cmd);
	TCLAP::ValueArg<int> index_of_tree("t", "treeIndex", "Defines tree index", false, 7, "integer", cmd);
	TCLAP::ValueArg<string> gamma_adj_folder_arg("g","gamma_folder_location","Defines location of folder containing gamma matrices",false,"./GammaAdjMatrices/GammaAdjMatrices_Full","string",cmd);
	TCLAP::ValueArg<string> clustered_mut_file_arg("i","clustered_mutations_file","Defines exact name of input file that contains clustered mutations",false,"./Clustered_mut.txt","string",cmd);
	TCLAP::ValueArg<string> objective_values_file_arg("v","objective_values_file","Defines exact location offile that objective values are printed to",false,"","string",cmd);	
	TCLAP::ValueArg<string> output_file_arg("o","output_file_location","Defines exact name of file where output is to be placed",false,"./Results.txt","string",cmd);

	TCLAP::ValueArg<double> cplex_time_limit_arg("d","time_limit","Time limit for cplex in hours",false,CPLEX_TIMELIM,"real",cmd);
	TCLAP::ValueArg<int> cplex_threads_arg("r", "cplex_Threads", "Defines number of threads for Cplex", false, CPLEX_THREADS, "integer", cmd);
	TCLAP::ValueArg<int> cplex_trelimit_arg("m", "cplex_Memory", "Defines max memory in Gb used by Cplex", false, CPLEX_TRELIM, "integer", cmd);

	TCLAP::SwitchArg quadratic_error_arg("q","qerr","Quadratic Error",cmd);
	TCLAP::SwitchArg node_has_mutation_arg("p","node_has_mutation", "True if every internal node has at least one mutation", cmd);
	
	
	cmd.parse(argc,argv);
	
	
	num_nodes			=	number_nodes.getValue();
	tree_index			=	index_of_tree.getValue();
	clustered_mut_file_location	=	clustered_mut_file_arg.getValue();
	gamma_adj_folder_location	=	gamma_adj_folder_arg.getValue();
	output_file_location		=	output_file_arg.getValue();
	cplex_time_limit 		= 	cplex_time_limit_arg.getValue();
	cplex_trelimit			=	cplex_trelimit_arg.getValue();
	cplex_threads			=	cplex_threads_arg.getValue();
	quadratic_error                 =       quadratic_error_arg.getValue();
	node_has_mutation		=	node_has_mutation_arg.getValue();
	objective_values_file_location  =	objective_values_file_arg.getValue();

	Read_clustered_mut_file(clustered_mut_file_location.c_str());


	ostringstream gamma_location;
	gamma_location << gamma_adj_folder_location << "/GammaMatrix" << num_nodes << ".txt";

	string gamma_matrix_file_location = gamma_location.str();	
	
	ifstream gamma_matrix_file(gamma_matrix_file_location.c_str());
	if (gamma_matrix_file == NULL)
		fprintf(stderr, "Problem with opening file %s \n", gamma_matrix_file_location.c_str());
	assert(gamma_matrix_file != NULL);


	gamma_matrix_file >> num_nodes >> num_trees;
	
	if( (tree_index < 0) || (tree_index >= num_trees) )
	{
		fprintf(stderr, "ERROR! Tree index out of range!\n");
		assert((tree_index >= 0) && (tree_index < num_trees));
	}
	
	for(int skip_tree = 0; skip_tree < tree_index; skip_tree++)
		for(int i = 0; i < num_nodes;i++)
			for(int j = 0; j < num_nodes; j++)
				gamma_matrix_file >> skip_int;
	
	
	for(int i = 0; i < num_nodes; i++)
		for(int j = 0; j < num_nodes; j++)
			gamma_matrix_file >> gamma_matrix[i][j];
	
	
	gamma_matrix_file.close();
	

	SolveInstance(tree_index, output_file_location.c_str());
	
	return 0;
	
}




int SolveInstance(int tree_index, const char *output_file_location)
{
	
	ofstream output_file(output_file_location);
	if(output_file == NULL)
		fprintf(stderr, "Problems with creating file %s \n.", output_file_location);
	assert(output_file != NULL);

		
	IloEnv env;
	IloModel model(env);
	
	IloArray<IloIntVarArray> delta(env, num_clusters);
	for(int i=0; i<num_clusters;i++)
	{
		delta[i] = IloIntVarArray(env, num_nodes, 0, 1);
	}
	

	/******
	
		This is constraint that every cluster is assigned to
	    	exactly one node.
	
	*******/
	
	for(int i=0; i<num_clusters; i++)
	{
		IloIntExpr first_deltaConstraint(env);
		for(int j=0; j<num_nodes; j++)
		{ 
			first_deltaConstraint += delta[i][j];
		}
		
		model.add(first_deltaConstraint == 1);
	}


	/******
	
		This is constraint that every node, except possibly root is
		assigned at least one cluster.
	
	*******/
	if(node_has_mutation)
	{
		for(int i=1; i<num_nodes; i++)
		{
			IloIntExpr second_deltaConstraint(env);
			for(int j=0; j<num_clusters; j++)
			{
				second_deltaConstraint += delta[j][i];
			}
		
			model.add(second_deltaConstraint >= 1);
		}
	}

	
	IloArray<IloFloatVarArray> alpha(env, num_nodes);
	for(int i=0; i<num_nodes; i++)
	{
		alpha[i] = IloFloatVarArray(env, num_samples, 0, 1);
	}
	
	for(int s=0; s<num_samples; s++)
	{
		IloExpr sum_nodes_frequencies_constraint(env);
		
		for(int i=0; i<num_nodes; i++)
		{	
				sum_nodes_frequencies_constraint += alpha[i][s];	
		}
	 	
		model.add(sum_nodes_frequencies_constraint == 1);
	}
	
/*	
	IloArray<IloArray<IloFloatVarArray> > x(env, num_clusters), y(env, num_clusters);
*
*/

	IloArray<IloArray<IloFloatVarArray> > x(env, num_mutations), y(env, num_mutations);
 
/*
 * for(int i=0;i<num_clusters;i++)
*/	
	for(int i=0; i<num_mutations; i++)
	{
		x[i] = IloArray<IloFloatVarArray>(env,num_nodes);
		y[i] = IloArray<IloFloatVarArray>(env,num_nodes);
	}
	
/*	
 *	for(int i=0; i<num_clusters;i++)
*/
	for(int i=0; i<num_mutations; i++)
	{
		for(int j=0; j<num_nodes; j++)
		{
			x[i][j] = IloFloatVarArray(env, num_samples, 0, 1);
			y[i][j] = IloFloatVarArray(env, num_samples, 0, 1);
		}
	}
	

	/*****
	
		Now the q[i][s] denotes the frequecy of all nodes in subtree rooted at i
	
	*****/

	IloArray<IloFloatVarArray> q(env, num_nodes);
	for(int i=0; i<num_nodes; i++)
	{
		q[i] = IloFloatVarArray(env, num_samples, 0, 1);
	}

	
	for(int s=0; s<num_samples; s++)
	{
		for(int j=0; j<num_nodes; j++)
		{
			
				IloExpr qConstraint(env);
	
				qConstraint -= q[j][s];
			
				for(int k=0; k<num_nodes; k++)
				{
					if(gamma_matrix[j][k] == 1)
					{
						qConstraint += alpha[k][s];
					}
				}

				model.add(qConstraint == 0);
		}
	}

	

 	IloExpr objective(env);

/*
	for(int i=0; i<num_clusters; i++)
	{
		for(int j=0; j<num_nodes; j++)
		{
			for(int s=0; s<num_samples; s++)
			{
			
				model.add(x[i][j][s] >= (f[i][s]-q[j][s]));
				model.add(x[i][j][s] >= (q[j][s]-f[i][s]));
				model.add(y[i][j][s] >= (1)*(delta[i][j]-1) + x[i][j][s]);

				//model.add(y[i][j][s]>=(max(f[i][s],1-f[i][s])*(delta[i][j]-1)+x[i][j][s]));
				
				
				if (quadratic_error)
				{
					objective += num_points_in_cluster[i] * y[i][j][s] * y[i][j][s];
				}
				else
				{
					objective += num_points_in_cluster[i] * y[i][j][s];
				}
			}
		}
	}
*/



	for(int i=0; i<num_mutations; i++)
	{
		for(int j=0; j<num_nodes; j++)
		{
			for(int s=0; s<num_samples; s++)
			{
		/*
				model.add(x[i][j][s] >= (f[i][s]-q[j][s]));
				model.add(x[i][j][s] >= (q[j][s]-f[i][s]));
		*/

				model.add(x[i][j][s] >= (original_f[i][s]-q[j][s]));
				model.add(x[i][j][s] >= (q[j][s]-original_f[i][s]));

				model.add(y[i][j][s] >= (1)*(delta[mutations_clusters[i]][j]-1) + x[i][j][s]);

				/*****************************************************

				--------Oktay's idea not implemented------------------ 

					model.add(y[i][j][s]>=(max(f[i][s],1-f[i][s])*(delta[i][j]-1)+x[i][j][s]));

				******************************************************/
				
				if (quadratic_error)
				{
					objective +=  y[i][j][s] * y[i][j][s];
				}
				else
				{
					objective +=  y[i][j][s];
				}
			}
		}
	}
	
	IloObjective objectiveExpression = IloMinimize(env, objective);
	
	model.add(objectiveExpression);
	
	IloCplex cplex(model);
	
	double hours_seconds = 3600.0;
	
	cplex.setParam(IloCplex::TiLim,(int)(cplex_time_limit * hours_seconds));
	cplex.setOut(env.getNullStream());
	cplex.setWarning(env.getNullStream());
	cplex.setParam(IloCplex::EpGap, 0.002);	
	// Set number of threads
	cplex.setParam(IloCplex::Threads, cplex_threads);
	
	// Set tree memory limit
	cplex.setParam(IloCplex::TreLim, cplex_trelimit*1024);
	
	// Deterministic when using multiple threads
	cplex.setParam(IloCplex::ParallelMode, 1);
	

	double startTime = cplex.getCplexTime();
	cplex.solve();
	double endTime = cplex.getCplexTime();
	
	double cplexHours = (endTime - startTime) / hours_seconds;
	
	
	IloCplex::CplexStatus status = cplex.getCplexStatus();
	
	int width = 20;
	output_file << setw(width) << "Objective_value " << cplex.getObjValue() << endl;
	output_file << setw(width) << "Num_clusters " << num_clusters << endl;
	output_file << setw(width) << "Tree_index " << tree_index << endl;
	output_file << setw(width) << "Num_nodes " << num_nodes << endl;
	output_file << setw(width) << "Num_samples " << num_samples << endl;
	output_file << setw(width) << "Num_mutations " << num_mutations << endl;
	output_file << endl;
	output_file << setw(width) << "Cplex_status " << status << endl;
	output_file << setw(width) << "Cplex_hours " << cplexHours << endl;
	output_file << setw(width) << "Cplex_memory_Gb " << cplex_trelimit << endl;
	output_file << setw(width) << "Cplex_threads " << cplex_threads << endl;
	output_file << setw(width) << "Quadratic_error ";
	if(quadratic_error)
		 output_file << "true" << endl;
	else 
		 output_file << "false" << endl;


	output_file << endl << endl;		


	output_file << "ALPHA_frequencies" << endl << endl;
	for(int sample=0; sample<num_samples; sample++)
	{
		for(int i=0; i<num_nodes; i++)
		{
			{
				output_file << setw(10) << cplex.getValue(alpha[i][sample]) << " ";
			}
		}
		
		output_file << endl;
	}


	output_file << endl << endl;

	
	output_file << "Q_frequencies" << endl << endl;
	for(int sample=0; sample<num_samples; sample++)
	{
		for(int i=0; i<num_nodes; i++)
		{
				output_file << setw(10) << cplex.getValue(q[i][sample]) << " ";	
		}
		
		output_file << endl;
	}


	output_file << endl << endl;



	int cluster_center_location[num_clusters];
	for(int i=0; i<num_clusters; i++)
	{
		for(int j=0; j<num_nodes; j++)
		{
			if(cplex.getValue(delta[i][j]) > 0.5)
			{
				cluster_center_location[i] = j;
				break;
			}
		}
	}



	for(int i= 0; i < num_mutations; i++)
		output_file << setw(3) << cluster_center_location[mutations_clusters[i]] << " ";
	

	output_file << endl << endl;

 	
	for (int i=0; i<num_mutations; i++)
	{
	
		output_file << setw(3) << cluster_center_location[mutations_clusters[i]] << " ";
	
		for (int s=0; s<num_samples; s++)
			output_file << setw(10) << original_f[i][s] << " ";
		
		output_file << endl;
	}
	
	// Print some of the most important information in
	// nicely formatted output file

	if(objective_values_file_location != "")
	{
	   
		ofstream scores_file;
	        scores_file.open(objective_values_file_location.c_str(), std::ofstream::out | std::ofstream::app);
	        assert(scores_file != NULL);
	       	scores_file << num_nodes << "_" << tree_index << " " << setw(10) << cplex.getObjValue() << " ";
	        scores_file <<  setw(10) << status << "      " << setw(10) << cplexHours << endl;
	        scores_file.close();
	                                                                                                                                   }	

}



void  Read_clustered_mut_file(const char *clustered_mut_file_location)
{
	
	ifstream clustered_mut_file(clustered_mut_file_location);
	if (clustered_mut_file == NULL)
		fprintf(stderr, "Problem with opening file %s \n", clustered_mut_file_location);
	assert(clustered_mut_file != NULL);

	string skip_string;
	
	clustered_mut_file >> skip_string >> num_mutations;
	clustered_mut_file >> skip_string >> num_clusters;
	clustered_mut_file >> skip_string >> num_samples;
	double error_frequency;
	clustered_mut_file >> skip_string >> error_frequency;
	
	for(int i=0; i<num_clusters; i++)
		for(int s=0; s<num_samples; s++)
			clustered_mut_file >> f[i][s];

	num_points_in_cluster =  new int[num_clusters+5];
	for(int i=0; i<num_clusters+5; i++)
		num_points_in_cluster[i]=0;

	mutations_clusters    =  new int[num_mutations+5];
	
	int read_int;
	for (int i=0; i<num_mutations; i++)
	{
		clustered_mut_file >> read_int;
		num_points_in_cluster[read_int]++;
		mutations_clusters[i] = read_int;
	}

	for(int i=0; i<num_mutations; i++)
		for(int s=0; s<num_samples; s++)		
			clustered_mut_file >> original_f[i][s];

	clustered_mut_file.close();
}

