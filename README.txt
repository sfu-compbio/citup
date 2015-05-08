
-- CITUP MANUAL --

CITUP (Clonality Inference in Tumors Using Phylogeny) is a bioinformatics tool that can be used to infer tumor heterogeneity using multiple samples from a single patient. Given mutational frequencies for each sample, CITUP uses an optimization based algorithm to find the evolutionary tree best explaining the data.

This manual contains all the necessary information to install and run CITUP. For questions and other inquiries please send an email to: smalikic@sfu.ca.

-- 1. SYSTEM REQUIREMENTS --

CITUP is mainly tested on 64-bit UNIX-based environments. It requires the following libraries and tools:
	
	make (version 3.81 or higher)
	g++ (GCC version 4.1.2 or higher)
	Boost library for C++ (version 1.53 or higher)
	Python (version 2.7)
	IBM ILOG CPLEX Optimization Studio
	
CITUP's clustering module also requires several Python libraries such as NumPy and ScikitLearn. These libraries can be downloaded from http://continuum.io/downloads. For visualization of the trees, you will also need Perl and Graphviz (http://www.graphviz.org/) .

-- 2. INSTALLATION --

To install CITUP, first go to the directory where the CITUP package is extracted. Change your directory to "./src" and modify the Boost and Cplex directories in "makefile" based on the absolute paths of these libraries in your system. Run make.

If the installation is successful, the following files will be placed under "./bin":
	
	AnalyzeOutputs.exe	
	CitupIter.exe
	CplexApplication.exe
	generatingInput.exe
	
In addition to these binaries, there should also be scripts named "citup_iter.py", "clustering.py", "pickBestTree.py" and "runCITUP.sh" in this directory.

-- 3. INPUT PREPARATION --

CITUP can be used with allelic frequencies of SNV and micro-indels calls as measured by deep sequencing.

The input file should start with the following three lines:

Num_mutations: <integer greater than or equal to 1>
Num_samples: <integer greater than or equal to 1>
Error_rate: <real value between (0.0-1.0)>

Followed by a space delimited <Num_mutations> x <Num_samples> size table of mutational frequencies. As an example, suppose you have the following mutations and the number of reads matching the reference and variant allele are depicted as "Reference" and "Variant" respectively for samples #1 and #2:
                               
Chr | Position | Mutation | Reference(#1) | Variant(#1) | Reference(#2) | Variant(#2)
1 | 10642787 | C->T | 436 | 351 | 721  | 513 
1 | 10696158 | C->A | 397 | 328 | 389  | 262 
2 | 10746481 | A->G | 250 | 199 | 1067 | 735 
9 | 27997115 | C->T | 574 | 423 | 1975 | 53  
9 | 28164689 | G->A | 780 | 150 | 643  | 104 

Then your input file should look like:

Num_mutations: 5
Num_samples: 2
Error_rate: 0.03

0.891995 0.831442
0.904828 0.804916
0.886414 0.81576
0.848546 0.0522682
0.322581 0.278447

Here, the mutational frequencies are calculated as: (2 x Variant)/(Reference + Variant). For instance, the first mutational frequency for sample 1 is 0.891995=(2 x 351)/(436 + 351) and for sample 2 is 0.831442=(2 x 513)/(721 + 513). Note that when converting read counts to mutational frequencies, it is advised that you should discard any mutation that is suspected or known to be homozygous (e.g. mutational frequency significantly greater than 1.00) or is located on a region with copy number not equal to 2 in any sample. Similarly, mutations on X and Y chromosomes should either be discarded or their frequencies should be calculated with respect to the gender of the patient.

Above, the error rate is used to determine how noisy the allelic frequencies are. For simulations, this is set to be the Gaussian noise that is used to generate the mutational frequencies from the given clonal frequencies. For real data, set this value low (0.03-0.05) for deep sequencing and higher (0.05-0.08) for lower coverage datasets.

-- 4. PRECOMPUTED TREES --

In addition to the mutational frequencies, CITUP requires the set of all possible tree topologies up to a certain number of nodes. These topologies are precomputed and included with the CITUP package under the directory ./GammaAdjMatrices. The location of this directory must be given as an argument to some of the executables (see below).

-- 5. RUNNING CITUP EXECUTABLES --

In this section we describe how to run the individual tools in the CITUP package and explain the most commonly used options. To get the full list of command line options for each *.exe file, run it with the argument "--help".

STEP1: Clustering

For large datasets, the input frequencies have to be clustered prior to processing with CPLEX in order to speed up the computation. We have provided a script named "clustering.py" under "./bin" to perform this step:

python ./clustering.py -i Frequencies.txt -o Clustered_mut.txt -c 5

where "Frequencies.txt" is the name of the input file containing the mutational frequencies as explained above. "Clustered_mut.txt" is the name for the output file. 5 is the number of clusters.

STEP2: CPLEX

In this step, we invoke the CPLEX optimizer to find the optimal tree topologies and clonal frequencies:

./CplexApplication.exe -n 4 -t 3 -i Clustered_mut.txt -o Results.txt -d 2 -r 20 -m 180 -p -q -g ../GammaAdjMatrices

The command above means that we are running CplexApplication.exe on the tree with index 3 (-t 3) from 4 node trees (-n 4). Maximum duration is 2 hours (-d 2) with up to 20 threads (-r 20) and at most 180 Gb of memory is used (-m 180). The option "-g ../GammaAdjMatrices" specifies the location where the precomputed trees are stored. "-p" option tells the program to require every node to be assigned at least one new mutation. "-q" option enables quadratic error (otherwise a linear error model is used).

-- 6. RUNNING CITUP IN BATCH MODE --

For your convenience, we have included a shell script named "runCITUP.sh" under "./bin" that can be used to run the programs described above on all trees up to a pre-determined number of nodes and report the best tree(s) given the input data. This is the preferred way to run CITUP. To run this script, first modify the parameters at the beginning of the script according to your system. Then run it as:

./runCITUP.sh ../running_directory 15

where "../running_directory" is a directory that contains a file named "Frequencies.txt" which is the file containing the mutational frequencies of the samples as described in the input preparation section. 15 is the number of clusters to be used in the clustering step. In the current distribution, maximum number of clusters allowed is 30. If this parameter is omitted or smaller than the maximum number of nodes allowed, then a separate clustering is performed for each tree size based on the number of nodes in the tree.

An example simulation dataset is included with this package under "./test". To demo CITUP, simply run:

./runCITUP.sh ../test

If this script runs ok, the best tree is reported to be "../test/Results_4_1.txt". This file contains the clonal frequencies (each row represents one sample) referred to as "ALPHA_frequencies" and Q frequencies which are the cumulative allele frequencies. At the end of this file, all input mutations are listed together with the first clone they emerge (i.e. the first column of each row contains the clone number). See the next section on how to visualize this data.

-- 7. RUNNING CITUP_ITER IN BATCH MODE --

Iterative citup can be run using the python script "citup_iter.py".  The script requires installation of the pypeliner module, downloadable at "http://bitbucket.org/dranew/pypeliner".

The following command will run iterative citup on the example data in ../test

python citup_iter.py --config citup_iter_config.ini ../test/Frequencies.txt BestTree.txt TreesTable.tsv

The file BestTree.txt will contain the optimal tree in the same format as for CITUP, and TreesTable.tsv will contain a table of objective values and BICs for all trees.  The configuration file citup_iter_config.ini contains the following configuration options:

min_nodes : minimum number of nodes for trees explored
max_nodes : maximum number of nodes for trees explored
max_children_per_node : maximum number of children per node for trees explored

-- 8. VISUALIZING THE RESULTS --

If you have Perl and Graphviz (http://www.graphviz.org/) installed, you can visualize the trees reported by CITUP using the script "visualizeResults.pl" located under "./bin". The following is an example usage:

perl ./visualizeResults.pl -r ../running_directory -g ../GammaAdjMatrices -o ../results_summary -n 4

where "../running_directory" is the path of the directory where the results are stored after running "runCitup.sh" and "../GammaAdjMatrices" is where the gamma and adjacency matrices are stored. The main output of this script will be a file called "index.html" under "../results_summary". The option "-n 4" tells the script to only visualize the best scoring trees that contain 4 nodes. If this option is omitted, all best scoring trees (from each tree size) will be visualized.

-- 9. EVALUATING CITUP --

SIMULATING INPUT: 

We also provide a tool to simulate a dataset in order to evaluate CITUP:

./generatingInput.exe -m 500 -s 5 -n 4 -e 0.05 -z 54783 -d 5.0 -i Input.txt -f Frequencies.txt -g ../GammaAdjMatrices/

This command will generate an input dataset with 500 mutations (-m 500) on 5 samples (-s 5) with allele frequency noise rate 0.05 (-e 0.05). The option "-n 4" determines the number of nodes in the tree (also the number of clones). The option "-d" can be used to set the Dirichlet parameter, which is the distribution used to produce the clonal frequencies. The option "-z" can be used to set the random number generator seed. 

The file "Frequencies.txt" is to be used in the clustering step (see above). File "Input.txt" contains information about the true tree topology and clonal frequencies. This file is later used for evaluation (see below).

ANALYZING THE RESULTS:

If the input is simulated as above, the output can be evaluated with the program AnalyzeOutputs.exe. The following command takes as input the file generated by generatingInput.exe (-i Input.txt) and also the output of CplexApplication.exe on the same dataset (-o Results.txt):

./AnalyzeOutputs.exe -i Input.txt -o Results.txt -g ../GammaAdjMatrices/

The output of this program is written on the standard output.

