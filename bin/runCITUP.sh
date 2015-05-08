#!/bin/sh

if [ $# -lt 1 ] 
then
	echo "Usage: $0 <working_directory> <num_clusters>"
	echo "NOTE: If <num_clusters> is omitted, tree-size clustering is enabled."
	echo "For more information see the README.txt file"
	exit
fi

# ------ change the following parameters according to your system -----
citup_directory=".."		# full path of the citup installation directory
cplex_time=1				# Max time allowed for CPLEX in hours
cplex_memory=10				# Max memory allowed for CPLEX in gb
cplex_threads=8				# Number of threads available to CPLEX
min_nodes=2					# Trees with less number of nodes will not be processed (min value is 2)
max_nodes=5					# Trees with more number of nodes will not be processed (max value is 10)
# ---------------------------------------------------------------------

# ------- do not modify the script beyond this line ----------

results_location=$1
num_clusters=0

if [ $# -gt 1 ]
then	
	num_clusters=$2
fi

if [ "$num_clusters" -lt "$max_nodes" ]
then
	num_clusters=0
fi

if [ "$num_clusters" -gt 30 ]
then
	num_clusters=30
fi

if [ "$min_nodes" -lt 2 ]
then
	echo "ERROR: min_nodes can not be set to a number less than 2!"
	exit
fi

if [ "$max_nodes" -gt 10 ]
then
	echo "ERROR: max_nodes can not be set to a number greater than 10!"
	exit
fi

if [ "$max_nodes" -lt "$min_nodes" ]
then
	echo "ERROR: max_nodes should be greater than min_nodes!"
	exit
fi

num_trees[2]=1
num_trees[3]=2
num_trees[4]=4
num_trees[5]=9
num_trees[6]=20
num_trees[7]=48
num_trees[8]=115
num_trees[9]=286
num_trees[10]=719

input_file="Clustered_mut.txt"
if [ "$num_clusters" -gt 1 ]; then
	python $citup_directory/bin/clustering.py -i $results_location/Frequencies.txt -o $results_location/Clustered_mut.txt -c $num_clusters	
fi

for num_nodes in $(seq $min_nodes $max_nodes)
do 			
	if [ "$num_clusters" -lt 1 ]; then
		python $citup_directory/bin/clustering.py -i $results_location/Frequencies.txt -o $results_location/Clustered_mut${num_nodes}.txt -c $num_nodes
		input_file="Clustered_mut${num_nodes}.txt"
	fi
		
	for (( tree_index = 0; tree_index < ${num_trees[$num_nodes]}; tree_index++ ))
	do
		$citup_directory/bin/CplexApplication.exe -i $results_location/$input_file -m $cplex_memory -d $cplex_time -q -p -r $cplex_threads -n $num_nodes -t $tree_index -o $results_location/Results_${num_nodes}_${tree_index}.txt -g $citup_directory/GammaAdjMatrices/
	done
done

# finally report the best tree(s)
python $citup_directory/bin/pickBestTree.py $results_location
