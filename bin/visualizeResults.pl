
#------------------------------------------------
#  Part of CITUP package
#------------------------------------------------

use strict;
use warnings;
use Getopt::Long;

sub usage
{
	print STDERR "\nUsage: perl ./visualizeResults.pl -r <results_directory> -g <gamma_adj_directory> -o <output_directory> [-n <num_of_nodes>]>\n\n";
	print STDERR "<results_directory>    : path to the directory that contains the CITUP results (result files should start with 'Results'\n";
	print STDERR "<gamma_adj_directory>  : path to the directory that contains the Gamma and Adjacency matrix files\n";
	print STDERR "<output_directory>     : path to a directory to be created for the output files\n";
	print STDERR "<num_of_nodes>         : (optional argument) If it is given, only the trees containing that many nodes are processed.\n\n";	
	exit(0);
}

if($#ARGV < 5)
{	usage(); }

my $resultsFolder;
my $treeFolder;
my $pngFolder;
my $nd = 0;
my $hlp = 0;
my $noPNG = 0;

my $result = GetOptions("results|r=s" => \$resultsFolder,
						"output|o=s" => \$pngFolder,
						"gamma|g=s" => \$treeFolder,
						"num|n=i" => \$nd,
						"nopng|p" => \$noPNG,
						"help|h" => \$hlp);

unless($result) {	usage(); }
if($hlp) {	usage(); }

my %bestScores;
opendir(DIR, $resultsFolder) || die "Can not open the directory $resultsFolder\n";
while(my $file = readdir(DIR))
{	
	if($file =~ /Results/)
	{	
		my $objective = 0.0;
		my $inputFile = $resultsFolder."/".$file;
		open(IFH, "<$inputFile") || die "can not open $inputFile\n";
		while(<IFH>)
		{
			my @parts = split;
			chomp @parts;
			if(/Objective_value/i)
			{	$objective = $parts[1];	}
			if(/Num_nodes/i)
			{	
				my $num = $parts[1];
				if(!defined $bestScores{$num} || $bestScores{$num} > $objective)
				{	$bestScores{$num} = $objective; }
				last;
			}
		}
		close(IFH);
	}
}
closedir(DIR);

mkdir($pngFolder, 0744) || die "can not create directory $pngFolder\n";
my $imagesFolder = $pngFolder."/images";
my $dotsFolder = $pngFolder."/dots";
mkdir($imagesFolder, 0744) || die "can not create directory $imagesFolder\n";
mkdir($dotsFolder, 0744) || die "can not create directory $dotsFolder\n";

my $htmlFile = $pngFolder."/index.html";
open(HFH, ">$htmlFile") || die "Can not open $htmlFile\n";
print HFH "<!DOCTYPE html>\n";
print HFH "<html><title>CITUP results</title><body>\n";

my $isInput = 0;
my $num = 1;
opendir(DIR, $resultsFolder) || die "Can not open the directory $resultsFolder\n";
while(my $file = readdir(DIR))
{	
	if($file =~ /Results/)
	{	$isInput = 0; }
	elsif($file =~ /Input/)
	{	$isInput = 1; }
	else
	{	next;	}

	my $numTree = 0;
	my $numSamples = 0;
	my $numNodes = 0;
	my $objective = 0.0;

	my @aFreq;
	my @qFreq;	
	my $inputFile = $resultsFolder."/".$file;
	open(IFH, "<$inputFile") || die "can not open $inputFile\n";
	while(<IFH>)
	{
		my @parts = split;
		chomp @parts;
		if(/Objective_value/i)
		{	$objective = $parts[1];	}
		elsif(/Tree_index/i)
		{	$numTree = $parts[1];	}
		elsif(/Num_nodes/i)
		{	$numNodes = $parts[1];	}
		elsif(/Num_samples/i)
		{	$numSamples = $parts[1];	}

		if(/ALPHA_freq/i)
		{
			my $dummy = <IFH>; 
			if($isInput) 
			{ $dummy = <IFH>; }
			
			for(my $i=0; $i<$numSamples; $i++)
			{
				my $str = <IFH>;
				chomp $str;
				push @aFreq, $str;
			}
		}

		if(/Q_freq/i && $objective < $bestScores{$numNodes} + 0.0001 && ($nd eq 0 || $numNodes eq $nd))
		{		
			my $dummy = <IFH>; 
			if($isInput) 
			{ 				
				print HFH "<h2>True tree</h2>\n";
				print HFH "<p>Nodes: $numNodes Index: $numTree</p>\n";
			}
			else
			{	
				print HFH "<h2>Solution ", $num, "</h2>\n";
				print HFH "<p>Nodes: $numNodes Index: $numTree Score: $objective</p>\n";
				$num+=1;
			}
			
			for(my $i=0; $i<$numSamples; $i++)
			{
				my $str = <IFH>;
				chomp $str;
				push @qFreq, $str;
			}
					
			my $treeFile = $treeFolder."/AdjacencyMatrix".$numNodes.".txt";
			open(TFH, "<$treeFile") || die "Can not open $treeFile\n";
			my @lines = <TFH>;
			chomp @lines;
			close(TFH);											

			for(my $i=0; $i<$numSamples; $i++)
			{				
				my $outputFile = "";
				my $pngFile = "";
				if($isInput) 
				{ 
					$pngFile = $imagesFolder."/true.".$i.".png"; 
					$outputFile = $dotsFolder."/true.".$i.".dot"; 
					print HFH "<img src=\"./images/true.", $i, ".png\" alt=\"./dots/true.", $i, ".dot\">\n";
				}
				else
				{
					$pngFile = $imagesFolder."/node".$numNodes."_tree".$numTree.".".$i.".png";
				 	$outputFile = $dotsFolder."/node".$numNodes."_tree".$numTree.".".$i.".dot";
					print HFH "<img src=\"./images/node",$numNodes,"_tree",$numTree,".",$i,".png\" alt=\"./dots/node",$numNodes,"_tree",$numTree,".",$i,".dot\">\n";
				}
				open(OFH, ">$outputFile") || die "can not open $outputFile\n";
				print OFH "digraph G{\n";
			
				my @fieldsA = split(' ', $aFreq[$i]);		
				my @fieldsQ = split(' ', $qFreq[$i]);
				for(my $j=0; $j<$numNodes; $j++)
				{
					my $color = "black";
					if($fieldsA[$j] < 0.01)
					{	$color = "gray"; }
					elsif($fieldsA[$j] > 0.5)
					{	$color = "tomato"; }
					
					print OFH $j, " [shape=box, color=", $color, ", fontsize=10, label=\"";						
					printf OFH "%1.2f\"];\n", $fieldsA[$j];
				} 	
				
				my @parts = split(' ', $lines[2+$numTree]);
				for(my $j=1; $j<$#parts; $j+=2)
				{
					print OFH $parts[$j], " -> ", $parts[$j+1], ";\n";
				}
				print OFH "}\n";
				close(OFH);
				
				unless($noPNG)
				{	system("dot -Tpng $outputFile -o $pngFile");	}
			}
		}
	}
	close(IFH);
}			
print HFH "</body></html>\n";
close(HFH);
closedir(DIR);
