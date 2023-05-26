#!/software/mambaforge/envs/Murat_scripts/bin/python

import argparse
import os
import sys
import subprocess
import textwrap
import time
import multiprocessing

parser = argparse.ArgumentParser(prog='python Phylo_2_genemap.py',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=textwrap.dedent('''\

Author: Murat Buyukyoruk

Associated lab: Wiedenheft lab

        Phylo_2_genemap help:

This script is developed to illustrate a phylogenetic tree (.newick) along side with a genemap that demonstrates 
the protein sequence alignment and domain prevelance.

R is required to generate plots with various CRAN packages (i.e., ggtree, gggenes, tidyverse, tidytree and their 
depencencies). The R script is generated to install CRAN packages if they were not available in the system, 

Syntax:

        python Phylo_2_genemap.py -t demo.nexus -o demo_tree_genemap -g genemap_dataframe.txt -a ANK_1 -c 
clade_info.txt

Example Dataframe for the genemap_dataframe.txt file (tab separated excel file is required):

        molecule	gene	start	end	strand	orientation
        NC_002967.9_2504	Length	1	934	forward	1
        NC_002967.9_2504	ANK_1	65	97	forward	1
        NC_002967.9_2504	ANK	247	312	forward	1
        NC_002967.9_2504	ANK	313	346	forward	1
        NC_002967.9_2504	ANK	347	380	forward	1
        NC_002967.9_2504	ANK	381	413	forward	1
        NC_002967.9_2504	ANK	543	575	forward	1
        NC_002967.9_2504	ANK	813	845	forward	1

Example Dataframe for the clade_info.txt file (tab separated excel file is required):

        Clade1	Clade2
        NZ_OX366336.1_932	NZ_CP025544.1_57
        NZ_OX366337.1_584	
        NZ_OX366343.1_632	
        NZ_OX366401.1_611	
        NZ_CP076228.1_128	
        NZ_JAATLD010000001.1_206	
        NZ_CP069053.1_339	
        NZ_VCEF01000099.1_2	
        NZ_VCEG01000064.1_7	
        NC_012416.1_91	
        NZ_CP042904.1_88	
        NZ_JAGKTH010000088.1_15	
        NZ_OX366369.1_350	

Phylo_2_genemap dependencies:

	R                                       refer to https://rstudio-education.github.io/hopr/starting.html

Input Paramaters (REQUIRED):
----------------------------
	-t/--tree		NEWICK		Specify a phylogenetic tree in .newick format.

	-o/--output		Output file	Specify a output file name for PNG and EPS files.

	-g/--genemap		Dataframe	Specify a genemap file that include the doamin name, start, stop, 
strand and orientation information for each domain occurs in protein accessions.

        -a/--anchor     	gene_name       Specify a gene name for aligning the proteins like an anchoring point.
	
	-c/--clade      	Dataframe       Specify a dataframe that includes the clades to color on phylogenetic 
tree. 

Basic Options:
--------------
	-h/--help		HELP		Shows this help text and exits the run.

      	'''))

parser.add_argument('-t', '--tree', required=True, type=str, dest='filename',
                    help='Specify a phylogenetic tree in .newick format.\n')
parser.add_argument('-o', '--output', required=True, dest='out',
                    help='Specify a output file name for PNG and EPS files.\n')
parser.add_argument('-g', '--genemap', required=True, dest='genemap',
                    help='Specify a genemap file.\n')
parser.add_argument('-a', '--anchor', required=True, dest='ank',
                    help='Specify a gene name for aligning the proteins.\n')
parser.add_argument('-c', '--clade', required=False, dest='clade', default= None,
                    help='Specify a dataframe that includes the clades to color on phylogenetic tree.\n')

results = parser.parse_args()
filename = results.filename
out = results.out
genemap = results.genemap
clade = results.clade
ank = results.ank

def spinning_cursor():
    while True:
        for cursor in '|/-\\':
            yield cursor

spinner = spinning_cursor()

def run():
    if clade == None:
        os.system("Rscript Phylo_2_genemap.R --tree " + filename + " --genemap " + genemap + " --anchor " + ank + 
" --output " + out + " > /dev/null")
    else:
        os.system("Rscript Phylo_2_genemap.R --tree " + filename + " --genemap " + genemap + " --anchor " + ank + 
" --clade " + clade + " --output " + out + " > /dev/null")


p = multiprocessing.Process(target=run)
p.start()
p.join(timeout=0)

while p.is_alive():
    sys.stdout.write("Generating Phylogeny and adding genemap panel " + next(spinner))
    sys.stdout.flush()
    time.sleep(0.1)
    sys.stdout.write('\r')

print("Raw plot is exported as PDF and SVG files and can be found in " + os.getcwd())

os.system("rm " + os.getcwd() + "/Rplots.pdf")

