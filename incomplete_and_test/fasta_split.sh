#!/bin/bash

#       fasta_split.sh
#       
#       2011 - Benjamin Tovar <scenesfromamemory4@gmail.com>
#       
#
# NAME OF THE PROGRAM: fasta_split.sh
# 
# DATE: 29/JUN/2011
# 
# AUTHOR: Benjamin Tovar
# 
# COMMENTS: This script will split a FASTA file that contains many sequences into
# Separated files containing one single sequence per file.
#
# I used this script to split a FASTA file containing ~16,500 sequences into
# ~16,500 independent files very quickly and simple.
#
################################################################################

# BEGINNING OF THE PROGRAM
echo
echo "This script relies on the program called \"csplit\" available in Linux/UNIX systems."

# Ask the user to type the name of the input FASTA file:
echo
echo -n "   Enter the name of the input file in FASTA format: "
	   read input_file

# Ask the user to type the name of the output FASTA files (note in line 73, 
# "--suffix="%02d.fa" means that every output file will have the extension "*.fa"
# if you like to use another extension (for example, you like that every output file have
# the extension *.fasta), just replace the ".fa" with ".fasta" this way -> "--suffix="%02d.fasta"

# The part of "02" in "--suffix="%02d.fa" means that every file will be named with two numbers 
# sorted by their occurrence in the original FASTA input file.

# For example: in an input file that contains 2 sequences and with "--suffix="%02d.fa"
# the output will be:
#
# outputfile-00.fa
# outputfile-01.fa
#
# For example: in an input file that contains 2 sequences and with "--suffix="%04d.fa"
# the output will be:
#
# outputfile-0000.fa
# outputfile-0001.fa
echo  	   
echo -n "   Enter the name of the output files: "
	   read output_file

# Ask the user to type the name of the output folder that will contain all the output files	  
echo 
echo -n "   Enter the name of the output folder: "
read dir
echo
echo "Creating directory called \"$dir\"."

# Create the output folder
mkdir $dir

# Copy the input FASTA file to the output folder
cp $input_file $dir

# Open to the output folder
cd $dir
echo
echo "   ...splitting FASTA file ...saving them in \"$dir\"."

# Splitting the input FASTA file (For some settings, read line 31 to 49)
csplit -z $input_file '/^>/' '{*}' --suffix="%02d.fa" --prefix=$output_file- -s

# Delete the input FASTA file in the output folder
rm $input_file
echo
echo "Printing output summary in a file called \"OUTPUT-SUMMARY.out\""

# Create the output summary
ls | sort | sed -e 's/OUTPUT-SUMMARY.out//g' > OUTPUT-SUMMARY.out 
echo
echo " ---- PROCESS DONE ---- "
echo
