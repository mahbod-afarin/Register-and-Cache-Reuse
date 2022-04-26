#!/bin/bash
#  
#SBATCH 
#SBATCH -o lab1.out
#SBATCH -e lab1.err
#SBATCH -J lab1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=11000
#SBATCH --time=8:0:0

##################################
# DO NOT CHANGE CODES ABOVE THIS #
##################################

##################################
# your own scripts		 #
# to generate standard output	 #
##################################

### code starts here ###
./data.sh mydata.txt

### code ends here ###
