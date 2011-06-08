#!/bin/bash
#
# Set the name of the job.
#$ -N gaussian
#
# Cluster queue to use
#$ -q lfj.q
#
# Make sure that the .e and .o file arrive in the 
#working directory
#$ -cwd
#
# Merge the standard out and standard error to one file
#$ -o screen.log -j y 
#
# My code is re-runnable
#$ -r y
#
# Mail options
#$ -M yan@mpi-hd.mpg.de -m bes 

./prop_1