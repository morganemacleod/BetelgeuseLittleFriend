#!/bin/bash
#SBATCH -J radvel           # Job name
#SBATCH -o radvel.o%j       # Name of stdout output file
#SBATCH -e radvel.e%j       # Name of stderr error file
#SBATCH -p icx      # Queue (partition) name
#SBATCH -N 1               # Total # of nodes 
#SBATCH -n 1             # Total # of mpi tasks
#SBATCH -t 12:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=morgan.macleod@cfa.harvard.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
####SBATCH -A myproject       # Allocation name (req'd if you have more than 1)

# Other commands must follow all #SBATCH directives...
module list
pwd
date

sh runscript.sh

