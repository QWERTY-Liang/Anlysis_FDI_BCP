#!/bin/bash
#SBATCH -n 8		# number of cores
#SBATCH -p compute	# compute queue
#SBATCH -t 4-00:00:00	# time (ddd-hh:mm:ss)
#SBATCH -J dn3v2d      # job name

# load up the correct modules

. /etc/profile.d/modules.sh
module load apps matlab

# call matlab non-interactively

matlab -nodisplay < dn_3v2d.m