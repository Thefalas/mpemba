#!/bin/sh
# usage: "./submit.sh 44" where 44 is the simulation number

#sbatch --exclusive ./Rough.sh 
sbatch -p normal --exclusive Rough.sh $1 
