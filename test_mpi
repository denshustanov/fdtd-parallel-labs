#!/bin/bash
mpic++ lab2_mpi.cpp -o lab2_mpi

iterations=$1
echo $iterations
for i in $(seq 1 $iterations)
do
	start=`date +%s%N`
	mpirun -np 4 ./lab2_mpi $2 $3 $4
	end=`date +%s%N`
	avg=$((avg+(end-start)/1000000))
done
echo 'avg = ' $((avg/iterations))
