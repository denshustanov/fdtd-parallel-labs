#!/bin/bash
g++ lab1_sequential.cpp -o lab1_sequential

iterations=$1
echo $iterations
for i in $(seq 1 $iterations)
do
	start=`date +%s%N`
	./lab1_sequential $2 $3 $4
	end=`date +%s%N`
	avg=$((avg+(end-start)/1000000))
done
echo 'avg = ' $((avg/iterations))
