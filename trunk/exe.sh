#!/bin/sh

for (( i=1; i<=10; i++ ))
do
    bin/bSSEIOMP $1 $2 $3 $4 $5 $6 
    echo 
done

