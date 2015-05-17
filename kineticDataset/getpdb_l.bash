#!/bin/bash

for i in `cat $1`
do
    getpdb $i
done
