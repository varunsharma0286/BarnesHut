#!/bin/bash

g++ -fopenmp -g -o Barnes Main.cpp -lm;

if [ ! -d output ]; then
	mkdir output
fi

