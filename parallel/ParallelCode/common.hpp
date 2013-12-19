#include<stdio.h>
#include "omp.h"
#include<math.h>
#include <time.h>
#include<stdlib.h>
#include <vector>
#include<sstream>
#include<fstream>
#include<iostream>
using namespace std;

#define THETA 0.5
#define XDIM 16384
#define YDIM 16384
#define ZDIM 16384
#define XMIN 0.0
#define YMIN 0.0
#define MASS 10.0
#define INITVELOCITY 0.0
#define START 0
#define END 1

#define NW 1
#define NE 2
#define SW 3
#define SE 4
#define INVALID -1
#define G 6.673e-11
#define TIMESTEP 10

