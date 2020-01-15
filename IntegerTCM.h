//
//  IntegerTCM.h
//  NewEntropyCoding_Enc
//
//  Created by Chang Sun on 2014-01-08.
//  Copyright (c) 2014 University of Waterloo. All rights reserved.
//
#pragma once
#ifndef _INTEGERTCM_H_
#define _INTEGERTCM_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "global.h"
#include "ImageDependent.h"

// Assume the max for DCT coefficients as 12bits. Need to be checked. this is the a in TCM paper 
#define MaxAmp 4096 
#define LambdaDelta 0.1
#define MinLikelyhood  -1.e30
//(-(1<<30))
#define StartPointProb 0.1

typedef struct {
	int count;
	double AcumAbsAmp;  // sum of the values befor this need for likelyhood 
	double AcumSampNum; // number of samples befor this yi
	double prob; // this should be 'b' in the TCM paper
	double lambda;
	double likelyhood;
} tBucket;

void TCMprocessOneSequence(int *C, int Len, int* peak, double*prob, double *lambda, double*Yc);
void  ComputeLikelyhood(int thePoint, int N, tBucket *buck, int peak);
int FindStartPoint(tBucket *buck, int peak, int N);
double ComputeLambdaGivenYc(double Yc, double sumYi, double totalNumYi); // this is algorithm 2 in TCM paper
void fnTCM(double ** DCTcoef); // DCTcoef is 512*512s

#endif
