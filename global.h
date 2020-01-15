#pragma once
#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#define FILTER_ECEB 1
#define NUMBER_OF_THREADS 1
//#define ROWS 512  //512//1024//512//720//1024//512//512//1024//512///512//20//
//#define COLS  512   //512//1024//512//1024//512//1280//1024//512//512//1024//512//1024//512//512//1024
//#define NUMBLK  4096  //16384//32400//14400//16384//16384//4096//16384///4096//4096//

#define NUM_CTX_PRED_DC 12

#define DC_RESCALE 1
#define PRINT_EXT_INFO 1
#define BIAS_CMP_LOCO 0//1:use the bias computation in LOCO to correct the prediction in DC coding

#define PRINT_IDX 1

#define pi  3.14159265358979323846f
#define sqr(x) ((x)*(x))  
#define sqterm(a,b) sqr(((double)(a))-((double)(b)))

#include <vector>


#define MAX_ROWS 4096
#define MAX_COLS 4096

extern int COLS; 
extern int ROWS; 
extern int NUMBLK; 
/*
#define num_blk 4096
#define BPPSCALE 5000
#define ESCALE 500
#define SER_LIM 8
#define ZERO 0
*/

#define EMP_MEAN_SUB 0

#define PRECISION 16

#define MAXIMUM(a,b) ((a>b)?(a):(b))
#define MINIMUM(a,b) ((a<b)?(a):(b))

#define SIGN(x) ((x>=0)?1:-1)

#define PERCENT 0.1
#define Lambdadelta 0.1

#define ACDEBUG 0
#define DEC_DEBUG 0
#define	ENC_DEBUG 0

#define UNIFORM_DIST 1

#define SUBMEAN 0
//#define TRUNLAP 1

//#define ImageDependet 1

//#define CONSCALE 0.9

//#define GROUPCOEFF 0

//#define LNUMGROUP 4
//#define RNUMGROUP 3

//#define UPDATEDIST 0

//#define SORTLTABLE 0

#define INFINITY 1000000000.0

#define ORIGLTABLE 0

//#define RLCODE 1
//#define LYCODE 0

//#define SCLPAR 0

#define DOUBLEL 0

typedef enum { FALSE, TRUE } boolean;

//#define SCALELAMBDA 4
//extern double dTestDist;
extern int rgiAlpSizeDC[4];
extern int EXT_SCH_NO;
extern std::vector<std::vector<int>> drgiQIdxSDQ;
extern std::vector<std::vector<int>> drgiQIdxHDQ;
extern double TAR_DIST;
extern int TAR_PSNR;
extern double SDQ_LAMBDA;
extern int SDQ, ADPT_DIST_PROFL;
extern char rgcIFName[256];
extern int drgnTrctLv[8][8];
extern double drgdEmpMean[8][8];
extern double drgdEmpScl[8][8];
extern int rgiRiceGolombOrder[6];
extern double dDCqstep;
extern int drgiBlk[8][8];
//extern int drgiDiffDC[ROWS / 8][COLS / 8];

extern int ZigZag[8][8];
extern double drgdZeroPercent[8][8];
extern double dMaxDC;
extern double dMinDC;
extern double dLScale;
extern double dScaleLambda;
//extern int EXTF[ROWS / 8][COLS / 8]; 
extern char rgcOFNameTQtbl[256];
extern char rgcOFNameAdpZZ[256];
extern char rgcIFNameLtbl[256];
extern char rgcIFNameQtbl[256];
extern double dAvgCostRD;
extern double drgdSumYmulC[8][8];
extern double drgdSumCmulC[8][8];
extern double drgdSumC[8][8];
extern double drgdSumY[8][8];
extern double dDistCurIterSDQ;
extern int iNumIter;
extern double dCostHDQ;
extern int TAR_PSNR_ARG;

extern int drgiLTblSDQ[8][8];
extern int drgiNumBlkClgt0[8][8];

extern int iOptMean;   // 0: Substract EmpMean 1: Substract Mean of this image
extern int iOptAdpQ;   // 0: Use Lambda Yc Aval Bval from 30 images 1: Use Lambda Yc Aval Bval from this image
extern int iOptSortL;  // 0: Use Zigzag Scan 1: Use adpative scan
extern int iOptAdpL;   // 0: Use default L table; 1: Use L calculated from Yc and Lambda
extern int iOptLYC;    // 0: Use Run Level Coding 1: Use Layer-based Coding
extern int iOptLapUnf; // 0: Use Truncated Laplacian Only 1: Use Truncated Laplacian with Uniform Dist.
extern int iOptLamScl;  // 0: Use tuned Scale 1: Use actual Lambda
extern int iRoseRD;

double PSNR2MSE(double psnr);
double MSE2PSNR(double mse);
int fnRound(double dRoundIn);

#endif