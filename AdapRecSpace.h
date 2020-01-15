/*
*  AdapRecSpace.h
*  mjpeg
*
*  Created by Chang Sun on 11-07-20.
*  Copyright 2011 Univ. of Waterloo All rights reserved.
*
*/
#pragma once
#define NUM_QTZ_GAMMA 64
#define NUM_QTZ_LAMBDA 256

#define MdlQ 1;//0: use real data to calculate distortion;1: use the TCM model to estimate distortion
#define MdlDZoneQ 1//1:use the TCM model with the deadzone to estimate distortion
#define MdlDZoneQForce 1//1:for subbands whose T=1, force L'=1 by using MdlDZone scheme
#define AdptQprime 0//1:for subbands whose T=1, allow diff Qprime for diff subband to achieve the defualt distortion using real data calculation
#define ExtRec 0//for subbands whose L>0 - 0: extend intervals for the uniform part; 1:extend reconstruction point for the uniform part
#define REC_DATA 1
#define EXT_Q_OFFSET 0
//#define TAR_DIST 10.0

/* Target distortion levels: 30 32 34 36 38 40 42 in PSNR*/
#define TAR_DIST_1 65.0250
#define TAR_DIST_2 41.0280
#define TAR_DIST_3 25.8869
#define TAR_DIST_4 16.3335 
#define TAR_DIST_5 10.3058
#define TAR_DIST_6  6.5025
#define TAR_DIST_7  4.1028

//#define sqr(x) ((x)*(x))  
#define cube(x) ((x)*(x)*(x)) 

typedef struct {
	int Gamma_idx[8][8]; // inedex of quantized Gamma for each freq 
	int Lambda_idx[8][8]; // inedex of quantized Lambda for each freq 
	int LTbl[8][8]; // for each freq you shoudl know number of levels of quantization L_k 
	double QTbl[8][8]; // q_k step size of quantization for each freq
	double UTbl[8][8]; //  u_k deadzone width of each freq 
	double DeltaTbl[8][8];
}OptDzoneTbl;

extern int iLLowDC, iLHighDC;
extern int irgdLowRSDQFlag[8][8];
extern int drgiQTbl[8][8];
extern int drgiLTbl[8][8]; // num of levels of inlier 
extern int drgiTTbl[8][8]; // T parameter : 1 skipp that freq, 0 : quantize that freq
extern int drgiLPrimeTbl[8][8]; // num of level outlier 
extern char RecordrgcIFName[256];
//extern int drgiExtF[ROWS / 8][COLS / 8], 
extern std::vector<std::vector<int>> drgiExtF;
extern std::vector<std::vector<int>> drgiExtFDec;
//extern int drgiExtFDec[ROWS / 8][COLS / 8];

//extern int drgiCBP[ROWS / 8][COLS / 8];
extern std::vector<std::vector<int>>drgiCBP;

extern int iMaxLprime;
extern int iMaxAbsDC;
extern int iMaxLDC;
extern double **Q_cache;
extern double dAdOnVrnceSDQ;
extern int iNumAdOnVrnceSDQ;
extern double drgdRealRate[8][8];
extern double drgdRealDist[8][8];
extern int drgiNumBlkWithinYc[8][8];
extern OptDzoneTbl OptDzoneQtzr;

void fnParaTCMQtz(int iRow, int iCol);
void fnOptMdlQ(double **drgdCoef);
int fnHDQ(double dCoef, int iBlkRow, int iBlkCol, int iFreqRow, int iFreqCol);
double fnDeQtz(int iLv, int iRowIdx, int iColIdx, int iRow, int iCol);
double fnDeQtzExt(int iLv, int iFreqRow, int iFreqCol);
double fnDeQtzSig(int iLv, int iFreqRow, int iFreqCol);
double fnCalDistDeriv(int iFreqRow, int iFreqCol, double dCurQ);
double fnCalDist(int iFreqRow, int iFreqCol, double dCurQ);
void fnUpdateQtbl();
void fnOptDistProfl(); 
