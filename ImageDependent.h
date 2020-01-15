/*
*  ImageDependent.h
*  mjpeg
*
*  Created by Jin Meng on 11-07-04.
*  Copyright 2011 __MyCompanyName__. All rights reserved.
*
*/
#pragma once
#ifndef _IMAGEDEPENDENT_H_
#define _IMAGEDEPENDENT_H_

struct YcNode_Alias {
	double m_dYc;
	double m_dLambda;
	double m_dYCum;
	int m_nY;
	double m_dLLD;

	struct YcNode_Alias* m_YNNext;
};

typedef struct YcNode_Alias YcNode;

typedef struct {
	double m_dYcBdy;
	double m_dYCum;
	int m_nY;
	double m_dMLD;

	YcNode* m_YNHead;
} YcItvl;

typedef struct {
	int m_iVal;
	int m_iIdx;
} IntScalar;

extern double drgdLambda[8][8], drgdGamma[8][8];
// the value for this will be found in integerTCM which will be called in encoder
extern double drgdYc[8][8];
extern double drgdBval[8][8];
extern double drgdAval[8][8];
extern double drgdMean[8][8];

double fnCalLambda(double dYc, double dCen, double dPcn);
void fnIntlYcItvl(YcItvl* YIpInst);
void fnIntlYcNode(YcNode* YNpInst, double dYc);
void fnFreeYcItvl(YcItvl* YIpInst);
void fnAddYcToItvl(YcItvl* YIpInst, double dYc);
YcNode* fnOptimalYc(YcItvl* YIpInst, double dTlYSum, double dAval, int iTlY);
void fnDistrParaTCM(double **drgdCoeff);
int fnGenRandom(int iLow, int iHigh);
int fnRandomPartition(double *rgdList, int iLow, int iHigh);
double fnRandomSelect(double *rgdList, int iLow, int iHigh, int i);
int fnCmpIntScalar(const void *a, const void *b);

#endif