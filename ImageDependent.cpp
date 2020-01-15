/*
*  ImageDependent.c
*  mjpeg
*
*  Created by Jin Meng on 11-07-04.
*  Copyright 2011 __MyCompanyName__. All rights reserved.
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include "global.h"
#include "ImageDependent.h"
#include <iostream>

/*
*  Variable: drgnTrctLv
*  Number of truncated levels of non-uniform quantizer
*  for each coefficient in 8x8 block
*/
int drgnTrctLv[8][8] =
{
	{ 671,127,63,31,15,7,7,5 },
{ 127, 63,31,15, 7,7,5,5 },
{ 63, 31,15, 7, 7,5,5,5 },
{ 31, 15, 7, 7, 5,5,5,5 },
{ 15,  7, 7, 5, 5,5,5,5 },
{ 7,  7, 5, 5, 5,5,5,5 },
{ 7,  5, 5, 5, 5,5,5,5 },
{ 5,  5, 5, 5, 5,5,5,5 }
};
/*
{
{405, 133, 111,  85,  93,  71,  69,  47},
{171, 137,  99,  67,  65,  43,  51,  35},
{139, 113,  63,  65,  55,  47,  41,  39},
{107,  55,  67,  77,  53,  37,  43,  33},
{ 93,  43,  45,  55,  51,  55,  35,  35},
{ 69,  37,  35,  35,  57,  63,  45,  35},
{ 47,  31,  31,  23,  33,  47,  51,  59},
{ 53,  35,  31,  25,  29,  43,  41,  35}
};
*/
double drgdBval[8][8] =
{
	{ 0.9872, 0.9780, 0.9861, 0.9869, 0.9910, 0.9908, 0.9934, 0.9925 },
{ 0.9879, 0.9942, 0.9925, 0.9910, 0.9922, 0.9899, 0.9946, 0.9943 },
{ 0.9905, 0.9940, 0.9893, 0.9920, 0.9908, 0.9920, 0.9934, 0.9955 },
{ 0.9897, 0.9871, 0.9920, 0.9942, 0.9930, 0.9919, 0.9953, 0.9962 },
{ 0.9905, 0.9866, 0.9892, 0.9931, 0.9941, 0.9963, 0.9961, 0.9976 },
{ 0.9896, 0.9881, 0.9889, 0.9922, 0.9964, 0.9977, 0.9981, 0.9985 },
{ 0.9867, 0.9889, 0.9899, 0.9896, 0.9953, 0.9980, 0.9987, 0.9993 },
{ 0.9904, 0.9900, 0.9907, 0.9922, 0.9956, 0.9983, 0.9989, 0.9984 }
};

double drgdLambda[8][8] =
{
	{ 394.27, 38.20, 21.76, 14.42, 10.86, 7.97, 6.14, 5.12 },
{ 41.07, 22.67, 15.58, 11.26, 8.68, 6.51, 5.24, 4.51 },
{ 22.96, 15.84, 11.90, 9.44, 7.43, 5.86, 4.72, 4.17 },
{ 15.07, 11.08, 9.41, 7.77, 6.30, 5.02, 4.31, 3.84 },
{ 10.94, 8.38, 7.27, 6.24, 5.22, 4.46, 3.83, 3.49 },
{ 8.04, 6.49, 5.67, 4.97, 4.41, 3.90, 3.66, 3.37 },
{ 6.05, 5.20, 4.66, 4.11, 3.88, 3.72, 3.65, 3.32 },
{ 5.39, 4.61, 4.18, 3.85, 3.61, 3.61, 3.43, 3.51 }
};

double drgdAval[8][8] =
{
	{ 1156.59, 669.10, 585.24, 523.29, 688.65, 564.65, 601.89, 325.93 },
{ 680.20, 477.25, 434.86, 310.65, 345.74, 245.71, 293.30, 158.47 },
{ 671.66, 464.78, 286.27, 315.93, 323.91, 275.44, 243.87, 186.96 },
{ 674.45, 275.95, 334.84, 471.16, 323.12, 202.94, 226.70, 132.31 },
{ 693.93, 228.19, 251.13, 351.99, 322.16, 323.34, 151.79, 109.90 },
{ 597.05, 202.05, 192.34, 188.80, 349.68, 390.59, 159.97, 87.80 },
{ 416.52, 161.08, 171.56, 111.62, 156.26, 168.01, 149.02, 149.13 },
{ 527.61, 228.75, 192.15, 117.93, 112.87, 138.09, 96.21, 88.53 }
};

double drgdYc[8][8] =
{
	{ 883.72, 240.38, 158.31, 111.13, 94.65, 70.33, 58.35, 45.54 },
{ 277.26, 175.72, 123.16, 86.90, 72.18, 52.03, 47.53, 38.09 },
{ 177.82, 129.69, 87.75, 76.23, 61.41, 49.68, 41.47, 37.33 },
{ 123.28, 80.25, 76.68, 70.56, 54.88, 41.62, 39.39, 33.90 },
{ 94.68, 61.37, 56.93, 55.13, 47.54, 43.39, 34.34, 31.86 },
{ 70.32, 49.42, 44.06, 41.07, 43.49, 41.30, 35.60, 31.38 },
{ 50.97, 39.93, 37.16, 31.12, 34.18, 36.34, 36.65, 35.55 },
{ 49.18, 38.40, 34.78, 30.72, 30.88, 35.06, 33.16, 32.40 }
};

double drgdMean[8][8] = { 0.0 };

int fnCmpIntScalar(const void *a, const void *b)
{
	if (((IntScalar *)a)->m_iVal < ((IntScalar *)b)->m_iVal)
	{
		return 1;
	}
	else if (((IntScalar *)a)->m_iVal >((IntScalar *)b)->m_iVal)
	{
		return -1;
	}
	else
	{
		return (((IntScalar *)a)->m_iIdx > ((IntScalar *)b)->m_iIdx) ? 1 : -1;
	}
}

/*
*  function: fnGenRandom
*  Generate a random number within a range
*/
int fnGenRandom(int iLow, int iHigh)
{
	return (iLow + (int)((iHigh - iLow + 1)*rand() / (double)(RAND_MAX + 1.0)));
}

/*
*  function: fnRandomPartition
*  Randomly select a value x and partition the list into two parts:
*  >=x and <x
*/
int fnRandomPartition(double *rgdList, int iLow, int iHigh)
{
	int iRIdx;
	int iCrIdx;
	int j;
	double x, temp;

	if ((iLow>iHigh) || (iLow<0) || (iHigh<0))
	{
		printf("error: out of range!\n");
		exit(0);
	}

	if (iLow == iHigh)
	{
		return iLow;
	}

	iRIdx = fnGenRandom(iLow, iHigh);
	x = rgdList[iRIdx];

	temp = rgdList[iRIdx];
	rgdList[iRIdx] = rgdList[iHigh];
	rgdList[iHigh] = temp;

	iCrIdx = iLow - 1;

	for (j = iLow; j <= iHigh - 1; j++)
	{
		if (rgdList[j] >= x)
		{
			iCrIdx++;
			temp = rgdList[j];
			rgdList[j] = rgdList[iCrIdx];
			rgdList[iCrIdx] = temp;
		}
	}

	rgdList[iHigh] = rgdList[iCrIdx + 1];
	rgdList[iCrIdx + 1] = x;

	return (iCrIdx + 1);
}

/*
*  function: fnRandomSelect
*  Choose a i-th largest element from a list
*/
double fnRandomSelect(double *rgdList, int iLow, int iHigh, int i)
{
	int q, k;

	if ((iLow>iHigh) || (iLow<0) || (iHigh<0))
	{
		printf("error: out of range!\n");
		exit(0);
	}

	if (iLow == iHigh)
	{
		return rgdList[iLow];
	}

	q = fnRandomPartition(rgdList, iLow, iHigh);
	k = q - iLow + 1;

	if (i == k)
	{
		return rgdList[q];
	}
	else if (i < k)
		return fnRandomSelect(rgdList, iLow, q - 1, i);
	else
		return fnRandomSelect(rgdList, q + 1, iHigh, i - k);
}

/*
*  function: fnCalLambda
*  Calculate Lambda by fix-point method
*/
double fnCalLambda(double dYc, double dCen, double dPcn)
{
	double doLambda = dCen;
	double dnLambda = dCen + dYc * exp(-dYc / doLambda) / (1 - exp(-dYc / doLambda));

	assert(dPcn>0);

	while (fabs(dnLambda - doLambda)>dPcn)
	{
		doLambda = dnLambda;
		dnLambda = dCen + dYc * exp(-dYc / dnLambda) / (1 - exp(-dYc / dnLambda));
	}

	return dnLambda;
}

/*
*  function: fnIntlYcItvl, fnIntlYcNode
*  Initialize instants of structures YcItvl and YcNode
*/
void fnIntlYcItvl(YcItvl* YIpInst)
{
	YIpInst->m_dYCum = 0.0;
	YIpInst->m_nY = 0;
	YIpInst->m_dMLD = 0.0;
	YIpInst->m_YNHead = NULL;
}

void fnIntlYcNode(YcNode* YNpInst, double dYc)
{
	YNpInst->m_dYc = dYc;
	YNpInst->m_dLLD = 0.0;
	YNpInst->m_dLambda = 0.0;
	YNpInst->m_nY = 1;
	YNpInst->m_dYCum = 0.0;
	YNpInst->m_YNNext = NULL;
}

/*
*  function: fnFreeYcItvl
*  Free the space of an instant of YcItvl
*/
void fnFreeYcItvl(YcItvl* YIpInst)
{
	YcNode *YNpCr, *YNpNx;

	YNpCr = YIpInst->m_YNHead;

	while (YNpCr != NULL)
	{
		YNpNx = YNpCr->m_YNNext;
		free(YNpCr);
		YNpCr = YNpNx;
	}

	YIpInst->m_YNHead = NULL;
}

/*
*  function: fnAddYcToItvl
*  Add a Yc into YcItvl
*/
void fnAddYcToItvl(YcItvl* YIpInst, double dYc)
{
	YcNode *YNpCr, *YNpPv;
	YcNode *YNpNew;

	if (YIpInst->m_YNHead == NULL)
	{
		YNpNew = (YcNode*)malloc(sizeof(YcNode));
		fnIntlYcNode(YNpNew, dYc);

		YIpInst->m_YNHead = YNpNew;
	}
	else
	{
		YNpCr = YIpInst->m_YNHead;
		if (YNpCr->m_dYc<dYc)
		{
			YNpNew = (YcNode*)malloc(sizeof(YcNode));
			fnIntlYcNode(YNpNew, dYc);

			YIpInst->m_YNHead = YNpNew;
			YNpNew->m_YNNext = YNpCr;
		}
		else
		{
			YNpPv = YNpCr;
			YNpCr = YNpCr->m_YNNext;

			while (YNpCr != NULL && YNpCr->m_dYc >= dYc)
			{
				YNpPv = YNpCr;
				YNpCr = YNpCr->m_YNNext;
			}

			if (YNpCr != NULL && YNpCr->m_dYc == dYc)
			{
				(YNpCr->m_nY)++;
			}
			else
			{
				YNpNew = (YcNode*)malloc(sizeof(YcNode));
				fnIntlYcNode(YNpNew, dYc);

				YNpPv->m_YNNext = YNpNew;
				YNpNew->m_YNNext = YNpCr;
			}
		}
	}

	YIpInst->m_dYCum += dYc;
	(YIpInst->m_nY)++;
}

/*
*  function: fnOptimalYc
*  Calculate Optimal Yc within a Yc interval
*/
YcNode* fnOptimalYc(YcItvl* YIpInst, double dTlYSum, double dAval, int iTlY)
{
	int nPvY;
	double dPvYSum;
	YcNode *YNCr, *YNOp;

	if (YIpInst->m_YNHead == NULL)
	{
		return NULL;
	}

	nPvY = YIpInst->m_nY;
	dPvYSum = YIpInst->m_dYCum;
	YNCr = YIpInst->m_YNHead;
	YNOp = YNCr;

	while (YNCr != NULL)
	{
		YNCr->m_dLambda = fnCalLambda(YNCr->m_dYc, (dTlYSum - dPvYSum) / (double)(iTlY - nPvY), Lambdadelta);

		if (nPvY)
		{
			YNCr->m_dLLD = nPvY * (log((double)(nPvY) / (double)(iTlY)) - log(2 * (dAval - YNCr->m_dYc)))
				+ (iTlY - nPvY)*(log(1 - (double)(nPvY) / (double)(iTlY)) - log(2 * YNCr->m_dLambda) - log(1 - exp(-YNCr->m_dYc / YNCr->m_dLambda)))
				- (dTlYSum - dPvYSum) / YNCr->m_dLambda;
		}
		else
		{
			YNCr->m_dLLD = iTlY * (-log(2 * YNCr->m_dLambda) - log(1 - exp(-YNCr->m_dYc / YNCr->m_dLambda))) - (dTlYSum - dPvYSum) / YNCr->m_dLambda;
		}

		dPvYSum += (double)(YNCr->m_nY)*YNCr->m_dYc;
		nPvY += YNCr->m_nY;

		YNCr->m_dYCum = dPvYSum - (double)(YNCr->m_nY)*YNCr->m_dYc;
		YNCr->m_nY = nPvY - YNCr->m_nY;

		if (YNCr->m_dLLD>YNOp->m_dLLD)
		{
			YNOp = YNCr;
		}

		YNCr = YNCr->m_YNNext;
	}

	return YNOp;
}

/*
*  function: fnDistrParaTCM
*  Update quantization parameters according to current image
*/
void fnDistrParaTCM(double** drgdCoeff)
{
	//std::cout << "tcm func " << std::endl; 
	int u, v, x, y, i, j;
	int iIdx, nPvY;
	double dPvYCum;
	double drgdYSum[8][8];
	double* rgdList = new double [NUMBLK];
	YcItvl *YIpQYc[8][8];
	int drgnItvl[8][8];
	YcNode *YNpCr{nullptr}, *YNpOp{ nullptr };

#if 0
	double dminLambda;
	double dmaxLambda;
	double dQstepLambda;
	double dminRatio;
	double dmaxRatio;
	double dQstepRatio;
	int drgiGroupMap[8][8];
	int trgiGroup[LNUMGROUP*RNUMGROUP][63][2];
	int rgiGroupSize[LNUMGROUP*RNUMGROUP];

	double *dpGList;
	int rgnGItvl[LNUMGROUP*RNUMGROUP];
	YcItvl *YIpGQYc[LNUMGROUP*RNUMGROUP];

	double dAval;
	double dYSum;
	double dYc;
	double dLambda;
#endif


	IntScalar rgISZZL[64];
	int iMaxL;
	int rgiPermute[64];


#if ZEROFORCE
	double dPerUnit;
#endif

	for (i = 0; i<8; i++)//freq 
	{
		for (j = 0; j<8; j++)
		{
			if ((i + j) == 0)
			{
				continue;
			}

			drgdAval[i][j] = 0.0;
			drgdYSum[i][j] = 0.0;

			y = 0;
			for (u = 0; u<ROWS / 8; u++)//block
			{
				for (v = 0; v<COLS / 8; v++)
				{
					rgdList[y] = fabs(drgdCoeff[8 * u + i][8 * v + j]);
					y++;
					drgdAval[i][j] = MAXIMUM(drgdAval[i][j], fabs(drgdCoeff[8 * u + i][8 * v + j]));
					drgdYSum[i][j] += fabs(drgdCoeff[8 * u + i][8 * v + j]);
				}
			}

			drgnItvl[i][j] = (int)(drgdAval[i][j] - fnRandomSelect(rgdList, 0, NUMBLK - 1, (int)(NUMBLK*PERCENT)));
			if (drgnItvl[i][j] == 0)
			{
				drgnItvl[i][j] = 1;
			}
			YIpQYc[i][j] = (YcItvl*)malloc(drgnItvl[i][j] * sizeof(YcItvl));

			for (x = 0; x<drgnItvl[i][j]; x++)
			{
				fnIntlYcItvl(&(YIpQYc[i][j][x]));
			}

			for (u = 0; u<ROWS / 8; u++)
			{
				for (v = 0; v<COLS / 8; v++)
				{
					iIdx = (int)(drgdAval[i][j] - fabs(drgdCoeff[8 * u + i][8 * v + j]));
					if (iIdx<drgnItvl[i][j])
					{
						fnAddYcToItvl(&(YIpQYc[i][j][iIdx]), fabs(drgdCoeff[8 * u + i][8 * v + j]));
					}
				}
			}

			//printf("AddYc,%d,%d\n",i,j);

			dPvYCum = 0.0;
			nPvY = 0;
			for (x = 0; x<drgnItvl[i][j]; x++)
			{
				dPvYCum += YIpQYc[i][j][x].m_dYCum;
				nPvY += YIpQYc[i][j][x].m_nY;

				YIpQYc[i][j][x].m_dYCum = dPvYCum - YIpQYc[i][j][x].m_dYCum;
				YIpQYc[i][j][x].m_nY = nPvY - YIpQYc[i][j][x].m_nY;
			}

			for (x = 0; x<drgnItvl[i][j]; x++)
			{
				YNpCr = fnOptimalYc(&(YIpQYc[i][j][x]), drgdYSum[i][j], drgdAval[i][j], NUMBLK);

				if ((x == 0) || ((YNpCr != NULL) && (YNpCr->m_dLLD>YNpOp->m_dLLD)))
				{
					YNpOp = YNpCr;
				}
			}

			//printf("OptimalYc,%d,%d\n",i,j);

			drgdYc[i][j] = YNpOp->m_dYc;
			drgdLambda[i][j] = YNpOp->m_dLambda;
			drgdBval[i][j] = 1.0 - (double)(YNpOp->m_nY) / (double)(NUMBLK);

			if (iOptAdpL)
			{
				drgnTrctLv[i][j] = (int)(dLScale*drgdLambda[i][j] * sqrt(pow(1 - exp(-drgdYc[i][j] / drgdLambda[i][j] / 3), 3.0) / (1 - exp(-drgdYc[i][j] / drgdLambda[i][j]))));
				//			drgnTrctLv[i][j] = (int)(CONSCALE*drgdLambda[i][j]);
				if (drgnTrctLv[i][j] % 2 == 0)
				{
					(drgnTrctLv[i][j])++;
				}
			}

			for (x = 0; x<drgnItvl[i][j]; x++)
			{
				fnFreeYcItvl(&(YIpQYc[i][j][x]));
			}

			free(YIpQYc[i][j]);

			//printf("Free,%d,%d\n",i,j);
		}
	}

#if 0
	dminLambda = INFINITY;
	dmaxLambda = 0.0;
	dminRatio = INFINITY;
	dmaxRatio = 0.0;
	for (i = 0; i<8; i++)
	{
		for (j = 0; j<8; j++)
		{
			if ((i + j) == 0)
			{
				continue;
			}

			dminLambda = MINIMUM(dminLambda, drgdLambda[i][j]);
			dmaxLambda = MAXIMUM(dmaxLambda, drgdLambda[i][j]);
			dminRatio = MINIMUM(dminRatio, drgdYc[i][j] / drgdLambda[i][j]);
			dmaxRatio = MAXIMUM(dmaxRatio, drgdYc[i][j] / drgdLambda[i][j]);
		}
	}

	dQstepLambda = (dmaxLambda - dminLambda) / LNUMGROUP;
	dQstepRatio = (dmaxRatio - dminRatio) / RNUMGROUP;
	for (i = 0; i<LNUMGROUP*RNUMGROUP; i++)
	{
		rgiGroupSize[i] = 0;
	}
	for (i = 0; i<8; i++)
	{
		for (j = 0; j<8; j++)
		{
			if ((i + j) == 0)
			{
				drgiGroupMap[i][j] = -1;
			}
			else
			{
				x = MINIMUM((LNUMGROUP - 1), ((int)((drgdLambda[i][j] - dminLambda) / dQstepLambda)));
				y = MINIMUM((RNUMGROUP - 1), ((int)((drgdYc[i][j] / drgdLambda[i][j] - dminRatio) / dQstepRatio)));
				drgiGroupMap[i][j] = y * LNUMGROUP + x;
				trgiGroup[y*LNUMGROUP + x][rgiGroupSize[y*LNUMGROUP + x]][0] = i;
				trgiGroup[y*LNUMGROUP + x][rgiGroupSize[y*LNUMGROUP + x]][1] = j;
				(rgiGroupSize[y*LNUMGROUP + x])++;
			}
		}
	}
	for (x = 0; x<LNUMGROUP*RNUMGROUP; x++)
	{
		if (rgiGroupSize[x])
		{
			dpGList = (double *)malloc(rgiGroupSize[x] * NUMBLK * sizeof(double));
			dAval = 0.0;
			dYSum = 0.0;

			for (y = 0; y<rgiGroupSize[x]; y++)
			{
				i = trgiGroup[x][y][0];
				j = trgiGroup[x][y][1];
				dAval = MAXIMUM(dAval, drgdAval[i][j]);
				dYSum += drgdYSum[i][j];

				for (u = 0; u<ROWS / 8; u++)
				{
					for (v = 0; v<COLS / 8; v++)
					{
						dpGList[y*NUMBLK + u * (ROWS / 8) + v] = fabs(drgdCoeff[8 * u + i][8 * v + j]);
					}
				}
			}

			rgnGItvl[x] = (int)(dAval - fnRandomSelect(dpGList, 0, rgiGroupSize[x] * NUMBLK - 1, (int)(rgiGroupSize[x] * NUMBLK*PERCENT)));
			if (rgnGItvl[x] == 0)
			{
				rgnGItvl[x] = 1;
			}
			YIpGQYc[x] = (YcItvl*)malloc(rgnGItvl[x] * sizeof(YcItvl));

			for (y = 0; y<rgnGItvl[x]; y++)
			{
				fnIntlYcItvl(&(YIpGQYc[x][y]));
			}

			for (y = 0; y<rgiGroupSize[x] * NUMBLK; y++)
			{
				iIdx = (int)(dAval - dpGList[y]);
				if (iIdx<rgnGItvl[x])
				{
					fnAddYcToItvl(&(YIpGQYc[x][iIdx]), dpGList[y]);
				}
			}

			//printf("AddYc,%d,%d\n",i,j);

			dPvYCum = 0.0;
			nPvY = 0;
			for (y = 0; y<rgnGItvl[x]; y++)
			{
				dPvYCum += YIpGQYc[x][y].m_dYCum;
				nPvY += YIpGQYc[x][y].m_nY;

				YIpGQYc[x][y].m_dYCum = dPvYCum - YIpGQYc[x][y].m_dYCum;
				YIpGQYc[x][y].m_nY = nPvY - YIpGQYc[x][y].m_nY;
			}

			for (y = 0; y<rgnGItvl[x]; y++)
			{
				YNpCr = fnOptimalYc(&(YIpGQYc[x][y]), dYSum, dAval, rgiGroupSize[x] * NUMBLK);

				if ((y == 0) || ((YNpCr != NULL) && (YNpCr->m_dLLD>YNpOp->m_dLLD)))
				{
					YNpOp = YNpCr;
				}
			}

			//printf("OptimalYc,%d,%d\n",i,j);

			for (y = 0; y<rgiGroupSize[x]; y++)
			{
				i = trgiGroup[x][y][0];
				j = trgiGroup[x][y][1];
#if UPDATEDIST
				drgdAval[i][j] = dAval;
				drgdYc[i][j] = YNpOp->m_dYc;
				drgdLambda[i][j] = YNpOp->m_dLambda;
				drgdBval[i][j] = 1.0 - (double)(YNpOp->m_nY) / (double)(rgiGroupSize[x] * NUMBLK);
#endif
				dYc = YNpOp->m_dYc;
				dLambda = YNpOp->m_dLambda;
				drgnTrctLv[i][j] = (int)(dLScale*dLambda*sqrt(pow(1 - exp(-dYc / dLambda / 3), 3.0) / (1 - exp(-dYc / dLambda))));
				if (drgnTrctLv[i][j] % 2 == 0)
				{
					(drgnTrctLv[i][j])++;
				}
			}

			for (y = 0; y<rgnGItvl[x]; y++)
			{
				fnFreeYcItvl(&(YIpGQYc[x][y]));
			}

			free(YIpGQYc[x]);
		}
	}

	for (i = 0; i<8; i++)
	{
		for (j = 0; j<8; j++)
		{
			printf("%2d ", drgiGroupMap[i][j]);
		}
		printf("\n");
	}

	printf("-----------------------------\n");
#endif


	if (iOptSortL)
	{
		iMaxL = 0;
		for (i = 0; i<8; i++)
		{
			for (j = 0; j<8; j++)
			{
				rgISZZL[ZigZag[i][j]].m_iVal = drgnTrctLv[i][j];
				rgISZZL[ZigZag[i][j]].m_iIdx = ZigZag[i][j];

				iMaxL = MAXIMUM(iMaxL, drgnTrctLv[i][j]);
			}
		}

		rgISZZL[0].m_iVal = iMaxL + 1;

		qsort(rgISZZL, 64, sizeof(IntScalar), fnCmpIntScalar);

		for (u = 0; u<64; u++)
		{
			rgiPermute[rgISZZL[u].m_iIdx] = u;
		}

		for (i = 0; i<8; i++)
		{
			for (j = 0; j<8; j++)
			{
				ZigZag[i][j] = rgiPermute[ZigZag[i][j]];
			}
		}

		/*for (i = 0; i<8; i++)
		{
			for (j = 0; j<8; j++)
			{
				printf("%2d ", ZigZag[i][j]);
			}
			printf("\n");
		}
		printf("-----------------------------\n");*/
	}


#if ZEROFORCE
	dPerUnit = 0.0;
	for (i = 0; i<8; i++)
	{
		for (j = 0; j<8; j++)
		{
			if ((i + j) == 0)
			{
				continue;
			}

			//			dPerUnit += 1.0/drgdLambda[i][j]/drgdLambda[i][j];
			//			dPerUnit += 1.0/(double)(drgnTrctLv[i][j]);
			dPerUnit += 1.0 / pow((double)(drgnTrctLv[i][j]), 1.0 / 2.0);
		}
	}

	dPerUnit = 63.0*ZEROPERCENT / dPerUnit;
	for (i = 0; i<8; i++)
	{
		for (j = 0; j<8; j++)
		{
			if ((i + j) == 0)
			{
				drgdZeroPercent[i][j] = 0.0;
				continue;
			}

			//			drgdZeroPercent[i][j] = dPerUnit/drgdLambda[i][j]/drgdLambda[i][j]; 
			//			drgdZeroPercent[i][j] = dPerUnit/(double)(drgnTrctLv[i][j]);
			drgdZeroPercent[i][j] = dPerUnit / pow((double)(drgnTrctLv[i][j]), 1.0 / 2.0);
		}
	}

	for (i = 0; i<8; i++)
	{
		for (j = 0; j<8; j++)
		{
			printf("%.4f ", drgdZeroPercent[i][j]);
		}
		printf("\n");
	}
	printf("-----------------------------\n");
#endif

	if (!iRoseRD)
	{
		/*printf("-----------------------------\n");
		for (i = 0; i<8; i++)
		{
			for (j = 0; j<8; j++)
			{
				printf("%3d ", drgnTrctLv[i][j]);
			}
			printf("\n");
		}*/
	}
	/*printf("now you see me");
	printf("---------------Yc--------------\n");
	for (i = 0; i<8; i++)
	{
		for (j = 0; j<8; j++)
		{
			printf("%6.2f ", drgdYc[i][j]);
		}
		printf("\n");
	}
*/
	/*printf("---------------A--------------\n");

	for (i = 0; i<8; i++)
	{
		for (j = 0; j<8; j++)
		{
			printf("%7.2f ", drgdAval[i][j]);
		}
		printf("\n");
	}
*/
	//printf("---------------Lambda--------------\n");

	//for (i = 0; i<8; i++)
	//{
	//	for (j = 0; j<8; j++)
	//	{
	//		printf("%6.2f ", drgdLambda[i][j]);
	//	}
	//	printf("\n");
	//}

	//printf("-------------B----------------\n");

	//for (i = 0; i<8; i++)
	//{
	//	for (j = 0; j<8; j++)
	//	{
	//		printf("%.6f ", drgdBval[i][j]);
	//	}
	//	printf("\n");
	//}
	delete[] rgdList;

}
