/*
*  Context.c
*  mjpeg
*
*  Created by Jin Meng on 11-06-08.
*  Copyright 2011 __MyCompanyName__. All rights reserved.
*
*/
#pragma warning(disable:4996)
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Context.h"
#include "global.h"
#include "HACACM.h"

#define F_BITS 8
Context *rgctxpCoeffLayerCtx[15];
Context drgulpNZCount[7][NZCTX];
Context drgulpLZCount[7][LZCTX];
Context *fnGetCoeffCtx(int ippCoeffBlk[8][8], int iRow, int iCol, int isVal, int nNZLayer);

unsigned short int Max_frequency = ((unsigned short int) 1 << F_BITS);  // maximum frequency of c0+c1

																		/*
																		*  variable: rgnCoeffLayerCtx
																		*  array to store number of contexts for coefficients (abs val) in each layer
																		*  nz: non-zero; no: non-one;
																		*                        layer: 1 2 3 4 5 6,7 8-15nz 8no 9no 10no 11no 12no 13no 14no 15no
																		*/
static int rgnCoeffLayerCtx[15] = { 1,5,6,4,6,  3,    10,  5,  5,   5,   5,   2,   2,   2,   2 };


//static int rgnAlphabetSize[]={};

/*
*  function: fnCreateCtx
*  Create an instance of context
*/
void fnCreateIniCtxEXTF(Context *ctxpInst, int CtxDstrMatrx[2][16], int iAlpSize, int iCtxIdx)
{
	int i;

	ctxpInst->m_iAlpSize = iAlpSize;
	ctxpInst->m_lpCumSum = (long *)malloc((iAlpSize + 1) * sizeof(long));
	ctxpInst->m_lpCumSum[0] = 0;

	for (i = 1; i <= iAlpSize; i++)
	{
		ctxpInst->m_lpCumSum[i] = ctxpInst->m_lpCumSum[i - 1] + CtxDstrMatrx[i - 1][iCtxIdx];
	}
}

void fnCreateIniCtxExtImg(Context *ctxpInst, int CtxDstrMatrx[6][64], int iAlpSize, int iCtxIdx, int iCoefIdx)
{
	int i;
	int iIniDistriCtx1[2] = { 960,40 };
	int iIniDistriCtx2[2] = { 880,120 };
	int iIniDistriCtx0Coef0[2] = { 1,1000 };

	ctxpInst->m_iAlpSize = iAlpSize;
	ctxpInst->m_lpCumSum = (long *)malloc((iAlpSize + 1) * sizeof(long));
	ctxpInst->m_lpCumSum[0] = 0;

	for (i = 1; i <= iAlpSize; i++)
	{
		if (iCtxIdx == 1 && iCoefIdx == 1)
			ctxpInst->m_lpCumSum[i] = ctxpInst->m_lpCumSum[i - 1] + 1;//iIniDistriCtx1[i-1];
																	  //else if (iCtxIdx==2)
																	  //ctxpInst->m_lpCumSum[i]=ctxpInst->m_lpCumSum[i-1]+iIniDistriCtx1[i-1];
		else
		{
			if (iCtxIdx == 0 && iCoefIdx == 0)
				ctxpInst->m_lpCumSum[i] = ctxpInst->m_lpCumSum[i - 1] + iIniDistriCtx0Coef0[i - 1];
			else
				ctxpInst->m_lpCumSum[i] = ctxpInst->m_lpCumSum[i - 1] + 1;
		}

	}
}
/*
*  function: fnCreateCtx
*  Create an instance of context
*/
void fnCreateCtx(Context *ctxpInst, int iAlpSize)
{
	int i;

	ctxpInst->m_iAlpSize = iAlpSize;
	ctxpInst->m_lpCumSum = (long *)malloc((iAlpSize + 1) * sizeof(long));

	ctxpInst->m_lpCumSum[0] = 0;
	for (i = 1; i <= iAlpSize; i++)
	{
		ctxpInst->m_lpCumSum[i] = ctxpInst->m_lpCumSum[i - 1] + 1;
	}

#if BAC_EN
	if (iAlpSize == 2)
	{
		ctxpInst->m_lpCumSum[0] = 1;
		ctxpInst->m_lpCumSum[1] = 1;
	}
#endif
}

/*
*  function: fnUpdateCtx
*  iSymbol>=0: Update a context when a symbol is encoded.
*  iSymbol=-1: Normalize counts to a certain level
*/
void fnUpdateCtx(Context *ctxpInst, int iSymbol)
{
	int i;
	long lpCounts, lpTmp;

	while (ctxpInst->m_lpCumSum[ctxpInst->m_iAlpSize] >= (1 << (PRECISION - 2)))
	{
		lpTmp = 0;
		for (i = 0; i<ctxpInst->m_iAlpSize; i++)
		{
			lpCounts = ctxpInst->m_lpCumSum[i + 1] - lpTmp;
			lpTmp = ctxpInst->m_lpCumSum[i + 1];
			lpCounts >>= 1;
			ctxpInst->m_lpCumSum[i + 1] = ctxpInst->m_lpCumSum[i] + ((lpCounts == 0) ? 1 : lpCounts);
		}
	}

	if (iSymbol >= 0)
	{
		for (i = iSymbol + 1; i <= ctxpInst->m_iAlpSize; i++)
		{
			ctxpInst->m_lpCumSum[i]++;
		}
	}
}

void fnUpdateBACCtx(Context *ctxpInst, int iSymbol)
{
	ctxpInst->m_lpCumSum[iSymbol]++;

	if ((ctxpInst->m_lpCumSum[0] + ctxpInst->m_lpCumSum[1]) > Max_frequency)
	{
		ctxpInst->m_lpCumSum[0] = (ctxpInst->m_lpCumSum[0] + 1) >> 1;
		ctxpInst->m_lpCumSum[1] = (ctxpInst->m_lpCumSum[1] + 1) >> 1;
	}
}

void fnFreeCtx(Context *ctxpInst)
{
	free(ctxpInst->m_lpCumSum);
}

/*
*  funciton: fnIntlCount
*  Initialize Counts for Coefficients nz and lz of each layer
*/

void fnIntlCount(char *cpCtxFilePath, int iMaxDCVal)
{
	FILE *fpCount;
	int i, j, k;
	int iAlpSize;
	long lSymbol;

	char rgcCoeff[256];
	char rgcNZ[256];
	char rgcLZ[256];

	char CoeffLayerCtx[] = "CoeffLayerCount.txt";
	char NZCtx[] = "NZStat.txt";
	char LZCtx[] = "LZStat.txt";

	strcpy(rgcCoeff, cpCtxFilePath);
	strcat(rgcCoeff, CoeffLayerCtx);
	strcpy(rgcNZ, cpCtxFilePath);
	strcat(rgcNZ, NZCtx);
	strcpy(rgcLZ, cpCtxFilePath);
	strcat(rgcLZ, LZCtx);

	/*
	*  Read Statistics from files CoeffLayerCount NZStat and LZStat
	*/
	//fpCount = fopen(rgcCoeff, "r");
	for (i = 0; i<15; i++)
	{
		rgctxpCoeffLayerCtx[i] = (Context *)malloc(rgnCoeffLayerCtx[i] * sizeof(Context));
		for (j = 0; j<rgnCoeffLayerCtx[i]; j++)
		{
			//fscanf(fpCount, "%d:,",&iAlpSize);
			if (i == 0)
			{
				fnCreateCtx(&(rgctxpCoeffLayerCtx[i][j]), (iMaxDCVal + 1));
				/*
				for (k=0; k<(iMaxDCVal+1); k++)
				{
				#if UNIFORM_DIST
				lSymbol = 1;
				#endif
				rgctxpCoeffLayerCtx[i][j].m_lpCumSum[k+1]=rgctxpCoeffLayerCtx[i][j].m_lpCumSum[k]+lSymbol;
				}

				for (k=0; k<iAlpSize; k++)
				{
				fscanf(fpCount,"%ld,",&lSymbol);
				}*/
			}
			else
			{
				iAlpSize = (i <= 5) ? (1 << (7 - i)) : 2;
				fnCreateCtx(&(rgctxpCoeffLayerCtx[i][j]), iAlpSize);
				/*
				for (k=0; k<iAlpSize; k++)
				{
				//fscanf(fpCount,"%ld,",&lSymbol);
				#if UNIFORM_DIST
				lSymbol = 1;
				#endif
				rgctxpCoeffLayerCtx[i][j].m_lpCumSum[k+1]=rgctxpCoeffLayerCtx[i][j].m_lpCumSum[k]+lSymbol;
				}*/
			}

			/*
			*  Normalize the counts to certain levels
			*/
			fnUpdateCtx(&(rgctxpCoeffLayerCtx[i][j]), -1);
			//fscanf(fpCount, "\n");
		}
	}
	//fclose(fpCount);

	//fpCount = fopen(rgcNZ, "r");
	for (i = 0; i<7; i++)
	{
		for (j = 0; j<5; j++)
		{
			fnCreateCtx(&(drgulpNZCount[i][j]), 2);
			/*
			for (k=0; k<2; k++)
			{
			//fscanf(fpCount, "%ld,",&lSymbol);
			#if UNIFORM_DIST
			lSymbol = 1;
			#endif
			drgulpNZCount[i][j].m_lpCumSum[k+1]=drgulpNZCount[i][j].m_lpCumSum[k]+lSymbol;
			}*/
			/*
			*  Normailize the counts to certain levels
			*/
			fnUpdateCtx(&(drgulpNZCount[i][j]), -1);
			//fscanf(fpCount, "\n");			
		}
	}
	//fclose(fpCount);

	//fpCount = fopen(rgcLZ, "r");
	for (i = 0; i<7; i++)
	{
		for (j = 0; j<5; j++)
		{
			fnCreateCtx(&(drgulpLZCount[i][j]), 2);
			/*
			for (k=0; k<2; k++)
			{
			//fscanf(fpCount, "%ld,",&lSymbol);
			#if UNIFORM_DIST
			lSymbol = 1;
			#endif
			drgulpLZCount[i][j].m_lpCumSum[k+1]=drgulpLZCount[i][j].m_lpCumSum[k]+lSymbol;
			}*/
			/*
			*  Normailize the counts to certain levels
			*/
			fnUpdateCtx(&(drgulpLZCount[i][j]), -1);
			//fscanf(fpCount, "\n");			
		}
	}
	//fclose(fpCount);
}

void fnFreeCount()
{
	int i, j;
	for (i = 0; i<15; i++)
	{
		for (j = 0; j<rgnCoeffLayerCtx[i]; j++)
		{
			fnFreeCtx(&(rgctxpCoeffLayerCtx[i][j]));
		}
		free(rgctxpCoeffLayerCtx[i]);
	}

	for (i = 0; i<7; i++)
	{
		for (j = 0; j<5; j++)
		{
			fnFreeCtx(&(drgulpNZCount[i][j]));
			fnFreeCtx(&(drgulpLZCount[i][j]));
		}
	}
}

/*
*  function: fnGetCoeffCtx
*  Get the context of a coefficient at (iRow,iCol).
*  Suppose DC coefficient has been differentiated.
*/
Context *fnGetCoeffCtx(int ippCoeffBlk[8][8], int iRow, int iCol, int isVal, int nNZLayer)
{
	int iMaxVal, nNZVal;
	int i;

	if (iRow == 0 && iCol == 0) /*DC*/
	{
		return &(rgctxpCoeffLayerCtx[0][0]);
	}

	switch (iRow + iCol)
	{
	case 1:
		if (abs(ippCoeffBlk[0][0])<4)
		{
			return &(rgctxpCoeffLayerCtx[1][0]);
		}
		else
		{
			if (abs(ippCoeffBlk[0][0])<12)
			{
				return &(rgctxpCoeffLayerCtx[1][1]);
			}
			else
			{
				if (abs(ippCoeffBlk[0][0])<28)
				{
					return &(rgctxpCoeffLayerCtx[1][2]);
				}
				else
				{
					if (abs(ippCoeffBlk[0][0])<64)
					{
						return &(rgctxpCoeffLayerCtx[1][3]);
					}
					else
					{
						return &(rgctxpCoeffLayerCtx[1][4]);
					}

				}

			}
		}
		break;
	case 2:
		iMaxVal = MAXIMUM(abs(ippCoeffBlk[0][1]), abs(ippCoeffBlk[1][0]));
		iMaxVal = MINIMUM(iMaxVal, ((drgnTrctLv[0][1] - 1) / 2));
		if (iMaxVal<2)
		{
			return &(rgctxpCoeffLayerCtx[2][0]);
		}
		else
		{
			return &(rgctxpCoeffLayerCtx[2][(int)ceil(log(iMaxVal + 1) / log(2)) - 1]);
			/*
			if (iMaxVal<4)
			{
			return &(rgctxpCoeffLayerCtx[2][1]);
			}
			else
			{
			if (iMaxVal<8)
			{
			return &(rgctxpCoeffLayerCtx[2][2]);
			}
			else
			{
			if (iMaxVal<16)
			{
			return &(rgctxpCoeffLayerCtx[2][3]);
			}
			else
			{
			if (iMaxVal<32)
			{
			return &(rgctxpCoeffLayerCtx[2][4]);
			}
			else
			{
			return &(rgctxpCoeffLayerCtx[2][5]);
			}
			}
			}
			}
			*/
		}
		break;
	case 3:
		iMaxVal = MAXIMUM(abs(ippCoeffBlk[0][1]), abs(ippCoeffBlk[1][0]));
		iMaxVal = MINIMUM(iMaxVal, ((drgnTrctLv[0][1] - 1) / 2));
		if (iMaxVal<2)
		{
			return &(rgctxpCoeffLayerCtx[2][0]);
		}
		else
		{
			if (iMaxVal<4)
			{
				return &(rgctxpCoeffLayerCtx[2][1]);
			}
			else
			{
				return &(rgctxpCoeffLayerCtx[3][(int)ceil(log(iMaxVal + 1) / log(2)) - 3]);
				/*
				if (iMaxVal<8)
				{
				return &(rgctxpCoeffLayerCtx[3][0]);
				}
				else
				{
				if (iMaxVal<16)
				{
				return &(rgctxpCoeffLayerCtx[3][1]);
				}
				else
				{
				if (iMaxVal<32)
				{
				return &(rgctxpCoeffLayerCtx[3][2]);
				}
				else
				{
				return &(rgctxpCoeffLayerCtx[3][3]);
				}
				}
				}
				*/
			}
		}
		break;
	case 4:
		iMaxVal = abs(ippCoeffBlk[0][3]);
		for (i = 1; i<4; i++)
		{
			iMaxVal = MAXIMUM(iMaxVal, abs(ippCoeffBlk[i][3 - i]));
		}
		iMaxVal = MINIMUM(iMaxVal, ((drgnTrctLv[0][3] - 1) / 2));
		for (i = 0; i<3; i++)
		{
			iMaxVal = MAXIMUM(iMaxVal, abs(ippCoeffBlk[i][2 - i]));
		}
		iMaxVal = MINIMUM(iMaxVal, ((drgnTrctLv[0][2] - 1) / 2));

		if (iMaxVal<1)
		{
			return &(rgctxpCoeffLayerCtx[4][5]);
		}
		else
		{
			return &(rgctxpCoeffLayerCtx[4][(int)ceil(log(iMaxVal + 1) / log(2)) - 1]);
			/*
			if (iMaxVal<2)
			{
			return &(rgctxpCoeffLayerCtx[4][0]);
			}
			else
			{
			if (iMaxVal<4)
			{
			return &(rgctxpCoeffLayerCtx[4][1]);
			}
			else
			{
			if (iMaxVal<8)
			{
			return &(rgctxpCoeffLayerCtx[4][2]);
			}
			else
			{
			if (iMaxVal<16)
			{
			return &(rgctxpCoeffLayerCtx[4][3]);
			}
			else
			{
			return &(rgctxpCoeffLayerCtx[4][4]);
			}
			}
			}
			}
			*/
		}
		break;
	case 5:
	case 6:
		iMaxVal = abs(ippCoeffBlk[0][4]);
		for (i = 1; i<5; i++)
		{
			iMaxVal = MAXIMUM(iMaxVal, abs(ippCoeffBlk[i][4 - i]));
		}
		iMaxVal = MINIMUM(iMaxVal, ((drgnTrctLv[0][4] - 1) / 2));

		if (iMaxVal<1)
		{
			return &(rgctxpCoeffLayerCtx[4][5]);
		}
		else
		{
			return &(rgctxpCoeffLayerCtx[5][(int)ceil(log(iMaxVal + 1) / log(2)) - 1]);
			/*
			if (iMaxVal<2)
			{
			return &(rgctxpCoeffLayerCtx[5][0]);
			}
			else
			{
			if (iMaxVal<4)
			{
			return &(rgctxpCoeffLayerCtx[5][1]);
			}
			else
			{
			return &(rgctxpCoeffLayerCtx[5][2]);
			}
			}
			*/
		}
		break;
	}

	assert(iRow + iCol>6);
	nNZVal = 0;
	iMaxVal = 0;

	switch (iRow + iCol)
	{
	case 7:
		if (iRow<2)
		{
			for (i = 3; i<7; i++)
			{
				nNZVal += (ippCoeffBlk[0][i] != 0);
				iMaxVal = MAXIMUM(iMaxVal, abs(ippCoeffBlk[0][i]));
			}
			for (i = 3; i<6; i++)
			{
				nNZVal += (ippCoeffBlk[1][i] != 0);
				iMaxVal = MAXIMUM(iMaxVal, abs(ippCoeffBlk[1][i]));
			}
			break;
		}
		if (iCol<2)
		{
			for (i = 3; i<7; i++)
			{
				nNZVal += (ippCoeffBlk[i][0] != 0);
				iMaxVal = MAXIMUM(iMaxVal, abs(ippCoeffBlk[i][0]));
			}
			for (i = 3; i<6; i++)
			{
				nNZVal += (ippCoeffBlk[i][1] != 0);
				iMaxVal = MAXIMUM(iMaxVal, abs(ippCoeffBlk[i][1]));
			}
			break;
		}
	case 8:
		if (iRow == 1)
		{
			for (i = 4; i<8; i++)
			{
				nNZVal += (ippCoeffBlk[0][i] != 0);
				iMaxVal = MAXIMUM(iMaxVal, abs(ippCoeffBlk[0][i]));
			}
			for (i = 4; i<7; i++)
			{
				nNZVal += (ippCoeffBlk[1][i] != 0);
				iMaxVal = MAXIMUM(iMaxVal, abs(ippCoeffBlk[1][i]));
			}
			break;
		}
		if (iCol == 1)
		{
			for (i = 4; i<8; i++)
			{
				nNZVal += (ippCoeffBlk[i][0] != 0);
				iMaxVal = MAXIMUM(iMaxVal, abs(ippCoeffBlk[i][0]));
			}
			for (i = 4; i<7; i++)
			{
				nNZVal += (ippCoeffBlk[i][1] != 0);
				iMaxVal = MAXIMUM(iMaxVal, abs(ippCoeffBlk[i][1]));
			}
			break;
		}
	default:
		if (iRow == 7 || iCol == 7)
		{
			nNZVal += (ippCoeffBlk[iRow - 1][iCol] != 0);
			iMaxVal = MAXIMUM(iMaxVal, abs(ippCoeffBlk[iRow - 1][iCol]));
			nNZVal += (ippCoeffBlk[iRow][iCol - 1] != 0);
			iMaxVal = MAXIMUM(iMaxVal, abs(ippCoeffBlk[iRow][iCol - 1]));

			nNZVal += (ippCoeffBlk[iRow - 2][iCol] != 0);
			iMaxVal = MAXIMUM(iMaxVal, abs(ippCoeffBlk[iRow - 2][iCol]));
			nNZVal += (ippCoeffBlk[iRow - 1][iCol - 1] != 0);
			iMaxVal = MAXIMUM(iMaxVal, abs(ippCoeffBlk[iRow - 1][iCol - 1]));
			nNZVal += (ippCoeffBlk[iRow][iCol - 2] != 0);
			iMaxVal = MAXIMUM(iMaxVal, abs(ippCoeffBlk[iRow][iCol - 2]));

			nNZVal += (ippCoeffBlk[iRow - 2][iCol - 1] != 0);
			iMaxVal = MAXIMUM(iMaxVal, abs(ippCoeffBlk[iRow - 2][iCol - 1]));
			nNZVal += (ippCoeffBlk[iRow - 1][iCol - 2] != 0);
			iMaxVal = MAXIMUM(iMaxVal, abs(ippCoeffBlk[iRow - 1][iCol - 2]));
		}
		else
		{
			for (i = 2; i >= -1; i--)
			{
				nNZVal += (ippCoeffBlk[iRow - i][iCol - 1 + i] != 0);
				iMaxVal = MAXIMUM(iMaxVal, abs(ippCoeffBlk[iRow - i][iCol - 1 + i]));
			}
			for (i = 2; i >= 0; i--)
			{
				nNZVal += (ippCoeffBlk[iRow - i][iCol - 2 + i] != 0);
				iMaxVal = MAXIMUM(iMaxVal, abs(ippCoeffBlk[iRow - i][iCol - 2 + i]));
			}
		}
		break;
	}

	if (isVal == 0)
	{
		if (nNZLayer<4)
		{
			return &(rgctxpCoeffLayerCtx[6][0]);
		}
		else
		{
			if (nNZLayer<8)
			{
				return &(rgctxpCoeffLayerCtx[6][1 + (int)ceil(log(nNZVal + 1) / log(2))]);
			}
			else
			{
				if (nNZVal<6)
				{
					return &(rgctxpCoeffLayerCtx[6][5 + (int)floor(nNZVal / 2)]);
				}
				else
				{
					return &(rgctxpCoeffLayerCtx[6][2 + nNZVal]);
				}
			}
		}
	}

	assert(isVal == 1);

	if (iRow + iCol<11)
	{
		if (nNZLayer<4)
		{
			return &(rgctxpCoeffLayerCtx[iRow + iCol][0]);
		}
		else
		{
			if (nNZLayer<8)
			{
				return &(rgctxpCoeffLayerCtx[iRow + iCol][1 + (int)(iMaxVal>1)]);
			}
			else
			{
				return &(rgctxpCoeffLayerCtx[iRow + iCol][3 + (int)(iMaxVal>1)]);
			}
		}
	}
	else
	{
		return &(rgctxpCoeffLayerCtx[iRow + iCol][(int)(iMaxVal>1)]);
	}
}

int fnCalNZLayer(int ippCoeffBlk[8][8])
{
	int nNZLayer = 0;
	int i;

	for (i = 0; i<6; i++)
	{
		nNZLayer += (ippCoeffBlk[i][5 - i] != 0);
	}
	for (i = 0; i<7; i++)
	{
		nNZLayer += (ippCoeffBlk[i][6 - i] != 0);
	}

	return nNZLayer;
}