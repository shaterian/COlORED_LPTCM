#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "global.h"
#include "AdapRecSpace.h"
#include "ImageDependent.h"
#include "HACACM.h"

#include<iostream>
ArithmeticCode* GACpEngine;
BitStream* GbspOutputBuffer;
Layer LyrInfo;
double dTestRate = 0.0;
double dTestRateHDQ = 0.0;
int iErrorNumSDQ = 0;
//int ctx_ctr[2][16] = {0};

static int rgiRowZZ[64] = { 0,0,1,2,1,0,0,1,2,3,4,3,2,1,0,0,1,
2,3,4,5,6,5,4,3,2,1,0,0,1,2,3,4,5,
6,7,7,6,5,4,3,2,1,2,3,4,5,6,7,7,6,
5,4,3,4,5,6,7,7,6,5,6,7,7 };
static int rgiColZZ[64] = { 0,1,0,0,1,2,3,2,1,0,0,1,2,3,4,5,4,
3,2,1,0,0,1,2,3,4,5,6,7,6,5,4,3,2,
1,0,1,2,3,4,5,6,7,7,6,5,4,3,2,3,4,
5,6,7,7,6,5,4,5,6,7,7,6,7 };

//int CtxDstrEXTF[2][16] = {
//	/*{965, 802, 452, 639, 770, 678, 444, 502, 532, 708, 378, 315, 397, 513, 228, 173},
//     {35,  198, 548, 361, 230, 322, 556, 498, 468, 292, 622, 685, 603, 487, 772, 827},*/
//	{965, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
//    {35,  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
//};
//
//int CtxDstrExt[6][64] = {1};

Context CtxDC[NUM_CTX_DC], *ExtCtx[NUM_CTX_EXT][BIT_LOC_EXT], *SigCtx[2][NUM_CTX_SIG][BIT_LOC_SIG], *AmpCtx[2][NUM_CTX_AMP], *ExtAmpCtx[NUM_CTX_EXTAMP], ExtFCtx[NUM_CTX_EXTF], AZFCtx[NUM_CTX_AZF];
Context *ExtAmpSeprtCtx[NUM_CTX_EXTAMP_SEPRT], SignCtx[2][NUM_CTX_SIGN][7];

int drgiCurIncCtx[4][3];
int *prgiSig[MAX_ROWS / 8][MAX_COLS / 8][BIT_LOC_SIG];
//int ipAZFCtxMainCtr[NUM_CTX_AZF],i2pAZFCtxCtr[NUM_CTX_AZF][2];
int **i2pAmpCtrEmp, i5pSigCtxCtr[2][NUM_CTX_SIG][2][7][2], *i3pAmpCtxCtr[2][7], **i3pAmpCtxElemtCtr[2], i5pSigCtxElemtCtr[2][NUM_CTX_SIG][2][63][2] = { 1 };
double d5pSigCtxEntpy[2][NUM_CTX_SIG][2][7][2], **d3pAmpCtxEntpy[2];//d2pAZFEntpy[NUM_CTX_AZF][2]
double dRateTmp = 0.0, drgdEmpEnt[8][8] = { 0.0 };
/***************************************************************/
///int vrgiSign[ROWS / 8][COLS / 8][8][8] = { 0 };
std::vector<std::vector<std::vector<std::vector<int>>>> vrgiSign;
//int vrgiExt[ROWS / 8][COLS / 8][8][8][BIT_LOC_EXT] = { 0 };
int vrgiExt[MAX_ROWS/8 ][MAX_COLS/8 ][8][8][BIT_LOC_EXT];

//int vrgiExtDec[ROWS / 8][COLS / 8][8][8][BIT_LOC_EXT] = { 0 };
int vrgiExtDec[MAX_ROWS / 8][MAX_COLS / 8][8][8][BIT_LOC_EXT];

//int vrgiHDQIdxDec[ROWS / 8][COLS / 8][8][8] = { 0 };
int vrgiHDQIdxDec[MAX_ROWS / 8][MAX_COLS / 8][8][8] = { 0 };

//int vrgiSig[ROWS / 8][COLS / 8][8][8][BIT_LOC_SIG] = { 0 };
int vrgiSig[MAX_ROWS / 8][MAX_COLS / 8][8][8][BIT_LOC_SIG] = { 0 };


/*
*  Function: fnElemtStup
*  Setup the layer information in preparing for the following entropy coding
*/
void fnLyrStup()
{
	int iFreqRow, iFreqCol, iFreqRow1, iFreqCol1, iScanLyrNum, iSwapFlag, iFreqRowTmp, iFreqColTmp, iRowZZplus1, iColZZplus1, iCurScnLyrIdx, iCurScnLyrNumElemt;
	int i, j, iTmp, iRowZZ, iColZZ, iZZIdx;
	int iMaxL = 0, iMaxLprime = 0, iLyrMaxL = 0, iLyrMaxLprime = 0, iMaxNumLyr = 0, iNumElemtLLr1 = 0;
	int iCtxLyrIdxMaxL[7], iCtxLyrIdxMaxLprime[7];

	LyrInfo.ExtNumElemt = 0;
	LyrInfo.SigNumElemt = 0;

	iZZIdx = 1;
	while (iZZIdx < 64)
	{
		iRowZZ = rgiRowZZ[iZZIdx];
		iColZZ = rgiColZZ[iZZIdx];
		if (drgiTTbl[iRowZZ][iColZZ] > 0)
		{
			LyrInfo.ExtScanRowIdx[LyrInfo.ExtNumElemt] = iRowZZ;
			LyrInfo.ExtScanColIdx[LyrInfo.ExtNumElemt] = iColZZ;
			LyrInfo.ExtNumElemt++;

			if (drgiLTbl[iRowZZ][iColZZ] > 0)
			{
				LyrInfo.SigScanRowIdx[LyrInfo.SigNumElemt] = iRowZZ;
				LyrInfo.SigScanColIdx[LyrInfo.SigNumElemt] = iColZZ;
				LyrInfo.SigNumElemt++;
			}

			//            if (drgiLTbl[iRowZZ][iColZZ] > 0)
			//            {
			//                LyrInfo.SigScanRowIdx[LyrInfo.SigNumElemt] = iRowZZ;
			//                LyrInfo.SigScanColIdx[LyrInfo.SigNumElemt] = iColZZ;
			//                LyrInfo.SigNumElemt++; 
			//            }
		}

		iZZIdx++;
	}

	//adptive scanning according to the alphabet size
#if ADP_SCN
	LyrInfo.ExtNumElemt = 0;

	iZZIdx = 1;
	iScanLyrNum = 1;
	while (iScanLyrNum < 14)
	{
		iScanLyrNum = rgiRowZZ[iZZIdx] + rgiColZZ[iZZIdx];
		iSwapFlag = 0;

		if (iScanLyrNum < 3 || iScanLyrNum > 12)
		{
			if (drgiLPrimeTbl[rgiColZZ[iZZIdx]][rgiRowZZ[iZZIdx]] > drgiLPrimeTbl[rgiRowZZ[iZZIdx]][rgiColZZ[iZZIdx]])
			{
				iSwapFlag = 1;
			}

			while (rgiRowZZ[iZZIdx] + rgiColZZ[iZZIdx] == iScanLyrNum)
			{
				if (iSwapFlag == 1)
				{
					iRowZZ = rgiColZZ[iZZIdx];
					iColZZ = rgiRowZZ[iZZIdx];
				}
				else
				{
					iRowZZ = rgiRowZZ[iZZIdx];
					iColZZ = rgiColZZ[iZZIdx];
				}

				if (drgiTTbl[iRowZZ][iColZZ] > 0)
				{
					LyrInfo.ExtScanRowIdx[LyrInfo.ExtNumElemt] = iRowZZ;
					LyrInfo.ExtScanColIdx[LyrInfo.ExtNumElemt] = iColZZ;
					LyrInfo.ExtNumElemt++;
				}

				iZZIdx++;
			}
		}

		else if (iScanLyrNum < 5 || iScanLyrNum > 10)
		{
			if (drgiLPrimeTbl[rgiColZZ[iZZIdx]][rgiRowZZ[iZZIdx]] + drgiLPrimeTbl[rgiColZZ[iZZIdx + 1]][rgiRowZZ[iZZIdx + 1]] > drgiLPrimeTbl[rgiRowZZ[iZZIdx]][rgiColZZ[iZZIdx]] + drgiLPrimeTbl[rgiRowZZ[iZZIdx + 1]][rgiColZZ[iZZIdx + 1]])
			{
				iSwapFlag = 1;
			}

			while (rgiRowZZ[iZZIdx] + rgiColZZ[iZZIdx] == iScanLyrNum)
			{
				if (iSwapFlag == 1)
				{
					iRowZZ = rgiColZZ[iZZIdx];
					iColZZ = rgiRowZZ[iZZIdx];
				}
				else
				{
					iRowZZ = rgiRowZZ[iZZIdx];
					iColZZ = rgiColZZ[iZZIdx];
				}

				if (drgiTTbl[iRowZZ][iColZZ] > 0)
				{
					LyrInfo.ExtScanRowIdx[LyrInfo.ExtNumElemt] = iRowZZ;
					LyrInfo.ExtScanColIdx[LyrInfo.ExtNumElemt] = iColZZ;
					LyrInfo.ExtNumElemt++;
				}

				iZZIdx++;
			}
		}


		else if (iScanLyrNum < 7 || iScanLyrNum > 8)
		{
			if (drgiLPrimeTbl[rgiColZZ[iZZIdx]][rgiRowZZ[iZZIdx]] + drgiLPrimeTbl[rgiColZZ[iZZIdx + 1]][rgiRowZZ[iZZIdx + 1]] + drgiLPrimeTbl[rgiColZZ[iZZIdx + 2]][rgiRowZZ[iZZIdx + 2]] > drgiLPrimeTbl[rgiRowZZ[iZZIdx]][rgiColZZ[iZZIdx]] + drgiLPrimeTbl[rgiRowZZ[iZZIdx + 1]][rgiColZZ[iZZIdx + 1]] + drgiLPrimeTbl[rgiRowZZ[iZZIdx + 2]][rgiColZZ[iZZIdx + 2]])
			{
				iSwapFlag = 1;
			}

			while (rgiRowZZ[iZZIdx] + rgiColZZ[iZZIdx] == iScanLyrNum)
			{
				if (iSwapFlag == 1)
				{
					iRowZZ = rgiColZZ[iZZIdx];
					iColZZ = rgiRowZZ[iZZIdx];
				}
				else
				{
					iRowZZ = rgiRowZZ[iZZIdx];
					iColZZ = rgiColZZ[iZZIdx];
				}

				if (drgiTTbl[iRowZZ][iColZZ] > 0)
				{
					LyrInfo.ExtScanRowIdx[LyrInfo.ExtNumElemt] = iRowZZ;
					LyrInfo.ExtScanColIdx[LyrInfo.ExtNumElemt] = iColZZ;
					LyrInfo.ExtNumElemt++;
				}

				iZZIdx++;
			}
		}

		else
		{
			if (drgiLPrimeTbl[rgiColZZ[iZZIdx]][rgiRowZZ[iZZIdx]] + drgiLPrimeTbl[rgiColZZ[iZZIdx + 1]][rgiRowZZ[iZZIdx + 1]] + drgiLPrimeTbl[rgiColZZ[iZZIdx + 2]][rgiRowZZ[iZZIdx + 2]] + drgiLPrimeTbl[rgiColZZ[iZZIdx + 3]][rgiRowZZ[iZZIdx + 3]] > drgiLPrimeTbl[rgiRowZZ[iZZIdx]][rgiColZZ[iZZIdx]] + drgiLPrimeTbl[rgiRowZZ[iZZIdx + 1]][rgiColZZ[iZZIdx + 1]] + drgiLPrimeTbl[rgiRowZZ[iZZIdx + 2]][rgiColZZ[iZZIdx + 2]] + drgiLPrimeTbl[rgiRowZZ[iZZIdx + 3]][rgiColZZ[iZZIdx + 3]])
			{
				iSwapFlag = 1;
			}

			while (rgiRowZZ[iZZIdx] + rgiColZZ[iZZIdx] == iScanLyrNum)
			{
				if (iSwapFlag == 1)
				{
					iRowZZ = rgiColZZ[iZZIdx];
					iColZZ = rgiRowZZ[iZZIdx];
				}
				else
				{
					iRowZZ = rgiRowZZ[iZZIdx];
					iColZZ = rgiColZZ[iZZIdx];
				}

				if (drgiTTbl[iRowZZ][iColZZ] > 0)
				{
					LyrInfo.ExtScanRowIdx[LyrInfo.ExtNumElemt] = iRowZZ;
					LyrInfo.ExtScanColIdx[LyrInfo.ExtNumElemt] = iColZZ;
					LyrInfo.ExtNumElemt++;
				}

				iZZIdx++;
			}
		}
	}

	LyrInfo.SigNumElemt = 0;

	iZZIdx = 1;
	iScanLyrNum = 1;
	while (iScanLyrNum < 14)
	{
		iScanLyrNum = rgiRowZZ[iZZIdx] + rgiColZZ[iZZIdx];
		LyrInfo.SwapFlag[iScanLyrNum] = 0;

		if (iScanLyrNum < 3 || iScanLyrNum > 12)
		{
			if (drgiLTbl[rgiColZZ[iZZIdx]][rgiRowZZ[iZZIdx]] > drgiLTbl[rgiRowZZ[iZZIdx]][rgiColZZ[iZZIdx]])
			{
				LyrInfo.SwapFlag[iScanLyrNum] = 1;
			}

			while (rgiRowZZ[iZZIdx] + rgiColZZ[iZZIdx] == iScanLyrNum)
			{
				if (LyrInfo.SwapFlag[iScanLyrNum] == 1)
				{
					iRowZZ = rgiColZZ[iZZIdx];
					iColZZ = rgiRowZZ[iZZIdx];
				}
				else
				{
					iRowZZ = rgiRowZZ[iZZIdx];
					iColZZ = rgiColZZ[iZZIdx];
				}

				if (drgiLTbl[iRowZZ][iColZZ] > 0)
				{
					LyrInfo.SigScanRowIdx[LyrInfo.SigNumElemt] = iRowZZ;
					LyrInfo.SigScanColIdx[LyrInfo.SigNumElemt] = iColZZ;
					LyrInfo.SigNumElemt++;
				}

				iZZIdx++;
			}
		}

		else if (iScanLyrNum < 5 || iScanLyrNum > 10)
		{
			if (drgiLTbl[rgiColZZ[iZZIdx]][rgiRowZZ[iZZIdx]] + drgiLTbl[rgiColZZ[iZZIdx + 1]][rgiRowZZ[iZZIdx + 1]] > drgiLTbl[rgiRowZZ[iZZIdx]][rgiColZZ[iZZIdx]] + drgiLTbl[rgiRowZZ[iZZIdx + 1]][rgiColZZ[iZZIdx + 1]])
			{
				LyrInfo.SwapFlag[iScanLyrNum] = 1;
			}

			while (rgiRowZZ[iZZIdx] + rgiColZZ[iZZIdx] == iScanLyrNum)
			{
				if (LyrInfo.SwapFlag[iScanLyrNum] == 1)
				{
					iRowZZ = rgiColZZ[iZZIdx];
					iColZZ = rgiRowZZ[iZZIdx];
				}
				else
				{
					iRowZZ = rgiRowZZ[iZZIdx];
					iColZZ = rgiColZZ[iZZIdx];
				}

				if (drgiLTbl[iRowZZ][iColZZ] > 0)
				{
					LyrInfo.SigScanRowIdx[LyrInfo.SigNumElemt] = iRowZZ;
					LyrInfo.SigScanColIdx[LyrInfo.SigNumElemt] = iColZZ;
					LyrInfo.SigNumElemt++;
				}

				iZZIdx++;
			}
		}


		else if (iScanLyrNum < 7 || iScanLyrNum > 8)
		{
			if (drgiLTbl[rgiColZZ[iZZIdx]][rgiRowZZ[iZZIdx]] + drgiLTbl[rgiColZZ[iZZIdx + 1]][rgiRowZZ[iZZIdx + 1]] + drgiLTbl[rgiColZZ[iZZIdx + 2]][rgiRowZZ[iZZIdx + 2]] > drgiLTbl[rgiRowZZ[iZZIdx]][rgiColZZ[iZZIdx]] + drgiLTbl[rgiRowZZ[iZZIdx + 1]][rgiColZZ[iZZIdx + 1]] + drgiLTbl[rgiRowZZ[iZZIdx + 2]][rgiColZZ[iZZIdx + 2]])
			{
				LyrInfo.SwapFlag[iScanLyrNum] = 1;
			}

			while (rgiRowZZ[iZZIdx] + rgiColZZ[iZZIdx] == iScanLyrNum)
			{
				if (LyrInfo.SwapFlag[iScanLyrNum] == 1)
				{
					iRowZZ = rgiColZZ[iZZIdx];
					iColZZ = rgiRowZZ[iZZIdx];
				}
				else
				{
					iRowZZ = rgiRowZZ[iZZIdx];
					iColZZ = rgiColZZ[iZZIdx];
				}

				if (drgiLTbl[iRowZZ][iColZZ] > 0)
				{
					LyrInfo.SigScanRowIdx[LyrInfo.SigNumElemt] = iRowZZ;
					LyrInfo.SigScanColIdx[LyrInfo.SigNumElemt] = iColZZ;
					LyrInfo.SigNumElemt++;
				}

				iZZIdx++;
			}
		}

		else
		{
			if (drgiLTbl[rgiColZZ[iZZIdx]][rgiRowZZ[iZZIdx]] + drgiLTbl[rgiColZZ[iZZIdx + 1]][rgiRowZZ[iZZIdx + 1]] + drgiLTbl[rgiColZZ[iZZIdx + 2]][rgiRowZZ[iZZIdx + 2]] + drgiLTbl[rgiColZZ[iZZIdx + 3]][rgiRowZZ[iZZIdx + 3]] > drgiLTbl[rgiRowZZ[iZZIdx]][rgiColZZ[iZZIdx]] + drgiLTbl[rgiRowZZ[iZZIdx + 1]][rgiColZZ[iZZIdx + 1]] + drgiLTbl[rgiRowZZ[iZZIdx + 2]][rgiColZZ[iZZIdx + 2]] + drgiLTbl[rgiRowZZ[iZZIdx + 3]][rgiColZZ[iZZIdx + 3]])
			{
				LyrInfo.SwapFlag[iScanLyrNum] = 1;
			}

			while (rgiRowZZ[iZZIdx] + rgiColZZ[iZZIdx] == iScanLyrNum)
			{
				if (LyrInfo.SwapFlag[iScanLyrNum] == 1)
				{
					iRowZZ = rgiColZZ[iZZIdx];
					iColZZ = rgiRowZZ[iZZIdx];
				}
				else
				{
					iRowZZ = rgiRowZZ[iZZIdx];
					iColZZ = rgiColZZ[iZZIdx];
				}

				if (drgiLTbl[iRowZZ][iColZZ] > 0)
				{
					LyrInfo.SigScanRowIdx[LyrInfo.SigNumElemt] = iRowZZ;
					LyrInfo.SigScanColIdx[LyrInfo.SigNumElemt] = iColZZ;
					LyrInfo.SigNumElemt++;
				}

				iZZIdx++;
			}
		}
	}

#endif

	/* setup the context layer info */
	for (i = 0; i<LyrInfo.ExtNumElemt; i++)
	{
		iFreqRow = LyrInfo.ExtScanRowIdx[i];
		iFreqCol = LyrInfo.ExtScanColIdx[i];

		iTmp = iFreqRow + iFreqCol;

		if (iTmp == 1)
		{
			LyrInfo.CtxLyr[i] = i;
			LyrInfo.CtxLyrSig[i] = i;
		}

		else if (iTmp < 5)
		{
			LyrInfo.CtxLyr[i] = iTmp;
			LyrInfo.CtxLyrSig[i] = 2;
		}

		else if (iTmp < 7)
		{
			LyrInfo.CtxLyr[i] = iTmp - 1;
			LyrInfo.CtxLyrSig[i] = 3;
		}

		else if (iTmp == 7 || iTmp == 8)
		{
			LyrInfo.CtxLyr[i] = 5;
			LyrInfo.CtxLyrSig[i] = 4;
		}

		else
		{
			LyrInfo.CtxLyr[i] = 6;
			LyrInfo.CtxLyrSig[i] = 3;
		}

		//		if (iTmp < 5)
		//		{
		//		    LyrInfo.CtxLyrSig[i] = i;
		//		}
		//        
		//		else
		//		{
		//			if (iTmp == 5)
		//		        LyrInfo.CtxLyrSig[i] = LyrInfo.ExtNumElemt-2;
		//			else
		//				LyrInfo.CtxLyrSig[i] = LyrInfo.ExtNumElemt-1;
		//		}
	}

	j = 0;
	for (i = 0; i<LyrInfo.ExtNumElemt; i++)
	{
		iFreqRow = LyrInfo.ExtScanRowIdx[i];
		iFreqCol = LyrInfo.ExtScanColIdx[i];

		if (iMaxLprime < drgiLPrimeTbl[iFreqRow][iFreqCol])
		{
			LyrInfo.LyrIdx_MaxLprime = i;
			iMaxLprime = drgiLPrimeTbl[iFreqRow][iFreqCol];
		}

		if (i > 0 && (LyrInfo.CtxLyr[i] - LyrInfo.CtxLyr[i - 1]) == 1)
		{
			iLyrMaxLprime = 0;
			j++;
		}

		if (iLyrMaxLprime < drgiLPrimeTbl[iFreqRow][iFreqCol])
		{
			iLyrMaxLprime = drgiLPrimeTbl[iFreqRow][iFreqCol];
			iCtxLyrIdxMaxLprime[j] = i;
		}

	}

	for (i = 0; i<LyrInfo.ExtNumElemt; i++)
	{
		LyrInfo.CtxLyrIdxMaxLprime[i] = iCtxLyrIdxMaxLprime[LyrInfo.CtxLyr[i]];
	}


	j = 0;
	iMaxNumLyr = 0;
	iCurScnLyrIdx = 1;
	iCurScnLyrNumElemt = 0;
	LyrInfo.SigNumLyr = 0;
	for (i = 0; i<LyrInfo.SigNumElemt; i++)
	{
		iFreqRow = LyrInfo.SigScanRowIdx[i];
		iFreqCol = LyrInfo.SigScanColIdx[i];

		LyrInfo.SigAphSize[i] = drgiLTbl[iFreqRow][iFreqCol];

		iTmp = iFreqRow + iFreqCol;
		//LyrInfo.ScanLyrIdxSig[i] = iTmp;

		if (iCurScnLyrIdx == iTmp)
		{
			iCurScnLyrNumElemt++;
		}
		else
		{
			LyrInfo.SigLyrNumElemt[iCurScnLyrIdx - 1] = iCurScnLyrNumElemt;
			LyrInfo.SigNumLyr++;
			iCurScnLyrNumElemt = 1;
			iCurScnLyrIdx++;
		}

		if (drgiLTbl[iFreqRow][iFreqCol] > 0)
		{
			iNumElemtLLr1++;
			iMaxNumLyr = MAXIMUM(iMaxNumLyr, iTmp);
		}

		if (iMaxL < drgiLTbl[iFreqRow][iFreqCol])
		{
			LyrInfo.LyrIdx_MaxL = i;
			iMaxL = drgiLTbl[iFreqRow][iFreqCol];
		}

		if (i > 0 && (LyrInfo.CtxLyr[i] - LyrInfo.CtxLyr[i - 1]) == 1)
		{
			iLyrMaxL = 0;
			j++;
		}

		if (iLyrMaxL < drgiLTbl[iFreqRow][iFreqCol])
		{
			iLyrMaxL = drgiLTbl[iFreqRow][iFreqCol];
			iCtxLyrIdxMaxL[j] = i;
		}

	}

	LyrInfo.SigLyrNumElemt[iCurScnLyrIdx - 1] = iCurScnLyrNumElemt;
	LyrInfo.SigNumLyr++;

	for (i = 0; i<LyrInfo.SigNumElemt; i++)
	{
		LyrInfo.CtxLyrIdxMaxL[i] = iCtxLyrIdxMaxL[LyrInfo.CtxLyr[i]];
	}

	if (iNumElemtLLr1 < 7)
		LyrInfo.Acum1Pos = 0;

	else
	{
		LyrInfo.Acum1Pos = INFINITY;
		for (i = 0; i<LyrInfo.SigNumElemt; i++)
		{
			iFreqRow = LyrInfo.SigScanRowIdx[i];
			iFreqCol = LyrInfo.SigScanColIdx[i];

			iTmp = iFreqRow + iFreqCol;
			if (iTmp == (iMaxNumLyr - 1))
			{
				LyrInfo.Acum1Pos = MINIMUM(LyrInfo.Acum1Pos, i);
			}
		}
	}

	/* setup the context template for Sig0-2 */
	for (i = 0; i<LyrInfo.SigNumElemt; i++)
	{
		iFreqRow = LyrInfo.SigScanRowIdx[i];
		iFreqCol = LyrInfo.SigScanColIdx[i];
		LyrInfo.SigNumCtx[i] = 0;

		if (iFreqRow == 7)
		{
			if (iFreqCol == 7)
				LyrInfo.SigNumCtx[i] = 0;

			else if (iFreqCol == 6)
			{
				if (LyrInfo.SigNumElemt == 63)
				{
					LyrInfo.SigNumCtx[i] = 1;
					LyrInfo.SigCtxMdl[i][0] = 63;
				}
			}

			else
			{
				for (j = i + 1; j<LyrInfo.SigNumElemt; j++)
				{
					iFreqRowTmp = LyrInfo.SigScanRowIdx[j];
					iFreqColTmp = LyrInfo.SigScanColIdx[j];

					if (iFreqRowTmp == iFreqRow && iFreqColTmp == (iFreqCol + 1))
					{
						LyrInfo.SigCtxMdl[i][LyrInfo.SigNumCtx[i]] = j;
						LyrInfo.SigNumCtx[i]++;
					}

					if (iFreqRowTmp == iFreqRow && iFreqColTmp == (iFreqCol + 2))
					{
						LyrInfo.SigCtxMdl[i][LyrInfo.SigNumCtx[i]] = j;
						LyrInfo.SigNumCtx[i]++;
					}
				}
			}
		}

		else if (iFreqRow == 6)
		{
			if (iFreqCol == 7)
			{
				if (LyrInfo.SigNumElemt == 63)
				{
					LyrInfo.SigNumCtx[i] = 1;
					LyrInfo.SigCtxMdl[i][0] = 63;
				}
			}

			else if (iFreqCol == 6)
			{
				if (LyrInfo.SigNumElemt == 63)
				{
					LyrInfo.SigCtxMdl[i][0] = 63;
					LyrInfo.SigNumCtx[i]++;
				}

				for (j = i + 1; j<LyrInfo.SigNumElemt; j++)
				{
					iFreqRowTmp = LyrInfo.SigScanRowIdx[j];
					iFreqColTmp = LyrInfo.SigScanColIdx[j];

					if (iFreqRowTmp == iFreqRow && iFreqColTmp == (iFreqCol + 1))
					{
						LyrInfo.SigCtxMdl[i][LyrInfo.SigNumCtx[i]] = j;
						LyrInfo.SigNumCtx[i]++;
					}

					if (iFreqRowTmp == (iFreqRow + 1) && iFreqColTmp == iFreqCol)
					{
						LyrInfo.SigCtxMdl[i][LyrInfo.SigNumCtx[i]] = j;
						LyrInfo.SigNumCtx[i]++;
					}
				}
			}

			else
			{
				for (j = i + 1; j<LyrInfo.SigNumElemt; j++)
				{
					iFreqRowTmp = LyrInfo.SigScanRowIdx[j];
					iFreqColTmp = LyrInfo.SigScanColIdx[j];

					if (iFreqRowTmp == iFreqRow && iFreqColTmp == (iFreqCol + 1))
					{
						LyrInfo.SigCtxMdl[i][LyrInfo.SigNumCtx[i]] = j;
						LyrInfo.SigNumCtx[i]++;
					}

					if (iFreqRowTmp == (iFreqRow + 1) && iFreqColTmp == iFreqCol)
					{
						LyrInfo.SigCtxMdl[i][LyrInfo.SigNumCtx[i]] = j;
						LyrInfo.SigNumCtx[i]++;
					}

					if (iFreqRowTmp == (iFreqRow + 1) && iFreqColTmp == (iFreqCol + 1))
					{
						LyrInfo.SigCtxMdl[i][LyrInfo.SigNumCtx[i]] = j;
						LyrInfo.SigNumCtx[i]++;
					}

					if (iFreqRowTmp == iFreqRow && iFreqColTmp == (iFreqCol + 2))
					{
						LyrInfo.SigCtxMdl[i][LyrInfo.SigNumCtx[i]] = j;
						LyrInfo.SigNumCtx[i]++;
					}
				}
			}
		}

		else if (iFreqCol == 7)
		{
			for (j = i + 1; j<LyrInfo.SigNumElemt; j++)
			{
				iFreqRowTmp = LyrInfo.SigScanRowIdx[j];
				iFreqColTmp = LyrInfo.SigScanColIdx[j];

				if (iFreqRowTmp == (iFreqRow + 1) && iFreqColTmp == iFreqCol)
				{
					LyrInfo.SigCtxMdl[i][LyrInfo.SigNumCtx[i]] = j;
					LyrInfo.SigNumCtx[i]++;
				}

				if (iFreqRowTmp == (iFreqRow + 2) && iFreqColTmp == iFreqCol)
				{
					LyrInfo.SigCtxMdl[i][LyrInfo.SigNumCtx[i]] = j;
					LyrInfo.SigNumCtx[i]++;
				}
			}
		}

		else if (iFreqCol == 6)
		{
			for (j = i + 1; j<LyrInfo.SigNumElemt; j++)
			{
				iFreqRowTmp = LyrInfo.SigScanRowIdx[j];
				iFreqColTmp = LyrInfo.SigScanColIdx[j];

				if (iFreqRowTmp == iFreqRow && iFreqColTmp == (iFreqCol + 1))
				{
					LyrInfo.SigCtxMdl[i][LyrInfo.SigNumCtx[i]] = j;
					LyrInfo.SigNumCtx[i]++;
				}

				if (iFreqRowTmp == (iFreqRow + 1) && iFreqColTmp == iFreqCol)
				{
					LyrInfo.SigCtxMdl[i][LyrInfo.SigNumCtx[i]] = j;
					LyrInfo.SigNumCtx[i]++;
				}

				if (iFreqRowTmp == (iFreqRow + 1) && iFreqColTmp == (iFreqCol + 1))
				{
					LyrInfo.SigCtxMdl[i][LyrInfo.SigNumCtx[i]] = j;
					LyrInfo.SigNumCtx[i]++;
				}

				if (iFreqRowTmp == (iFreqRow + 2) && iFreqColTmp == iFreqCol)
				{
					LyrInfo.SigCtxMdl[i][LyrInfo.SigNumCtx[i]] = j;
					LyrInfo.SigNumCtx[i]++;
				}
			}
		}

		else
		{
			for (j = i + 1; j<LyrInfo.SigNumElemt; j++)
			{
				iFreqRowTmp = LyrInfo.SigScanRowIdx[j];
				iFreqColTmp = LyrInfo.SigScanColIdx[j];

				if (iFreqRowTmp == iFreqRow && iFreqColTmp == (iFreqCol + 1))
				{
					LyrInfo.SigCtxMdl[i][LyrInfo.SigNumCtx[i]] = j;
					LyrInfo.SigNumCtx[i]++;
				}

				if (iFreqRowTmp == (iFreqRow + 1) && iFreqColTmp == iFreqCol)
				{
					LyrInfo.SigCtxMdl[i][LyrInfo.SigNumCtx[i]] = j;
					LyrInfo.SigNumCtx[i]++;
				}

				if (iFreqRowTmp == (iFreqRow + 1) && iFreqColTmp == (iFreqCol + 1))
				{
					LyrInfo.SigCtxMdl[i][LyrInfo.SigNumCtx[i]] = j;
					LyrInfo.SigNumCtx[i]++;
				}

				if (iFreqRowTmp == (iFreqRow + 2) && iFreqColTmp == iFreqCol)
				{
					LyrInfo.SigCtxMdl[i][LyrInfo.SigNumCtx[i]] = j;
					LyrInfo.SigNumCtx[i]++;
				}

				if (iFreqRowTmp == iFreqRow && iFreqColTmp == (iFreqCol + 2))
				{
					LyrInfo.SigCtxMdl[i][LyrInfo.SigNumCtx[i]] = j;
					LyrInfo.SigNumCtx[i]++;
				}
			}
		}

		//        if (i < (LyrInfo.SigNumElemt-1))
		//        {    
		//            iRowZZplus1 = LyrInfo.SigScanRowIdx[i+1];
		//            iColZZplus1 = LyrInfo.SigScanColIdx[i+1];
		//            
		//            if ((iRowZZplus1+iColZZplus1) == (iFreqRow+iFreqCol))
		//            {
		//#if IMMED_NEIBR
		//                if (iRowZZplus1 == iFreqRow+1 || iRowZZplus1 == iFreqCol-1) 
		//                {
		//                    LyrInfo.SigCtxMdl[i][LyrInfo.SigNumCtx[i]] = i+1;
		//                    LyrInfo.SigNumCtx[i]++;
		//                }
		//#else
		//                LyrInfo.SigCtxMdl[i][LyrInfo.SigNumCtx[i]] = i+1;
		//                LyrInfo.SigNumCtx[i]++;
		//#endif
		//            }
		//        }
	}
}


/*
* FUNCTION: fnMemAllocCtr
* allocate memory for the counters and initialize AC for each context
*/
void fnMemAllocCtr()
{
	int i, j, k, l, iFreqRow, iFreqCol;

	for (i = 0; i<NUM_CTX_DC; i++)
	{
#if DC_SEP_CTX
		fnCreateCtx(&(CtxDC[i]), rgiAlpSizeDC[i]);//
#else
		fnCreateCtx(&(CtxDC[i]), 1 + iMaxLDC);
#endif
	}

	for (i = 0; i<NUM_CTX_AZF; i++)
	{
		fnCreateCtx(&AZFCtx[i], 2);
	}

	for (i = 0; i<NUM_CTX_EXTF; i++)
	{
		fnCreateCtx(&ExtFCtx[i], 2);
	}

	for (i = 0; i<NUM_CTX_EXT; i++)
	{
		for (j = 0; j<BIT_LOC_EXT; j++)
		{
			ExtCtx[i][j] = (Context *)malloc(sizeof(Context)*LyrInfo.ExtNumElemt);
		}
	}

	for (i = 0; i<NUM_CTX_EXTAMP_SEPRT; i++)
	{
		ExtAmpSeprtCtx[i] = (Context *)malloc(sizeof(Context)*LyrInfo.ExtNumElemt);
	}

	for (i = 0; i<2; i++)
	{
		for (j = 0; j<NUM_CTX_SIGN; j++)
		{
			for (k = 0; k<7; k++)
			{
				fnCreateCtx(&(SignCtx[i][j][k]), 2);
			}
		}
	}

	for (i = 0; i<NUM_CTX_EXT; i++)
	{
		for (j = 0; j<BIT_LOC_EXT; j++)
		{
			for (k = 0; k<LyrInfo.ExtNumElemt; k++)
			{
				fnCreateCtx(&(ExtCtx[i][j][k]), 2);
			}
		}
	}

	for (i = 0; i<NUM_CTX_EXTAMP_SEPRT; i++)
	{
		for (j = 0; j<LyrInfo.ExtNumElemt; j++)
		{
			iFreqRow = LyrInfo.ExtScanRowIdx[j];
			iFreqCol = LyrInfo.ExtScanColIdx[j];
			fnCreateCtx(&(ExtAmpSeprtCtx[i][j]), MAXIMUM(1, drgiLPrimeTbl[iFreqRow][iFreqCol] - BIT_LOC_EXT + 1));
		}
	}

	for (i = 0; i<2; i++)
	{
		for (j = 0; j<NUM_CTX_SIG; j++)
		{
			for (k = 0; k<BIT_LOC_SIG; k++)
			{
				SigCtx[i][j][k] = (Context *)malloc(sizeof(Context)*LyrInfo.SigNumElemt);
			}

		}

		for (j = 0; j<NUM_CTX_AMP; j++)
		{
			AmpCtx[i][j] = (Context *)malloc(sizeof(Context)*LyrInfo.SigNumElemt);
		}

	}

	for (i = 0; i<2; i++)
	{
		for (j = 0; j<NUM_CTX_SIG; j++)
		{
			for (k = 0; k<BIT_LOC_SIG; k++)
			{
				for (l = 0; l<LyrInfo.SigNumElemt; l++)
				{
					fnCreateCtx(&(SigCtx[i][j][k][l]), 2);
				}
			}
		}
	}

	for (i = 0; i<2; i++)
	{
		for (j = 0; j<NUM_CTX_AMP; j++)
		{
			for (k = 0; k<LyrInfo.SigNumElemt; k++)
			{
				fnCreateCtx(&(AmpCtx[i][j][k]), MAXIMUM(1, LyrInfo.SigAphSize[k] - BIT_LOC_SIG + 1));
			}
		}
	}

	for (i = 0; i<NUM_CTX_EXTAMP; i++)
	{
		ExtAmpCtx[i] = (Context *)malloc(sizeof(Context)*LyrInfo.ExtNumElemt);
	}

	for (i = 0; i<NUM_CTX_EXTAMP; i++)
	{
		for (j = 0; j<LyrInfo.ExtNumElemt; j++)
		{
			fnCreateCtx(&(ExtAmpCtx[i][j]), MAXIMUM(1, iMaxLprime - BIT_LOC_EXT + 1));
		}
	}

	if (SDQ)
	{
		for (i = 0; i<4; i++)
		{
			for (j = 0; j<3; j++)
			{
				drgiCurIncCtx[i][j] = i>j ? 1 : 0;
			}
		}

		for (i = 0; i<ROWS / 8; i++)
		{
			for (j = 0; j<COLS / 8; j++)
			{
				for (k = 0; k<BIT_LOC_SIG; k++)
				{
					prgiSig[i][j][k] = (int *)malloc(sizeof(int)*LyrInfo.SigNumElemt);
				}
			}
		}

		for (i = 0; i<2; i++)
		{
			for (j = 0; j<7; j++)
			{
				i3pAmpCtxCtr[i][j] = (int *)malloc(sizeof(int)*(MAXIMUM(1, LyrInfo.SigAphSize[LyrInfo.CtxLyrIdxMaxL[j]] - BIT_LOC_SIG + 1)));
			}

			d3pAmpCtxEntpy[i] = (double **)malloc(sizeof(double *)*LyrInfo.SigNumElemt);
		}

		for (i = 0; i<2; i++)
		{
			for (j = 0; j<LyrInfo.SigNumElemt; j++)
			{
				d3pAmpCtxEntpy[i][j] = (double *)malloc(sizeof(double)*(LyrInfo.SigAphSize[j] - BIT_LOC_SIG + 1));
			}
		}

#if TEST_RD
		i2pAmpCtrEmp = (int **)malloc(sizeof(int *)*LyrInfo.SigNumElemt);
		for (j = 0; j<LyrInfo.SigNumElemt; j++)
		{
			i2pAmpCtrEmp[j] = (int *)malloc(sizeof(int)*(LyrInfo.SigAphSize[j] + 1));
		}
		for (j = 0; j<LyrInfo.SigNumElemt; j++)
		{
			for (k = 0; k <= LyrInfo.SigAphSize[j]; k++)
			{
				i2pAmpCtrEmp[j][k] = 0;
			}
		}

		for (i = 0; i<2; i++)
			i3pAmpCtxElemtCtr[i] = (int **)malloc(sizeof(int *)*LyrInfo.SigNumElemt);
		for (i = 0; i<2; i++)
			for (j = 0; j<LyrInfo.SigNumElemt; j++)
				i3pAmpCtxElemtCtr[i][j] = (int *)malloc(sizeof(int)*(LyrInfo.SigAphSize[j] - BIT_LOC_SIG + 1));
		for (i = 0; i<2; i++)
			for (j = 0; j<LyrInfo.SigNumElemt; j++)
				for (k = 0; k<(LyrInfo.SigAphSize[j] - BIT_LOC_SIG + 1); k++)
					i3pAmpCtxElemtCtr[i][j][k] = 0;
#endif
	}
}


/*
* FUNCTION: fnCtrIni
* initialize the counters for SDQ algorithm
*/
void fnCtrIni()
{
	int i, j, k, l, m, iFreqRow, iFreqCol;

	//    for (i=0; i<NUM_CTX_AZF; i++)
	//    {
	//        ipAZFCtxMainCtr[i] = 2;
	//        
	//        for (j=0; j<2; j++)
	//        {
	//            i2pAZFCtxCtr[i][j] = 1;
	//        }
	//    }

	for (i = 0; i<2; i++)
	{
		for (j = 0; j<NUM_CTX_SIG; j++)
		{
			for (k = 0; k<2; k++)
			{
				for (l = 0; l<7; l++)
				{
					for (m = 0; m<2; m++)
					{
						i5pSigCtxCtr[i][j][k][l][m] = 1;
					}
				}
			}
		}

		for (j = 0; j<7; j++)
		{
			for (k = 0; k<MAXIMUM(1, LyrInfo.SigAphSize[LyrInfo.CtxLyrIdxMaxL[j]] - BIT_LOC_SIG + 1); k++)
			{
				i3pAmpCtxCtr[i][j][k] = 1;
			}
		}
	}

}


/*
* FUNCTION fnACItrfc:
* arithmetic coding for the symbols under corresponding context models
*/
void fnACItrfc(Context *input_ctx, int input)
{
#if BAC_EN
	fnBACEncode(GACpEngine, input_ctx, input, GbspOutputBuffer);
#else
	fnACEncode(GACpEngine, input_ctx, input, GbspOutputBuffer);
#endif
}

void fnACItrfcCoeffDC(Context *input_ctx, int input)
{
	fnACEncode(GACpEngine, input_ctx, input, GbspOutputBuffer);
}

void fnACItrfcVoid(Context *input_ctx, int input)
{
	fnACEncodeVoid(GACpEngine, input_ctx, input, GbspOutputBuffer);
}

void fnACmAlpItrfc(Context *input_ctx, int input_sym, int Alp)
{
	fnACEncodeMAlp(GACpEngine, input_ctx, input_sym, GbspOutputBuffer, Alp);
}

/*
* FUNCTION fnACItrfcByps:
* arithmetic coding using by-pass mode for symbols without context models
*/
void fnACItrfcByps(int input)
{
	fnACEncodebyPass(GACpEngine, input, GbspOutputBuffer);
}


/*
*  Function: fnGetCtxAZF
*  Context modeling for AZF based on the neighboring blocks (top-left, top, top-right, left;
*  i.e. Block index: 0-[i-1][j-1], 1-[i-1][j], 2-[i-1][j+1], 3-[i][j-1])
*/
Context *fnGetCtxAZF(int iRow, int iCol, int iCurAZF)
{
	int iCtxNum = 0, iExtCtxNum;

	//using already encoded 4 AZF neighbors
	if (TAR_PSNR <= 32)
	{
		if (iRow == 0 && iCol == 0)
			iCtxNum = 0;

		else if (iRow == 0 && iCol != 0)
		{
			iCtxNum = drgiCBP[iRow][iCol - 1];
		}

		else if (iRow != 0 && iCol == 0)
		{
			iCtxNum = drgiCBP[iRow - 1][iCol] + drgiCBP[iRow - 1][iCol + 1];
		}

		else if (iRow != 0 && iCol == (COLS / 8 - 1))
		{
			iCtxNum = drgiCBP[iRow - 1][iCol - 1] + drgiCBP[iRow - 1][iCol] + drgiCBP[iRow][iCol - 1];
		}

		else
		{
			iCtxNum = drgiCBP[iRow - 1][iCol - 1] + drgiCBP[iRow - 1][iCol] + drgiCBP[iRow - 1][iCol + 1] + drgiCBP[iRow][iCol - 1];
		}
	}

	else
	{
		if (iRow == 0 && iCol == 0)
			iCtxNum = 0;

		else if (iRow == 0 && iCol != 0)
		{
			iCtxNum = drgiCBP[iRow][iCol - 1];
			if (iCol>1)
				iCtxNum += drgiCBP[iRow][iCol - 2];
		}

		else if (iRow != 0 && iCol == 0)
		{
			iCtxNum = drgiCBP[iRow - 1][iCol] + drgiCBP[iRow - 1][iCol + 1];
			if (iRow>1)
				iCtxNum += drgiCBP[iRow - 2][iCol];
		}

		else if (iRow != 0 && iCol == (COLS / 8 - 1))
		{
			iCtxNum = drgiCBP[iRow - 1][iCol - 1] + drgiCBP[iRow - 1][iCol] + drgiCBP[iRow][iCol - 1];
			if (iRow>1)
				iCtxNum += drgiCBP[iRow - 2][iCol];
		}

		else
		{
			iCtxNum = drgiCBP[iRow - 1][iCol - 1] + drgiCBP[iRow - 1][iCol] + drgiCBP[iRow - 1][iCol + 1] + drgiCBP[iRow][iCol - 1];
			if (iRow>1)
				iCtxNum += drgiCBP[iRow - 2][iCol];
			if (iCol>1)
				iCtxNum += drgiCBP[iRow][iCol - 2];
		}
	}

	//using already encoded EXTF info of 2 neighbors
	if (iRow<(ROWS / 8 - 1))
	{
		if (iCol == (COLS / 8 - 1))
			iExtCtxNum = drgiExtF[iRow + 1][iCol];

		else
			iExtCtxNum = drgiExtF[iRow + 1][iCol] + drgiExtF[iRow][iCol + 1];//+ drgiExtF[iRow+1][iCol+1];
	}

	else
	{
		if (iCol<(COLS / 8 - 1))
			iExtCtxNum = drgiExtF[iRow][iCol + 1];
		else
			iExtCtxNum = 0;
	}

	iCtxNum += iExtCtxNum;
	iCtxNum = MINIMUM(NUM_CTX_AZF - 1, iCtxNum);

	return &(AZFCtx[iCtxNum]);
}

/*
*  Function: fnCalCtxAZF
*  Collect statistics of AZF for SDQ algorithm
*/
//int fnCalCtxAZF(int iRow, int iCol)
//{
//	int iCtxNum=0, iExtCtxNum;
//    
//    //using already encoded 4 AZF neighbors
//    if (TAR_PSNR <= 32)
//    {
//        if (iRow==0 && iCol==0)
//            iCtxNum = 0;
//        
//        else if (iRow==0 && iCol!=0)
//        {
//            iCtxNum = drgiCBP[iRow][iCol-1];
//        }
//        
//        else if (iRow!=0 && iCol==0)
//        {
//            iCtxNum = drgiCBP[iRow-1][iCol] + drgiCBP[iRow-1][iCol+1];
//        }   
//        
//        else if (iRow!=0 && iCol==(COLS/8-1))
//        {
//            iCtxNum = drgiCBP[iRow-1][iCol-1] + drgiCBP[iRow-1][iCol] + drgiCBP[iRow][iCol-1];
//        }
//        
//        else
//        {
//            iCtxNum = drgiCBP[iRow-1][iCol-1] + drgiCBP[iRow-1][iCol] + drgiCBP[iRow-1][iCol+1] + drgiCBP[iRow][iCol-1];
//        }
//    }
//    
//    else
//    {
//        if (iRow==0 && iCol==0)
//            iCtxNum = 0;
//        
//        else if (iRow==0 && iCol!=0)
//        {
//            iCtxNum = drgiCBP[iRow][iCol-1];
//            if (iCol>1)
//                iCtxNum += drgiCBP[iRow][iCol-2];
//        }
//        
//        else if (iRow!=0 && iCol==0)
//        {
//            iCtxNum = drgiCBP[iRow-1][iCol] + drgiCBP[iRow-1][iCol+1];
//            if (iRow>1)
//                iCtxNum += drgiCBP[iRow-2][iCol];
//        }   
//        
//        else if (iRow!=0 && iCol==(COLS/8-1))
//        {
//            iCtxNum = drgiCBP[iRow-1][iCol-1] + drgiCBP[iRow-1][iCol] + drgiCBP[iRow][iCol-1];
//            if (iRow>1)
//                iCtxNum += drgiCBP[iRow-2][iCol];
//        }
//        
//        else
//        {
//            iCtxNum = drgiCBP[iRow-1][iCol-1] + drgiCBP[iRow-1][iCol] + drgiCBP[iRow-1][iCol+1] + drgiCBP[iRow][iCol-1];
//            if (iRow>1)
//                iCtxNum += drgiCBP[iRow-2][iCol];
//            if (iCol>1)
//                iCtxNum += drgiCBP[iRow][iCol-2];
//        }
//    }
//    
//    //using already encoded EXTF info of 2 neighbors
//    if (iRow<(ROWS/8-1))
//    {
//        if (iCol==(COLS/8-1))
//            iExtCtxNum = drgiExtF[iRow+1][iCol];
//        
//        else
//            iExtCtxNum = drgiExtF[iRow+1][iCol] + drgiExtF[iRow][iCol+1];//+ drgiExtF[iRow+1][iCol+1];
//    }
//    
//    else
//    {
//        if (iCol<(COLS/8-1))
//            iExtCtxNum = drgiExtF[iRow][iCol+1];
//        else
//            iExtCtxNum = 0;
//    }
//    
//    iCtxNum += iExtCtxNum;
//    iCtxNum = MINIMUM(NUM_CTX_AZF-1, iCtxNum);
//        
//    return iCtxNum;
//}

/*
*  Function: fnGetCtxExtF
*  Context modeling for ExtF based on the neighboring blocks (top-left, top, top-right, left;
*  i.e. Block index: 0-[i-1][j-1], 1-[i-1][j], 2-[i-1][j+1], 3-[i][j-1])
*/
Context *fnGetCtxExtF(int iRow, int iCol)
{
	int iCtxNum = 0;
	int SchmExtF = 1;

	if (SchmExtF == 0)//6 context models using 6 neighboring blocks (with context quantization and context merging)
	{
		if (iRow == 0 && iCol == 0)
			iCtxNum = 0;

		else if (iRow == 0 && iCol>0)
		{
			iCtxNum = drgiExtF[iRow][iCol - 1];

			if (iCol > 1)
				iCtxNum = iCtxNum + drgiExtF[iRow][iCol - 2];
		}

		else if (iRow>0 && iCol == 0)
		{
			iCtxNum = drgiExtF[iRow - 1][iCol] + drgiExtF[iRow - 1][iCol + 1];

			if (iRow > 1)
				iCtxNum = iCtxNum + drgiExtF[iRow - 2][iCol];
		}

		else if (iRow>0 && iCol == (COLS / 8 - 1))
		{
			iCtxNum = drgiExtF[iRow - 1][iCol - 1] + drgiExtF[iRow - 1][iCol] + drgiExtF[iRow][iCol - 1] + drgiExtF[iRow][iCol - 2];

			if (iRow > 1)
				iCtxNum = iCtxNum + drgiExtF[iRow - 2][iCol];
		}

		else
		{
			iCtxNum = drgiExtF[iRow - 1][iCol] + drgiExtF[iRow][iCol - 1] + drgiExtF[iRow - 1][iCol - 1] + drgiExtF[iRow - 1][iCol + 1];

			if (iRow > 1)
				iCtxNum = iCtxNum + drgiExtF[iRow - 2][iCol];

			if (iCol > 1)
				iCtxNum = iCtxNum + drgiExtF[iRow][iCol - 2];
		}

		if (iCtxNum == 5 || iCtxNum == 6)
			iCtxNum = 5;
	}

	else//16 context models using 4 neighboring blocks (no context quantization)
	{
		if (iRow == 0 && iCol == 0)
			iCtxNum = 0;

		else if (iRow == 0 && iCol != 0)
		{
			iCtxNum = 8 * drgiExtF[iRow][iCol - 1];
		}

		else if (iRow>0 && iCol == 0)
		{
			iCtxNum = 2 * drgiExtF[iRow - 1][iCol] + 4 * drgiExtF[iRow - 1][iCol + 1];
		}

		else if (iRow>0 && iCol == (COLS / 8 - 1))
		{
			iCtxNum = drgiExtF[iRow - 1][iCol - 1] + 2 * drgiExtF[iRow - 1][iCol] + 8 * drgiExtF[iRow][iCol - 1];
		}

		else
		{
			iCtxNum = drgiExtF[iRow - 1][iCol - 1] + 2 * drgiExtF[iRow - 1][iCol] + 4 * drgiExtF[iRow - 1][iCol + 1] + 8 * drgiExtF[iRow][iCol - 1];
		}
	}

	return &(ExtFCtx[iCtxNum]);
}

Context *fnGetCtxSigExt(int iAcumSig, int iElemtIdx, int iRow, int iCol, int iBlkRow, int iBlkCol, int iBitPos, int iCurSigExt)
{
	int iCtxNumBlkExt = 0, iCtxNumBlkNeibr = 0, iCtxNumSubndNeibr = 0, iCtxNum, i, j, iRowZZplus1, iColZZplus1, iRowZZMinus1, iColZZMinus1, iElemtIdxTmp;

	//explore the correlation among 5 block neighbors of the current encoding coeff
	i = iBlkRow;
	j = iBlkCol;

#if SIG_SCH
	if (i == 7 && j<7)
		iCtxNumBlkNeibr = vrgiSig[iRow][iCol][i][j + 1][iBitPos];
	else if (i == 7 && j == 7)
		iCtxNumBlkNeibr = 0;
	else if (j == 7)
		iCtxNumBlkNeibr = vrgiSig[iRow][iCol][i + 1][j][iBitPos];
	else
		iCtxNumBlkNeibr = vrgiSig[iRow][iCol][i + 1][j + 1][iBitPos] + vrgiSig[iRow][iCol][i][j + 1][iBitPos] + vrgiSig[iRow][iCol][i + 1][j][iBitPos];

	if (iElemtIdx < (LyrInfo.SigNumElemt - 1))
	{
		iRowZZplus1 = LyrInfo.SigScanRowIdx[iElemtIdx + 1];
		iColZZplus1 = LyrInfo.SigScanColIdx[iElemtIdx + 1];

		if ((iRowZZplus1 + iColZZplus1) == (i + j))
		{
#if IMMED_NEIBR
			if (iRowZZplus1 == i + 1 || iRowZZplus1 == i - 1)
				iCtxNumBlkNeibr += vrgiSig[iRow][iCol][iRowZZplus1][iColZZplus1][iBitPos];
#else
			iCtxNumBlkNeibr += vrgiSig[iRow][iCol][iRowZZplus1][iColZZplus1][iBitPos];
#endif
		}
	}

#else
	if (i == 7)
	{
		if (j == 7)
			iCtxNumBlkNeibr = 0;

		else if (j == 6)
			iCtxNumBlkNeibr = vrgiSig[iRow][iCol][7][7][iBitPos];

		else
			iCtxNumBlkNeibr = vrgiSig[iRow][iCol][i][j + 1][iBitPos] + vrgiSig[iRow][iCol][i][j + 2][iBitPos];
	}

	else if (i == 6)
	{
		if (j == 7)
			iCtxNumBlkNeibr = vrgiSig[iRow][iCol][7][7][iBitPos];

		else if (j == 6)
			iCtxNumBlkNeibr = vrgiSig[iRow][iCol][i][j + 1][iBitPos] + vrgiSig[iRow][iCol][i + 1][j][iBitPos] + vrgiSig[iRow][iCol][7][7][iBitPos];

		else
			iCtxNumBlkNeibr = vrgiSig[iRow][iCol][i][j + 1][iBitPos] + vrgiSig[iRow][iCol][i + 1][j][iBitPos] + vrgiSig[iRow][iCol][i + 1][j + 1][iBitPos] + vrgiSig[iRow][iCol][i][j + 2][iBitPos];
	}

	else if (j == 7)
	{
		iCtxNumBlkNeibr = vrgiSig[iRow][iCol][i + 1][j][iBitPos] + vrgiSig[iRow][iCol][i + 2][j][iBitPos];
	}

	else if (j == 6)
	{
		iCtxNumBlkNeibr = vrgiSig[iRow][iCol][i + 1][j][iBitPos] + vrgiSig[iRow][iCol][i + 2][j][iBitPos] + vrgiSig[iRow][iCol][i][j + 1][iBitPos] + vrgiSig[iRow][iCol][i + 1][j + 1][iBitPos];
	}

	else
	{
		iCtxNumBlkNeibr = vrgiSig[iRow][iCol][i + 1][j][iBitPos] + vrgiSig[iRow][iCol][i + 2][j][iBitPos] + vrgiSig[iRow][iCol][i][j + 1][iBitPos] + vrgiSig[iRow][iCol][i + 1][j + 1][iBitPos] + vrgiSig[iRow][iCol][i][j + 2][iBitPos];
	}

#if JUNE26_SIG_LYR_NEIBER
	iElemtIdxTmp = iElemtIdx;
	while (iElemtIdxTmp < (LyrInfo.SigNumElemt - 1))
	{
		iRowZZplus1 = LyrInfo.SigScanRowIdx[iElemtIdxTmp + 1];
		iColZZplus1 = LyrInfo.SigScanColIdx[iElemtIdxTmp + 1];

		if ((iRowZZplus1 + iColZZplus1) == (i + j))
		{
			if (vrgiExt[iRow][iCol][iRowZZplus1][iColZZplus1][0] == 0)
			{
				iCtxNumBlkNeibr += vrgiSig[iRow][iCol][iRowZZplus1][iColZZplus1][iBitPos];
				break;
			}
		}
		else
		{
			break;
		}
		iElemtIdxTmp++;
	}
#else
	if (iElemtIdx < (LyrInfo.SigNumElemt - 1))
	{
		iRowZZplus1 = LyrInfo.SigScanRowIdx[iElemtIdx + 1];
		iColZZplus1 = LyrInfo.SigScanColIdx[iElemtIdx + 1];

		if ((iRowZZplus1 + iColZZplus1) == (i + j))
		{
#if IMMED_NEIBR
			if (iRowZZplus1 == i + 1 || iRowZZplus1 == i - 1)
				iCtxNumBlkNeibr += vrgiSig[iRow][iCol][iRowZZplus1][iColZZplus1][iBitPos];
#else
			iCtxNumBlkNeibr += vrgiSig[iRow][iCol][iRowZZplus1][iColZZplus1][iBitPos];
#endif
		}
	}
#endif
#endif

	//    if (iBitPos>=0)
	//    {
	//explore the correlation among 4 future neighbors of the current encoding coeff using Ext info
	if (i == 0)
	{
		if (j>0)
			iCtxNumBlkExt = vrgiExt[iRow][iCol][i][j - 1][0];
		else
			iCtxNumBlkExt = 0;
	}

	else if (j == 1)
		iCtxNumBlkExt = vrgiExt[iRow][iCol][i - 1][j][0];

	else
		iCtxNumBlkExt = vrgiExt[iRow][iCol][i - 1][j][0] + vrgiExt[iRow][iCol][i][j - 1][0] + vrgiExt[iRow][iCol][i - 1][j - 1][0];

	if (iElemtIdx > 0)
	{
		iRowZZMinus1 = LyrInfo.SigScanRowIdx[iElemtIdx - 1];
		iColZZMinus1 = LyrInfo.SigScanColIdx[iElemtIdx - 1];

		if ((iRowZZMinus1 + iColZZMinus1) == (i + j))
		{
#if IMMED_NEIBR
			if (iRowZZMinus1 == i + 1 || iRowZZMinus1 == i - 1)
				iCtxNumBlkExt += vrgiExt[iRow][iCol][iRowZZMinus1][iColZZMinus1][0];
#else
			iCtxNumBlkExt += vrgiExt[iRow][iCol][iRowZZMinus1][iColZZMinus1][0];
#endif
		}
	}
	//    }
	//    else
	//    {
	//        //explore the correlation among 4 future neighbors of the current encoding coeff using the previous Sig bit-plane info
	//        if (i==0)
	//        {
	//            if (j>0)
	//                iCtxNumBlkExt = vrgiSig[iRow][iCol][i][j-1][iBitPos-1];
	//            else
	//                iCtxNumBlkExt = 0;
	//        }
	//        
	//        else if (j==1)
	//            iCtxNumBlkExt = vrgiSig[iRow][iCol][i-1][j][iBitPos-1];
	//        
	//        else
	//            iCtxNumBlkExt = vrgiSig[iRow][iCol][i-1][j][iBitPos-1]+vrgiSig[iRow][iCol][i][j-1][iBitPos-1]+vrgiSig[iRow][iCol][i-1][j-1][iBitPos-1];
	//        
	//        if (iElemtIdx > 0)
	//        {    
	//            iRowZZMinus1 = LyrInfo.SigScanRowIdx[iElemtIdx-1];
	//            iColZZMinus1 = LyrInfo.SigScanColIdx[iElemtIdx-1];
	//            
	//            if ((iRowZZMinus1+iColZZMinus1) == (i+j))
	//            {
	//#if IMMED_NEIBR
	//                if (iRowZZMinus1 == i+1 || iRowZZMinus1 == i-1)   
	//                    iCtxNumBlkExt += vrgiSig[iRow][iCol][iRowZZMinus1][iColZZMinus1][iBitPos-1];
	//#else
	//                iCtxNumBlkExt += vrgiSig[iRow][iCol][iRowZZMinus1][iColZZMinus1][iBitPos-1];
	//#endif
	//            }
	//            
	//        }
	//    }

	iCtxNum = iCtxNumBlkExt + iCtxNumBlkNeibr;

#if INTER_BLK_SIG_HELP
	//explore the correlation among 6 subband neighbors of the current encoding coeff
	if (iRow == 0)
	{
		if (iCol == 0)
			iCtxNumSubndNeibr = 0;

		else
		{
			iCtxNumSubndNeibr = vrgiSig[iRow][iCol - 1][iBlkRow][iBlkCol][iBitPos];
			if (iCol>1)
			{
				iCtxNumSubndNeibr += vrgiSig[iRow][iCol - 2][iBlkRow][iBlkCol][iBitPos];
			}
		}
	}

	else if (iCol == 0)
	{
		iCtxNumSubndNeibr = vrgiSig[iRow - 1][iCol][iBlkRow][iBlkCol][iBitPos] + vrgiSig[iRow - 1][iCol + 1][iBlkRow][iBlkCol][iBitPos];
		if (iRow>1)
		{
			iCtxNumSubndNeibr += vrgiSig[iRow - 2][iCol][iBlkRow][iBlkCol][iBitPos];
		}
	}

	else if (iCol == (COLS / 8 - 1))
	{
		iCtxNumSubndNeibr = vrgiSig[iRow - 1][iCol][iBlkRow][iBlkCol][iBitPos] + vrgiSig[iRow][iCol - 1][iBlkRow][iBlkCol][iBitPos] + vrgiSig[iRow - 1][iCol - 1][iBlkRow][iBlkCol][iBitPos] + vrgiSig[iRow][iCol - 2][iBlkRow][iBlkCol][iBitPos];
		if (iRow>1)
			iCtxNumSubndNeibr += vrgiSig[iRow - 2][iCol][iBlkRow][iBlkCol][iBitPos];
	}

	else
	{
		iCtxNumSubndNeibr = vrgiSig[iRow - 1][iCol][iBlkRow][iBlkCol][iBitPos] + vrgiSig[iRow][iCol - 1][iBlkRow][iBlkCol][iBitPos] + vrgiSig[iRow - 1][iCol - 1][iBlkRow][iBlkCol][iBitPos] + vrgiSig[iRow - 1][iCol + 1][iBlkRow][iBlkCol][iBitPos];
		if (iRow>1)
			iCtxNumSubndNeibr += vrgiSig[iRow - 2][iCol][iBlkRow][iBlkCol][iBitPos];

		if (iCol>1)
			iCtxNumSubndNeibr += vrgiSig[iRow][iCol - 2][iBlkRow][iBlkCol][iBitPos];
	}

	//    iCtxNumSubndNeibr=iCtxNumSubndNeibr>2?1:0;
	//    iCtxNumBlkNeibr=MINIMUM(iCtxNumBlkNeibr,3);
	iCtxNum += iCtxNumSubndNeibr;
	if (iBitPos == 0)
		iCtxNum = MINIMUM(iCtxNum, 10);
	else
		iCtxNum = MINIMUM(iCtxNum, 8);
#endif

#if BIN_SCH
	if (iBitPos < 8)
	{
#if SIG_MERG
		return &(SigCtx[1][iCtxNum][iBitPos>1 ? 1 : iBitPos][LyrInfo.CtxLyr[iElemtIdx]]);
#else
#if SIG_SEPRT
		if (iElemtIdx < LyrInfo.Acum1Pos)
			return &(SigCtx[1][iCtxNum][iBitPos][iElemtIdx]);
		else
			return &(SigCtx[1][MINIMUM(iAcumSig, 3)][iBitPos][iElemtIdx]);
#else
#if LST2LYR_MERG
		return &(SigCtx[1][iCtxNum][iBitPos>1 ? 1 : iBitPos][LyrInfo.CtxLyrSig[iElemtIdx]]);
#else
		return &(SigCtx[1][iCtxNum][iBitPos>1 ? 1 : iBitPos][iElemtIdx]);
#endif 
#endif
#endif
	}

	else
	{
		if (BIT_LOC_SIG > 2)
		{
#if SIG_MERG
			return &(SigCtx[1][iCtxNum][1][LyrInfo.CtxLyr[iElemtIdx]]);
#else
#if SIG_SEPRT
			if (iElemtIdx < LyrInfo.Acum1Pos)
				return &(SigCtx[1][iCtxNum][1][iElemtIdx]);
			else
				return &(SigCtx[1][MINIMUM(iAcumSig, 3)][1][iElemtIdx]);
#else
			return &(SigCtx[1][iCtxNum][iBitPos][iElemtIdx]);
#endif
#endif
		}
	}

#else
	return &(SigCtx[1][iCtxNum][0][iElemtIdx]);
#endif

}

Context *fnGetCtxSig(int iAcumSig, int iElemtIdx, int iRow, int iCol, int iBlkRow, int iBlkCol, int iBitPos, int iCurSig)
{
	int iCtxNumBlkNeibr = 0, iCtxNumSubndNeibr = 0, iCtxNum, i, j, iRowZZplus1, iColZZplus1, iRowZZMinus1, iColZZMinus1, iCtxNumBlkExt, iElemtIdxTmp;

	//explore the correlation among 5 block neighbors of the current encoding coeff
	i = iBlkRow;
	j = iBlkCol;

#if SIG_SCH
	if (i == 7 && j<7)
		iCtxNumBlkNeibr = vrgiSig[iRow][iCol][i][j + 1][iBitPos];
	else if (i == 7 && j == 7)
		iCtxNumBlkNeibr = 0;
	else if (j == 7)
		iCtxNumBlkNeibr = vrgiSig[iRow][iCol][i + 1][j][iBitPos];
	else
		iCtxNumBlkNeibr = vrgiSig[iRow][iCol][i + 1][j + 1][iBitPos] + vrgiSig[iRow][iCol][i][j + 1][iBitPos] + vrgiSig[iRow][iCol][i + 1][j][iBitPos];

	if (iElemtIdx < (LyrInfo.SigNumElemt - 1))
	{
		iRowZZplus1 = LyrInfo.SigScanRowIdx[iElemtIdx + 1];
		iColZZplus1 = LyrInfo.SigScanColIdx[iElemtIdx + 1];

		if ((iRowZZplus1 + iColZZplus1) == (i + j))
		{
#if IMMED_NEIBR
			if (iRowZZplus1 == i + 1 || iRowZZplus1 == i - 1)
				iCtxNumBlkNeibr += vrgiSig[iRow][iCol][iRowZZplus1][iColZZplus1][iBitPos];
#else
			iCtxNumBlkNeibr += vrgiSig[iRow][iCol][iRowZZplus1][iColZZplus1][iBitPos];
#endif
		}
	}

#else

	if (i == 7)
	{
		if (j == 7)
			iCtxNumBlkNeibr = 0;

		else if (j == 6)
			iCtxNumBlkNeibr = vrgiSig[iRow][iCol][7][7][iBitPos];

		else
			iCtxNumBlkNeibr = vrgiSig[iRow][iCol][i][j + 1][iBitPos] + vrgiSig[iRow][iCol][i][j + 2][iBitPos];
	}

	else if (i == 6)
	{
		if (j == 7)
			iCtxNumBlkNeibr = vrgiSig[iRow][iCol][7][7][iBitPos];

		else if (j == 6)
			iCtxNumBlkNeibr = vrgiSig[iRow][iCol][i][j + 1][iBitPos] + vrgiSig[iRow][iCol][i + 1][j][iBitPos] + vrgiSig[iRow][iCol][7][7][iBitPos];

		else
			iCtxNumBlkNeibr = vrgiSig[iRow][iCol][i][j + 1][iBitPos] + vrgiSig[iRow][iCol][i + 1][j][iBitPos] + vrgiSig[iRow][iCol][i + 1][j + 1][iBitPos] + vrgiSig[iRow][iCol][i][j + 2][iBitPos];
	}

	else if (j == 7)
	{
		iCtxNumBlkNeibr = vrgiSig[iRow][iCol][i + 1][j][iBitPos] + vrgiSig[iRow][iCol][i + 2][j][iBitPos];
	}

	else if (j == 6)
	{
		iCtxNumBlkNeibr = vrgiSig[iRow][iCol][i + 1][j][iBitPos] + vrgiSig[iRow][iCol][i + 2][j][iBitPos] + vrgiSig[iRow][iCol][i][j + 1][iBitPos] + vrgiSig[iRow][iCol][i + 1][j + 1][iBitPos];
	}

	else
	{
		iCtxNumBlkNeibr = vrgiSig[iRow][iCol][i + 1][j][iBitPos] + vrgiSig[iRow][iCol][i + 2][j][iBitPos] + vrgiSig[iRow][iCol][i][j + 1][iBitPos] + vrgiSig[iRow][iCol][i + 1][j + 1][iBitPos] + vrgiSig[iRow][iCol][i][j + 2][iBitPos];
	}

#if JUNE26_SIG_LYR_NEIBER
	iElemtIdxTmp = iElemtIdx;
	while (iElemtIdxTmp < (LyrInfo.SigNumElemt - 1))
	{
		iRowZZplus1 = LyrInfo.SigScanRowIdx[iElemtIdxTmp + 1];
		iColZZplus1 = LyrInfo.SigScanColIdx[iElemtIdxTmp + 1];

		if ((iRowZZplus1 + iColZZplus1) == (i + j))
		{
			if (vrgiExt[iRow][iCol][iRowZZplus1][iColZZplus1][0] == 0)
			{
				iCtxNumBlkNeibr += vrgiSig[iRow][iCol][iRowZZplus1][iColZZplus1][iBitPos];
				break;
			}
		}
		else
		{
			break;
		}
		iElemtIdxTmp++;
	}
#else
	if (iElemtIdx < (LyrInfo.SigNumElemt - 1))
	{
		iRowZZplus1 = LyrInfo.SigScanRowIdx[iElemtIdx + 1];
		iColZZplus1 = LyrInfo.SigScanColIdx[iElemtIdx + 1];

		if ((iRowZZplus1 + iColZZplus1) == (i + j))
		{
#if IMMED_NEIBR
			if (iRowZZplus1 == i + 1 || iRowZZplus1 == i - 1)
				iCtxNumBlkNeibr += vrgiSig[iRow][iCol][iRowZZplus1][iColZZplus1][iBitPos];
#else
			iCtxNumBlkNeibr += vrgiSig[iRow][iCol][iRowZZplus1][iColZZplus1][iBitPos];
#endif
		}
	}
#endif
#endif

	//    if (iBitPos>0)
	//    {
	//        //explore the correlation among 4 future neighbors of the current encoding coeff using the previous Sig bit-plane info
	//        if (i==0)
	//        {
	//            if (j>0)
	//                iCtxNumBlkExt = vrgiSig[iRow][iCol][i][j-1][iBitPos-1];
	//            else
	//                iCtxNumBlkExt = 0;
	//        }
	//        
	//        else if (j==1)
	//            iCtxNumBlkExt = vrgiSig[iRow][iCol][i-1][j][iBitPos-1];
	//        
	//        else
	//            iCtxNumBlkExt = vrgiSig[iRow][iCol][i-1][j][iBitPos-1]+vrgiSig[iRow][iCol][i][j-1][iBitPos-1]+vrgiSig[iRow][iCol][i-1][j-1][iBitPos-1];
	//        
	//        if (iElemtIdx > 0)
	//        {    
	//            iRowZZMinus1 = LyrInfo.SigScanRowIdx[iElemtIdx-1];
	//            iColZZMinus1 = LyrInfo.SigScanColIdx[iElemtIdx-1];
	//            
	//            if ((iRowZZMinus1+iColZZMinus1) == (i+j))
	//            {
	//#if IMMED_NEIBR
	//                if (iRowZZMinus1 == i+1 || iRowZZMinus1 == i-1)   
	//                    iCtxNumBlkExt += vrgiSig[iRow][iCol][iRowZZMinus1][iColZZMinus1][iBitPos-1];
	//#else
	//                iCtxNumBlkExt += vrgiSig[iRow][iCol][iRowZZMinus1][iColZZMinus1][iBitPos-1];
	//#endif
	//            }
	//            
	//        }
	//    }
	//    
	//    iCtxNum = MINIMUM(iCtxNumBlkExt+iCtxNumBlkNeibr, 4);


#if 0
	//explore the correlation among 6 subband neighbors of the current encoding coeff
	if (iRow == 0)
	{
		if (iCol == 0)
			iCtxNumSubndNeibr = 0;

		else
		{
			iCtxNumSubndNeibr = vrgiExt[iRow][iCol - 1][iBlkRow][iBlkCol][0];
			if (iCol>1)
			{
				iCtxNumSubndNeibr += vrgiExt[iRow][iCol - 2][iBlkRow][iBlkCol][0];
			}
		}
	}

	else if (iCol == 0)
	{
		iCtxNumSubndNeibr = vrgiExt[iRow - 1][iCol][iBlkRow][iBlkCol][0] + vrgiExt[iRow - 1][iCol + 1][iBlkRow][iBlkCol][0];
		if (iRow>1)
		{
			iCtxNumSubndNeibr += vrgiExt[iRow - 2][iCol][iBlkRow][iBlkCol][0];
		}
	}

	else if (iCol == (COLS / 8 - 1))
	{
		iCtxNumSubndNeibr = vrgiExt[iRow - 1][iCol][iBlkRow][iBlkCol][0] + vrgiExt[iRow][iCol - 1][iBlkRow][iBlkCol][0] + vrgiExt[iRow - 1][iCol - 1][iBlkRow][iBlkCol][0] + vrgiExt[iRow][iCol - 2][iBlkRow][iBlkCol][0];
		if (iRow>1)
			iCtxNumSubndNeibr += vrgiExt[iRow - 2][iCol][iBlkRow][iBlkCol][0];
	}

	else
	{
		iCtxNumSubndNeibr = vrgiExt[iRow - 1][iCol][iBlkRow][iBlkCol][0] + vrgiExt[iRow][iCol - 1][iBlkRow][iBlkCol][0] + vrgiExt[iRow - 1][iCol - 1][iBlkRow][iBlkCol][0] + vrgiExt[iRow - 1][iCol + 1][iBlkRow][iBlkCol][0];
		if (iRow>1)
			iCtxNumSubndNeibr += vrgiExt[iRow - 2][iCol][iBlkRow][iBlkCol][0];

		if (iCol>1)
			iCtxNumSubndNeibr += vrgiExt[iRow][iCol - 2][iBlkRow][iBlkCol][0];
	}

	//    iCtxNumSubndNeibr=iCtxNumSubndNeibr>2?1:0;
	//    iCtxNumBlkNeibr=MINIMUM(iCtxNumBlkNeibr,3);
	//    iCtxNum = iCtxNumBlkNeibr*2+iCtxNumSubndNeibr;
	iCtxNum = iCtxNumBlkNeibr + iCtxNumSubndNeibr;
	iCtxNum = MINIMUM(iCtxNum, 4);
	//#else
	//    iCtxNum = MINIMUM(iCtxNumBlkNeibr,4);
#endif

#if INTER_BLK_SIG_HELP
	//explore the correlation among 6 subband neighbors of the current encoding coeff
	if (iRow == 0)
	{
		if (iCol == 0)
			iCtxNumSubndNeibr = 0;

		else
		{
			iCtxNumSubndNeibr = vrgiSig[iRow][iCol - 1][iBlkRow][iBlkCol][iBitPos];
			if (iCol>1)
			{
				iCtxNumSubndNeibr += vrgiSig[iRow][iCol - 2][iBlkRow][iBlkCol][iBitPos];
			}
		}
	}

	else if (iCol == 0)
	{
		iCtxNumSubndNeibr = vrgiSig[iRow - 1][iCol][iBlkRow][iBlkCol][iBitPos] + vrgiSig[iRow - 1][iCol + 1][iBlkRow][iBlkCol][iBitPos];
		if (iRow>1)
		{
			iCtxNumSubndNeibr += vrgiSig[iRow - 2][iCol][iBlkRow][iBlkCol][iBitPos];
		}
	}

	else if (iCol == (COLS / 8 - 1))
	{
		iCtxNumSubndNeibr = vrgiSig[iRow - 1][iCol][iBlkRow][iBlkCol][iBitPos] + vrgiSig[iRow][iCol - 1][iBlkRow][iBlkCol][iBitPos] + vrgiSig[iRow - 1][iCol - 1][iBlkRow][iBlkCol][iBitPos] + vrgiSig[iRow][iCol - 2][iBlkRow][iBlkCol][iBitPos];
		if (iRow>1)
			iCtxNumSubndNeibr += vrgiSig[iRow - 2][iCol][iBlkRow][iBlkCol][iBitPos];
	}

	else
	{
		iCtxNumSubndNeibr = vrgiSig[iRow - 1][iCol][iBlkRow][iBlkCol][iBitPos] + vrgiSig[iRow][iCol - 1][iBlkRow][iBlkCol][iBitPos] + vrgiSig[iRow - 1][iCol - 1][iBlkRow][iBlkCol][iBitPos] + vrgiSig[iRow - 1][iCol + 1][iBlkRow][iBlkCol][iBitPos];
		if (iRow>1)
			iCtxNumSubndNeibr += vrgiSig[iRow - 2][iCol][iBlkRow][iBlkCol][iBitPos];

		if (iCol>1)
			iCtxNumSubndNeibr += vrgiSig[iRow][iCol - 2][iBlkRow][iBlkCol][iBitPos];
	}

	//    iCtxNumSubndNeibr=iCtxNumSubndNeibr>2?1:0;
	//    iCtxNumBlkNeibr=MINIMUM(iCtxNumBlkNeibr,3);
	iCtxNum = iCtxNumBlkNeibr + iCtxNumSubndNeibr;//+iCtxNumBlkNeibr;
	if (iBitPos == 0)
		iCtxNum = MINIMUM(iCtxNum, 9);
	else
		iCtxNum = MINIMUM(iCtxNum, 7);

#else
	iCtxNum = MINIMUM(0 + iCtxNumBlkNeibr, 4);
#endif

#if BIN_SCH
	if (iBitPos < 8)
	{
#if SIG_MERG
		return &(SigCtx[0][iCtxNum][iBitPos>1 ? 1 : iBitPos][LyrInfo.CtxLyr[iElemtIdx]]);
#else
#if SIG_SEPRT
		if (iElemtIdx < LyrInfo.Acum1Pos)
			return &(SigCtx[0][iCtxNum][iBitPos][iElemtIdx]);
		else
			return &(SigCtx[0][MINIMUM(iAcumSig, 3)][iBitPos][iElemtIdx]);
#else

#if LST2LYR_MERG
		return &(SigCtx[0][iCtxNum][iBitPos>1 ? 1 : iBitPos][LyrInfo.CtxLyrSig[iElemtIdx]]);
#else
		return &(SigCtx[0][iCtxNum][iBitPos>1 ? 1 : iBitPos][iElemtIdx]);
#endif

#endif
#endif
	}

	else
	{
		if (BIT_LOC_SIG > 2)
		{
#if SIG_MERG
			return &(SigCtx[0][iCtxNum][1][LyrInfo.CtxLyr[iElemtIdx]]);
#else
#if SIG_SEPRT
			if (iElemtIdx < LyrInfo.Acum1Pos)
				return &(SigCtx[0][iCtxNum][1][iElemtIdx]);
			else
				return &(SigCtx[0][MINIMUM(iAcumSig, 3)][1][iElemtIdx]);
#else
			return &(SigCtx[0][iCtxNum][1][iElemtIdx]);
#endif
#endif
		}
	}

#else
	return &(SigCtx[0][iCtxNum][0][iElemtIdx]);
#endif

}

Context *fnGetCtxExt(int iElemtIdx, int iRow, int iCol, int iBlkRow, int iBlkCol, int iBitPos)
{
	int iCtxNumSubndNeibr = 0, iCtxNumBlkNeibr = 0, iCtxNumBlkNeibrPrevBitPos = 0, iCtxNum, i, j, iRowZZplus1, iColZZplus1;

	//explore the correlation among 5 block neighbors of the current encoding coeff
	i = iBlkRow;
	j = iBlkCol;
	if (i == 7)
	{
		if (j == 7)
			iCtxNumBlkNeibr = 0;

		else if (j == 6)
			iCtxNumBlkNeibr = vrgiExt[iRow][iCol][7][7][iBitPos];
		else
			iCtxNumBlkNeibr = vrgiExt[iRow][iCol][i][j + 1][iBitPos] + vrgiExt[iRow][iCol][i][j + 2][iBitPos];
	}

	else if (i == 6)
	{
		if (j == 7)
			iCtxNumBlkNeibr = vrgiExt[iRow][iCol][7][7][iBitPos];

		else if (j == 6)
			iCtxNumBlkNeibr = vrgiExt[iRow][iCol][i][j + 1][iBitPos] + vrgiExt[iRow][iCol][i + 1][j][iBitPos] + vrgiExt[iRow][iCol][7][7][iBitPos];

		else
			iCtxNumBlkNeibr = vrgiExt[iRow][iCol][i][j + 1][iBitPos] + vrgiExt[iRow][iCol][i + 1][j][iBitPos] + vrgiExt[iRow][iCol][i + 1][j + 1][iBitPos] + vrgiExt[iRow][iCol][i][j + 2][iBitPos];
	}

	else if (j == 7)
	{
		iCtxNumBlkNeibr = vrgiExt[iRow][iCol][i + 1][j][iBitPos] + vrgiExt[iRow][iCol][i + 2][j][iBitPos];
	}

	else if (j == 6)
	{
		iCtxNumBlkNeibr = vrgiExt[iRow][iCol][i + 1][j][iBitPos] + vrgiExt[iRow][iCol][i + 2][j][iBitPos] + vrgiExt[iRow][iCol][i][j + 1][iBitPos] + vrgiExt[iRow][iCol][i + 1][j + 1][iBitPos];
	}

	else
	{
		iCtxNumBlkNeibr = vrgiExt[iRow][iCol][i + 1][j][iBitPos] + vrgiExt[iRow][iCol][i + 2][j][iBitPos] + vrgiExt[iRow][iCol][i][j + 1][iBitPos] + vrgiExt[iRow][iCol][i + 1][j + 1][iBitPos] + vrgiExt[iRow][iCol][i][j + 2][iBitPos];
	}

	if (iElemtIdx < (LyrInfo.ExtNumElemt - 1))
	{
		iRowZZplus1 = LyrInfo.ExtScanRowIdx[iElemtIdx + 1];
		iColZZplus1 = LyrInfo.ExtScanColIdx[iElemtIdx + 1];

		if ((iRowZZplus1 + iColZZplus1) == (i + j))
		{
#if IMMED_NEIBR
			if (iRowZZplus1 == i + 1 || iRowZZplus1 == i - 1)
				iCtxNumBlkNeibr += vrgiExt[iRow][iCol][iRowZZplus1][iColZZplus1][iBitPos];
#else
			iCtxNumBlkNeibr += vrgiExt[iRow][iCol][iRowZZplus1][iColZZplus1][iBitPos];
#endif
		}

	}

	/* explore the correlation among 6 subband neighbours of the current encoding coeff */
	if (iRow == 0)
	{
		if (iCol == 0)
			iCtxNumSubndNeibr = 0;

		else
		{
			iCtxNumSubndNeibr = vrgiExt[iRow][iCol - 1][iBlkRow][iBlkCol][iBitPos];
			if (iCol>1)
			{
				iCtxNumSubndNeibr += vrgiExt[iRow][iCol - 2][iBlkRow][iBlkCol][iBitPos];
			}
		}
	}

	else if (iCol == 0)
	{
		iCtxNumSubndNeibr = vrgiExt[iRow - 1][iCol][iBlkRow][iBlkCol][iBitPos] + vrgiExt[iRow - 1][iCol + 1][iBlkRow][iBlkCol][iBitPos];
		if (iRow>1)
		{
			iCtxNumSubndNeibr += vrgiExt[iRow - 2][iCol][iBlkRow][iBlkCol][iBitPos];
		}
	}

	else if (iCol == (COLS / 8 - 1))
	{
		iCtxNumSubndNeibr = vrgiExt[iRow - 1][iCol][iBlkRow][iBlkCol][iBitPos] + vrgiExt[iRow][iCol - 1][iBlkRow][iBlkCol][iBitPos] + vrgiExt[iRow - 1][iCol - 1][iBlkRow][iBlkCol][iBitPos] + vrgiExt[iRow][iCol - 2][iBlkRow][iBlkCol][iBitPos];
		if (iRow>1)
			iCtxNumSubndNeibr += vrgiExt[iRow - 2][iCol][iBlkRow][iBlkCol][iBitPos];
	}

	else
	{
		iCtxNumSubndNeibr = vrgiExt[iRow - 1][iCol][iBlkRow][iBlkCol][iBitPos] + vrgiExt[iRow][iCol - 1][iBlkRow][iBlkCol][iBitPos] + vrgiExt[iRow - 1][iCol - 1][iBlkRow][iBlkCol][iBitPos] + vrgiExt[iRow - 1][iCol + 1][iBlkRow][iBlkCol][iBitPos];
		if (iRow>1)
			iCtxNumSubndNeibr += vrgiExt[iRow - 2][iCol][iBlkRow][iBlkCol][iBitPos];

		if (iCol>1)
			iCtxNumSubndNeibr += vrgiExt[iRow][iCol - 2][iBlkRow][iBlkCol][iBitPos];
	}

#if EXT_BLK_NEIBR_PREVBIT
	if (iBitPos > 0)
	{
		if (i == 0)
		{
			if (j == 1)
				iCtxNumBlkNeibrPrevBitPos = vrgiExt[iRow][iCol][i][j - 1][iBitPos - 1];
			else
				iCtxNumBlkNeibrPrevBitPos = vrgiExt[iRow][iCol][i][j - 1][iBitPos - 1] + vrgiExt[iRow][iCol][i][j - 2][iBitPos - 1];
		}

		else if (i == 1)
		{
			if (j == 0)
				iCtxNumBlkNeibrPrevBitPos = vrgiExt[iRow][iCol][i - 1][j][iBitPos - 1];
			else if (j == 1)
				iCtxNumBlkNeibrPrevBitPos = vrgiExt[iRow][iCol][i - 1][j][iBitPos - 1] + vrgiExt[iRow][iCol][i][j - 1][iBitPos - 1] + vrgiExt[iRow][iCol][i - 1][j - 1][iBitPos - 1];
			else
				iCtxNumBlkNeibrPrevBitPos = vrgiExt[iRow][iCol][i - 1][j][iBitPos - 1] + vrgiExt[iRow][iCol][i][j - 1][iBitPos - 1] + vrgiExt[iRow][iCol][i - 1][j - 1][iBitPos - 1] + vrgiExt[iRow][iCol][i][j - 2][iBitPos - 1];
		}

		else
		{
			if (j == 0)
				iCtxNumBlkNeibrPrevBitPos = vrgiExt[iRow][iCol][i - 1][j][iBitPos - 1] + vrgiExt[iRow][iCol][i - 2][j][iBitPos - 1];
			else if (j == 1)
				iCtxNumBlkNeibrPrevBitPos = vrgiExt[iRow][iCol][i - 1][j][iBitPos - 1] + vrgiExt[iRow][iCol][i][j - 1][iBitPos - 1] + vrgiExt[iRow][iCol][i - 1][j - 1][iBitPos - 1] + vrgiExt[iRow][iCol][i - 2][j][iBitPos - 1];
			else
				iCtxNumBlkNeibrPrevBitPos = vrgiExt[iRow][iCol][i - 1][j][iBitPos - 1] + vrgiExt[iRow][iCol][i][j - 1][iBitPos - 1] + vrgiExt[iRow][iCol][i - 1][j - 1][iBitPos - 1] + vrgiExt[iRow][iCol][i - 2][j][iBitPos - 1] + vrgiExt[iRow][iCol][i][j - 2][iBitPos - 1];
		}

		if (iElemtIdx > 0)
		{
			iRowZZsub1 = LyrInfo.ExtScanRowIdx[iElemtIdx - 1];
			iColZZsub1 = LyrInfo.ExtScanColIdx[iElemtIdx - 1];

			if ((iRowZZsub1 + iColZZsub1) == (i + j))
			{
#if IMMED_NEIBR
				if (iRowZZsub1 == i + 1 || iRowZZsub1 == i - 1)
					iCtxNumBlkNeibrPrevBitPos += vrgiExt[iRow][iCol][iRowZZsub1][iColZZsub1][iBitPos - 1];
#else
				iCtxNumBlkNeibrPrevBitPos += vrgiExt[iRow][iCol][iRowZZsub1][iColZZsub1][iBitPos - 1];
#endif
			}

		}
	}
#endif

#if BIN_SCH
	if (iBitPos < 8)
	{
		if (iBitPos == 0)
		{
			//iCtxNum = 3*MINIMUM(iCtxNumBlkNeibr,2)+MINIMUM(iCtxNumSubndNeibr,2);
			iCtxNum = MINIMUM(iCtxNumBlkNeibr + iCtxNumSubndNeibr, 7);
			return &(ExtCtx[iCtxNum][0][LyrInfo.CtxLyr[iElemtIdx]]);
		}
		else
		{
#if EXT_BLK_NEIBR_PREVBIT
			iCtxNum = MINIMUM(iCtxNumBlkNeibr + iCtxNumBlkNeibrPrevBitPos, 5);
#else
			iCtxNum = MINIMUM(iCtxNumBlkNeibr + iCtxNumSubndNeibr, 5);
#endif
			return &(ExtCtx[iCtxNum][1][LyrInfo.CtxLyr[iElemtIdx]]);
		}
	}
	else
	{
		if (BIT_LOC_EXT > 2)
		{
			iCtxNum = iCtxNumBlkNeibr;
			return &(ExtCtx[iCtxNum][1][0]);
		}
	}

#else
	return &(ExtCtx[iCtxNum][0][LyrInfo.CtxLyr[iElemtIdx]]);
#endif

}

/*
*  Function: fnGetCtxSign
*  Context modeling for the sign based on the neighboring blocks (top-left, top, top-right, left;
*  i.e. Block index: 0-[i-1][j-1], 1-[i-1][j], 2-[i-1][j+1], 3-[i][j-1])
*/
Context *fnGetCtxSign(int iExt, int iRow, int iCol, int iBlkRowIdx, int iBlkColIdx)
{
	//std::cout << "ROWS is " << iRow << "cols is " << iCol << std::endl;
	
	int iCtxNum = 0;

	if (iRow == 0 && iCol == 0)  
		iCtxNum = 0;

	else if (iRow == 0 && iCol != 0)
	{

		iCtxNum = vrgiSign[iRow][iCol - 1][iBlkRowIdx][iBlkColIdx];
	}

	else if (iRow != 0 && iCol == 0)
	{
		iCtxNum = vrgiSign[iRow - 1][iCol][iBlkRowIdx][iBlkColIdx] + vrgiSign[iRow - 1][iCol + 1][iBlkRowIdx][iBlkColIdx];
	}

	else if (iRow != 0 && iCol == (COLS/8 - 1))
	{
		 
		iCtxNum = vrgiSign[iRow - 1][iCol - 1][iBlkRowIdx][iBlkColIdx] + vrgiSign[iRow - 1][iCol][iBlkRowIdx][iBlkColIdx] + vrgiSign[iRow][iCol - 1][iBlkRowIdx][iBlkColIdx];
	}

	//else  if (iRow != 0 && iCol == (COLS/8 -1 ) )
	//{

	//	iCtxNum = vrgiSign[iRow - 1][iCol - 1][iBlkRowIdx][iBlkColIdx] + vrgiSign[iRow - 1][iCol][iBlkRowIdx][iBlkColIdx] + vrgiSign[iRow ][0][iBlkRowIdx][iBlkColIdx] + vrgiSign[iRow][iCol - 1][iBlkRowIdx][iBlkColIdx];
	//}

	else
	{
		//std::cout << "second  if " << std::endl;
		iCtxNum = vrgiSign[iRow - 1][iCol - 1][iBlkRowIdx][iBlkColIdx] + vrgiSign[iRow - 1][iCol][iBlkRowIdx][iBlkColIdx] + vrgiSign[iRow - 1][iCol + 1][iBlkRowIdx][iBlkColIdx] + vrgiSign[iRow][iCol - 1][iBlkRowIdx][iBlkColIdx];
	}

	if (iCtxNum < 0)
		iCtxNum = MAXIMUM(iCtxNum, -3) + 3;
	else
		iCtxNum = MINIMUM(iCtxNum, 3) + 3;

	//if (iCtxNum < 0)
	//	iCtxNum = 0;
	//else if (iCtxNum > 0)
	//	iCtxNum = 1;
	//else
	//	iCtxNum = 2;

	if (iBlkRowIdx == 0 && iBlkColIdx == 0)
		return &(SignCtx[iExt][iCtxNum][0]);

	else if (iBlkRowIdx == 0 && iBlkColIdx == 1)
		return &(SignCtx[iExt][iCtxNum][1]);

	else if (iBlkRowIdx == 1 && iBlkColIdx == 0)
		return &(SignCtx[iExt][iCtxNum][2]);

	else if (iBlkRowIdx == 0 && iBlkColIdx < 5)
		return &(SignCtx[iExt][iCtxNum][3]);

	else if (iBlkColIdx == 0 && iBlkRowIdx < 5)
		return &(SignCtx[iExt][iCtxNum][4]);

	else if ((iBlkRowIdx == 1 && iBlkColIdx < 4) || (iBlkRowIdx == 2 && iBlkColIdx < 4) || (iBlkRowIdx == 3 && iBlkColIdx < 3))
		return &(SignCtx[iExt][iCtxNum][5]);

	else
		return &(SignCtx[iExt][iCtxNum][6]);

}

/*
*  Function: fnEncSign
*  Encode the sign information. Context modeling for the sign based on the neighboring blocks (top-left, top, top-right, left;
*  i.e. Block index: 0-[i-1][j-1], 1-[i-1][j], 2-[i-1][j+1], 3-[i][j-1])
*/
void fnEncSign(int iRowIdx, int iColIdx, int iBlkRowIdx, int iBlkColIdx, int iCoef, int iExt)
{

#if SIGN_CTX_MODL
	fnACItrfc(fnGetCtxSign(iExt, iCoefIdx, iRowIdx, iColIdx, iBlkRowIdx, iBlkColIdx), iCoef<0 ? 0 : 1);
#else
#if SIGN_AC2AC3_CTX_MODL
	if (iBlkRowIdx + iBlkColIdx == 1)
	{
		fnACItrfc(fnGetCtxSign(iExt, iRowIdx, iColIdx, iBlkRowIdx, iBlkColIdx), iCoef<0 ? 0 : 1);
	}
	else
	{
		fnACItrfcByps((iCoef<0 ? 0 : 1));
	}
#else
	fnACItrfcByps((iCoef<0 ? 0 : 1));
#endif
#endif

}

/*
*  Function: fnEncBlkBiLvImg
*  Encode the bi-level image
*/
void fnEncBlkBiLvImg(int *rgdiQntSeq, int iBlkRow, int iBlkCol)
{
	int iCurBitLoc, iElemtIdx, iCurExt, iCurExtAmp, iFreqRow, iFreqCol;

	fnACItrfc(fnGetCtxExtF(iBlkRow, iBlkCol), drgiExtF[iBlkRow][iBlkCol]);
	drgiExtFDec[iBlkRow][iBlkCol] = drgiExtF[iBlkRow][iBlkCol];

	if (drgiExtF[iBlkRow][iBlkCol] == 1)
	{
		for (iElemtIdx = (LyrInfo.ExtNumElemt - 1); iElemtIdx >= 0; iElemtIdx--)
		{
			iFreqRow = LyrInfo.ExtScanRowIdx[iElemtIdx];
			iFreqCol = LyrInfo.ExtScanColIdx[iElemtIdx];
			iCurExt = vrgiExt[iBlkRow][iBlkCol][iFreqRow][iFreqCol][0];
			vrgiExtDec[iBlkRow][iBlkCol][iFreqRow][iFreqCol][0] = iCurExt;

			/* encode Ext */
			fnACItrfc(fnGetCtxExt(iElemtIdx, iBlkRow, iBlkCol, iFreqRow, iFreqCol, 0), iCurExt);

			if (iCurExt == 1)
			{
				/* encode Sign */
				fnEncSign(iBlkRow, iBlkCol, iFreqRow, iFreqCol, rgdiQntSeq[iElemtIdx], 0);
				vrgiSign[iBlkRow][iBlkCol][iFreqRow][iFreqCol] = rgdiQntSeq[iElemtIdx]<0 ? -1 : 1;

				/* encode ExtAmp */
				iCurExtAmp = abs(rgdiQntSeq[iElemtIdx]);
#if BIN_SCH     
				iCurBitLoc = 1;
				while (iCurBitLoc < BIT_LOC_EXT)
				{
					if (drgiLPrimeTbl[iFreqRow][iFreqCol] > iCurBitLoc  && iCurExt > 0)
					{
						iCurExt = iCurExtAmp>iCurBitLoc ? 1 : 0;
						fnACItrfc(fnGetCtxExt(iElemtIdx, iBlkRow, iBlkCol, iFreqRow, iFreqCol, iCurBitLoc), iCurExt);
						vrgiExt[iBlkRow][iBlkCol][iFreqRow][iFreqCol][iCurBitLoc] = iCurExt;
						vrgiExtDec[iBlkRow][iBlkCol][iFreqRow][iFreqCol][iCurBitLoc] = iCurExt;
						iCurBitLoc++;
					}

					else
						break;
				}

				if (drgiLPrimeTbl[iFreqRow][iFreqCol] > BIT_LOC_EXT && iCurExt > 0)
				{
#if ENC_EXTAMP_SEPRT
					//fnACItrfc(&(ExtAmpSeprtCtx[0][LyrInfo.LyrIdx_MaxLprime]), iCurExtAmp-BIT_LOC_EXT);
					fnACmAlpItrfc(&(ExtAmpSeprtCtx[0][LyrInfo.LyrIdx_MaxLprime]), iCurExtAmp - BIT_LOC_EXT, drgiLPrimeTbl[iFreqRow][iFreqCol] - BIT_LOC_EXT + 1);
					//fnACmAlpItrfc(&(ExtAmpSeprtCtx[0][LyrInfo.CtxLyrIdxMaxLprime[iElemtIdx]]), iCurExtAmp-BIT_LOC_EXT, drgiLPrimeTbl[iFreqRow][iFreqCol]-BIT_LOC_EXT+1);
					vrgiHDQIdxDec[iBlkRow][iBlkCol][iFreqRow][iFreqCol] = iCurExtAmp - BIT_LOC_EXT;
#else
					fnACItrfc(&(ExtAmpCtx[0][iElemtIdx]), iCurExtAmp - BIT_LOC_EXT);
#endif
				}

				vrgiHDQIdxDec[iBlkRow][iBlkCol][iFreqRow][iFreqCol] = iCurExtAmp - BIT_LOC_EXT;

#else                

#if ENC_EXTAMP_SEPRT
				fnACItrfc(&(ExtAmpSeprtCtx[0][iElemtIdx]), iCurExtAmp - 1);
#else
				fnACItrfc(&(ExtAmpCtx[0][iElemtIdx]), iCurExtAmp - 1);
#endif
#endif
			}
		}
	}
}

/*
*  Function: fnEncOneBlk
*  Encode the truncated laplacian part
*/
void fnEncOneBlk(int *rgdiQntSeq, int iBlkRow, int iBlkCol)
{
	int iFreqRow, iFreqCol, iCurBitLoc, iElemtIdx, iAcumSigCur, iAcumSigNxt, iCurAmpExt, iCurSigExt, iCurAmp, iCurSig;

	if (drgiExtF[iBlkRow][iBlkCol] == 1)
	{
		iAcumSigCur = 0;
		for (iElemtIdx = (LyrInfo.SigNumElemt - 1); iElemtIdx >= 0; iElemtIdx--)
		{
			iFreqRow = LyrInfo.SigScanRowIdx[iElemtIdx];
			iFreqCol = LyrInfo.SigScanColIdx[iElemtIdx];

			if (vrgiExt[iBlkRow][iBlkCol][iFreqRow][iFreqCol][0] == 0)
			{
				iCurAmpExt = abs(rgdiQntSeq[iElemtIdx]);
				iCurSigExt = iCurAmpExt>0 ? 1 : 0;
				iAcumSigNxt = iAcumSigCur + iCurSigExt;

				if (drgiLTbl[iFreqRow][iFreqCol] > 0)
				{
					/* encode SigExt */
					fnACItrfc(fnGetCtxSigExt(iAcumSigCur, iElemtIdx, iBlkRow, iBlkCol, iFreqRow, iFreqCol, 0, iCurSigExt), iCurSigExt);

					/* encode the level info */
					if (iCurSigExt == 1)
					{
						vrgiSig[iBlkRow][iBlkCol][iFreqRow][iFreqCol][0] = 1;

						/* encode Sign */
						fnEncSign(iBlkRow, iBlkCol, iFreqRow, iFreqCol, rgdiQntSeq[iElemtIdx], 0);
						vrgiSign[iBlkRow][iBlkCol][iFreqRow][iFreqCol] = rgdiQntSeq[iElemtIdx]<0 ? -1 : 1;

						if (drgiLTbl[iFreqRow][iFreqCol] > 1)
						{
#if BIN_SCH
							iCurBitLoc = 1;
							while (iCurBitLoc < BIT_LOC_SIG)
							{
								if (drgiLTbl[iFreqRow][iFreqCol] > iCurBitLoc && iCurSigExt == 1)
								{
									iCurSigExt = iCurAmpExt>iCurBitLoc ? 1 : 0;
									fnACItrfc(fnGetCtxSigExt(iAcumSigCur, iElemtIdx, iBlkRow, iBlkCol, iFreqRow, iFreqCol, iCurBitLoc, iCurSigExt), iCurSigExt);
									vrgiSig[iBlkRow][iBlkCol][iFreqRow][iFreqCol][iCurBitLoc] = iCurSigExt;
									iCurBitLoc++;
								}

								else
									break;
							}

							/* encode AmpExt */
							if (drgiLTbl[iFreqRow][iFreqCol] > BIT_LOC_SIG && iCurSigExt == 1)
							{
								fnACmAlpItrfc(&(AmpCtx[1][0][LyrInfo.CtxLyrIdxMaxL[iElemtIdx]]), iCurAmpExt - BIT_LOC_SIG, drgiLTbl[iFreqRow][iFreqCol] - BIT_LOC_SIG + 1);
							}

							vrgiHDQIdxDec[iBlkRow][iBlkCol][iFreqRow][iFreqCol] = iCurAmpExt - BIT_LOC_SIG;

#else
							fnACItrfc(&(AmpCtx[1][0][iElemtIdx]), iCurAmpExt - 1);
							//fnACItrfc(&(AmpCtx[0][iCtxSigExtTmp][iElemtIdx]), iCurAmpExt-1);
#endif
						}
					}
				}

				iAcumSigCur = iAcumSigNxt;
			}

			else
			{
				for (iCurBitLoc = 0; iCurBitLoc<BIT_LOC_SIG; iCurBitLoc++)
				{
					vrgiSig[iBlkRow][iBlkCol][iFreqRow][iFreqCol][iCurBitLoc] = 1;
				}
			}
		}
	}

	else
	{
		fnACItrfc(fnGetCtxAZF(iBlkRow, iBlkCol, drgiCBP[iBlkRow][iBlkCol]), drgiCBP[iBlkRow][iBlkCol]);

		if (drgiCBP[iBlkRow][iBlkCol] == 1)
		{
			iAcumSigCur = 0;
			for (iElemtIdx = (LyrInfo.SigNumElemt - 1); iElemtIdx >= 0; iElemtIdx--)
			{
				iFreqRow = LyrInfo.SigScanRowIdx[iElemtIdx];
				iFreqCol = LyrInfo.SigScanColIdx[iElemtIdx];

				iCurAmp = abs(rgdiQntSeq[iElemtIdx]);
				iCurSig = iCurAmp>0 ? 1 : 0;
				iAcumSigNxt = iAcumSigCur + iCurSig;

				/* encode Sig */
				if (drgiLTbl[iFreqRow][iFreqCol] > 0)
				{
					fnACItrfc(fnGetCtxSig(iAcumSigCur, iElemtIdx, iBlkRow, iBlkCol, iFreqRow, iFreqCol, 0, iCurSig), iCurSig);

					/* encode the level info */
					if (iCurSig == 1)
					{
						vrgiSig[iBlkRow][iBlkCol][iFreqRow][iFreqCol][0] = 1;

						/* encode Sign */
						fnEncSign(iBlkRow, iBlkCol, iFreqRow, iFreqCol, rgdiQntSeq[iElemtIdx], 0);
						vrgiSign[iBlkRow][iBlkCol][iFreqRow][iFreqCol] = rgdiQntSeq[iElemtIdx]<0 ? -1 : 1;

						if (drgiLTbl[iFreqRow][iFreqCol] > 1)
						{
#if BIN_SCH
							iCurBitLoc = 1;
							while (iCurBitLoc < BIT_LOC_SIG)
							{
								if (drgiLTbl[iFreqRow][iFreqCol] > iCurBitLoc && iCurSig > 0)
								{
									iCurSig = iCurAmp>iCurBitLoc ? 1 : 0;
									fnACItrfc(fnGetCtxSig(iAcumSigCur, iElemtIdx, iBlkRow, iBlkCol, iFreqRow, iFreqCol, iCurBitLoc, iCurSig), iCurSig);
									vrgiSig[iBlkRow][iBlkCol][iFreqRow][iFreqCol][iCurBitLoc] = iCurSig;
									iCurBitLoc++;
								}

								else
									break;
							}

							/* encode Amp */
							if (drgiLTbl[iFreqRow][iFreqCol] > BIT_LOC_SIG && iCurSig == 1)
							{
								fnACmAlpItrfc(&(AmpCtx[0][0][LyrInfo.CtxLyrIdxMaxL[iElemtIdx]]), iCurAmp - BIT_LOC_SIG, drgiLTbl[iFreqRow][iFreqCol] - BIT_LOC_SIG + 1);
								//fnACItrfc(&(AmpCtx[0][0][iElemtIdx]), iCurAmp-BIT_LOC_SIG);
							}

							vrgiHDQIdxDec[iBlkRow][iBlkCol][iFreqRow][iFreqCol] = iCurAmp - BIT_LOC_SIG;

#else
							fnACItrfc(&(AmpCtx[0][0][iElemtIdx]), iCurAmp - 1);
							//fnACItrfc(&(AmpCtx[0][iCtxSigTmp][iElemtIdx]), iCurAmp-1);
#endif
						}
					}
				}

				iAcumSigCur = iAcumSigNxt;
			}
		}
	}
}

void fnEncDC(int iOptCtxDC, int iIdxDC)
{
	int iCtxIncDC, iSymblNxtCtxDC;

#if DC_SEP_CTX
	if (iIdxDC > (rgiAlpSizeDC[iOptCtxDC] - 2))
	{
		fnACItrfcCoeffDC(&(CtxDC[iOptCtxDC]), rgiAlpSizeDC[iOptCtxDC] - 1);
		iSymblNxtCtxDC = iIdxDC - rgiAlpSizeDC[iOptCtxDC] + 1;
		if (iSymblNxtCtxDC == 0 && iOptCtxDC == NUM_CTX_DC - 1);
		else
		{
			iCtxIncDC = 1;
			while (iSymblNxtCtxDC > (rgiAlpSizeDC[iOptCtxDC + iCtxIncDC] - 2) && iOptCtxDC < NUM_CTX_DC - 1)
			{
				fnACItrfcCoeffDC(&(CtxDC[iOptCtxDC + iCtxIncDC]), rgiAlpSizeDC[iOptCtxDC + iCtxIncDC] - 1);
				iSymblNxtCtxDC = iSymblNxtCtxDC - rgiAlpSizeDC[iOptCtxDC + iCtxIncDC] + 1;
				iCtxIncDC++;
			}

			fnACItrfcCoeffDC(&(CtxDC[iOptCtxDC + iCtxIncDC]), iSymblNxtCtxDC);
		}
	}

	else
		fnACItrfcCoeffDC(&(CtxDC[iOptCtxDC]), iIdxDC);
#else
	fnACItrfcCoeffDC(&(CtxDC[iOptCtxDC]), iIdxDC);
#endif
}

/*
*  Function: fnEncOneBlk
*  Encode the truncated laplacian part
*/
void fnColecStatisElemt(int iQntIdx, int iElemtIdx, int iBlkRow, int iBlkCol, int iFreqRow, int iFreqCol)
{
	int iCurBitLoc = 0, iCurAmp, iCurSig;

	iCurAmp = abs(iQntIdx);
	iCurSig = iCurAmp>0 ? 1 : 0;

	fnColecStatisSig(drgiExtF[iBlkRow][iBlkCol], iBlkRow, iBlkCol, iCurBitLoc, iElemtIdx, iCurSig);

	if (iCurSig == 1)
	{
		iCurBitLoc = 1;
		while (iCurBitLoc < BIT_LOC_SIG)
		{
			if (drgiLTbl[iFreqRow][iFreqCol] > iCurBitLoc && iCurSig == 1)
			{
				iCurSig = iCurAmp>iCurBitLoc ? 1 : 0;
				fnColecStatisSig(drgiExtF[iBlkRow][iBlkCol], iBlkRow, iBlkCol, iCurBitLoc, iElemtIdx, iCurSig);
				iCurBitLoc++;
			}

			else
			{
				while (iCurBitLoc < BIT_LOC_SIG)
				{
					prgiSig[iBlkRow][iBlkCol][iCurBitLoc++][iElemtIdx] = 0;
				}

				break;
			}
		}

		/* encode Amp */
		if (drgiLTbl[iFreqRow][iFreqCol] > BIT_LOC_SIG && iCurSig == 1)
		{
			i3pAmpCtxCtr[drgiExtF[iBlkRow][iBlkCol]][LyrInfo.CtxLyr[iElemtIdx]][iCurAmp - BIT_LOC_SIG]++;
#if TEST_RD
			i3pAmpCtxElemtCtr[drgiExtF[iBlkRow][iBlkCol]][iElemtIdx][iCurAmp - BIT_LOC_SIG]++;
#endif
		}
	}

	else
	{
		for (iCurBitLoc = 1; iCurBitLoc < BIT_LOC_SIG; iCurBitLoc++)
		{
			prgiSig[iBlkRow][iBlkCol][iCurBitLoc][iElemtIdx] = 0;
		}
	}
}

void fnColecStatisSig(int iCurEXTF, int iBlkRow, int iBlkCol, int iCurBitLoc, int iElemtIdx, int iCurSig)
{
	int iCtxIdx, iElemtIdxTmp, iRowZZplus1, iColZZplus1, iCtxNum = 0;
	int i, j, iCtxNumBlkExt = 0, iRowZZMinus1, iColZZMinus1, iCtxNumSubndNeibr = 0;

	i = LyrInfo.SigScanRowIdx[iElemtIdx];
	j = LyrInfo.SigScanColIdx[iElemtIdx];

	for (iCtxIdx = 0; iCtxIdx<LyrInfo.SigNumCtx[iElemtIdx]; iCtxIdx++)
	{
		iCtxNum += prgiSig[iBlkRow][iBlkCol][iCurBitLoc][LyrInfo.SigCtxMdl[iElemtIdx][iCtxIdx]];
	}

#if JUNE26_SIG_LYR_NEIBER
	iElemtIdxTmp = iElemtIdx;
	while (iElemtIdxTmp < (LyrInfo.SigNumElemt - 1))
	{
		iRowZZplus1 = LyrInfo.SigScanRowIdx[iElemtIdxTmp + 1];
		iColZZplus1 = LyrInfo.SigScanColIdx[iElemtIdxTmp + 1];

		if ((iRowZZplus1 + iColZZplus1) == (i + j))
		{
			if (vrgiExt[iBlkRow][iBlkCol][iRowZZplus1][iColZZplus1][0] == 0)
			{
				iCtxNum += vrgiSig[iBlkRow][iBlkCol][iRowZZplus1][iColZZplus1][iCurBitLoc];
				break;
			}
		}
		else
		{
			break;
		}
		iElemtIdxTmp++;
	}
#else
	if (iElemtIdx < (LyrInfo.SigNumElemt - 1))
	{
		iRowZZplus1 = LyrInfo.SigScanRowIdx[iElemtIdx + 1];
		iColZZplus1 = LyrInfo.SigScanColIdx[iElemtIdx + 1];

		if ((iRowZZplus1 + iColZZplus1) == (i + j))
		{
#if IMMED_NEIBR
			if (iRowZZplus1 == i + 1 || iRowZZplus1 == i - 1)
				iCtxNum += vrgiSig[iBlkRow][iBlkCol][iRowZZplus1][iColZZplus1][iCurBitLoc];
#else
			iCtxNum += vrgiSig[iBlkRow][iBlkCol][iRowZZplus1][iColZZplus1][iCurBitLoc];
#endif
		}
	}
#endif

	if (drgiExtF[iBlkRow][iBlkCol] == 1)
	{
		if (i == 0)
		{
			if (j>0)
				iCtxNumBlkExt = vrgiExt[iBlkRow][iBlkCol][i][j - 1][0];
			else
				iCtxNumBlkExt = 0;
		}

		else if (j == 1)
			iCtxNumBlkExt = vrgiExt[iBlkRow][iBlkCol][i - 1][j][0];

		else
			iCtxNumBlkExt = vrgiExt[iBlkRow][iBlkCol][i - 1][j][0] + vrgiExt[iBlkRow][iBlkCol][i][j - 1][0] + vrgiExt[iBlkRow][iBlkCol][i - 1][j - 1][0];

		if (iElemtIdx > 0)
		{
			iRowZZMinus1 = LyrInfo.SigScanRowIdx[iElemtIdx - 1];
			iColZZMinus1 = LyrInfo.SigScanColIdx[iElemtIdx - 1];

			if ((iRowZZMinus1 + iColZZMinus1) == (i + j))
			{
#if IMMED_NEIBR
				if (iRowZZMinus1 == i + 1 || iRowZZMinus1 == i - 1)
					iCtxNumBlkExt += vrgiExt[iBlkRow][iBlkCol][iRowZZMinus1][iColZZMinus1][0];
#else
				iCtxNumBlkExt += vrgiExt[iBlkRow][iBlkCol][iRowZZMinus1][iColZZMinus1][0];
#endif
			}
		}

		iCtxNum += iCtxNumBlkExt;
	}

#if INTER_BLK_SIG_HELP
	//explore the correlation among 6 subband neighbors of the current encoding coeff
	if (iBlkRow == 0)
	{
		if (iBlkCol == 0)
			iCtxNumSubndNeibr = 0;

		else
		{
			iCtxNumSubndNeibr = prgiSig[iBlkRow][iBlkCol - 1][iCurBitLoc][iElemtIdx];
			if (iBlkCol>1)
			{
				iCtxNumSubndNeibr += prgiSig[iBlkRow][iBlkCol - 2][iCurBitLoc][iElemtIdx];
			}
		}
	}

	else if (iBlkCol == 0)
	{
		iCtxNumSubndNeibr = prgiSig[iBlkRow - 1][iBlkCol][iCurBitLoc][iElemtIdx] + prgiSig[iBlkRow - 1][iBlkCol + 1][iCurBitLoc][iElemtIdx];
		if (iBlkRow>1)
		{
			iCtxNumSubndNeibr += prgiSig[iBlkRow - 2][iBlkCol][iCurBitLoc][iElemtIdx];
		}
	}

	else if (iBlkCol == (COLS / 8 - 1))
	{
		iCtxNumSubndNeibr = prgiSig[iBlkRow - 1][iBlkCol][iCurBitLoc][iElemtIdx] + prgiSig[iBlkRow][iBlkCol - 1][iCurBitLoc][iElemtIdx] + prgiSig[iBlkRow - 1][iBlkCol - 1][iCurBitLoc][iElemtIdx] + prgiSig[iBlkRow][iBlkCol - 2][iCurBitLoc][iElemtIdx];
		if (iBlkRow>1)
			iCtxNumSubndNeibr += prgiSig[iBlkRow - 2][iBlkCol][iCurBitLoc][iElemtIdx];
	}

	else
	{
		iCtxNumSubndNeibr = prgiSig[iBlkRow - 1][iBlkCol][iCurBitLoc][iElemtIdx] + prgiSig[iBlkRow][iBlkCol - 1][iCurBitLoc][iElemtIdx] + prgiSig[iBlkRow - 1][iBlkCol - 1][iCurBitLoc][iElemtIdx] + prgiSig[iBlkRow - 1][iBlkCol + 1][iCurBitLoc][iElemtIdx];
		if (iBlkRow>1)
			iCtxNumSubndNeibr += prgiSig[iBlkRow - 2][iBlkCol][iCurBitLoc][iElemtIdx];

		if (iBlkCol>1)
			iCtxNumSubndNeibr += prgiSig[iBlkRow][iBlkCol - 2][iCurBitLoc][iElemtIdx];
	}

	iCtxNum += iCtxNumSubndNeibr;

	if (drgiExtF[iBlkRow][iBlkCol] == 1)
	{
		if (iCurBitLoc == 0)
			iCtxNum = MINIMUM(iCtxNum, 10);
		else
			iCtxNum = MINIMUM(iCtxNum, 8);
	}

	else
	{
		if (iCurBitLoc == 0)
			iCtxNum = MINIMUM(iCtxNum, 9);
		else
			iCtxNum = MINIMUM(iCtxNum, 7);
	}
#endif

	i5pSigCtxCtr[iCurEXTF][iCtxNum][iCurBitLoc>0 ? 1 : 0][LyrInfo.CtxLyr[iElemtIdx]][iCurSig]++;

#if TEST_RD
	i5pSigCtxElemtCtr[iCurEXTF][iCtxNum][iCurBitLoc>0 ? 1 : 0][iElemtIdx][iCurSig]++;
#endif

	prgiSig[iBlkRow][iBlkCol][iCurBitLoc][iElemtIdx] = iCurSig;
}

void fnColecStatisBlk(std::vector<std::vector<int>> drgiQtzIdx, int iBlkRow, int iBlkCol)
{
	int iElemtIdx, iFreqRow, iFreqCol, iCurBitLoc;//, iCtxNumAZF;

	if (drgiExtF[iBlkRow][iBlkCol] == 0)
	{
		//        iCtxNumAZF = fnCalCtxAZF(iBlkRow, iBlkCol);
		//        ipAZFCtxMainCtr[iCtxNumAZF]++;
		//        i2pAZFCtxCtr[iCtxNumAZF][drgiCBP[iBlkRow][iBlkCol]]++; 

		if (drgiCBP[iBlkRow][iBlkCol] == 1)
		{
			for (iElemtIdx = (LyrInfo.SigNumElemt - 1); iElemtIdx >= 0; iElemtIdx--)
			{
				iFreqRow = LyrInfo.SigScanRowIdx[iElemtIdx];
				iFreqCol = LyrInfo.SigScanColIdx[iElemtIdx];

				fnColecStatisElemt(drgiQtzIdx[iBlkRow * 8 + iFreqRow][iBlkCol * 8 + iFreqCol], iElemtIdx, iBlkRow, iBlkCol, iFreqRow, iFreqCol);
#if TEST_RD
				i2pAmpCtrEmp[iElemtIdx][abs(drgiQtzIdx[iBlkRow * 8 + iFreqRow][iBlkCol * 8 + iFreqCol])]++;
#endif
			}
		}
		else
		{
			for (iElemtIdx = (LyrInfo.SigNumElemt - 1); iElemtIdx >= 0; iElemtIdx--)
			{
				iFreqRow = LyrInfo.SigScanRowIdx[iElemtIdx];
				iFreqCol = LyrInfo.SigScanColIdx[iElemtIdx];
#if TEST_RD
				//fnColecStatisElemt(drgiQtzIdx[iBlkRow*8+iFreqRow][iBlkCol*8+iFreqCol], iElemtIdx, iBlkRow, iBlkCol, iFreqRow, iFreqCol);
				i2pAmpCtrEmp[iElemtIdx][0]++;
#else
				for (iCurBitLoc = 0; iCurBitLoc<BIT_LOC_SIG; iCurBitLoc++)
				{
					prgiSig[iBlkRow][iBlkCol][iCurBitLoc][iElemtIdx] = 0;
				}
#endif
			}
		}
	}

	else
	{
		for (iElemtIdx = (LyrInfo.SigNumElemt - 1); iElemtIdx >= 0; iElemtIdx--)
		{
			iFreqRow = LyrInfo.SigScanRowIdx[iElemtIdx];
			iFreqCol = LyrInfo.SigScanColIdx[iElemtIdx];

			if (vrgiExt[iBlkRow][iBlkCol][iFreqRow][iFreqCol][0] == 0)
			{
				fnColecStatisElemt(drgiQtzIdx[iBlkRow * 8 + iFreqRow][iBlkCol * 8 + iFreqCol], iElemtIdx, iBlkRow, iBlkCol, iFreqRow, iFreqCol);
#if TEST_RD
				i2pAmpCtrEmp[iElemtIdx][abs(drgiQtzIdx[iBlkRow * 8 + iFreqRow][iBlkCol * 8 + iFreqCol])]++;
#endif
			}

			else
			{
				for (iCurBitLoc = 0; iCurBitLoc<BIT_LOC_SIG; iCurBitLoc++)
				{
					prgiSig[iBlkRow][iBlkCol][iCurBitLoc][iElemtIdx] = 1;
				}
			}
		}
	}
}

void fnBlkSDQ(int iBlkRow, int iBlkCol, double** drgdOrigCoef)
{
	int iElemtIdx, iLyrElemtIdx, iLyrElemtIdxCache[8], iLyrIdx, iFreqRow, iFreqCol, iNumLyrElemt, iFlagCBP = 0, iFlagFrstElemt;
	double dAbsOrigCoef;
	StateSDQ TrellisSDQ;

	iElemtIdx = LyrInfo.SigNumElemt - 1;
	for (iLyrIdx = LyrInfo.SigNumLyr - 1; iLyrIdx >= 0; iLyrIdx--)
	{
		iFlagFrstElemt = 0;
		iNumLyrElemt = 0;
		for (iLyrElemtIdx = 0; iLyrElemtIdx<LyrInfo.SigLyrNumElemt[iLyrIdx]; iLyrElemtIdx++)
		{
			iFreqRow = LyrInfo.SigScanRowIdx[iElemtIdx];
			iFreqCol = LyrInfo.SigScanColIdx[iElemtIdx];
			dAbsOrigCoef = fabs(drgdOrigCoef[iBlkRow * 8 + iFreqRow][iBlkCol * 8 + iFreqCol]);

			if (iLyrElemtIdx == 0 || iFlagFrstElemt == 1)
			{
				if (vrgiExt[iBlkRow][iBlkCol][iFreqRow][iFreqCol][0] == 0)
				{
					iLyrElemtIdxCache[iNumLyrElemt++] = iElemtIdx;
					fnInitTrellis(iBlkRow, iBlkCol, iFreqRow, iFreqCol, iElemtIdx, dAbsOrigCoef, &TrellisSDQ);
					iFlagFrstElemt = 0;
				}
				else
				{
					iFlagFrstElemt = 1;
				}
			}
			else
			{
				if (vrgiExt[iBlkRow][iBlkCol][iFreqRow][iFreqCol][0] == 0)
				{
					iLyrElemtIdxCache[iNumLyrElemt] = iElemtIdx;
					fnContTrellis(iBlkRow, iBlkCol, iFreqRow, iFreqCol, iLyrElemtIdxCache, iNumLyrElemt, dAbsOrigCoef, &TrellisSDQ);
					iNumLyrElemt++;
				}
			}

			iElemtIdx--;
		}

		if (iNumLyrElemt > 0)
		{
			fnEndTrellis(iBlkRow, iBlkCol, drgdOrigCoef, iLyrElemtIdxCache, iNumLyrElemt - 1, &iFlagCBP, &TrellisSDQ);
		}
	}

	if (iFlagCBP == 0)
	{
		drgiCBP[iBlkRow][iBlkCol] = 0;
	}
	else
	{
		drgiCBP[iBlkRow][iBlkCol] = 1;
	}
}

void fnEndTrellis(int iBlkRow, int iBlkCol, double ** drgdOrigCoef, int iLyrElemtIdxCache[8], int iCurLyrElemtIdx, int *iCBP, StateSDQ* CurLyrTrellis)
{
	int iCurBitLoc, iIdxSDQ, iFreqRow, iFreqCol, iSign, iPrvStat, iOptIdxAbs;
	double dMinLyrCost;// dCurLyrDist=0.0;

	dMinLyrCost = CurLyrTrellis->s_cost[iCurLyrElemtIdx][0];
	iPrvStat = 0;
	for (iIdxSDQ = 1; iIdxSDQ <= MINIMUM(3, LyrInfo.SigAphSize[iLyrElemtIdxCache[iCurLyrElemtIdx]]); iIdxSDQ++)
	{
		if (dMinLyrCost > CurLyrTrellis->s_cost[iCurLyrElemtIdx][iIdxSDQ])
		{
			dMinLyrCost = CurLyrTrellis->s_cost[iCurLyrElemtIdx][iIdxSDQ];
			iPrvStat = iIdxSDQ;
		}
	}

	dAvgCostRD += CurLyrTrellis->s_cost[iCurLyrElemtIdx][iPrvStat];

	//    if (CurLyrTrellis->s_cost[iCurLyrElemtIdx][iPrvStat]<0)//>CurLyrTrellis->min_cost_SDQ)
	//    {
	//        iErrorNumSDQ++;
	//        printf("error: SDQ cost is higer than HDQ cost (HDQ: %4.4f/ SDQ: %4.4f)!\n", CurLyrTrellis->min_cost_SDQ, CurLyrTrellis->s_cost[iCurLyrElemtIdx][iPrvStat]);
	//        
	//        while(iCurLyrElemtIdx >= 0)
	//        {
	//            iFreqRow = LyrInfo.SigScanRowIdx[iLyrElemtIdxCache[iCurLyrElemtIdx]];
	//            iFreqCol = LyrInfo.SigScanColIdx[iLyrElemtIdxCache[iCurLyrElemtIdx]];
	//            iSign = drgdOrigCoef[iBlkRow*8+iFreqRow][iBlkCol*8+iFreqCol]<0?-1:1;
	//            drgiQIdxSDQ[iBlkRow*8+iFreqRow][iBlkCol*8+iFreqCol] = iSign*CurLyrTrellis->opt_amp[iCurLyrElemtIdx];
	//            if (CurLyrTrellis->opt_amp[iCurLyrElemtIdx] > 2)
	//            {
	//                iCurBitLoc = 0;
	//                while (iCurBitLoc < BIT_LOC_SIG)
	//                {
	//                    prgiSig[iBlkRow][iBlkCol][iCurBitLoc++][iLyrElemtIdxCache[iCurLyrElemtIdx]] = 1;
	//                }
	//            }
	//            else 
	//            {
	//                iCurBitLoc = 0;
	//                while (iCurBitLoc < BIT_LOC_SIG)
	//                {
	//                    prgiSig[iBlkRow][iBlkCol][iCurBitLoc][iLyrElemtIdxCache[iCurLyrElemtIdx]] = iPrvStat>iCurBitLoc?1:0;
	//                    iCurBitLoc++;
	//                }
	//            }
	//            
	//            //dCurLyrDist += sqterm(fnDeQtzSig(drgiQIdxSDQ[iBlkRow*8+iFreqRow][iBlkCol*8+iFreqCol], iFreqRow, iFreqCol), drgdOrigCoef[iBlkRow*8+iFreqRow][iBlkCol*8+iFreqCol]);
	//            
	//            if (CurLyrTrellis->opt_amp[iCurLyrElemtIdx] > 0)
	//            {
	//                *iCBP = 1;
	//            }
	//            
	//            iCurLyrElemtIdx--;
	//            
	//        }
	//    }

	if (CurLyrTrellis->s_cost[iCurLyrElemtIdx][iPrvStat] - CurLyrTrellis->total_cost_HDQ > 0.00000000001)
	{
		printf("Error: SDQ coding cost is larger than HDQ! HDQ=%4.4f, SDQ=%4.4f\n", CurLyrTrellis->total_cost_HDQ, CurLyrTrellis->s_cost[iCurLyrElemtIdx][iPrvStat]);

		//        printf("The coding cost for the current layer has been reduced by at least 90 percent! (HDQ: %4.4f/ SDQ: %4.4f)!\n", CurLyrTrellis->total_cost_HDQ, CurLyrTrellis->s_cost[iCurLyrElemtIdx][iPrvStat]);
		iErrorNumSDQ++;
	}

	//else{
	/* record the optimal quantized index sequence */
	while (iCurLyrElemtIdx >= 0)
	{
		iFreqRow = LyrInfo.SigScanRowIdx[iLyrElemtIdxCache[iCurLyrElemtIdx]];
		iFreqCol = LyrInfo.SigScanColIdx[iLyrElemtIdxCache[iCurLyrElemtIdx]];
		iSign = drgdOrigCoef[iBlkRow * 8 + iFreqRow][iBlkCol * 8 + iFreqCol]<0 ? -1 : 1;
		if (iPrvStat > 2)
		{
			iOptIdxAbs = CurLyrTrellis->s_opt_amp[iCurLyrElemtIdx];
			iCurBitLoc = 0;
			while (iCurBitLoc < BIT_LOC_SIG)
			{
				prgiSig[iBlkRow][iBlkCol][iCurBitLoc++][iLyrElemtIdxCache[iCurLyrElemtIdx]] = 1;
			}
		}
		else
		{
			iOptIdxAbs = iPrvStat;
			iCurBitLoc = 0;
			while (iCurBitLoc < BIT_LOC_SIG)
			{
				prgiSig[iBlkRow][iBlkCol][iCurBitLoc][iLyrElemtIdxCache[iCurLyrElemtIdx]] = iPrvStat>iCurBitLoc ? 1 : 0;
				iCurBitLoc++;
			}
		}
		drgiQIdxSDQ[iBlkRow * 8 + iFreqRow][iBlkCol * 8 + iFreqCol] = iSign * iOptIdxAbs;

		//dCurLyrDist += sqterm(fnDeQtzSig(drgiQIdxSDQ[iBlkRow*8+iFreqRow][iBlkCol*8+iFreqCol], iFreqRow, iFreqCol), drgdOrigCoef[iBlkRow*8+iFreqRow][iBlkCol*8+iFreqCol]);

		if (iPrvStat > 0)
		{
			*iCBP = 1;
			drgdSumY[iFreqRow][iFreqCol] += fabs(drgdOrigCoef[iBlkRow * 8 + iFreqRow][iBlkCol * 8 + iFreqCol]);
			drgdSumYmulC[iFreqRow][iFreqCol] += fabs(drgdOrigCoef[iBlkRow * 8 + iFreqRow][iBlkCol * 8 + iFreqCol] * drgiQIdxSDQ[iBlkRow * 8 + iFreqRow][iBlkCol * 8 + iFreqCol]);
			drgdSumCmulC[iFreqRow][iFreqCol] += (double)abs(drgiQIdxSDQ[iBlkRow * 8 + iFreqRow][iBlkCol * 8 + iFreqCol] * drgiQIdxSDQ[iBlkRow * 8 + iFreqRow][iBlkCol * 8 + iFreqCol]);
			drgdSumC[iFreqRow][iFreqCol] += (double)abs(drgiQIdxSDQ[iBlkRow * 8 + iFreqRow][iBlkCol * 8 + iFreqCol]);

			drgiNumBlkClgt0[iFreqRow][iFreqCol]++;

			if (iOptIdxAbs > drgiLTblSDQ[iFreqRow][iFreqCol])
				drgiLTblSDQ[iFreqRow][iFreqCol] = iOptIdxAbs;
		}

		dTestRate += CurLyrTrellis->s_rate[iCurLyrElemtIdx][iPrvStat];
		iPrvStat = CurLyrTrellis->s_prvste[iCurLyrElemtIdx--][iPrvStat];
		//}
	}
	dTestRateHDQ += CurLyrTrellis->total_rate_HDQ;
	dCostHDQ += CurLyrTrellis->total_cost_HDQ;
	//dTestDist += dCurLyrDist;
}

void fnContTrellis(int iBlkRow, int iBlkCol, int iFreqRow, int iFreqCol, int iLyrElemtIdxCache[8], int iCurLyrElemtIdx, double dOrigCoef, StateSDQ* CurLyrTrellis)
{
	int iIdxSDQ, iIdxHDQ, iPrvIdxHDQ, iIdxPrvSDQ, iCurFixCtx[4], iCurCtx, iCurElemt, iPrvElemt, iCurBitLoc = 0;
	double dCostTmp, dCurRateTmp, dCostMin = INFINITY;
	double dCurRate, dMinStateRate, dCurDist, dCurStateCost, dMinStateCost, dCurRateSign;

	iCurElemt = iLyrElemtIdxCache[iCurLyrElemtIdx];
	iPrvElemt = iLyrElemtIdxCache[iCurLyrElemtIdx - 1];

	iIdxHDQ = abs(drgiQIdxHDQ[iBlkRow * 8 + iFreqRow][iBlkCol * 8 + iFreqCol]);
	iPrvIdxHDQ = abs(drgiQIdxHDQ[iBlkRow * 8 + LyrInfo.SigScanRowIdx[iPrvElemt]][iBlkCol * 8 + LyrInfo.SigScanColIdx[iPrvElemt]]);

	for (iCurBitLoc = 0; iCurBitLoc<MINIMUM(3, LyrInfo.SigAphSize[iCurElemt]); iCurBitLoc++)
	{
		iCurFixCtx[iCurBitLoc] = fnCalFixCtxSig(iBlkRow, iBlkCol, iCurBitLoc, iCurElemt);
		if (iCurFixCtx[iCurBitLoc] < 0)
			printf("Error: context modeling!\n");
	}

	for (iIdxSDQ = 0; iIdxSDQ<3; iIdxSDQ++)
	{
		if (iIdxSDQ <= LyrInfo.SigAphSize[iCurElemt])
		{
			dCurRateSign = iIdxSDQ == 0 ? 0.0 : SDQ_LAMBDA;//including the rate for encoding the Sign 
			dMinStateCost = INFINITY;

			for (iIdxPrvSDQ = 0; iIdxPrvSDQ <= MINIMUM(3, LyrInfo.SigAphSize[iPrvElemt]); iIdxPrvSDQ++)
			{
				dCurRate = 0.0;
				for (iCurBitLoc = 0; iCurBitLoc<iIdxSDQ; iCurBitLoc++)
				{
					iCurCtx = drgiCurIncCtx[iIdxPrvSDQ][iCurBitLoc] + iCurFixCtx[iCurBitLoc];
					if (iCurCtx > 7)
					{
						iCurCtx = fnTrucCtx(iBlkRow, iBlkCol, iCurCtx, iCurBitLoc);
					}
					dCurRate += d5pSigCtxEntpy[drgiExtF[iBlkRow][iBlkCol]][iCurCtx][iCurBitLoc>0 ? 1 : 0][LyrInfo.CtxLyr[iCurElemt]][1];
				}

				dCurRate += dCurRateSign;
				if (LyrInfo.SigAphSize[iCurElemt] > iIdxSDQ)
				{
					iCurCtx = drgiCurIncCtx[iIdxPrvSDQ][iCurBitLoc] + iCurFixCtx[iCurBitLoc];
					if (iCurCtx > 7)
					{
						iCurCtx = fnTrucCtx(iBlkRow, iBlkCol, iCurCtx, iCurBitLoc);
					}

					dCurRate += d5pSigCtxEntpy[drgiExtF[iBlkRow][iBlkCol]][iCurCtx][iIdxSDQ>0 ? 1 : 0][LyrInfo.CtxLyr[iCurElemt]][0];
				}
				dCurStateCost = dCurRate + CurLyrTrellis->s_cost[iCurLyrElemtIdx - 1][iIdxPrvSDQ];

				if (dCurStateCost < dMinStateCost)
				{
					dMinStateCost = dCurStateCost;
					CurLyrTrellis->s_prvste[iCurLyrElemtIdx][iIdxSDQ] = iIdxPrvSDQ;
					CurLyrTrellis->s_rate[iCurLyrElemtIdx][iIdxSDQ] = dCurRate;
				}

				if (iPrvIdxHDQ < 3)
				{
					if (iIdxPrvSDQ == iPrvIdxHDQ && iIdxHDQ == iIdxSDQ)
					{
						CurLyrTrellis->total_rate_HDQ += dCurRate;
						CurLyrTrellis->total_cost_HDQ += (dCurRate + sqterm(dOrigCoef, fnDeQtzSig(iIdxHDQ, iFreqRow, iFreqCol)));
					}
				}
				else
				{
					if (iIdxPrvSDQ == 3 && iIdxHDQ == iIdxSDQ)
					{
						CurLyrTrellis->total_rate_HDQ += dCurRate;
						CurLyrTrellis->total_cost_HDQ += (dCurRate + sqterm(dOrigCoef, fnDeQtzSig(iIdxHDQ, iFreqRow, iFreqCol)));
					}
				}
			}

			dCurDist = sqterm(dOrigCoef, fnDeQtzSig(iIdxSDQ, iFreqRow, iFreqCol));
			if (iIdxSDQ == 0)
			{
				CurLyrTrellis->s_dist[iCurLyrElemtIdx] = dCurDist + CurLyrTrellis->s_dist[iCurLyrElemtIdx - 1];
			}

			CurLyrTrellis->s_cost[iCurLyrElemtIdx][iIdxSDQ] = dMinStateCost + dCurDist;
		}

		else
		{
			break;
		}
	}

	if (LyrInfo.SigAphSize[iCurElemt] > 2)
	{
		CurLyrTrellis->s_cost[iCurLyrElemtIdx][3] = INFINITY;
		for (iIdxSDQ = 3; iIdxSDQ <= LyrInfo.SigAphSize[iCurElemt]; iIdxSDQ++)
		{
			if (LyrInfo.SigAphSize[iCurElemt] > 3)
			{
				dCurRate = d3pAmpCtxEntpy[drgiExtF[iBlkRow][iBlkCol]][iCurElemt][iIdxSDQ - BIT_LOC_SIG];
			}
			else
			{
				dCurRate = 0.0;
			}
			dCurDist = sqterm(dOrigCoef, fnDeQtzSig(iIdxSDQ, iFreqRow, iFreqCol));
			dCurStateCost = dCurRate + dCurDist;

			if (dCurStateCost < CurLyrTrellis->s_cost[iCurLyrElemtIdx][3])
			{
				CurLyrTrellis->s_cost[iCurLyrElemtIdx][3] = dCurStateCost;
				CurLyrTrellis->s_rate[iCurLyrElemtIdx][3] = dCurRate;
				CurLyrTrellis->s_opt_amp[iCurLyrElemtIdx] = iIdxSDQ;
			}

			if (iIdxHDQ == iIdxSDQ)
			{
				CurLyrTrellis->total_rate_HDQ += dCurRate;//dCurStateCost;
				CurLyrTrellis->total_cost_HDQ += dCurStateCost;
			}
		}

		dMinStateCost = INFINITY;
		for (iIdxPrvSDQ = 0; iIdxPrvSDQ <= MINIMUM(3, LyrInfo.SigAphSize[iPrvElemt]); iIdxPrvSDQ++)
		{
			dCurRate = SDQ_LAMBDA;
			for (iCurBitLoc = 0; iCurBitLoc<3; iCurBitLoc++)
			{
				iCurCtx = drgiCurIncCtx[iIdxPrvSDQ][iCurBitLoc] + iCurFixCtx[iCurBitLoc];
				if (iCurCtx > 7)
				{
					iCurCtx = fnTrucCtx(iBlkRow, iBlkCol, iCurCtx, iCurBitLoc);
				}
				dCurRate += d5pSigCtxEntpy[drgiExtF[iBlkRow][iBlkCol]][iCurCtx][iCurBitLoc>0 ? 1 : 0][LyrInfo.CtxLyr[iCurElemt]][1];
			}

			dCurStateCost = dCurRate + CurLyrTrellis->s_cost[iCurLyrElemtIdx][3] + CurLyrTrellis->s_cost[iCurLyrElemtIdx - 1][iIdxPrvSDQ];

			if (dCurStateCost < dMinStateCost)
			{
				dMinStateCost = dCurStateCost;
				CurLyrTrellis->s_prvste[iCurLyrElemtIdx][3] = iIdxPrvSDQ;
				dMinStateRate = dCurRate;
			}

			if (iPrvIdxHDQ < 3)
			{
				if (iIdxPrvSDQ == iPrvIdxHDQ && iIdxHDQ > 2)
				{
					CurLyrTrellis->total_rate_HDQ += dCurRate;
					CurLyrTrellis->total_cost_HDQ += dCurRate;
				}
			}
			else
			{
				if (iIdxPrvSDQ == 3 && iIdxHDQ > 2)
				{
					CurLyrTrellis->total_rate_HDQ += dCurRate;
					CurLyrTrellis->total_cost_HDQ += dCurRate;
				}
			}
		}

		CurLyrTrellis->s_cost[iCurLyrElemtIdx][3] = dMinStateCost;
		CurLyrTrellis->s_rate[iCurLyrElemtIdx][3] += dMinStateRate;
	}
}

void fnInitTrellis(int iBlkRow, int iBlkCol, int iFreqRow, int iFreqCol, int iCurElemtIdx, double dOrigCoef, StateSDQ* CurLyrTrellis)
{
	int iIdxSDQ, iIdxHDQ, iCurFixCtx, iCurBitLoc = 0;
	double dCurRate, dCurRateTmp, dCurDist, dCurStateCost;

	CurLyrTrellis->total_rate_HDQ = 0.0;
	//    if (iBlkRow==3 && iBlkCol == 43 && iFreqRow == 0 && iFreqCol == 1)
	//    {
	//        iIdxSDQ=0;
	//    }

	iIdxHDQ = abs(drgiQIdxHDQ[iBlkRow * 8 + iFreqRow][iBlkCol * 8 + iFreqCol]);
	for (iIdxSDQ = 0; iIdxSDQ<3; iIdxSDQ++)
	{
		if (iIdxSDQ <= LyrInfo.SigAphSize[iCurElemtIdx])
		{
			dCurRate = iIdxSDQ == 0 ? 0.0 : SDQ_LAMBDA;//including the rate for encoding the Sign 
			for (iCurBitLoc = 0; iCurBitLoc<iIdxSDQ; iCurBitLoc++)
			{
				iCurFixCtx = fnCalFixCtxSig(iBlkRow, iBlkCol, iCurBitLoc, iCurElemtIdx);
				if (iCurFixCtx > 7)
				{
					iCurFixCtx = fnTrucCtx(iBlkRow, iBlkCol, iCurFixCtx, iCurBitLoc);
				}
				dCurRate += d5pSigCtxEntpy[drgiExtF[iBlkRow][iBlkCol]][iCurFixCtx][iCurBitLoc>0 ? 1 : 0][LyrInfo.CtxLyr[iCurElemtIdx]][1];
			}

			dCurDist = sqterm(dOrigCoef, fnDeQtzSig(iIdxSDQ, iFreqRow, iFreqCol));
			//            if (iIdxSDQ == 0)
			//            {
			//                CurLyrTrellis->s_dist[0] = dCurDist;
			//            }

			CurLyrTrellis->s_cost[0][iIdxSDQ] = dCurDist + dCurRate;
			CurLyrTrellis->s_rate[0][iIdxSDQ] = dCurRate;
			if (LyrInfo.SigAphSize[iCurElemtIdx] > iIdxSDQ)
			{
				iCurFixCtx = fnCalFixCtxSig(iBlkRow, iBlkCol, iCurBitLoc, iCurElemtIdx);
				if (iCurFixCtx > 7)
				{
					iCurFixCtx = fnTrucCtx(iBlkRow, iBlkCol, iCurFixCtx, iCurBitLoc);
				}
				CurLyrTrellis->s_cost[0][iIdxSDQ] += d5pSigCtxEntpy[drgiExtF[iBlkRow][iBlkCol]][iCurFixCtx][iIdxSDQ>0 ? 1 : 0][LyrInfo.CtxLyr[iCurElemtIdx]][0];
				CurLyrTrellis->s_rate[0][iIdxSDQ] += d5pSigCtxEntpy[drgiExtF[iBlkRow][iBlkCol]][iCurFixCtx][iIdxSDQ>0 ? 1 : 0][LyrInfo.CtxLyr[iCurElemtIdx]][0];
			}

			if (iIdxSDQ == iIdxHDQ)
			{
				CurLyrTrellis->total_rate_HDQ = CurLyrTrellis->s_rate[0][iIdxHDQ];//CurLyrTrellis->s_cost[0][iIdxHDQ];
				CurLyrTrellis->total_cost_HDQ = CurLyrTrellis->s_cost[0][iIdxHDQ];
			}
		}

		else
		{
			break;
		}
	}

	if (LyrInfo.SigAphSize[iCurElemtIdx] > 2)
	{
		iCurFixCtx = fnCalFixCtxSig(iBlkRow, iBlkCol, iCurBitLoc, iCurElemtIdx);
		if (iCurFixCtx > 7)
		{
			iCurFixCtx = fnTrucCtx(iBlkRow, iBlkCol, iCurFixCtx, iCurBitLoc);
		}
		dCurRate += d5pSigCtxEntpy[drgiExtF[iBlkRow][iBlkCol]][iCurFixCtx][1][LyrInfo.CtxLyr[iCurElemtIdx]][1];
		dCurRateTmp = dCurRate;

		CurLyrTrellis->s_cost[0][3] = INFINITY;
		for (iIdxSDQ = 3; iIdxSDQ <= LyrInfo.SigAphSize[iCurElemtIdx]; iIdxSDQ++)
		{
			if (LyrInfo.SigAphSize[iCurElemtIdx] > 3)
			{
				dCurRate = dCurRateTmp + d3pAmpCtxEntpy[drgiExtF[iBlkRow][iBlkCol]][iCurElemtIdx][iIdxSDQ - BIT_LOC_SIG];
			}
			dCurDist = sqterm(dOrigCoef, fnDeQtzSig(iIdxSDQ, iFreqRow, iFreqCol));
			dCurStateCost = dCurRate + dCurDist;

			if (dCurStateCost < CurLyrTrellis->s_cost[0][3])
			{
				CurLyrTrellis->s_cost[0][3] = dCurStateCost;
				CurLyrTrellis->s_rate[0][3] = dCurRate;
				CurLyrTrellis->s_opt_amp[0] = iIdxSDQ;
			}

			if (iIdxSDQ == iIdxHDQ)
			{
				CurLyrTrellis->total_rate_HDQ = dCurRate;// dCurStateCost;
				CurLyrTrellis->total_cost_HDQ = dCurStateCost;
			}
		}
	}
}

int fnTrucCtx(int iBlkRow, int iBlkCol, int iFixCtx, int iBitLoc)
{
	if (drgiExtF[iBlkRow][iBlkCol] == 0)
	{
		if (iBitLoc == 0)
		{
			iFixCtx = 7;
		}
		else
		{
			iFixCtx = MINIMUM(9, iFixCtx);
		}
	}
	else
	{
		if (iBitLoc == 0)
		{
			iFixCtx = 8;
		}
		else
		{
			iFixCtx = MINIMUM(10, iFixCtx);
		}
	}

	return iFixCtx;
}

void fnCalEntrpy()
{
	int i, j, k, l, m, iSumCtr, iFreqRow, iFreqCol, iIdx, iCtrTmp[8][8] = { 0 };
	double dProb, dCurEmpEnt;

	for (i = 0; i<2; i++)
	{
		for (j = 0; j<NUM_CTX_SIG; j++)
		{
			for (k = 0; k<2; k++)
			{
				for (l = 0; l<7; l++)
				{
					iSumCtr = i5pSigCtxCtr[i][j][k][l][0] + i5pSigCtxCtr[i][j][k][l][1];

					for (m = 0; m<2; m++)
					{
						if (i5pSigCtxCtr[i][j][k][l][m] == 1)
							d5pSigCtxEntpy[i][j][k][l][m] = 0.0;
						else
							d5pSigCtxEntpy[i][j][k][l][m] = SDQ_LAMBDA * (log10((double)iSumCtr) / log10(2.0) - log10((double)i5pSigCtxCtr[i][j][k][l][m]) / log10(2.0));

						if (d5pSigCtxEntpy[i][j][k][l][m] < 0.00000000001)
							d5pSigCtxEntpy[i][j][k][l][m] = SDQ_LAMBDA * 0.00000000001;
						else if (d5pSigCtxEntpy[i][j][k][l][m] < 0)
						{
							printf("Sig self-information is nagative (%2.6f)!\n", d5pSigCtxEntpy[i][j][k][l][m]);
						}

						//dRateTmp += d5pSigCtxEntpy[i][j][k][l][m]*(double)i5pSigCtxCtr[i][j][k][l][m];
						//printf("Sig self-information is %2.6f.\n", d5pSigCtxEntpy[i][j][k][l][m]);
					}
				}
			}
		}

		for (j = 0; j<LyrInfo.SigNumElemt; j++)
		{
			if (LyrInfo.SigAphSize[j] > BIT_LOC_SIG)
			{
				iSumCtr = 0;
				for (k = 0; k<LyrInfo.SigAphSize[j] - BIT_LOC_SIG + 1; k++)
				{
					iSumCtr += i3pAmpCtxCtr[i][LyrInfo.CtxLyr[j]][k];
				}

				for (k = 0; k<LyrInfo.SigAphSize[j] - BIT_LOC_SIG + 1; k++)
				{
					if (i3pAmpCtxCtr[i][LyrInfo.CtxLyr[j]][k] == 1)
						d3pAmpCtxEntpy[i][j][k] = 0.0;
					else
						d3pAmpCtxEntpy[i][j][k] = SDQ_LAMBDA * (log10((double)iSumCtr) / log10(2.0) - log10((double)i3pAmpCtxCtr[i][LyrInfo.CtxLyr[j]][k]) / log10(2.0));

					if (d3pAmpCtxEntpy[i][j][k] < 0.00000000001)
						d3pAmpCtxEntpy[i][j][k] = SDQ_LAMBDA * 0.00000000001;
					else if (d3pAmpCtxEntpy[i][j][k] < 0)
					{
						printf("Amp self-information is nagative (%2.6f)!\n", d3pAmpCtxEntpy[i][j][k]);
					}

					//dRateTmp += d3pAmpCtxEntpy[i][j][k]*(double)i3pAmpCtxCtr[i][j][k];
					//printf("Amp self-information is %2.6f.\n", d3pAmpCtxEntpy[i][j][k]);
				}
			}
		}
	}

#if TEST_RD
	if (iNumIter == 0)
	{
		for (iIdx = 0; iIdx<LyrInfo.SigNumElemt; iIdx++)
		{
			iFreqRow = LyrInfo.SigScanRowIdx[iIdx];
			iFreqCol = LyrInfo.SigScanColIdx[iIdx];
			l = LyrInfo.CtxLyr[iIdx];

			iSumCtr = 0;
			for (k = 0; k <= LyrInfo.SigAphSize[iIdx]; k++)
			{
				iSumCtr += i2pAmpCtrEmp[iIdx][k];
			}
			iCtrTmp[iFreqRow][iFreqCol] = iSumCtr;

			for (k = 0; k <= LyrInfo.SigAphSize[iIdx]; k++)
			{
				if (i2pAmpCtrEmp[iIdx][k] > 0)
				{
					dProb = (double)i2pAmpCtrEmp[iIdx][k] / (double)iSumCtr;
					dCurEmpEnt = -dProb * log10(dProb) / log10(2.0);
					drgdEmpEnt[iFreqRow][iFreqCol] += dCurEmpEnt;
				}

				if (drgdEmpEnt[iFreqRow][iFreqCol] < 0)
				{
					printf("Amp emprical entropy is nagative (%2.6f)!\n", drgdEmpEnt[iFreqRow][iFreqCol]);
				}
			}

			//iCtrTmp[iFreqRow][iFreqCol] = 0;
			for (i = 0; i<2; i++)
			{
				for (j = 0; j<NUM_CTX_SIG; j++)
				{
					for (k = 0; k<2; k++)
					{
						for (m = 0; m<2; m++)
						{
							drgdRealRate[iFreqRow][iFreqCol] += d5pSigCtxEntpy[i][j][k][l][m] * (double)i5pSigCtxElemtCtr[i][j][k][iIdx][m] / SDQ_LAMBDA;
						}
					}
				}

				for (k = 0; k<LyrInfo.SigAphSize[iIdx] - BIT_LOC_SIG + 1; k++)
				{
					//            dProb = (double)i3pAmpCtxCtr[i][LyrInfo.CtxLyr[j]][k]/(double)iSumCtr;
					//            dCurEmpEnt = -dProb*log10(dProb)/log10(2.0);
					drgdRealRate[iFreqRow][iFreqCol] += d3pAmpCtxEntpy[i][iIdx][k] * (double)i3pAmpCtxElemtCtr[i][iIdx][k] / SDQ_LAMBDA;
				}
			}
		}

		printf("---------------Emperical Entropy (zero order)--------------\n");
		for (iFreqRow = 0; iFreqRow<8; iFreqRow++)
		{
			for (iFreqCol = 0; iFreqCol<8; iFreqCol++)
			{
				if (drgiLTbl[iFreqRow][iFreqCol] > 0 && iFreqRow + iFreqCol>0)
					printf("%3.4f ", drgdEmpEnt[iFreqRow][iFreqCol]);
				else
					printf("%3.4f ", 0.0);
			}

			printf("\n");
		}

		printf("---------------Emperical Entropy (with context)--------------\n");
		for (iFreqRow = 0; iFreqRow<8; iFreqRow++)
		{
			for (iFreqCol = 0; iFreqCol<8; iFreqCol++)
			{
				if (drgiLTbl[iFreqRow][iFreqCol] > 0 && (iFreqRow + iFreqCol) > 0)
				{
					//                    if (iCtrTmp[iFreqRow][iFreqCol] != drgiNumBlkWithinYc[iFreqRow][iFreqCol])
					//                        printf("(%d,%d): %d, %d\n", iFreqRow, iFreqCol, iCtrTmp[iFreqRow][iFreqCol], drgiNumBlkWithinYc[iFreqRow][iFreqCol]);
					printf("%3.4f ", drgdRealRate[iFreqRow][iFreqCol] / drgiNumBlkWithinYc[iFreqRow][iFreqCol]);
				}

				else
					printf("%3.4f ", 0.0);
			}

			printf("\n");
		}
	}
#endif
}

int fnCalFixCtxSig(int iBlkRow, int iBlkCol, int iCurBitLoc, int iElemtIdx)
{
	int iCtxIdx, iCtxNum = 0;
	int i, j, iCtxNumBlkExt = 0, iRowZZMinus1, iColZZMinus1, iCtxNumSubndNeibr = 0;

	if (LyrInfo.SigNumCtx[iElemtIdx] > 0)
	{
		for (iCtxIdx = 0; iCtxIdx<LyrInfo.SigNumCtx[iElemtIdx]; iCtxIdx++)
		{
			iCtxNum += prgiSig[iBlkRow][iBlkCol][iCurBitLoc][LyrInfo.SigCtxMdl[iElemtIdx][iCtxIdx]];
		}
	}

	if (drgiExtF[iBlkRow][iBlkCol] == 1)
	{
		i = LyrInfo.SigScanRowIdx[iElemtIdx];
		j = LyrInfo.SigScanColIdx[iElemtIdx];

		if (i == 0)
		{
			if (j>0)
				iCtxNumBlkExt = vrgiExt[iBlkRow][iBlkCol][i][j - 1][0];
			else
				iCtxNumBlkExt = 0;
		}

		else if (j == 1)
			iCtxNumBlkExt = vrgiExt[iBlkRow][iBlkCol][i - 1][j][0];

		else
			iCtxNumBlkExt = vrgiExt[iBlkRow][iBlkCol][i - 1][j][0] + vrgiExt[iBlkRow][iBlkCol][i][j - 1][0] + vrgiExt[iBlkRow][iBlkCol][i - 1][j - 1][0];

		if (iElemtIdx > 0)
		{
			iRowZZMinus1 = LyrInfo.SigScanRowIdx[iElemtIdx - 1];
			iColZZMinus1 = LyrInfo.SigScanColIdx[iElemtIdx - 1];

			if ((iRowZZMinus1 + iColZZMinus1) == (i + j))
			{
#if IMMED_NEIBR
				if (iRowZZMinus1 == i + 1 || iRowZZMinus1 == i - 1)
					iCtxNumBlkExt += vrgiExt[iBlkRow][iBlkCol][iRowZZMinus1][iColZZMinus1][0];
#else
				iCtxNumBlkExt += vrgiExt[iBlkRow][iBlkCol][iRowZZMinus1][iColZZMinus1][0];
#endif
			}
		}

		iCtxNum += iCtxNumBlkExt;
	}

#if INTER_BLK_SIG_HELP
	//explore the correlation among 6 subband neighbors of the current encoding coeff
	if (iBlkRow == 0)
	{
		if (iBlkCol == 0)
			iCtxNumSubndNeibr = 0;

		else
		{
			iCtxNumSubndNeibr = prgiSig[iBlkRow][iBlkCol - 1][iCurBitLoc][iElemtIdx];
			if (iBlkCol>1)
			{
				iCtxNumSubndNeibr += prgiSig[iBlkRow][iBlkCol - 2][iCurBitLoc][iElemtIdx];
			}
		}
	}

	else if (iBlkCol == 0)
	{
		iCtxNumSubndNeibr = prgiSig[iBlkRow - 1][iBlkCol][iCurBitLoc][iElemtIdx] + prgiSig[iBlkRow - 1][iBlkCol + 1][iCurBitLoc][iElemtIdx];
		if (iBlkRow>1)
		{
			iCtxNumSubndNeibr += prgiSig[iBlkRow - 2][iBlkCol][iCurBitLoc][iElemtIdx];
		}
	}

	else if (iBlkCol == (COLS / 8 - 1))
	{
		iCtxNumSubndNeibr = prgiSig[iBlkRow - 1][iBlkCol][iCurBitLoc][iElemtIdx] + prgiSig[iBlkRow][iBlkCol - 1][iCurBitLoc][iElemtIdx] + prgiSig[iBlkRow - 1][iBlkCol - 1][iCurBitLoc][iElemtIdx] + prgiSig[iBlkRow][iBlkCol - 2][iCurBitLoc][iElemtIdx];
		if (iBlkRow>1)
			iCtxNumSubndNeibr += prgiSig[iBlkRow - 2][iBlkCol][iCurBitLoc][iElemtIdx];
	}

	else
	{
		iCtxNumSubndNeibr = prgiSig[iBlkRow - 1][iBlkCol][iCurBitLoc][iElemtIdx] + prgiSig[iBlkRow][iBlkCol - 1][iCurBitLoc][iElemtIdx] + prgiSig[iBlkRow - 1][iBlkCol - 1][iCurBitLoc][iElemtIdx] + prgiSig[iBlkRow - 1][iBlkCol + 1][iCurBitLoc][iElemtIdx];
		if (iBlkRow>1)
			iCtxNumSubndNeibr += prgiSig[iBlkRow - 2][iBlkCol][iCurBitLoc][iElemtIdx];

		if (iBlkCol>1)
			iCtxNumSubndNeibr += prgiSig[iBlkRow][iBlkCol - 2][iCurBitLoc][iElemtIdx];
	}

	iCtxNum += iCtxNumSubndNeibr;
#endif

	return iCtxNum;
}


void fnFreeMemHACACM()
{
	int i, j, k;

	for (i = 0; i<NUM_CTX_EXTAMP_SEPRT; i++)
	{
		free(ExtAmpSeprtCtx[i]);
	}

	for (i = 0; i<NUM_CTX_EXT; i++)
	{
		for (j = 0; j<BIT_LOC_EXT; j++)
		{
			free(ExtCtx[i][j]);
		}
	}

	for (i = 0; i<2; i++)
	{
		for (j = 0; j<NUM_CTX_SIG; j++)
		{
			for (k = 0; k<BIT_LOC_SIG; k++)
			{
				free(SigCtx[i][j][k]);
			}
		}
	}

	for (i = 0; i<2; i++)
	{
		for (j = 0; j<NUM_CTX_AMP; j++)
		{
			free(AmpCtx[i][j]);
		}
	}

	for (i = 0; i<NUM_CTX_EXTAMP; i++)
	{
		free(ExtAmpCtx[i]);
	}

	if (SDQ)
	{
		for (i = 0; i<2; i++)
		{
			for (j = 0; j<7; j++)
			{
				free(i3pAmpCtxCtr[i][j]);
			}

			for (j = 0; j<LyrInfo.SigNumElemt; j++)
			{
				free(d3pAmpCtxEntpy[i][j]);
			}

			free(d3pAmpCtxEntpy[i]);
		}

		for (i = 0; i<ROWS / 8; i++)
		{
			for (j = 0; j<COLS / 8; j++)
			{
				for (k = 0; k<BIT_LOC_SIG; k++)
				{
					free(prgiSig[i][j][k]);
				}
			}
		}
#if TEST_RD
		for (j = 0; j<LyrInfo.SigNumElemt; j++)
		{
			free(i2pAmpCtrEmp[j]);
		}
		free(i2pAmpCtrEmp);

		for (i = 0; i<2; i++)
		{
			for (j = 0; j<LyrInfo.SigNumElemt; j++)
			{
				free(i3pAmpCtxElemtCtr[i][j]);
			}
			free(i3pAmpCtxElemtCtr[i]);
		}
#endif
	}
}
