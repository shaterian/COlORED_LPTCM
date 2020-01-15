/*
* This file is responsible for the quantization 
*  AdapRecSpace.c
*  mjpeg
*	
*  Created by Chang Sun on 11-07-20.
*  Copyright 2011 Univ. of Waterloo. All rights reserved.
*
*/
#pragma warning(disable:4996)
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "ImageDependent.h"
#include "AdapRecSpace.h"
#include "HACACM.h"
#include <iostream>
double drgdGamma[8][8];
double **Q_cache;
OptDzoneTbl OptDzoneQtzr;

int iLLowDC, iLHighDC;
double drgdVar[8][8];
int drgdQprimeTbl[8][8];
int drgiLTbl[8][8];
int drgiLPrimeTbl[8][8];
int drgiTTbl[8][8] = { 0 };
std::vector<std::vector<int>> drgiExtF; // [ROWS / 8] [COLS / 8] = { 0 }
//int drgiExtFDec[ROWS / 8][COLS / 8] = { 0 };
std::vector<std::vector<int>> drgiExtFDec; 
//int drgiCBP[ROWS / 8][COLS / 8] = { 0 };

std::vector<std::vector<int>>drgiCBP;

int drgiNumBlkWithinYc[8][8] = { 0 };
double drgdSumYmulY[8][8] = { 0.0 };
double dAdOnVrnceSDQ = 0.0;
int iNumAdOnVrnceSDQ = 0;

char rgcRecFName[256];
int iMaxLprime, iMaxLDC, iMaxAbsDC;

int iOptQPrime, iOptQ;
double rgdGammaCache[NUM_QTZ_GAMMA];
double rgdLambdaCache[NUM_QTZ_LAMBDA];

double drgdTarRate[8][8] = { 0.0 };
double drgdTarDist[8][8] = { 0.0 };
double drgdRealRate[8][8] = { 0.0 };
double drgdRealDist[8][8] = { 0.0 };
double drgdOptMSE[8][8] = { 0.0 };
int drgiOptPSNR[8][8] = { 0 };

double drgdRealVSModlRate[8][8] =
{
	{ 1.0000,0.8020,0.8020,0.8034,0.8298,0.8916,0.8920,0.9038 },
{ 0.7986,0.8054,0.8060,0.8208,0.9167,0.9583,1.0067,1.0907 },
{ 0.8125,0.8088,0.8285,0.8877,0.8833,0.9455,1.1062,1.1452 },
{ 0.8842,0.8880,0.9036,0.9152,0.9742,1.3037,1.5145,1.0000 },
{ 0.9446,0.8615,0.9218,0.9926,1.2390,1.4402,1.0000,1.0000 },
{ 1.0061,1.0620,1.2394,1.2727,1.8805,1.0000,1.0000,1.0000 },
{ 1.0000,1.2158,1.5212,1.5559,1.0000,1.0000,1.0000,1.0000 },
{ 1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 }
};

double drgdRealVSModlDist[8][8] =
{
	{ 1.0000,1.0412,0.9994,0.9929,0.9341,0.9144,0.8463,0.8411 },
{ 1.0405,1.0146,0.9766,0.9154,0.8802,0.8181,0.8012,0.8539 },
{ 0.9717,0.9540,0.9055,0.9388,0.8335,0.8057,0.8830,0.8598 },
{ 0.8398,0.9247,0.9042,0.8199,0.7648,0.8117,0.8475,0.8004 },
{ 0.8506,0.8170,0.7839,0.7510,0.8725,0.8072,0.7552,0.7532 },
{ 0.8758,0.8568,0.8243,0.7976,0.7788,0.7632,0.7453,0.7540 },
{ 0.7663,0.8203,0.7725,0.7969,0.7783,0.7823,1.0000,1.0000 },
{ 0.7671,0.7689,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 }
};

// this function gets a freq location and for that locatio 
// returns the quantized version of Gamma , Lambda , Yc
void fnParaTCMQtz(int iRow, int iCol)
{
	int i{ 0 };
	int iIdxGamma{ 0 };
	int iIdxLambda{0};
	double dQtzDist, dDiff;
	FILE *fpIF;

	drgdGamma[iRow][iCol] = drgdYc[iRow][iCol] / drgdLambda[iRow][iCol];

	// finding the nearest value from GammaChache
	dDiff = 1000.0;
	for (i = 0; i<NUM_QTZ_GAMMA; i++)
	{
		dQtzDist = fabs(drgdGamma[iRow][iCol] - rgdGammaCache[i]);
		if (dQtzDist < dDiff)
		{
			dDiff = dQtzDist;
			iIdxGamma = i;
		}
	}
	// finding the nearest value from GammaChache
	dDiff = 1000.0;
	for (i = 0; i<NUM_QTZ_LAMBDA; i++)
	{
		dQtzDist = fabs(drgdLambda[iRow][iCol] - rgdLambdaCache[i]);
		if (dQtzDist < dDiff)
		{
			dDiff = dQtzDist;
			iIdxLambda = i;
		}
	}

	OptDzoneQtzr.Gamma_idx[iRow][iCol] = iIdxGamma;
	OptDzoneQtzr.Lambda_idx[iRow][iCol] = iIdxLambda;

	/* quantized Gamma, Lamnda and Yc */
	drgdGamma[iRow][iCol] = rgdGammaCache[iIdxGamma];
	drgdLambda[iRow][iCol] = rgdLambdaCache[iIdxLambda];
	drgdYc[iRow][iCol] = drgdGamma[iRow][iCol] * drgdLambda[iRow][iCol];
}
// optimum distortion profile
void fnOptDistProfl()
{
	int i, j, k, iCurTarPSNR, iPrvTarPSNR, iNumDist, rgiTarPSNR[NUM_DIST_TRELIS];
	double dCurCost, dMinCost, rgdTarMSE[NUM_DIST_TRELIS], **R_cache, **D_cache, rgdRateTCM[8][8][NUM_DIST_TRELIS], rgdDistTCM[8][8][NUM_DIST_TRELIS];
	char rgcIFNameRtbl[256] = { 0 };
	char rgcIFNameDtbl[256] = { 0 };
	FILE *fpIF;

	R_cache = (double **)malloc(sizeof(double *)*NUM_QTZ_GAMMA);
	for (i = 0; i<NUM_QTZ_GAMMA; i++)
	{
		R_cache[i] = (double *)malloc(sizeof(double)*NUM_QTZ_LAMBDA);
	}

	D_cache = (double **)malloc(sizeof(double *)*NUM_QTZ_GAMMA);
	for (i = 0; i<NUM_QTZ_GAMMA; i++)
	{
		D_cache[i] = (double *)malloc(sizeof(double)*NUM_QTZ_LAMBDA);
	}

	iPrvTarPSNR = TAR_PSNR_ARG - ((int)ceil((double)NUM_DIST_TRELIS / 2.0) - 1) * 25;

	iNumDist = 0;
	while (iNumDist < NUM_DIST_TRELIS)
	{
		iCurTarPSNR = iPrvTarPSNR;
		sprintf(rgcIFNameRtbl, "/Users/changsun/Dropbox/Work/MATLAB/training_images/TCM_parameters/CnstrDzone_Rate_PSNR%d.txt", iCurTarPSNR);
		sprintf(rgcIFNameDtbl, "/Users/changsun/Dropbox/Work/MATLAB/training_images/TCM_parameters/CnstrDzone_Dist_PSNR%d.txt", iCurTarPSNR);
		/* read the R and D tables from the files correspond to the target distortion */
		fpIF = fopen(rgcIFNameRtbl, "r");
		for (i = 0; i<NUM_QTZ_GAMMA; i++)
		{
			for (j = 0; j<NUM_QTZ_LAMBDA; j++)
			{
				fscanf(fpIF, "%lf ", &R_cache[i][j]);
			}
		}
		fclose(fpIF);

		fpIF = fopen(rgcIFNameDtbl, "r");
		for (i = 0; i<NUM_QTZ_GAMMA; i++)
		{
			for (j = 0; j<NUM_QTZ_LAMBDA; j++)
			{
				fscanf(fpIF, "%lf ", &D_cache[i][j]);
			}
		}
		fclose(fpIF);

		for (i = 0; i<8; i++)
		{
			for (j = 0; j<8; j++)
			{
				if (i + j > 0)
				{
					rgdRateTCM[i][j][iNumDist] = R_cache[OptDzoneQtzr.Gamma_idx[i][j]][OptDzoneQtzr.Lambda_idx[i][j]];
					rgdDistTCM[i][j][iNumDist] = D_cache[OptDzoneQtzr.Gamma_idx[i][j]][OptDzoneQtzr.Lambda_idx[i][j]];
					//                    if (i==1 && j==0)
					//                        printf("Gamma=%4.4f, Lambda=%4.4f, R=%2.4f, D=%4.4f\n", drgdGamma[i][j], drgdLambda[i][j], rgdRateTCM[i][j][iNumDist], rgdDistTCM[i][j][iNumDist]);
				}
			}
		}

		rgdTarMSE[iNumDist] = PSNR2MSE((double)iCurTarPSNR / 100.0);
		rgiTarPSNR[iNumDist] = iCurTarPSNR;

		iPrvTarPSNR += 25;
		iNumDist++;
	}

	for (i = 0; i<8; i++)
	{
		for (j = 0; j<8; j++)
		{
			if (i + j > 0)
			{
				dMinCost = INFINITY;
				for (k = 0; k<NUM_DIST_TRELIS; k++)
				{
					//                    dCurCost = SDQ_LAMBDA*drgdRealVSModlRate[i][j]*rgdRateTCM[i][j][k] + drgdRealVSModlDist[i][j]*rgdDistTCM[i][j][k];
					dCurCost = SDQ_LAMBDA * rgdRateTCM[i][j][k] + rgdDistTCM[i][j][k];

					if (dMinCost > dCurCost)
					{
						dMinCost = dCurCost;
						drgdOptMSE[i][j] = rgdTarMSE[k];
						drgiOptPSNR[i][j] = rgiTarPSNR[k];
						drgdTarRate[i][j] = rgdRateTCM[i][j][k];
						drgdTarDist[i][j] = rgdDistTCM[i][j][k];
					}
				}
			}
		}
	}

	for (i = 0; i<NUM_QTZ_GAMMA; i++)
	{
		free(R_cache[i]);
		free(D_cache[i]);
	}
	free(R_cache);
	free(D_cache);
}


/*
*  Function: fnOptMdlQ
*  Find the optimal quantization step size achieving the target distortion
*/
void fnOptMdlQ(double **drgdCoef)
{
	int i, j, iRow, iCol, iLprime, iCurTarPSNR;
	double dQprime, dVar, dQStep, dVarTmp, dTAR_DIST;

	int **L_cache;
	FILE *fpDF, *fpIF;

#if DIFFDIST
	TAR_PSNR_ARG = 3000;
#endif
	// a 2D matrix of 64*256
	L_cache = (int **)malloc(sizeof(int *)*NUM_QTZ_GAMMA);
	for (i = 0; i<NUM_QTZ_GAMMA; i++)
	{
		L_cache[i] = (int *)malloc(sizeof(int)*NUM_QTZ_LAMBDA);
	}
	// a 2D matrix of 64*256
	Q_cache = (double **)malloc(sizeof(double *)*NUM_QTZ_GAMMA);
	for (i = 0; i<NUM_QTZ_GAMMA; i++)
	{
		Q_cache[i] = (double *)malloc(sizeof(double)*NUM_QTZ_LAMBDA);
	}

	/* read the Gamma and Lambda tables from the files */
	fpIF = fopen("D:\\Research\\new_codec_test\\images\\GammaCache.txt", "r");
	for (i = 0; i<NUM_QTZ_GAMMA; i++)
	{
		fscanf(fpIF, "%lf", &rgdGammaCache[i]);
	}
	fclose(fpIF);

	fpIF = fopen("D:\\Research\\new_codec_test\\images\\LambdaCache.txt", "r");
	for (i = 0; i<NUM_QTZ_LAMBDA; i++)
	{
		fscanf(fpIF, "%lf", &rgdLambdaCache[i]);
	}
	fclose(fpIF);

	iMaxLprime = 0;
	//dTAR_DIST = TAR_DIST;
	for (iRow = 0; iRow<8; iRow++)
	{
		for (iCol = 0; iCol<8; iCol++)
		{
			fnParaTCMQtz(iRow, iCol);
		}
	}
	/*printf("---------------Quantized Yc--------------\n");
	for (i = 0; i<8; i++)
	{
		for (j = 0; j<8; j++)
		{
			printf("%6.2f ", drgdYc[i][j]);
		}
		printf("\n");
	}
*/
	if (ADPT_DIST_PROFL)
		fnOptDistProfl();

	else
	{
#if HDQ_EMP_DIST_PROFL
#else

		/*sprintf(rgcIFNameLtbl, "D:\\Research\\new_codec_test\\my_data\\CnstrDzone_OptL_PSNR%d.txt" , TAR_PSNR_ARG);
		sprintf(rgcIFNameQtbl, "D:\\Research\\new_codec_test\\my_data\\CnstrDzone_OptQ_PSNR%d.txt" , TAR_PSNR_ARG);
		*/
		sprintf(rgcIFNameLtbl, "D:\\Research\\New_Optimization_Qtables\\files\\Ltables\\CnstrDzone_OptL_PSNR%d.txt", TAR_PSNR_ARG);
		sprintf(rgcIFNameQtbl, "D:\\Research\\New_Optimization_Qtables\\files\\Qtables\\CnstrDzone_OptQ_PSNR%d.txt", TAR_PSNR_ARG);
		std::cout << rgcIFNameQtbl << std::endl;
		/* read the L and Q tables from the files correspond to the target distortion */
		fpIF = fopen(rgcIFNameLtbl, "r");
		for (i = 0; i<NUM_QTZ_GAMMA; i++)
		{
			for (j = 0; j<NUM_QTZ_LAMBDA; j++)
			{
				fscanf(fpIF, "%d ", &L_cache[i][j]);
			}
		}
		fclose(fpIF);

		fpIF = fopen(rgcIFNameQtbl, "r");
		for (i = 0; i<NUM_QTZ_GAMMA; i++)
		{
			for (j = 0; j<NUM_QTZ_LAMBDA; j++)
			{
				fscanf(fpIF, "%lf ", &Q_cache[i][j]);
			}
		}
		fclose(fpIF);
#endif
	}

	for (iRow = 0; iRow<8; iRow++)
	{
		for (iCol = 0; iCol<8; iCol++)
		{
			if ((iRow + iCol) == 0) // DC coefficients
			{
				OptDzoneQtzr.QTbl[iRow][iCol] = sqrt(TAR_DIST*12.0);

				iLHighDC = (int)floor(dMaxDC / OptDzoneQtzr.QTbl[iRow][iCol] + 0.5);
				iLLowDC = (int)floor(fabs(dMinDC) / OptDzoneQtzr.QTbl[iRow][iCol] + 0.5);

				iMaxAbsDC = MAXIMUM(iLHighDC, iLLowDC);

				iMaxLDC = iLHighDC + iLLowDC;//MAXIMUM(2*iLHighDC+iLLowDC,iLHighDC+2*iLLowDC);

				dQprime = OptDzoneQtzr.QTbl[iRow][iCol];//the quantization step size for uniform part accross all frequencies
														//is selected to be equal to that of DC's

				drgiLTbl[iRow][iCol] = iMaxLDC;
			}

			else
			{
				//                if (i+j<6)
				//                    TAR_DIST = dTAR_DIST;
				//                else 
				//                    TAR_DIST = (255*255/pow(10.0,42/10.0));

				dVar = 0.0;

				for (i = 0; i < ROWS / 8; i++)
				{
					for (j = 0; j < COLS / 8; j++)
					{
						dVar += drgdCoef[i * 8 + iRow][j * 8 + iCol] * drgdCoef[i * 8 + iRow][j * 8 + iCol];
					}
				}

				dVar /= NUMBLK;
				drgdVar[iRow][iCol] = dVar;

				if (ADPT_DIST_PROFL)
				{
					iCurTarPSNR = drgiOptPSNR[iRow][iCol];
					dTAR_DIST = drgdOptMSE[iRow][iCol];
				}

				else
				{
#if HDQ_EMP_DIST_PROFL
					if (TAR_PSNR_ARG < 4000)
					{
						if (iRow + iCol == 1)
						{
							iCurTarPSNR = TAR_PSNR_ARG;
						}
						else if (iRow + iCol == 2 || iRow + iCol == 3)
						{
							iCurTarPSNR = TAR_PSNR_ARG + 25;//- 25;
						}
						else if (iRow + iCol == 4 || iRow + iCol == 5)
						{
							iCurTarPSNR = TAR_PSNR_ARG + 50;// - 50;
						}
						else if (iRow + iCol == 6 || iRow + iCol == 7)
						{
							iCurTarPSNR = TAR_PSNR_ARG + 75;// - 50;
						}

						else if (iRow + iCol == 8 || iRow + iCol == 9)
						{
							iCurTarPSNR = TAR_PSNR_ARG + 100;// - 50;
						}

						else if (iRow + iCol == 10 || iRow + iCol == 11)
						{
							iCurTarPSNR = TAR_PSNR_ARG + 125;// - 50;
						}

						else
						{
							iCurTarPSNR = TAR_PSNR_ARG + 150;
						}
					}

					else
					{
						iCurTarPSNR = TAR_PSNR_ARG;
					}

					dTAR_DIST = PSNR2MSE((double)iCurTarPSNR / 100.0);
#else
					dTAR_DIST = TAR_DIST;
#endif
				}

				if (dTAR_DIST >= dVar)
				{
					OptDzoneQtzr.QTbl[iRow][iCol] = 0;
					drgiLTbl[iRow][iCol] = 0;
					drgiTTbl[iRow][iCol] = 0;
					//#if ADPT_DIST_PROFL
					//                    dAdOnVrnceSDQ += dVar;
					//                    iNumAdOnVrnceSDQ++;
					//#endif
				}

				else
				{
					if (ADPT_DIST_PROFL || HDQ_EMP_DIST_PROFL)
					{
						///* read the L and Q tables from the files correspond to the target distortion */
						//sprintf(rgcIFNameLtbl, "D:\\Research\\New_Codec_chung_code\\encoder_proj\\images\\Dropbox data\\CnstrDzone_OptL_PSNR4275.txt");
						//sprintf(rgcIFNameQtbl, "D:\\Research\\New_Codec_chung_code\\encoder_proj\\images\\Qtable.txt");
						//sprintf(rgcIFNameLtbl, "D:\\Users\\sshateri\\Desktop\\new_codec_test\\my_data\\CnstrDzone_OptL_PSNR%d.txt" , TAR_PSNR_ARG);
						//sprintf(rgcIFNameQtbl, "D:\\Users\\sshateri\\Desktop\\new_codec_test\\my_data\\CnstrDzone_OptQ_PSNR%d.txt" , TAR_PSNR_ARG);
						sprintf(rgcIFNameLtbl, "D:\\Research\\New_Optimization_Qtables\\files\\Ltable\\CnstrDzone_OptL_PSNR%d.txt", iCurTarPSNR); //TAR_PSNR_ARG);
						sprintf(rgcIFNameQtbl, "D:\\Research\\New_Optimization_Qtables\\files\\Qtable\\CnstrDzone_OptQ_PSNR%d.txt", iCurTarPSNR);// TAR_PSNR_ARG);
						//std::cout << rgcIFNameQtbl << std::endl; 
						fpIF = fopen(rgcIFNameLtbl, "r");
						for (i = 0; i<NUM_QTZ_GAMMA; i++)
						{
							for (j = 0; j<NUM_QTZ_LAMBDA; j++)
							{
								fscanf(fpIF, "%d ", &L_cache[i][j] );
							}
						}
						fclose(fpIF);

						fpIF = fopen(rgcIFNameQtbl, "r");
						for (i = 0; i<NUM_QTZ_GAMMA; i++)
						{
							for (j = 0; j<NUM_QTZ_LAMBDA; j++)
							{
								fscanf(fpIF, "%lf ", &Q_cache[i][j]);
							}
						}
						fclose(fpIF);
					
					}
					if (   TAR_PSNR_ARG < 4400 && TAR_PSNR_ARG > 2800 )
					{
						OptDzoneQtzr.LTbl[iRow][iCol] = MAXIMUM(1, L_cache[OptDzoneQtzr.Gamma_idx[iRow][iCol]][OptDzoneQtzr.Lambda_idx[iRow][iCol]]);
						OptDzoneQtzr.QTbl[iRow][iCol] = Q_cache[OptDzoneQtzr.Gamma_idx[iRow][iCol]][OptDzoneQtzr.Lambda_idx[iRow][iCol]];
						OptDzoneQtzr.UTbl[iRow][iCol] = drgdGamma[iRow][iCol] - (OptDzoneQtzr.LTbl[iRow][iCol] - 1) * OptDzoneQtzr.QTbl[iRow][iCol];
						OptDzoneQtzr.DeltaTbl[iRow][iCol] = OptDzoneQtzr.UTbl[iRow][iCol] + 1 + OptDzoneQtzr.QTbl[iRow][iCol] / (1 - exp(OptDzoneQtzr.QTbl[iRow][iCol]));

						drgiTTbl[iRow][iCol] = 1;
						drgiLPrimeTbl[iRow][iCol] = (int)ceil((drgdAval[iRow][iCol] - drgdYc[iRow][iCol]) / dQprime);
						iMaxLprime = MAXIMUM(iMaxLprime, drgiLPrimeTbl[iRow][iCol]);
						drgiLTbl[iRow][iCol] = OptDzoneQtzr.LTbl[iRow][iCol] - 1;

						/* scale the Q, U and Delta according to Lambda */
						OptDzoneQtzr.QTbl[iRow][iCol] *= drgdLambda[iRow][iCol];
						OptDzoneQtzr.UTbl[iRow][iCol] *= drgdLambda[iRow][iCol];
						OptDzoneQtzr.DeltaTbl[iRow][iCol] *= drgdLambda[iRow][iCol];
					}
					else
					{
						OptDzoneQtzr.LTbl[iRow][iCol] = MAXIMUM(1, L_cache[OptDzoneQtzr.Gamma_idx[iRow][iCol]][OptDzoneQtzr.Lambda_idx[iRow][iCol]]);
						OptDzoneQtzr.QTbl[iRow][iCol] = Q_cache[OptDzoneQtzr.Gamma_idx[iRow][iCol]][OptDzoneQtzr.Lambda_idx[iRow][iCol]];
						OptDzoneQtzr.UTbl[iRow][iCol] = drgdGamma[iRow][iCol] - (OptDzoneQtzr.LTbl[iRow][iCol] ) * OptDzoneQtzr.QTbl[iRow][iCol];
						OptDzoneQtzr.DeltaTbl[iRow][iCol] = OptDzoneQtzr.UTbl[iRow][iCol] + 1 + OptDzoneQtzr.QTbl[iRow][iCol] / (1 - exp(OptDzoneQtzr.QTbl[iRow][iCol]));

						drgiTTbl[iRow][iCol] = 1;
						drgiLPrimeTbl[iRow][iCol] = (int)ceil((drgdAval[iRow][iCol] - drgdYc[iRow][iCol]) / dQprime);
						iMaxLprime = MAXIMUM(iMaxLprime, drgiLPrimeTbl[iRow][iCol]);
						drgiLTbl[iRow][iCol] = OptDzoneQtzr.LTbl[iRow][iCol];

						/* scale the Q, U and Delta according to Lambda */
						OptDzoneQtzr.QTbl[iRow][iCol] *= drgdLambda[iRow][iCol];
						OptDzoneQtzr.UTbl[iRow][iCol] *= drgdLambda[iRow][iCol];
						OptDzoneQtzr.DeltaTbl[iRow][iCol] *= drgdLambda[iRow][iCol];
					}
					//drgiNumBlkWithinYc[iRow][iCol] = (int)(drgdBval[iRow][iCol]*NUMBLK+.5);//this value is approximated by B table which will be used in Q table updating. Note this value doesn't have to be accurate when calculating the distortion.
				}
			}
		}
	}

	//#if ADPT_DIST_PROFL  
	//    if (iNumAdOnVrnceSDQ > 0)
	//        dAdOnVrnceSDQ /= (double)(63-iNumAdOnVrnceSDQ);
	//    //dAdOnVrnceSDQ = 0;
	//#endif

	printf("hell ya");
	printf("---------------T table--------------\n");
	for (iRow = 0; iRow<8; iRow++)
	{
		for (iCol = 0; iCol<8; iCol++)
		{
			printf("%1d ", drgiTTbl[iRow][iCol]);
		}

		printf("\n");
	}

	printf("---------------L table--------------\n");
	for (iRow = 0; iRow<8; iRow++)
	{
		for (iCol = 0; iCol<8; iCol++)
		{
			printf("%d ", drgiLTbl[iRow][iCol]);
		}

		printf("\n");
	}

	printf("---------------L' table--------------\n");
	for (iRow = 0; iRow<8; iRow++)
	{
		for (iCol = 0; iCol<8; iCol++)
		{
			printf("%d ", drgiLPrimeTbl[iRow][iCol]);
		}

		printf("\n");
	}
	printf("oh ya ");
	printf("---------------Q table--------------\n");
	for (iRow = 0; iRow<8; iRow++)
	{
		for (iCol = 0; iCol<8; iCol++)
		{
			printf("%3.4f ", OptDzoneQtzr.QTbl[iRow][iCol]);
		}

		printf("\n");
	}

	printf("---------------U table--------------\n");
	for (iRow = 0; iRow<8; iRow++)
	{
		for (iCol = 0; iCol<8; iCol++)
		{
			if (drgiLTbl[iRow][iCol] > 0)
				printf("%3.4f ", OptDzoneQtzr.UTbl[iRow][iCol]);
			else
				printf("%3.4f ", 0.0);
		}

		printf("\n");
	}
	printf("var table\n");
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			printf("%3.4f ", drgdVar[i][j]);
		
		}
		printf("\n");
	}

#if TEST_RD
	if (ADPT_DIST_PROFL)
	{
		printf("---------------Target Rate--------------\n");
		for (iRow = 0; iRow<8; iRow++)
		{
			for (iCol = 0; iCol<8; iCol++)
			{
				if (drgiLTbl[iRow][iCol] > 0)
					printf("%3.4f ", drgdTarRate[iRow][iCol]);
				else
					printf("%3.4f ", 0.0);
			}

			printf("\n");
		}

		printf("---------------Target Distortion--------------\n");
		for (iRow = 0; iRow<8; iRow++)
		{
			for (iCol = 0; iCol<8; iCol++)
			{
				if (iRow + iCol == 0)
					printf("%3.4f ", PSNR2MSE((double)TAR_PSNR_ARG / 100.0));
				else
				{
					if (drgiLTbl[iRow][iCol] > 0)
						printf("%3.4f ", drgdTarDist[iRow][iCol]);
					else
						printf("%3.4f ", PSNR2MSE((double)TAR_PSNR_ARG / 100.0));
				}
			}

			printf("\n");
		}
	}
#endif

	for (i = 0; i<NUM_QTZ_GAMMA; i++)
	{
		free(L_cache[i]);
	}
	free(L_cache);
}

/*
*  Function: fnDeQtz
*  De-quantize a 8x8 block of indices quantized by the
*  hirarchical uniform quantizer with deadzone (HUQD)
*/
double fnDeQtz(int iLv, int iBlkRow, int iBlkCol, int iFreqRow, int iFreqCol)
{
	int iLvSign, iLvAbs;

	if ((iFreqRow + iFreqCol) == 0)
	{
		if (iLv == 0)
			return 0;

		else
		{
			iLvSign = (iLv >= 0) ? 1 : -1;
			iLvAbs = abs(iLv);

			return iLvSign * ((double)iLvAbs*OptDzoneQtzr.QTbl[0][0]);
		}
	}

	else
	{
		if (vrgiExt[iBlkRow][iBlkCol][iFreqRow][iFreqCol][0] == 0)
		{
			if (iLv == 0)
				return 0;
			else
			{
				iLvSign = (iLv >= 0) ? 1 : -1;
				iLvAbs = abs(iLv);
				return iLvSign * (OptDzoneQtzr.DeltaTbl[iFreqRow][iFreqCol] + OptDzoneQtzr.QTbl[iFreqRow][iFreqCol] * (iLvAbs - 1));
			}
		}

		else
		{
			iLvSign = (iLv >= 0) ? 1 : -1;
			iLvAbs = abs(iLv);
			return iLvSign * (drgdYc[iFreqRow][iFreqCol] + ((double)iLvAbs - 0.5)*OptDzoneQtzr.QTbl[0][0]);
		}
	}
}


/*
*  Function: fnHDQ
*  Quantize a DCT Coefficient by the deadzone plus uniform quantization
*/
int fnHDQ(double dCoef, int iBlkRow, int iBlkCol, int iFreqRow, int iFreqCol)
{

	int iLv, iSign, iQtzIdx, iBitLocIdx;
	double dDctCoeffAbs;

	iSign = (dCoef >= 0) ? 1 : -1;
	dDctCoeffAbs = fabs(dCoef);
	  
	if ((iFreqRow + iFreqCol) == 0)
	{
		iLv = (int)floor(dDctCoeffAbs / OptDzoneQtzr.QTbl[iFreqRow][iFreqCol] + 0.5);
		iQtzIdx = iSign * iLv;
	}

	else
	{
		if (dDctCoeffAbs <= drgdYc[iFreqRow][iFreqCol])
		{

			drgdSumYmulY[iFreqRow][iFreqCol] += dDctCoeffAbs * dDctCoeffAbs;
			drgiNumBlkWithinYc[iFreqRow][iFreqCol]++;

			if (dDctCoeffAbs <= OptDzoneQtzr.UTbl[iFreqRow][iFreqCol])
			{
				iQtzIdx = 0;
			}
			else
			{
				iLv = (int)ceil((dDctCoeffAbs - OptDzoneQtzr.UTbl[iFreqRow][iFreqCol]) / OptDzoneQtzr.QTbl[iFreqRow][iFreqCol]);
				iQtzIdx = iSign * iLv;
			
				drgiCBP[iBlkRow][iBlkCol] = 1;
				//vrgiSig[iBlkRow][iBlkCol][iFreqRow][iFreqCol][0] = 1;
			}
		}

		else
		{

			iLv = (int)ceil((dDctCoeffAbs - drgdYc[iFreqRow][iFreqCol]) / OptDzoneQtzr.QTbl[0][0]);
			iQtzIdx = iSign * iLv;

			//drgiCBP[iBlkRow][iBlkCol] = 1; 
			drgiExtF[iBlkRow][iBlkCol] = 1;
			vrgiExt[iBlkRow][iBlkCol][iFreqRow][iFreqCol][0] = 1;
			//            for (iBitLocIdx=0;iBitLocIdx<BIT_LOC_SIG;iBitLocIdx++)
			//                vrgiSig[iBlkRow][iBlkCol][iFreqRow][iFreqCol][iBitLocIdx] = 1;
		}

	}

	return iQtzIdx;
}

/*
*  Function: fnDeQtz
*  De-quantize a 8x8 block of indices quantized by the
*  hirarchical uniform quantizer with deadzone (HUQD)
*/
double fnDeQtzExt(int iLv, int iFreqRow, int iFreqCol)
{
	int iLvSign, iLvAbs;

	iLvSign = (iLv >= 0) ? 1 : -1;
	iLvAbs = abs(iLv);
	return iLvSign * (drgdYc[iFreqRow][iFreqCol] + ((double)iLvAbs - 0.5)*OptDzoneQtzr.QTbl[0][0]);
}


/*
*  Function: fnDeQtz
*  De-quantize a 8x8 block of indices quantized by the
*  hirarchical uniform quantizer with deadzone (HUQD)
*/
double fnDeQtzSig(int iLv, int iFreqRow, int iFreqCol)
{
	int iLvSign, iLvAbs;

	if (iLv == 0)
		return 0;
	else
	{
		iLvSign = (iLv >= 0) ? 1 : -1;
		iLvAbs = abs(iLv);
		return iLvSign * (OptDzoneQtzr.DeltaTbl[iFreqRow][iFreqCol] + OptDzoneQtzr.QTbl[iFreqRow][iFreqCol] * (iLvAbs - 1));
	}
}

/*
*  Function: fnUpdateQtbl
*  Update the quantization parameters in the iterative algorithm
*  perform right after SDQ when the index sequence is fixed
*/
void fnUpdateQtbl()
{
	int i, j, k, m, iOptGammaIdx, iOptLambdaIdx, iOptLambdaIdxDeriv, iOptGammaIdxDeriv;
	double dCurDistDeriv, dCurDist, dMinDistDerivAbs, dMinDist;

	for (i = 0; i<8; i++)
	{
		for (j = 0; j<8; j++)
		{
			if (drgiLTblSDQ[i][j] > 0)
			{
				if (drgiLTblSDQ[i][j] != drgiLTbl[i][j])
					printf("L is changed during SDQ! L_HDQ(%d,%d)=%d, L_SDQ(%d,%d)=%d. \n", i, j, drgiLTbl[i][j], i, j, drgiLTblSDQ[i][j]);

				dMinDist = INFINITY;
				dMinDistDerivAbs = INFINITY;
				//for(k=1;k<256;k++)
				//for(k=0;k<NUM_QTZ_LAMBDA;k++)
				for (k = OptDzoneQtzr.Lambda_idx[i][j]; k <= OptDzoneQtzr.Lambda_idx[i][j]; k++)
				{
					//for (m=OptDzoneQtzr.Gamma_idx[i][j];m<=OptDzoneQtzr.Gamma_idx[i][j];m++)
					for (m = 0; m<NUM_QTZ_GAMMA; m++)
					{
						dCurDist = fnCalDist(i, j, Q_cache[m][k]);
						//dCurDistDeriv = fabs(fnCalDistDeriv(i,j,k/drgdLambda[i][j]));//Q_cache[m][k]));

						if (dCurDist < 0.0)
						{
							printf("error: distortion=%4.4f is nagative!\n", dCurDist);
							break;
						}

						//                        if (dMinDistDerivAbs > dCurDistDeriv)
						//                        {
						//                            dMinDistDerivAbs = dCurDistDeriv;
						//                            iOptGammaIdx = m;
						//                            iOptLambdaIdx = k;
						//                        }

						if (dMinDist > dCurDist)
						{
							dMinDist = dCurDist;
							iOptGammaIdx = m;
							iOptLambdaIdx = k;
						}
					}
				}

				//                if (iOptGammaIdxDeriv != iOptGammaIdx || iOptLambdaIdxDeriv != iOptLambdaIdx)
				//                    printf("error: when deciding the optimal Lambda index! iOptLambdaIdxDeriv=%d, iOptLambdaIdx=%d\n", iOptLambdaIdxDeriv, iOptLambdaIdx);

				//printf("FreqRow=%d, FreqCol=%d, PrvDist=%4.4f, CurDist=%4.4f\n", i,j,fnCalDist(i,j,Q_cache[OptDzoneQtzr.Gamma_idx[i][j]][OptDzoneQtzr.Lambda_idx[i][j]]), dMinDist);
				if (fnCalDist(i, j, Q_cache[OptDzoneQtzr.Gamma_idx[i][j]][OptDzoneQtzr.Lambda_idx[i][j]]) < dMinDist)
					printf("Error: when updating the Q table!\n");

				OptDzoneQtzr.QTbl[i][j] = Q_cache[iOptGammaIdx][iOptLambdaIdx];
				OptDzoneQtzr.UTbl[i][j] = drgdGamma[i][j] - drgiLTblSDQ[i][j] * OptDzoneQtzr.QTbl[i][j];
				OptDzoneQtzr.DeltaTbl[i][j] = OptDzoneQtzr.UTbl[i][j] + 1 + OptDzoneQtzr.QTbl[i][j] / (1 - exp(OptDzoneQtzr.QTbl[i][j]));

				/* scale the Q, U and Delta according to Lambda */
				OptDzoneQtzr.QTbl[i][j] *= drgdLambda[i][j];
				OptDzoneQtzr.UTbl[i][j] *= drgdLambda[i][j];
				OptDzoneQtzr.DeltaTbl[i][j] *= drgdLambda[i][j];

				dDistCurIterSDQ += (double)drgiNumBlkWithinYc[i][j] * dMinDist;
			}
		}
	}
}

double fnCalDistDeriv(int iFreqRow, int iFreqCol, double dCurQ)
{
	double dFnQ, dFnDerivQ, dTemp, dDistDeriv;

	dTemp = dCurQ / (1 - exp(dCurQ));
	dFnQ = drgdYc[iFreqRow][iFreqCol] + drgdLambda[iFreqRow][iFreqCol] * (1 + dTemp - dCurQ * (drgiLTblSDQ[iFreqRow][iFreqCol] + 1));
	dFnDerivQ = drgdLambda[iFreqRow][iFreqCol] * (1 / (1 - exp(dCurQ)) + dTemp * exp(dCurQ) - drgiLTblSDQ[iFreqRow][iFreqCol] - 1);

	dDistDeriv = dFnDerivQ * ((double)drgiNumBlkClgt0[iFreqRow][iFreqCol] - drgdSumY[iFreqRow][iFreqCol]) + (drgdSumCmulC[iFreqRow][iFreqCol] * drgdLambda[iFreqRow][iFreqCol] * dCurQ - drgdSumYmulC[iFreqRow][iFreqCol] + dFnQ * drgdSumC[iFreqRow][iFreqCol])*drgdLambda[iFreqRow][iFreqCol];

	dDistDeriv = 2.0*dDistDeriv / (double)drgiNumBlkClgt0[iFreqRow][iFreqCol];

	return dDistDeriv;
}

double fnCalDist(int iFreqRow, int iFreqCol, double dCurQ)
{
	double dFnQ, dTemp, dDist;

	dTemp = dCurQ / (1 - exp(dCurQ));
	dFnQ = drgdYc[iFreqRow][iFreqCol] + drgdLambda[iFreqRow][iFreqCol] * (1.0 + dTemp - dCurQ * ((double)drgiLTblSDQ[iFreqRow][iFreqCol] + 1.0));
	dDist = drgdSumYmulY[iFreqRow][iFreqCol] + dFnQ * dFnQ*(double)drgiNumBlkClgt0[iFreqRow][iFreqCol] - 2.0*dFnQ*drgdSumY[iFreqRow][iFreqCol] + (drgdSumCmulC[iFreqRow][iFreqCol] * drgdLambda[iFreqRow][iFreqCol] * dCurQ - 2.0*drgdSumYmulC[iFreqRow][iFreqCol] + 2.0*dFnQ*drgdSumC[iFreqRow][iFreqCol])*drgdLambda[iFreqRow][iFreqCol] * dCurQ;

	dDist = (double)drgiNumBlkWithinYc[iFreqRow][iFreqCol];

	return dDist;
}
