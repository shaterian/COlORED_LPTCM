#pragma warning(disable:4996)
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include<fstream>
#include <vector>
#include <regex>
#define HAVE_STRUCT_TIMESPEC
#include <pthread.h>
#include <opencv2/opencv.hpp>
using namespace std;
#include "global.h"
#include "HACACM.h"
#include "ImageDependent.h"
#include "AdapRecSpace.h"
#include "FastDCT.h"
#include "IntegerTCM.h"

int iFilterECEB = 0;

#if MULTIRESOLUTION
long int iOutlierRate;
long int iDcRate;
unsigned char OBF_IMG[ROWS / 8][COLS / 8] = { 0 };
unsigned char OTF_IMG[ROWS][COLS] = { 0 };
unsigned char OTFndDC_IMG[ROWS][COLS] = { 0 };
unsigned char NoOutlier_IMG[ROWS][COLS] = { 0 };
#endif

#if ENC_DEBUG
#include "debug.h"
char rgcDebugrgcIFName[] = "/Users/j4meng/jpeg_debug";
char rgcBSDebug[] = "/Users/j4meng/jpeg_bs_debug";
#endif
/*making size */

int ROWS;
int COLS;
char inpath[1024] = { 0 }; 
int NUMBLK;
// 3 channel 
unsigned char** Y_picture_buff;
unsigned char** Cb_picture_buff;
unsigned char** Cr_picture_buff;
std::vector<float> psnr_vec; 
std::vector<float> bpp_vec;
// ************ encode function *******************
void encode(unsigned char** drgucInput_Y , char  *input , int ch , int psnr );
/* Options */
int iOptMean = 0;   // 0: Substract EmpMean 1: Substract Mean of this image
int iOptAdpQ = 0;   // 0: Use Lambda Yc Aval Bval from 30 images 1: Use Lambda Yc Aval Bval from this image it is not used 
int iOptSortL = 0;  // 0: Use Zigzag Scan 1: Use adpative scan
int iOptAdpL = 1;   // 0: Use default L table; 1: Use L calculated from Yc and Lambda
int iOptLYC = 0;    // 0: Use Run Level Coding 1: Use Layer-based Coding
int iOptLapUnf = 1; // 0: Use Truncated Laplacian Only 1: Use Truncated Laplacian with Uniform Dist.
int iOptLamScl = 0;  // 0: Use tuned Scale 1: Use actual Lambda
int iRoseRD = 0;     // 1: Use recon space cal by Rose Alogirthm
int EXT_SCH_NO, iSDQOutPutPer = 0;
double dDeltaDC, dAvgCostRD, dDistCurIterSDQ;
int iBetaDC, iNumIter = 0;
int iImgIdxNum;
//double dTestDist=0.0;

double dScaleLambda = 3.0;
double dLScale = 1.0;
double dRate, dDist;
double drgdSumYmulC[8][8];
double drgdSumCmulC[8][8];
double drgdSumC[8][8];
double drgdSumY[8][8];
int drgiLTblSDQ[8][8];
int drgiNumBlkClgt0[8][8];

/* Global variables */
double dDctSclF[8][8];
double TAR_DIST, SDQ_LAMBDA, dCostHDQ;
int SDQ, ADPT_DIST_PROFL = 0;
int TAR_PSNR;
int TAR_PSNR_ARG;
std::vector<std::vector<int>> drgiQIdxHDQ;
char input[256] = { 0 };
char outpath[2048] = { 0 }; 
char textpath[2048] = {}; 
char jpg_path[2048] = {}; 
char recpath[2048] = { 0 };
char rgcIFName[256] = { 0 };
char rgcOFName[256] = { 0 };
char rgcRecOFName[256] = { 0 };
char rgctextOFName[256] = { 0 };
char rgcOFNameTQtbl[256] = { 0 };
char rgcIFNameLtbl[256] = { 0 };
char rgcIFNameQtbl[256] = { 0 };
char rgcOFNameDCIdxAjst[256] = { 0 };

//std::vector<std::vector<int>> drgiQIdxHDQAjstDC;//[ROWS / 8] [COLS / 8] ;
//int drgiQIdxHDQAjstDC[ROWS / 8] [COLS / 8] ;
std::vector<std::vector<int>> drgiQIdxHDQAjstDC; 
std::vector<std::vector<int>> drgiQIdxHDQAjstDCDec;
std::vector<std::vector<double>> drgdQIdxHDQDCPred;
//int drgiQIdxSDQ[ROWS][COLS] = { 0 };
std::vector<std::vector<int>> drgiQIdxSDQ;

double dAcumDiffAvg[NUM_CTX_PRED_DC] = { 0.0 };
int rgiAlpSizeDC[4] = { 26,50,82,256 };//{18,40,72,256};//86,

double drgdEmpMean[8][8] =
{
	{ -140.5908,    0.1368,   -0.3032,   -0.2300,   -0.3505,   -0.2354,   -0.1299,   -0.0804 },
{ 0.8442,   -0.0381,    0.0610,    0.0746,   -0.0336,    0.0263,    0.0173,    0.0715 },
{ -0.5942,   -0.1179,    0.1848,    0.0553,    0.0460,    0.0659,   -0.0047,   -0.0246 },
{ 0.1403,    0.0119,    0.0283,   -0.0063,    0.0220,    0.0292,   -0.0287,    0.0569 },
{ -0.3204,    0.0066,   -0.0448,    0.0224,    0.0384,   -0.0277,    0.0201,    0.0122 },
{ 0.1264,    0.0284,    0.0036,    0.0110,   -0.0215,   -0.0069,   -0.0103,    0.1085 },
{ -0.1338,    0.0260,   -0.0067,    0.0138,   -0.0004,   -0.0161,    0.0009,    0.0052 },
{ 0.0269,    0.0559,    0.0091,    0.0622,    0.0193,    0.0834,   -0.0146,    0.2621 }
};

/*
*  Variable: drgdEmpScl
*  Empirical scale for each coefficient in 8x8 block
*  calculated by 30 test images
*/

double drgdEmpScl[8][8] =
{
	{ 377.536,44.264,24.226,16.065,11.791,8.741,6.616,5.515 },
{ 44.805,23.840,16.633,12.131, 9.352,7.185,5.589,4.765 },
{ 24.803,16.728,12.973,10.165, 8.149,6.448,5.124,4.384 },
{ 16.471,12.275,10.134, 8.271, 6.805,5.536,4.595,4.000 },
{ 11.947, 9.375, 8.013, 6.721, 5.600,4.688,4.026,3.584 },
{ 8.936, 7.2583, 6.291, 5.364, 4.628,4.057,3.762,3.425 },
{ 6.978, 5.783, 5.156, 4.540, 4.096,3.828,3.717,3.353 },
{ 6.045, 5.108, 4.596, 4.156, 3.796,3.701,3.478,3.537 }
};


/*
*  variable: rgiRiceGolombOrder
*  order of Rice-Golomb Code used for different layers
*                      layer: 1,2,3,4,5,6+
*/
int rgiRiceGolombOrder[6] = { 6,5,4,3,2,1 };
int drgiBlk[8][8];

double dMaxDC = 0;
double dMinDC = 1024;

double dDCqstep;
#if ZEROFORCE
double drgdZeroPercent[8][8];
#endif

int ZigZag[8][8] = { 0, 1, 5, 6,14,15,27,28,
2, 4, 7,13,16,26,29,42,
3, 8,12,17,25,30,41,43,
9,11,18,24,31,40,44,53,
10,19,23,32,39,45,52,54,
20,22,33,38,46,51,55,60,
21,34,37,47,50,56,59,61,
35,36,48,49,57,58,62,63 };

int ZigZag_OneDir[8][8] =
{ 0, 1, 3, 6,14,15,27,28,
2, 4, 7,13,16,26,29,42,
5, 8,12,17,25,30,41,43,
9,11,18,24,31,40,44,53,
10,19,23,32,39,45,52,54,
20,22,33,38,46,51,55,60,
21,34,37,47,50,56,59,61,
35,36,48,49,57,58,62,63 };

static double alpha(int n)
{
	if (n == 0) return 1.0f / (double)sqrt(8);
	else        return 1.0f / 2.0f;
}

double MSE2PSNR(double MSE)
{
	double PSNR;
	PSNR = 48.1308 - 10.0*log10(MSE);
	return PSNR;
}

int cmp(const void *a, const void *b)
{
	return *(double *)a > *(double *)b ? 1 : -1;
}

void fnCalDCqstep(double** drgdCoeff, double** drgdRecCoeff)
{
	double drgdMSE[8][8];
	int i, j, u, v;
	double dminMSE, davgMSE;

	dminMSE = INFINITY;
	davgMSE = 0.0;

	for (u = 0; u<8; u++)
		for (v = 0; v<8; v++)
		{
			if (u + v)
			{
				drgdMSE[u][v] = 0.0;
				for (i = 0; i<ROWS / 8; i++)
					for (j = 0; j<COLS / 8; j++)
					{
						drgdMSE[u][v] += sqterm(drgdCoeff[i * 8 + u][j * 8 + v], drgdRecCoeff[i * 8 + u][j * 8 + v]);
					}

				dminMSE = MINIMUM((dminMSE), (drgdMSE[u][v]));
				davgMSE += drgdMSE[u][v];
			}
		}

	davgMSE /= 63.0;
	// 	dDCqstep = (sqrt(12*dminMSE/NUMBLK));
	dDCqstep = (sqrt(12 * davgMSE / NUMBLK));

	printf("DCqstep: %.2f\n", dDCqstep);
}

void _prepost(double *p, int w[8], int m)
{
	p[0] += p[w[7]], p[w[7]] = p[0] * 0.5F - p[w[7]];
	p[w[1]] += p[w[6]], p[w[6]] = p[w[1]] * 0.5F - p[w[6]];
	p[w[2]] += p[w[5]], p[w[5]] = p[w[2]] * 0.5F - p[w[5]];
	p[w[3]] += p[w[4]], p[w[4]] = p[w[3]] * 0.5F - p[w[4]];

	if (m)
	{
		p[w[4]] *= 1.40F, p[w[5]] *= 1.12F, p[w[6]] *= 1.14F, p[w[7]] *= 1.19F;
		p[w[7]] -= p[w[6]] * 0.125F;
		p[w[6]] += (p[w[7]] - p[w[5]])*0.375F;
		p[w[5]] += p[w[6]] * 0.5F - p[w[4]] * 0.375F;
		p[w[4]] += p[w[5]] * 0.75F;
	}
	else
	{
		p[w[4]] -= p[w[5]] * 0.75F;
		p[w[5]] -= p[w[6]] * 0.5F - p[w[4]] * 0.375F;
		p[w[6]] -= (p[w[7]] - p[w[5]])*0.375F;
		p[w[7]] += p[w[6]] * 0.125F;
		p[w[4]] /= 1.40F, p[w[5]] /= 1.12F, p[w[6]] /= 1.14F, p[w[7]] /= 1.19F;
	}

	p[w[7]] = p[0] * 0.5F - p[w[7]], p[0] -= p[w[7]];
	p[w[6]] = p[w[1]] * 0.5F - p[w[6]], p[w[1]] -= p[w[6]];
	p[w[5]] = p[w[2]] * 0.5F - p[w[5]], p[w[2]] -= p[w[5]];
	p[w[4]] = p[w[3]] * 0.5F - p[w[4]], p[w[3]] -= p[w[4]];
}

// calculate the distortion
float cal_dist(double** drgdCoeff, double** drgdRecCoeff, unsigned char** drgucInput_Y)
{
	//std::cout << "in cal dist " << ROWS << std::endl; 
	float image_psnr; 
	double dMSE = 0, tmpd;
	int i, j, u, v, x, y;
	//unsigned char drgucRecImage[ROWS][COLS];
	
	unsigned char** drgucRecImage;
	drgucRecImage = (unsigned char**)malloc(ROWS * sizeof(unsigned char*));
	for (unsigned int i = 0; i < ROWS; i++) {
		drgucRecImage[i] = (unsigned char*)malloc(sizeof(unsigned char) * COLS);
	}

	//double drgdRecImage[ROWS][COLS];
	double dIDCTblk[8][8];
	double drgdMSE[8][8];
	FILE *fpRF;

	double dPI = atan(1.0) * 4.0;
	/*

	FILE *fpMSE;*/

	for (u = 0; u<8; u++)
	{
		for (v = 0; v<8; v++)
		{
			drgdMSE[u][v] = 0.0;
		}
	}

	for (i = 0; i<ROWS / 8; i++)
		for (j = 0; j<COLS / 8; j++)
		{
			for (u = 0; u<8; u++)
				for (v = 0; v<8; v++)
				{
					dMSE += sqterm(drgdCoeff[i * 8 + u][j * 8 + v], drgdRecCoeff[i * 8 + u][j * 8 + v]);
					drgdMSE[u][v] += sqterm(drgdCoeff[i * 8 + u][j * 8 + v], drgdRecCoeff[i * 8 + u][j * 8 + v]);
#if TEST_RD
					if (drgdCoeff[i * 8 + u][j * 8 + v] <= drgdYc[u][v])
					{
						drgdRealDist[u][v] += sqterm(drgdCoeff[i * 8 + u][j * 8 + v], drgdRecCoeff[i * 8 + u][j * 8 + v]);
					}
#endif
				}
		}

	dMSE; // /= (COLS * ROWS);
	//printf("Average distortion in MSE in DCT domain is: %5.4f\n", dMSE);
	printf("Average distortion in PSNR in DCT domain is: %5.4f\n", MSE2PSNR(dMSE));

	//printf("------------------------------\n");
	printf("The target distortion (PSNR) is %2.2fdB\n\n", (double)TAR_PSNR_ARG / 100.0);
#if TEST_RD
	printf("---------------Real Distortion--------------\n");
	for (u = 0; u<8; u++)
	{
		for (v = 0; v<8; v++)
		{
			if (u + v == 0)
				printf("%3.4f ", PSNR2MSE((double)TAR_PSNR_ARG / 100.0));
			else
			{
				if (drgiLTbl[u][v] > 0)
					printf("%3.4f ", drgdRealDist[u][v] / (double)drgiNumBlkWithinYc[u][v]);
				else
					printf("%3.4f ", PSNR2MSE((double)TAR_PSNR_ARG / 100.0));
			}
		}
		printf("\n");
	}
#endif

	/*printf("------------------------------\n");

	printf("Distortion (PSNR) profile in DCT domain:\n");
	for (u = 0; u<8; u++)
	{
		for (v = 0; v<8; v++)
		{
			printf("%3.3f ", MSE2PSNR(drgdMSE[u][v] / NUMBLK));
		}
		printf("\n");
	}
	printf("------------------------------\n");

	*/// Distortion in pixel domain
	//if (iFilterECEB)
	//{
	//	//#if FILTER_ECEB
	//	int _l[2][8];
	//	//double dIDCT[ROWS][COLS];
	//	std::vector<std::vector<double>>dIDCT(ROWS, std::vector<double>(COLS)); 
	//	for (i = 0; i<ROWS / 8; i++)
	//		for (j = 0; j<COLS / 8; j++)
	//		{
	//			for (u = 0; u < 8; u++)
	//				for (v = 0; v < 8; v++)
	//					dIDCTblk[u][v] = drgdRecCoeff[i * 8 + u][j * 8 + v] / dDctSclF[u][v];

	//			IDCT(dIDCTblk);

	//			for (u = 0; u < 8; u++)
	//			{
	//				for (v = 0; v < 8; v++)
	//				{
	//					dIDCT[i * 8 + u][j * 8 + v] = dIDCTblk[u][v];
	//				}
	//			}
	//		}

	//	double begin = clock();

	//	for (i = 0; i<8; i++)
	//		_l[0][i] = i, _l[1][i] = COLS * i;

	//	for (int x = 4; x<COLS - 8; x += 8)
	//		for (y = 0; y<ROWS; y++)
	//			_prepost(dIDCT[y] + x, _l[0], 0);
	//	for (int y = 4; y<ROWS - 8; y += 8)
	//		for (x = 0; x<COLS; x++)
	//			_prepost(dIDCT[y] + x, _l[1], 0);
	//	double end = clock();
	//	double time1 = ((double)end - begin) / CLOCKS_PER_SEC;
	//	printf("CPU time for Post-Filtering algorithm: %2.6f seconds \n", time1);

	//	for (i = 0; i<ROWS / 8; i++)
	//		for (j = 0; j<COLS / 8; j++)
	//		{
	//			for (u = 0; u < 8; u++)
	//			{
	//				for (v = 0; v < 8; v++)
	//				{
	//					dIDCTblk[u][v] = dIDCT[i * 8 + u][j * 8 + v];
	//					if (dIDCTblk[u][v] >= 0)
	//						drgiBlk[u][v] = (int)(dIDCTblk[u][v] + .5) + 128;
	//					else
	//						drgiBlk[u][v] = (int)(dIDCTblk[u][v] - .5) + 128;
	//					if (drgiBlk[u][v] < 0)
	//						drgiBlk[u][v] = 0;
	//					if (drgiBlk[u][v] > 255)
	//						drgiBlk[u][v] = 255;
	//					drgucRecImage[i * 8 + u][j * 8 + v] = (unsigned char)drgiBlk[u][v];
	//				}
	//			}
	//		}
	//}

	//else
	//{
		for (i = 0; i<ROWS / 8; i++)
			for (j = 0; j<COLS / 8; j++)
			{
				for (u = 0; u < 8; u++)
					for (v = 0; v < 8; v++)
						dIDCTblk[u][v] = drgdRecCoeff[i * 8 + u][j * 8 + v] / dDctSclF[u][v];

				IDCT(dIDCTblk);

				for (u = 0; u < 8; u++)
				{
					for (v = 0; v < 8; v++)
					{
						if (dIDCTblk[u][v] >= 0)
							drgiBlk[u][v] = (int)(dIDCTblk[u][v] + .5) + 128;
						else
							drgiBlk[u][v] = (int)(dIDCTblk[u][v] - .5) + 128;
						if (drgiBlk[u][v] < 0)
							drgiBlk[u][v] = 0;
						if (drgiBlk[u][v] > 255)
							drgiBlk[u][v] = 255;
						drgucRecImage[i * 8 + u][j * 8 + v] = (unsigned char)drgiBlk[u][v];
					}
				}
			}
	//}

	dMSE = 0.0;
	for (i = 0; i<ROWS; i++)
		for (j = 0; j<COLS; j++)
			dMSE += sqterm(drgucInput_Y[i][j], drgucRecImage[i][j]);
	dMSE /= (COLS*ROWS);
	//printf("Average distortion in MSE in image domain is: %5.3f\n", mse_image);
	printf("Average distortion in PSNR in image domain is: %5.4f\n", MSE2PSNR(dMSE));
	image_psnr = MSE2PSNR(dMSE); 
	dDist = MSE2PSNR(dMSE);

	/* output the reconstructed image */
#if OUTPUT_REC_IMG
	if (SDQ == 1)
	{
		if (iSDQOutPutPer == 1)
		{
			sprintf(rgcRecOFName, "D:\\Research\\New_Codec_chung_code\\encoder_proj\\images\\rec.pgm", iImgIdxNum, dRate, dDist);

			fpRF = fopen(rgcRecOFName, "wb");
			//fprintf(fpRF, "P5\n");
			//fprintf(fpRF, "# CREATOR: XV Version 3.00  Rev: 3/30/93\n");
			//fprintf(fpRF, "%d %d\n", COLS, ROWS);
			//fprintf(fpRF, "255\n");
			fwrite(drgucRecImage, sizeof(unsigned char), COLS*ROWS, fpRF);
			fclose(fpRF);
		}
	}
	else
	{
		//std::cout << "here  is rec " << std::endl; 
		//sprintf(rgcRecOFName, "D:\\Research\\new_codec_test\\images\\generated_images\\Lena\\rec.yuv");//, iImgIdxNum, dRate, dDist);
		std::cout << rgcRecOFName << std::endl;
		fpRF = fopen(rgcRecOFName, "wb");
		std::cout << rgcRecOFName << std::endl;
		//fprintf(fpRF, "P5\n");
		//fprintf(fpRF, "# CREATOR: XV Version 3.00  Rev: 3/30/93\n");
		//fprintf(fpRF, "%d %d\n", COLS, ROWS);
		//fprintf(fpRF, "255\n");
		
		for (int r = 0; r < ROWS; r++) {
			fwrite(drgucRecImage[r], sizeof(unsigned char), COLS, fpRF);
		}
		fclose(fpRF);
	}


#endif

	double dTmpCoeff;
#if DIFFDIST
	fpRF = fopen("DiffD_Rec.pgm", "wb");
	//fprintf(fpRF, "P5\n");
	//fprintf(fpRF, "# CREATOR: XV Version 3.00  Rev: 3/30/93\n");
	//fprintf(fpRF, "%d %d\n", COLS, ROWS);
	//fprintf(fpRF, "255\n");
	fwrite(drgucRecImage, sizeof(char), COLS*ROWS, fpRF);
	fclose(fpRF);

	for (i = 0; i<ROWS / 8; i++)
		for (j = 0; j<COLS / 8; j++)
		{
			for (u = 0; u < 8; u++)
				for (v = 0; v < 8; v++)
				{
					dTmpCoeff = drgdRecCoeff[i * 8 + u][j * 8 + v];
					if (OTF_IMG[i * 8 + u][j * 8 + v] == 1 && u + v>0)
						dTmpCoeff = dTmpCoeff;
					else if (OTF_IMG[i * 8 + u][j * 8 + v] == 0 && u + v>0)
						dTmpCoeff = 0;

					dIDCTblk[u][v] = dTmpCoeff / dDctSclF[u][v];
				}

			IDCT(dIDCTblk);

			for (u = 0; u < 8; u++)
			{
				for (v = 0; v < 8; v++)
				{
					if (dIDCTblk[u][v] >= 0)
						drgiBlk[u][v] = (int)(dIDCTblk[u][v] + .5) + 128;
					else
						drgiBlk[u][v] = (int)(dIDCTblk[u][v] - .5) + 128;
					if (drgiBlk[u][v] < 0)
						drgiBlk[u][v] = 0;
					if (drgiBlk[u][v] > 255)
						drgiBlk[u][v] = 255;
					OTFndDC_IMG[i * 8 + u][j * 8 + v] = (unsigned char)drgiBlk[u][v];
				}
			}
		}

	fpRF = fopen("DiffD_Outlier_Rec.pgm", "wb");
	fprintf(fpRF, "P5\n");
	fprintf(fpRF, "# CREATOR: XV Version 3.00  Rev: 3/30/93\n");
	fprintf(fpRF, "%d %d\n", COLS, ROWS);
	fprintf(fpRF, "255\n");
	fwrite(OTFndDC_IMG, sizeof(char), COLS*ROWS, fpRF);
	fclose(fpRF);

	for (i = 0; i<ROWS / 8; i++)
		for (j = 0; j<COLS / 8; j++)
		{
			for (u = 0; u < 8; u++)
				for (v = 0; v < 8; v++)
				{
					dTmpCoeff = drgdRecCoeff[i * 8 + u][j * 8 + v];

					if (OTF_IMG[i * 8 + u][j * 8 + v] == 1 && u + v>0)
						dTmpCoeff = (dTmpCoeff<0 ? -1 : 1)*drgdYc[u][v];

					dIDCTblk[u][v] = dTmpCoeff / dDctSclF[u][v];
				}

			IDCT(dIDCTblk);

			for (u = 0; u < 8; u++)
			{
				for (v = 0; v < 8; v++)
				{
					if (dIDCTblk[u][v] >= 0)
						drgiBlk[u][v] = (int)(dIDCTblk[u][v] + .5) + 128;
					else
						drgiBlk[u][v] = (int)(dIDCTblk[u][v] - .5) + 128;
					if (drgiBlk[u][v] < 0)
						drgiBlk[u][v] = 0;
					if (drgiBlk[u][v] > 255)
						drgiBlk[u][v] = 255;
					NoOutlier_IMG[i * 8 + u][j * 8 + v] = (unsigned char)drgiBlk[u][v];
				}
			}
		}

	fpRF = fopen("DiffD_Inlier_Rec.pgm", "wb");
	fprintf(fpRF, "P5\n");
	fprintf(fpRF, "# CREATOR: XV Version 3.00  Rev: 3/30/93\n");
	fprintf(fpRF, "%d %d\n", COLS, ROWS);
	fprintf(fpRF, "255\n");
	fwrite(NoOutlier_IMG, sizeof(char), COLS*ROWS, fpRF);
	fclose(fpRF);
#else
#if MULTIRESOLUTION
	printf("Rate for OBF is: 0.007 (bpp)\n");
	printf("Rate for Outlier flag is: 0.081 (bpp)\n");
	printf("Rate for Outlier is: %1.4f (bpp)\n", (double)iOutlierRate / ROWS / COLS);
	printf("Rate for DC is: %1.4f (bpp)\n", (double)(iDcRate) / ROWS / COLS);
	printf("Rate for Outlier+DC is: %1.4f (bpp)\n", (double)(iDcRate + iOutlierRate) / ROWS / COLS);
	printf("Rate for No Outlier is: %1.4f (bpp)\n", dRate - (double)iOutlierRate / ROWS / COLS);

	fpRF = fopen("OBF_Rec.pgm", "wb");
	fprintf(fpRF, "P5\n");
	fprintf(fpRF, "# CREATOR: XV Version 3.00  Rev: 3/30/93\n");
	fprintf(fpRF, "%d %d\n", COLS / 8, ROWS / 8);
	fprintf(fpRF, "255\n");
	fwrite(OBF_IMG, sizeof(char), (COLS / 8)*(ROWS / 8), fpRF);
	fclose(fpRF);

	fpRF = fopen("Rec.pgm", "wb");
	fprintf(fpRF, "P5\n");
	fprintf(fpRF, "# CREATOR: XV Version 3.00  Rev: 3/30/93\n");
	fprintf(fpRF, "%d %d\n", COLS, ROWS);
	fprintf(fpRF, "255\n");
	fwrite(drgucRecImage, sizeof(char), COLS*ROWS, fpRF);
	fclose(fpRF);

#if MULTIRESOLUTION_DC
	for (i = 0; i<ROWS / 8; i++)
		for (j = 0; j<COLS / 8; j++)
		{
			for (u = 0; u < 8; u++)
				for (v = 0; v < 8; v++)
				{
					dTmpCoeff = drgdRecCoeff[i * 8 + u][j * 8 + v];
					if (u + v>0)
						dTmpCoeff = 0;

					dIDCTblk[u][v] = dTmpCoeff / dDctSclF[u][v];
				}

			IDCT(dIDCTblk);

			for (u = 0; u < 8; u++)
			{
				for (v = 0; v < 8; v++)
				{
					if (dIDCTblk[u][v] >= 0)
						drgiBlk[u][v] = (int)(dIDCTblk[u][v] + .5) + 128;
					else
						drgiBlk[u][v] = (int)(dIDCTblk[u][v] - .5) + 128;
					if (drgiBlk[u][v] < 0)
						drgiBlk[u][v] = 0;
					if (drgiBlk[u][v] > 255)
						drgiBlk[u][v] = 255;
					OTFndDC_IMG[i * 8 + u][j * 8 + v] = (unsigned char)drgiBlk[u][v];
				}
			}
		}

	fpRF = fopen("DC_Rec.pgm", "wb");
	fprintf(fpRF, "P5\n");
	fprintf(fpRF, "# CREATOR: XV Version 3.00  Rev: 3/30/93\n");
	fprintf(fpRF, "%d %d\n", COLS, ROWS);
	fprintf(fpRF, "255\n");
	fwrite(OTFndDC_IMG, sizeof(char), COLS*ROWS, fpRF);
	fclose(fpRF);

	for (i = 0; i<ROWS / 8; i++)
		for (j = 0; j<COLS / 8; j++)
		{
			for (u = 0; u < 8; u++)
				for (v = 0; v < 8; v++)
				{
					dTmpCoeff = drgdRecCoeff[i * 8 + u][j * 8 + v];
					if (OTF_IMG[i * 8 + u][j * 8 + v] == 1 && u + v>0)
						dTmpCoeff = (dTmpCoeff<0 ? -1 : 1)*drgdYc[u][v];
					else if (OTF_IMG[i * 8 + u][j * 8 + v] == 0 && u + v>0)
						dTmpCoeff = 0;

					dIDCTblk[u][v] = dTmpCoeff / dDctSclF[u][v];
				}

			IDCT(dIDCTblk);

			for (u = 0; u < 8; u++)
			{
				for (v = 0; v < 8; v++)
				{
					if (dIDCTblk[u][v] >= 0)
						drgiBlk[u][v] = (int)(dIDCTblk[u][v] + .5) + 128;
					else
						drgiBlk[u][v] = (int)(dIDCTblk[u][v] - .5) + 128;
					if (drgiBlk[u][v] < 0)
						drgiBlk[u][v] = 0;
					if (drgiBlk[u][v] > 255)
						drgiBlk[u][v] = 255;
					OTFndDC_IMG[i * 8 + u][j * 8 + v] = (unsigned char)drgiBlk[u][v];
				}
			}
		}

	fpRF = fopen("OutlierFlag_DC_Rec.pgm", "wb");
	fprintf(fpRF, "P5\n");
	fprintf(fpRF, "# CREATOR: XV Version 3.00  Rev: 3/30/93\n");
	fprintf(fpRF, "%d %d\n", COLS, ROWS);
	fprintf(fpRF, "255\n");
	fwrite(OTFndDC_IMG, sizeof(char), COLS*ROWS, fpRF);
	fclose(fpRF);

	for (i = 0; i<ROWS / 8; i++)
		for (j = 0; j<COLS / 8; j++)
		{
			for (u = 0; u < 8; u++)
				for (v = 0; v < 8; v++)
				{
					dTmpCoeff = drgdRecCoeff[i * 8 + u][j * 8 + v];
					if (OTF_IMG[i * 8 + u][j * 8 + v] == 1 && u + v>0)
						dTmpCoeff = dTmpCoeff;
					else if (OTF_IMG[i * 8 + u][j * 8 + v] == 0 && u + v>0)
						dTmpCoeff = 0;

					dIDCTblk[u][v] = dTmpCoeff / dDctSclF[u][v];
				}

			IDCT(dIDCTblk);

			for (u = 0; u < 8; u++)
			{
				for (v = 0; v < 8; v++)
				{
					if (dIDCTblk[u][v] >= 0)
						drgiBlk[u][v] = (int)(dIDCTblk[u][v] + .5) + 128;
					else
						drgiBlk[u][v] = (int)(dIDCTblk[u][v] - .5) + 128;
					if (drgiBlk[u][v] < 0)
						drgiBlk[u][v] = 0;
					if (drgiBlk[u][v] > 255)
						drgiBlk[u][v] = 255;
					OTFndDC_IMG[i * 8 + u][j * 8 + v] = (unsigned char)drgiBlk[u][v];
				}
			}
		}

	fpRF = fopen("Outlier_DC_Rec.pgm", "wb");
	fprintf(fpRF, "P5\n");
	fprintf(fpRF, "# CREATOR: XV Version 3.00  Rev: 3/30/93\n");
	fprintf(fpRF, "%d %d\n", COLS, ROWS);
	fprintf(fpRF, "255\n");
	fwrite(OTFndDC_IMG, sizeof(char), COLS*ROWS, fpRF);
	fclose(fpRF);

	for (i = 0; i<ROWS / 8; i++)
		for (j = 0; j<COLS / 8; j++)
		{
			for (u = 0; u < 8; u++)
				for (v = 0; v < 8; v++)
				{
					dTmpCoeff = drgdRecCoeff[i * 8 + u][j * 8 + v];

					if (OTF_IMG[i * 8 + u][j * 8 + v] == 1 && u + v>0)
						dTmpCoeff = (dTmpCoeff<0 ? -1 : 1)*drgdYc[u][v];

					dIDCTblk[u][v] = dTmpCoeff / dDctSclF[u][v];
				}

			IDCT(dIDCTblk);

			for (u = 0; u < 8; u++)
			{
				for (v = 0; v < 8; v++)
				{
					if (dIDCTblk[u][v] >= 0)
						drgiBlk[u][v] = (int)(dIDCTblk[u][v] + .5) + 128;
					else
						drgiBlk[u][v] = (int)(dIDCTblk[u][v] - .5) + 128;
					if (drgiBlk[u][v] < 0)
						drgiBlk[u][v] = 0;
					if (drgiBlk[u][v] > 255)
						drgiBlk[u][v] = 255;
					NoOutlier_IMG[i * 8 + u][j * 8 + v] = (unsigned char)drgiBlk[u][v];
				}
			}
		}

	fpRF = fopen("Inlier_DC_Rec.pgm", "wb");
	fprintf(fpRF, "P5\n");
	fprintf(fpRF, "# CREATOR: XV Version 3.00  Rev: 3/30/93\n");
	fprintf(fpRF, "%d %d\n", COLS, ROWS);
	fprintf(fpRF, "255\n");
	fwrite(NoOutlier_IMG, sizeof(char), COLS*ROWS, fpRF);
	fclose(fpRF);
#else
	for (i = 0; i<ROWS / 8; i++)
		for (j = 0; j<COLS / 8; j++)
		{
			for (u = 0; u < 8; u++)
				for (v = 0; v < 8; v++)
				{
					dTmpCoeff = drgdRecCoeff[i * 8 + u][j * 8 + v];
					if (OTF_IMG[i * 8 + u][j * 8 + v] == 1 && u + v>0)
						dTmpCoeff = drgdYc[u][v];
					else if (OTF_IMG[i * 8 + u][j * 8 + v] == 0 && u + v>0)
						dTmpCoeff = 0;
					else if (u + v == 0)
						dTmpCoeff = 0;

					dIDCTblk[u][v] = dTmpCoeff / dDctSclF[u][v];
				}

			IDCT(dIDCTblk);

			for (u = 0; u < 8; u++)
			{
				for (v = 0; v < 8; v++)
				{
					if (dIDCTblk[u][v] >= 0)
						drgiBlk[u][v] = (int)(dIDCTblk[u][v] + .5) + 128;
					else
						drgiBlk[u][v] = (int)(dIDCTblk[u][v] - .5) + 128;
					if (drgiBlk[u][v] < 0)
						drgiBlk[u][v] = 0;
					if (drgiBlk[u][v] > 255)
						drgiBlk[u][v] = 255;
					OTFndDC_IMG[i * 8 + u][j * 8 + v] = (unsigned char)drgiBlk[u][v];
				}
			}
		}

	fpRF = fopen("OutlierFlag_Rec.pgm", "wb");
	fprintf(fpRF, "P5\n");
	fprintf(fpRF, "# CREATOR: XV Version 3.00  Rev: 3/30/93\n");
	fprintf(fpRF, "%d %d\n", COLS, ROWS);
	fprintf(fpRF, "255\n");
	fwrite(OTFndDC_IMG, sizeof(char), COLS*ROWS, fpRF);
	fclose(fpRF);

	for (i = 0; i<ROWS / 8; i++)
		for (j = 0; j<COLS / 8; j++)
		{
			for (u = 0; u < 8; u++)
				for (v = 0; v < 8; v++)
				{
					dTmpCoeff = drgdRecCoeff[i * 8 + u][j * 8 + v];
					if (OTF_IMG[i * 8 + u][j * 8 + v] == 1 && u + v>0)
						dTmpCoeff = dTmpCoeff;
					else if (OTF_IMG[i * 8 + u][j * 8 + v] == 0 && u + v>0)
						dTmpCoeff = 0;
					else if (u + v == 0)
						dTmpCoeff = 0;

					dIDCTblk[u][v] = dTmpCoeff / dDctSclF[u][v];
				}

			IDCT(dIDCTblk);

			for (u = 0; u < 8; u++)
			{
				for (v = 0; v < 8; v++)
				{
					if (dIDCTblk[u][v] >= 0)
						drgiBlk[u][v] = (int)(dIDCTblk[u][v] + .5) + 128;
					else
						drgiBlk[u][v] = (int)(dIDCTblk[u][v] - .5) + 128;
					if (drgiBlk[u][v] < 0)
						drgiBlk[u][v] = 0;
					if (drgiBlk[u][v] > 255)
						drgiBlk[u][v] = 255;
					OTFndDC_IMG[i * 8 + u][j * 8 + v] = (unsigned char)drgiBlk[u][v];
				}
			}
		}

	fpRF = fopen("Outlier_Rec.pgm", "wb");
	fprintf(fpRF, "P5\n");
	fprintf(fpRF, "# CREATOR: XV Version 3.00  Rev: 3/30/93\n");
	fprintf(fpRF, "%d %d\n", COLS, ROWS);
	fprintf(fpRF, "255\n");
	fwrite(OTFndDC_IMG, sizeof(char), COLS*ROWS, fpRF);
	fclose(fpRF);

	for (i = 0; i<ROWS / 8; i++)
		for (j = 0; j<COLS / 8; j++)
		{
			for (u = 0; u < 8; u++)
				for (v = 0; v < 8; v++)
				{
					dTmpCoeff = drgdRecCoeff[i * 8 + u][j * 8 + v];

					if (OTF_IMG[i * 8 + u][j * 8 + v] == 1 && u + v>0)
						dTmpCoeff = (dTmpCoeff<0 ? -1 : 1)*drgdYc[u][v];
					else if (u + v == 0)
						dTmpCoeff = 0;

					dIDCTblk[u][v] = dTmpCoeff / dDctSclF[u][v];
				}

			IDCT(dIDCTblk);

			for (u = 0; u < 8; u++)
			{
				for (v = 0; v < 8; v++)
				{
					if (dIDCTblk[u][v] >= 0)
						drgiBlk[u][v] = (int)(dIDCTblk[u][v] + .5) + 128;
					else
						drgiBlk[u][v] = (int)(dIDCTblk[u][v] - .5) + 128;
					if (drgiBlk[u][v] < 0)
						drgiBlk[u][v] = 0;
					if (drgiBlk[u][v] > 255)
						drgiBlk[u][v] = 255;
					NoOutlier_IMG[i * 8 + u][j * 8 + v] = (unsigned char)drgiBlk[u][v];
				}
			}
		}

	fpRF = fopen("Inlier_Rec.pgm", "wb");
	fprintf(fpRF, "P5\n");
	fprintf(fpRF, "# CREATOR: XV Version 3.00  Rev: 3/30/93\n");
	fprintf(fpRF, "%d %d\n", COLS, ROWS);
	fprintf(fpRF, "255\n");
	fwrite(NoOutlier_IMG, sizeof(char), COLS*ROWS, fpRF);
	fclose(fpRF);
#endif
#endif
#endif

	
	for (i = 0; i < ROWS; i++)
	{
		free(drgucRecImage[i]);
	}
	free(drgucRecImage);

	return image_psnr; 
}

float write_buff(double** drgdCoeff, double** drgdRecCoeff, unsigned char** drgucInput_ch , int component )
{
	//std::cout << "in cal dist " << ROWS << std::endl; 
	float image_psnr;
	double dMSE = 0, tmpd;
	int i, j, u, v, x, y;
	
	unsigned char** drgucRecImage;
	drgucRecImage = (unsigned char**)malloc(ROWS * sizeof(unsigned char*));
	for (unsigned int i = 0; i < ROWS; i++) {
		drgucRecImage[i] = (unsigned char*)malloc(sizeof(unsigned char) * COLS);
	}

	//double drgdRecImage[ROWS][COLS];
	double dIDCTblk[8][8];
	double drgdMSE[8][8];
	FILE* fpRF;

	double dPI = atan(1.0) * 4.0;
	/*

	FILE *fpMSE;*/

	for (u = 0; u < 8; u++)
	{
		for (v = 0; v < 8; v++)
		{
			drgdMSE[u][v] = 0.0;
		}
	}

	for (i = 0; i < ROWS / 8; i++)
		for (j = 0; j < COLS / 8; j++)
		{
			for (u = 0; u < 8; u++)
				for (v = 0; v < 8; v++)
				{
					dMSE += sqterm(drgdCoeff[i * 8 + u][j * 8 + v], drgdRecCoeff[i * 8 + u][j * 8 + v]);
					drgdMSE[u][v] += sqterm(drgdCoeff[i * 8 + u][j * 8 + v], drgdRecCoeff[i * 8 + u][j * 8 + v]);
				}
		}

	dMSE /= (COLS * ROWS);
	std::cout << "dMSE " << dMSE << std::endl;
	//printf("Average distortion in MSE in DCT domain is: %5.4f\n", dMSE);
	printf("Average distortion in PSNR in DCT domain is: %5.4f\n", MSE2PSNR(dMSE));

	//printf("------------------------------\n");
	printf("The target distortion (PSNR) is %2.2fdB\n\n", (double)TAR_PSNR_ARG / 100.0);
	for (i = 0; i < ROWS / 8; i++)
		for (j = 0; j < COLS / 8; j++)
		{
			for (u = 0; u < 8; u++)
				for (v = 0; v < 8; v++)
					dIDCTblk[u][v] = drgdRecCoeff[i * 8 + u][j * 8 + v] / dDctSclF[u][v];

			IDCT(dIDCTblk);

			for (u = 0; u < 8; u++)
			{
				for (v = 0; v < 8; v++)
				{
					if (dIDCTblk[u][v] >= 0)
						drgiBlk[u][v] = (int)(dIDCTblk[u][v] + .5) + 128;
					else
						drgiBlk[u][v] = (int)(dIDCTblk[u][v] - .5) + 128;
					if (drgiBlk[u][v] < 0)
						drgiBlk[u][v] = 0;
					if (drgiBlk[u][v] > 255)
						drgiBlk[u][v] = 255;
					drgucRecImage[i * 8 + u][j * 8 + v] = (unsigned char)drgiBlk[u][v];
				}
			}
		}
	dMSE = 0.0;
	for (i = 0; i < ROWS; i++)
		for (j = 0; j < COLS; j++)
			dMSE += sqterm(drgucInput_ch[i][j], drgucRecImage[i][j]);
	dMSE /= (COLS * ROWS);
	printf("Average distortion in PSNR in image domain is: %5.4f\n", MSE2PSNR(dMSE));
	image_psnr = MSE2PSNR(dMSE);
	dDist = MSE2PSNR(dMSE);

	switch (component)
	{
	case 1 :
		for (int r= 0 ; r < ROWS; r++)
			for (int c = 0 ; c < COLS; c++)
				Y_picture_buff[r][c] = drgucRecImage[r][c]; 
		break;

	case 2 :
		for (int r= 0 ; r < ROWS; r++)
			for (int c = 0 ; c < COLS; c++)
				Cb_picture_buff[r][c] = drgucRecImage[r][c];
		break;

	case 3 :

		for (int r = 0; r < ROWS; r++)
			for (int c = 0 ; c < COLS; c++)
				Cr_picture_buff[r][c] = drgucRecImage[r][c];
		break;
	}
	std::cout << "rec is done " << std::endl;
	return image_psnr;
}


double PSNR2MSE(double psnr)
{
	return (255 * 255 / pow(10.0, psnr / 10.0));
}

double fnPredCoding(int iRow, int iCol, double dScl)
{
	int iW, iNW, iN, iNN, iNE, iNNE, iWW;
	double dPred, dDv, dDh, dD; 
	 
	if (iRow == 0 && iCol == 0)
	{
		 
		dDeltaDC = 255.0;
		dPred = 0.0;
	}

	else if (iRow == 0 && iCol != 0)
	{
		dDeltaDC = 255.0;
		dPred = drgiQIdxHDQAjstDC[iRow][iCol - 1];
	}

	else if (iRow != 0 && iCol == 0)
	{
		dDeltaDC = 255.0;
		dPred = drgiQIdxHDQAjstDC[iRow - 1][iCol];
	}

	else if (iRow == 1 || iCol == 7  )
	{
		
		if (drgiQIdxHDQAjstDC[iRow - 1][iCol - 1] >= MAXIMUM(drgiQIdxHDQAjstDC[iRow - 1][iCol], drgiQIdxHDQAjstDC[iRow][iCol - 1]))
		{
			dDeltaDC = 255.0;
			dPred = MINIMUM(drgiQIdxHDQAjstDC[iRow - 1][iCol], drgiQIdxHDQAjstDC[iRow][iCol - 1]);
		}

		else if (drgiQIdxHDQAjstDC[iRow - 1][iCol - 1] <= MINIMUM(drgiQIdxHDQAjstDC[iRow - 1][iCol], drgiQIdxHDQAjstDC[iRow][iCol - 1]))
		{
			dDeltaDC = 255.0;
			dPred = MAXIMUM(drgiQIdxHDQAjstDC[iRow - 1][iCol], drgiQIdxHDQAjstDC[iRow][iCol - 1]);
		}

		else
		{
			dDeltaDC = 0.0;
			dPred = drgiQIdxHDQAjstDC[iRow - 1][iCol] - drgiQIdxHDQAjstDC[iRow - 1][iCol - 1] + drgiQIdxHDQAjstDC[iRow][iCol - 1];
		}
	}

	else
	{
		iN = drgiQIdxHDQAjstDC[iRow - 1][iCol];
		iW = drgiQIdxHDQAjstDC[iRow][iCol - 1];
		
		iNW = drgiQIdxHDQAjstDC[iRow - 1][iCol - 1];
		
		iNN = drgiQIdxHDQAjstDC[iRow - 2][iCol];
		
		if (iCol == (COLS / 8 - 1) )
			iNE = drgiQIdxHDQAjstDC[iRow ][0];
		else
		{
			iNE = drgiQIdxHDQAjstDC[iRow - 1][iCol + 1];
		}


		if (iCol == (COLS / 8 - 1))
		{
			
			iNNE = drgiQIdxHDQAjstDC[iRow - 1][0];
		}
		else
			iNNE = drgiQIdxHDQAjstDC[iRow - 2][iCol + 1];
		
		if (iCol == 1 )
			iWW = drgiQIdxHDQAjstDC[iRow - 1][COLS / 8 - 1];
		else
			iWW = drgiQIdxHDQAjstDC[iRow][iCol - 2];

		dDv = (double)(abs(iW - iNW) + abs(iN - iNN) + abs(iNE - iNNE));
		dDh = (double)(abs(iW - iWW) + abs(iN - iNW) + abs(iN - iNE));
		dD = dScl * (dDv - dDh);
		dDeltaDC = dScl * (dDv + dDh) + 2 * fabs(drgdQIdxHDQDCPred[iRow][iCol - 1]);//;+drgdQIdxHDQDCPred[iRow-1][iCol]);

		if (dD>80.0)
			dPred = (double)iW;
		else if (dD<-80.0)
			dPred = (double)iN;
		else
		{
			dPred = (double)(iW + iN) / 2.0 + (double)(iNE - iNW) / 4.0;
			if (dD > 32.0)
				dPred = (dPred + (double)iW) / 2.0;
			else if (dD > 8.0)
				dPred = (3.0*dPred + (double)iW) / 4.0;
			else if (dD < -32.0)
				dPred = (dPred + (double)iN) / 2.0;
			else if (dD < -8.0)
				dPred = (3.0*dPred + (double)iN) / 4.0;
		}
	}

	//iBetaDC = iN<dPred?1:0 + 2*(iW<dPred?1:0);// + 4*(iNW<dPred?1:0);

	return dPred;
}

int fnRound(double dRoundIn)
{
	int iSign;

	iSign = dRoundIn<0 ? -1 : 1;

	if (iSign == 1)
		return ((int)(.5 + dRoundIn));
	else
		return (-(int)(.5 - dRoundIn));
}
void write_picture(unsigned char ** Y_pic, unsigned char** Cb_pic, unsigned char** Cr_pic) {
	//strcat(outpath, "\\%d\\%s%d.yuv");
	strcat(recpath, "\\%d\\%s%d.yuv");
	strcat(textpath, "\\%d\\%s%d.txt");
	std::cout << TAR_PSNR_ARG << std::endl;

	char* outrec = new char[strlen(input) + 8];
	memcpy(outrec, input, strlen(input) - 4);
	outrec[strlen(input) - 4] = '\0';
	strcat(outrec, "_");
	sprintf(rgcRecOFName, recpath, TAR_PSNR_ARG / 100, outrec, TAR_PSNR_ARG / 100);
	sprintf(rgctextOFName, textpath, TAR_PSNR_ARG / 100, outrec, TAR_PSNR_ARG / 100);

	std::cout << "textfile is " << rgctextOFName << std::endl;

	std::cout << rgcRecOFName << std::endl;
	FILE* fpRF; 
	fpRF = fopen(rgcRecOFName, "wb");
	std::cout << rgcRecOFName << std::endl;
	for (int r = 0; r < ROWS; r++) {
		fwrite(Y_pic[r], sizeof(unsigned char), COLS, fpRF);
	}
	for (int r = 0; r < ROWS; r++) {
		fwrite(Cb_pic[r], sizeof(unsigned char), COLS, fpRF);
	}
	for (int r = 0; r < ROWS; r++) {
		fwrite(Cr_pic[r], sizeof(unsigned char), COLS, fpRF);
	}
	fclose(fpRF);

	
	float psnr_total;
	float mse_total;
	psnr_total = (6*psnr_vec[0] + psnr_vec[1] + psnr_vec[2])/8;
	std::cout << psnr_total << std::endl; 
	//psnr_total = MSE2PSNR(mse_total); 
	
	float total_bits; 
	total_bits = bpp_vec[0] + bpp_vec[1] + bpp_vec[2];
	std::cout << "total_bits  " << total_bits << std::endl;

	FILE* fptxt;
	fptxt = fopen(rgctextOFName, "wb");
	fprintf(fptxt, "PSNR\t%f\n", psnr_total);
	fprintf(fptxt, "BPP\t%f\n", total_bits);
	
	fclose(fptxt);
}



/* main function */
int main(int argc, char *argv[])
{
	sprintf(inpath, argv[1]);
	sprintf(input, argv[2]);
	//sprintf(outpath, argv[3]);
	sprintf(recpath, argv[3]);
	sprintf(textpath, argv[4]);
	sprintf(jpg_path, argv[5]); 
	TAR_PSNR_ARG = atoi(argv[6]);
	TAR_DIST = PSNR2MSE((double)TAR_PSNR_ARG / 100.0);
	TAR_PSNR = (int)((double)TAR_PSNR_ARG / 100.0 + .5);

	//char regx_str[] = "[0-9]+(x)[0-9]+";
	//regex rg(regx_str);
	//std::cmatch match;
	//regex_search(input, match, rg);
	//std::string subj = match.str();
	//std::regex n("[0-9]+");
	//std::sregex_iterator next(subj.begin(), subj.end(), n);
	//std::sregex_iterator end_str;
	//std::vector<double> size;
	//while (next != end_str) {

	//	std::smatch m = *next;
	//	size.push_back(stod(m.str()));
	//	next++;
	//}
	//double row = size[1];
	//double col = size[0];

	cv::Mat jpg_image = cv::imread(jpg_path);
	double row = jpg_image.rows;
	double col = jpg_image.cols;

	std::cout << "row" << row << "col" << col << std::endl;
	ROWS = 8 * ceil(row / 8);
	COLS = 8 * ceil(col / 8);

	NUMBLK = ROWS * COLS / 64;

	Y_picture_buff = (unsigned char**)malloc(ROWS * sizeof(unsigned char*));
	for (unsigned int i = 0; i < ROWS; i++) {
		Y_picture_buff[i] = (unsigned char*)malloc(sizeof(unsigned char) * COLS);
	}

	Cb_picture_buff = (unsigned char**)malloc(ROWS * sizeof(unsigned char*));
	for (unsigned int i = 0; i < ROWS; i++) {
		Cb_picture_buff[i] = (unsigned char*)malloc(sizeof(unsigned char) * COLS);
	}

	Cr_picture_buff = (unsigned char**)malloc(ROWS * sizeof(unsigned char*));
	for (unsigned int i = 0; i < ROWS; i++) {
		Cr_picture_buff[i] = (unsigned char*)malloc(sizeof(unsigned char) * COLS);
	}

	unsigned char** drgucInput_Y;
	drgucInput_Y = (unsigned char**)malloc(ROWS * sizeof(unsigned char*));
	for (unsigned int i = 0; i < ROWS; i++) {
		drgucInput_Y[i] = (unsigned char*)malloc(sizeof(unsigned char) * COLS);
	}

	unsigned char** image;
	image = (unsigned char**)malloc( 3 * (int)row * sizeof(unsigned char*));
	for (unsigned int i = 0; i < 3 * (int)row; i++) {
		image[i] = (unsigned char*)malloc(sizeof(unsigned char) *  (int)col);
	}

	unsigned char** drgucInput_U;
	drgucInput_U = (unsigned char**)malloc(ROWS * sizeof(unsigned char*));
	for (unsigned int i = 0; i < ROWS; i++) {
		drgucInput_U[i] = (unsigned char*)malloc(sizeof(unsigned char) * COLS);
	}


	unsigned char** drgucInput_V;
	drgucInput_V = (unsigned char**)malloc(ROWS * sizeof(unsigned char*));
	for (unsigned int i = 0; i < ROWS; i++) {
		drgucInput_V[i] = (unsigned char*)malloc(sizeof(unsigned char) * COLS);
	}


	/********************************************************************/

	/************Globals**********/

	drgiQIdxHDQ.resize(ROWS, std::vector<int>(COLS));
	
	drgiQIdxSDQ.resize(ROWS, std::vector<int>(COLS));
	drgiQIdxHDQAjstDCDec.resize(ROWS / 8, std::vector<int>(COLS / 8, 0));
	drgdQIdxHDQDCPred.resize(ROWS / 8, std::vector<double>(COLS / 8, 0));
	drgiQIdxHDQAjstDC.resize(ROWS / 8, std::vector<int>(COLS / 8, 0));
	/*********class variables******************/
	drgiExtF.resize(ROWS / 8, std::vector<int>(COLS / 8, 0));
	drgiExtFDec.resize(ROWS / 8, std::vector<int>(COLS / 8, 0));
	drgiCBP.resize(ROWS / 8, std::vector<int>(COLS / 8, 0));
	vrgiSign.resize(ROWS / 8);
	for (int i = 0; i < ROWS / 8; i++)
		vrgiSign[i].resize(COLS / 8, std::vector< std::vector<int>>(8, std::vector<int>(8)));

	FILE* fpIF;
	double TAR_RATE, dThreRate = 0.0001, dPI = atan(1.0) * 4.0;
	if (argc < 5)
	{
		printf("Usage: NewEntropyCoding_Enc Input_raw_file Output_file Target_PSNR_multiple100(integer num) SDQ-on/off SDQ_lambda(optional) ADPT_DIST_PROFL-on/off(optional) Target_rate_SDQ(optional/need further operation) Turn_on_ECEBfilter NOTE: image size and num of DCT blocks should be set in global.h \n");
		exit(1);

	}

	strcat(inpath, "\\");

	strcat(inpath, input);

	sprintf(rgcIFName, inpath);

	//std::cout << "input is" << rgcIFName << std::endl;

	
	SDQ = atoi(argv[7]);
	if (SDQ == 0)
		ADPT_DIST_PROFL = 0;
	fpIF = fopen(rgcIFName, "rb");

	for (int i = 0; i < 3 * row; i++) {
		fread(image[i], sizeof(unsigned char), col, fpIF);
	}
	fclose(fpIF);
	
	for (int i = 0; i < row; i++)
		for (int j = 0; j < col; j++)
			drgucInput_Y[i][j] = image[i][j]; 
	
	int r = 0; 
	int c = 0; 

	for (int i = row ; i < 2*row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			drgucInput_U[r][c] = image[i][j];
			c++;
		}
		c = 0; 
		r++; 
	}
	r = 0; 
	c = 0; 

	for (int i = 2*row; i < 3*row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			drgucInput_V[r][c] = image[i][j] ;
			c++;
		}
		c = 0;
		r++;    
	} 
	double time, begin, end;
	begin = clock(); 
	
	encode(drgucInput_Y, input, 1, TAR_PSNR_ARG);
	std::cout << "first is done" << std::endl; 
	encode(drgucInput_U, input, 2, TAR_PSNR_ARG);
	std::cout << "second is done" << std::endl;
	encode(drgucInput_V, input, 3, TAR_PSNR_ARG);
	std::cout << "third is done" << std::endl;
	write_picture(Y_picture_buff , Cb_picture_buff, Cr_picture_buff);

	end = clock(); 
	time = (end - begin) / CLOCKS_PER_SEC; 
	cout << "time of encoding is" << time << endl; 
	
	
	for (int i = 0; i < ROWS; i++)
	{
		free(drgucInput_Y[i]);
	}
	free(drgucInput_Y);


	for (int i = 0; i < ROWS; i++)
	{
		free(drgucInput_U[i]);
	}
	free(drgucInput_U);


	for (int i = 0; i < ROWS; i++)
	{
		free(drgucInput_V[i]);
	}
	free(drgucInput_V);

	for (int i = 0; i < ROWS; i++)
	{
		free(Y_picture_buff[i]);
	}
	free(Y_picture_buff);


	for (int i = 0; i < ROWS; i++)
	{
		free(Cb_picture_buff[i]);
	}
	free(Cb_picture_buff);


	for (int i = 0; i < ROWS; i++)
	{
		free(Cr_picture_buff[i]);
	}
	free(Cr_picture_buff);
	//std::cin.get();
	return 0;
}


void encode(unsigned char** drgucInput_Y ,char * input , int ch , int psnr   ) {


	TAR_PSNR_ARG = psnr;
	
	
	/***************  Globals  *****************/
	initialize_arithmetic_encoder();
	for (int i = 0; i < ROWS; i++)
		std::fill(drgiQIdxHDQ[i].begin(), drgiQIdxHDQ[i].end() , 0 );
	
	for (int i = 0; i < ROWS; i++)
		std::fill(drgiQIdxSDQ[i].begin(), drgiQIdxSDQ[i].end(), 0);


	for (int i = 0; i < ROWS/8; i++)
		std::fill(drgiQIdxHDQAjstDCDec[i].begin(), drgiQIdxHDQAjstDCDec[i].end(), 0);

	for (int i = 0; i < ROWS/8; i++)
		std::fill(drgdQIdxHDQDCPred[i].begin(), drgdQIdxHDQDCPred[i].end(), 0);

	for (int i = 0; i < ROWS/8; i++)
		std::fill(drgiQIdxHDQAjstDC[i].begin(), drgiQIdxHDQAjstDC[i].end(), 0);
	/****************************   class variable    *************************/
	for (int i = 0; i < ROWS / 8; i++)
		std::fill(drgiExtF[i].begin(), drgiExtF[i].end(), 0);

	for (int i = 0; i < ROWS / 8; i++)
		std::fill(drgiExtFDec[i].begin(), drgiExtFDec[i].end(), 0);

	for (int i = 0; i < ROWS / 8; i++)
		std::fill(drgiCBP[i].begin(), drgiCBP[i].end(), 0);


	for (int i = 0; i < ROWS / 8; i++)
	{
		for (int j = 0; j < COLS/8; j++)
		{
			for (int k = 0; k < 8; k++)
				for (int m = 0; m < 8; m++)

					vrgiSign[i][j][k][m] = 0; 
		}

	}
	// HACAM  VARS

	for (int i = 0; i < MAX_ROWS / 8; i++)
	{
		for (int j = 0; j < MAX_COLS / 8; j++)
		{
			for (int k = 0; k < 8; k++)
				for (int m = 0; m < 8; m++)
					for (int n = 0 ; n < BIT_LOC_EXT ; n++)
						vrgiExt[i][j][k][m][n] = 0;
		}

	}


	for (int i = 0; i < MAX_ROWS / 8; i++)
	{
		for (int j = 0; j < MAX_COLS / 8; j++)
		{
			for (int k = 0; k < 8; k++)
				for (int m = 0; m < 8; m++)
					for (int n = 0; n < BIT_LOC_EXT; n++)
						vrgiExtDec[i][j][k][m][n] = 0;
		}

	}



	for (int i = 0; i < MAX_ROWS / 8; i++)
	{
		for (int j = 0; j < MAX_COLS / 8; j++)
		{
			for (int k = 0; k < 8; k++)
				for (int m = 0; m < 8; m++)
					for (int n = 0; n < BIT_LOC_EXT; n++)
						vrgiSig[i][j][k][m][n] = 0;
		}

	}

	for (int i = 0; i < MAX_ROWS / 8; i++)
	{
		for (int j = 0; j < MAX_COLS / 8; j++)
		{
			for (int k = 0; k < 8; k++)
				for (int m = 0; m < 8; m++)

					vrgiHDQIdxDec[i][j][k][m] = 0;
		}

	}

	dRateTmp = 0; 
	i2pAmpCtrEmp = 0; 
	dTestRate = 0; 
	dTestRateHDQ = 0; 
	iErrorNumSDQ = 0 ; 
	 
	/************************************************/


	//double drgdCoeff[ROWS][COLS];
	double** drgdCoeff = new double* [ROWS];
	for (int c = 0; c < ROWS; c++)
		drgdCoeff[c] = new double[COLS] {0.0};

	//double drgdRecCoeff[ROWS][COLS] = { 0.0 }; 
	double** drgdRecCoeff = new double* [ROWS];
	for (int c{ 0 }; c < ROWS; c++)
		drgdRecCoeff[c] = new double[COLS] {0.0};
	int** DiffIdxDC;
	std::vector<std::vector<int>>drgiQIdxDC(ROWS / 8, std::vector<int>(COLS / 8));
	int* RefIdxDC, ** TmpIdxDC, * rgiExtQntSeq, * rgiSigQntSeq, iOptCtxDC, iCurQntIdx;
	double temp1, dAvgCostRDPrev;// dDistOutYcPlusDC=0.0;
	double drgdBlk[8][8];
	double rgdDeltaDC[3] = { 90.0,180.0,265.0 };//{90.0,180.0,265.0};
	double dAvgDist, dDistSDQ, dRateSDQ, dMinCost, dScl, dPredDC, dPredDCAjst[NUM_CTX_PRED_DC], dPredDCAcum, dCurDiff;
	int iCurDiff, iPredAC, iCurSign, iCurAbs, iA, iB, iC, iCorrect, iFreqRow, iFreqCol;
	int iRealEncIdx[MAX_ROWS / 8][MAX_COLS / 8] = { 0 };
	double begin, end, time0, time1, time2, dTotlBegn, dTotlEnd, dTotlTime;
	int iRateOvrhd;
	int iCounter, iCoeffIdx, iPredCtxDC, iCtxDCIdx, iCtrDiff[NUM_CTX_PRED_DC] = { 0 }, iCurPredDCAjst;
	double dLambdaL, dLambdaH, dAcumDiff[NUM_CTX_PRED_DC] = { 0.0 };
	int rgiZigZagIdx[64];
	int rgiZigZagExt[64];
	int rgiZigZagIdxPred[64], iCurBitLoc, iCurDecIdx;
	register int i, j, u, v, iIdx;
	for (i = 0; i < 8; i++)
		for (j = 0; j < 8; j++)
			dDctSclF[i][j] = _C[i] * _C[j] * 16.0F;


	/* allocate memory */
	GACpEngine = (ArithmeticCode*)malloc(sizeof(ArithmeticCode));
	// initialize arithmatic coding
	fnIntlAC(GACpEngine);

	GbspOutputBuffer = (BitStream*)malloc(sizeof(BitStream));
	fnIntlBitStream(GbspOutputBuffer, IntlBfSz);

	/* forward DCT */
	dTotlBegn = clock();
	begin = clock();
	for (i = 0; i < ROWS / 8; i++)
	{
		for (j = 0; j < COLS / 8; j++)
		{
			/* change the format of unsigned char to integer before computation and remove the mean */
			for (u = 0; u < 8; u++)
				for (v = 0; v < 8; v++)
					drgdBlk[u][v] = (double)((int)drgucInput_Y[i * 8 + u][j * 8 + v] - 128);

			FDCT(drgdBlk);

			for (u = 0; u < 8; u++)
			{
				for (v = 0; v < 8; v++)
				{
#if EMP_MEAN_SUB
					drgdCoeff[i * 8 + u][j * 8 + v] = drgdBlk[u][v] / dDctSclF[u][v] - drgdEmpMean[u][v];
#else
					drgdCoeff[i * 8 + u][j * 8 + v] = drgdBlk[u][v] / dDctSclF[u][v];
#endif
					dMaxDC = MAXIMUM(dMaxDC, drgdCoeff[i * 8][j * 8]);
					dMinDC = MINIMUM(dMinDC, drgdCoeff[i * 8][j * 8]);
				}
			}
		}
	}

	end = clock();
	time0 = ((double)end - begin) / CLOCKS_PER_SEC;


	begin = clock();
#if INTEGER_TCM
	std::cout << "TCM " << std::endl;
	fnTCM(drgdCoeff);
#else
	fnDistrParaTCM(drgdCoeff);
#endif
	end = clock();
	time1 = ((double)end - begin) / CLOCKS_PER_SEC;
	//printf("CPU time for TCM algorithm: %2.6f seconds \n", time1);

	begin = clock();
	/* setup the reconstruction space */
	fnOptMdlQ(drgdCoeff);

	dScl = 255.0f / (double)iMaxLDC;
	DiffIdxDC = (int**)malloc((iMaxLDC + 1) * sizeof(int*));
	for (j = 0; j <= iMaxLDC; j++)
		DiffIdxDC[j] = (int*)malloc((iMaxLDC + 1) * sizeof(int));

	RefIdxDC = (int*)malloc((1 + 2 * iMaxLDC) * sizeof(int));

	TmpIdxDC = (int**)malloc((iMaxLDC + 1) * sizeof(int*));
	for (j = 0; j <= iMaxLDC; j++)
		TmpIdxDC[j] = (int*)malloc((iMaxLDC + 1) * sizeof(int));

	for (i = 0; i < NUM_CTX_DC; i++)
	{
		rgiAlpSizeDC[i] = fnRound((double)rgiAlpSizeDC[i] / dScl);
	}
	rgiAlpSizeDC[NUM_CTX_DC - 1] = iMaxLDC + 1;

	/* setup the Diff index information for DC */
	RefIdxDC[0] = 0;
	iIdx = 1;
	for (i = 1; i <= iMaxLDC; i++)
	{
		RefIdxDC[iIdx++] = i;
		RefIdxDC[iIdx++] = -i;
	}

	for (j = 0; j <= iMaxLDC; j++)
	{
		for (i = 0; i <= iMaxLDC; i++)
			TmpIdxDC[j][i] = i - j;

		iIdx = 0;
		for (v = 0; v <= 2 * iMaxLDC; v++)
		{
			for (i = 0; i <= iMaxLDC; i++)
			{
				if (TmpIdxDC[j][i] == RefIdxDC[v])
				{
					DiffIdxDC[j][i] = iIdx++;
					break;
				}
			}
		}

		if (iIdx != iMaxLDC + 1)
			printf("Error: Diff Index Matrix is wrong!\n");
	}


	fnLyrStup();
	fnMemAllocCtr();


	rgiExtQntSeq = (int*)malloc(sizeof(int) * LyrInfo.ExtNumElemt);
	rgiSigQntSeq = (int*)malloc(sizeof(int) * LyrInfo.SigNumElemt);
	if (SDQ)
		fnCtrIni();

	for (i = 0; i < ROWS / 8; i++)
	{
		for (j = 0; j < COLS / 8; j++)
		{

			for (iIdx = 0; iIdx < LyrInfo.ExtNumElemt; iIdx++)
			{
				iFreqRow = LyrInfo.ExtScanRowIdx[iIdx];
				iFreqCol = LyrInfo.ExtScanColIdx[iIdx];

				/* Quantization */
				rgiExtQntSeq[iIdx] = fnHDQ(drgdCoeff[i * 8 + iFreqRow][j * 8 + iFreqCol], i, j, iFreqRow, iFreqCol);
				drgiQIdxHDQ[i * 8 + iFreqRow][j * 8 + iFreqCol] = rgiExtQntSeq[iIdx];
			}

			/* Bi-level image coding */
			fnEncBlkBiLvImg(rgiExtQntSeq, i, j);
		}
	}
///	std::cout << "multiresolution" << std::endl;
#if MULTIRESOLUTION
	iOutlierRate = GbspOutputBuffer->m_nBits;
#endif

	for (i = 0; i < ROWS / 8; i++)
	{
		for (j = 0; j < COLS / 8; j++)
		{
			iCurQntIdx = fnHDQ(drgdCoeff[i * 8][j * 8], i, j, 0, 0);// returns back the quantized DC coeffs 

			drgiQIdxHDQAjstDC[i][j] = iCurQntIdx + iLLowDC;
			dPredDC = fnPredCoding(i, j, dScl);
			iCtxDCIdx = 0;
			iOptCtxDC = NUM_CTX_DC - 1;
			while (iCtxDCIdx < (NUM_CTX_DC - 1))
			{
				if (dDeltaDC < rgdDeltaDC[0])
				{
					iOptCtxDC = iCtxDCIdx;
					break;
				}
				iCtxDCIdx++;
			}
			iPredCtxDC = iOptCtxDC;//+drgiExtF[i][j]*NUM_CTX_DC;

			dPredDCAjst[iPredCtxDC] = dPredDC + dAcumDiffAvg[iPredCtxDC];
			if (dPredDCAjst[iPredCtxDC] > iMaxLDC)
				dPredDCAjst[iPredCtxDC] = iMaxLDC;
			if (dPredDCAjst[iPredCtxDC] < 0)
				dPredDCAjst[iPredCtxDC] = 0;

			//dCurDiff = (double)drgiQIdxHDQAjstDC[i][j] - (double)fnRound(dPredDCAjst[iPredCtxDC]);
			dCurDiff = (double)iCurQntIdx + (double)iLLowDC - (double)fnRound(dPredDCAjst[iPredCtxDC]);
			iCurDiff = fnRound(dCurDiff);
			drgdQIdxHDQDCPred[i][j] = iCurDiff;
			dAcumDiff[iPredCtxDC] += iCurDiff;
			iCtrDiff[iPredCtxDC]++;

			//iCurDiff = dAcumDiffAvg[iPredCtxDC] < 0? -iCurDiff:iCurDiff;
			iCurPredDCAjst = fnRound(dPredDCAjst[iPredCtxDC]);
			//if (iCurDiff+iCurPredDCAjst < 0 || iCurDiff+iCurPredDCAjst > iMaxLDC)
			// iCurDiff = -iCurDiff;
			iRealEncIdx[i][j] = DiffIdxDC[iCurPredDCAjst][iCurDiff + iCurPredDCAjst];
			//            if (drgiQIdxHDQAjstDC[i][j]!=iCurDiff+iCurPredDCAjst)
			//                printf("error:encoding DC\n");

#if 0//DC_RESCALE
			if (iCtrDiff[iPredCtxDC] == 128)
			{
				iCtrDiff[iPredCtxDC] = 64;
				dAcumDiff[iPredCtxDC] /= 2.0;
			}
#endif
			dAcumDiffAvg[iPredCtxDC] = dAcumDiff[iPredCtxDC] / (double)iCtrDiff[iPredCtxDC];

			fnEncDC(iOptCtxDC, iRealEncIdx[i][j]);
		}
	}

#if MULTIRESOLUTION
	iDcRate = GbspOutputBuffer->m_nBits;
#endif

	double dTimeACBegin, dTimeACEnd;
	dTimeACBegin = clock();
	for (i = 0; i < ROWS / 8; i++)
	{
		for (j = 0; j < COLS / 8; j++)
		{
			/* encode the quantized coefficients within Yc */
			//begin = clock();
			if (SDQ)
			{
				fnColecStatisBlk(drgiQIdxHDQ, i, j);
			}
			else
			{
				for (iIdx = 0; iIdx < LyrInfo.SigNumElemt; iIdx++)
				{
					iFreqRow = LyrInfo.SigScanRowIdx[iIdx];
					iFreqCol = LyrInfo.SigScanColIdx[iIdx];
					rgiSigQntSeq[iIdx] = drgiQIdxHDQ[i * 8 + iFreqRow][j * 8 + iFreqCol];
				}
				fnEncOneBlk(rgiSigQntSeq, i, j);
			}

#if EMP_MEAN_SUB
			/* Add back Mean */
			for (u = 0; u < 8; u++)
			{
				for (v = 0; v < 8; v++)
				{
					drgdCoeff[i * 8 + u][j * 8 + v] += drgdEmpMean[u][v];
					drgdRecCoeff[i * 8 + u][j * 8 + v] += drgdEmpMean[u][v];
				}
			}
#endif
#if ENC_DEBUG
			fnPrintDataArray(drgiBlk, rgcDebugrgcIFName);
#endif

		}
	}

	dTimeACEnd = clock();
	//printf("\n AC running time = %f\n", ((double)dTimeACEnd - dTimeACBegin) / CLOCKS_PER_SEC);
	if (SDQ)
	{
		//iNumIter = 0;
		dAvgCostRDPrev = INFINITY;
		for (;;)
		{
			fnCalEntrpy();

			dAvgCostRD = 0.0;
			dDistCurIterSDQ = 0.0;
			for (i = 0; i < 8; i++)
			{
				for (j = 0; j < 8; j++)
				{
					drgdSumYmulC[i][j] = 0.0;
					drgdSumCmulC[i][j] = 0.0;
					drgdSumC[i][j] = 0.0;
					drgdSumY[i][j] = 0.0;
					drgiLTblSDQ[i][j] = 0;
					drgiNumBlkClgt0[i][j] = 0;
				}
			}

			for (i = 0; i < ROWS / 8; i++)
			{
				for (j = 0; j < COLS / 8; j++)
				{
					if (drgiCBP[i][j] == 1)
						fnBlkSDQ(i, j, drgdCoeff);
				}
			}

			dAvgCostRD /= COLS * ROWS;

			//fnUpdateQtbl();

			if (iNumIter > -1)
				break;
			// check whether the iteration converges
			if (dAvgCostRDPrev - dAvgCostRD < 0.001)
				break;

			if (iNumIter == 0)
			{
				printf("The average HDQ RD cost is: %4.4f\n", (dCostHDQ / ((double)(COLS * ROWS))));
				printf("The average SDQ RD cost for the %dth iteration is: %4.4f\n", iNumIter + 1, dAvgCostRD);
			}
			else
				printf("The average SDQ RD cost for the %dth iteration is: %4.4f\n", iNumIter + 1, dAvgCostRD);


			fnCtrIni();
			for (i = 0; i < ROWS / 8; i++)
			{
				for (j = 0; j < COLS / 8; j++)
				{
					fnColecStatisBlk(drgiQIdxSDQ, i, j);
				}
			}

			dAvgCostRDPrev = dAvgCostRD;
			iNumIter++;
		}
		//dAvgDist = MSE2PSNR((dDistCurIterSDQ + dDistOutYcPlusDC)/ROWS/COLS);

		for (i = 0; i < ROWS / 8; i++)
		{
			for (j = 0; j < COLS / 8; j++)
			{
				for (iIdx = 0; iIdx < LyrInfo.SigNumElemt; iIdx++)
				{
					iFreqRow = LyrInfo.SigScanRowIdx[iIdx];
					iFreqCol = LyrInfo.SigScanColIdx[iIdx];
					rgiSigQntSeq[iIdx] = drgiQIdxSDQ[i * 8 + iFreqRow][j * 8 + iFreqCol];
				}

				/* encode the quantized coefficients within Yc */
				//begin = clock();
				fnEncOneBlk(rgiSigQntSeq, i, j);
				//end = clock();
				//time = ((double)end-begin)/CLOCKS_PER_SEC;
				//printf("%2.6f seconds has been used in coding block at the %dth ROWS and %dth colum\n", time,i,j);

#if EMP_MEAN_SUB
				/* Add back Mean */
				for (u = 0; u < 8; u++)
				{
					for (v = 0; v < 8; v++)
					{
						drgdCoeff[i * 8 + u][j * 8 + v] += drgdEmpMean[u][v];
						drgdRecCoeff[i * 8 + u][j * 8 + v] += drgdEmpMean[u][v];
					}
				}
#endif
			}
		}
	}
	end = clock();
	time2 = ((double)end - begin) / CLOCKS_PER_SEC;
	dTotlEnd = clock();
	dTotlTime = (dTotlEnd - dTotlBegn) / CLOCKS_PER_SEC;
	/* calculate distortion (see the decoder project for real decoding) */
	for (i = 0; i < NUM_CTX_PRED_DC; i++)
	{
		dPredDCAjst[i] = 0.0;
		iCtrDiff[i] = 0;
		dAcumDiff[i] = 0.0;
	}

	/*for (i = 0; i < ROWS / 8; i++)
	{
		for (j = 0; j < COLS / 8; j++)
		{
			fprintf(fpIF, "%d ", iRealEncIdx[i][j]);
		}
		fprintf(fpIF, "\n");
	}
	fclose(fpIF);
*/
	for (i = 0; i < ROWS / 8; i++)
	{
		for (j = 0; j < COLS / 8; j++)
		{
			drgiQIdxDC[i][j] = iRealEncIdx[i][j];
		}
	}

	//fclose(fpIF);

	for (i = 0; i < ROWS / 8; i++)
	{
		for (j = 0; j < COLS / 8; j++)
		{
			/* Predictive Coding for DC */
			dPredDC = fnPredCoding(i, j, dScl);
			iCtxDCIdx = 0;
			iOptCtxDC = NUM_CTX_DC - 1;
			while (iCtxDCIdx < (NUM_CTX_DC - 1))
			{
				if (dDeltaDC < rgdDeltaDC[0])
				{
					iOptCtxDC = iCtxDCIdx;
					break;
				}
				iCtxDCIdx++;
			}
			iPredCtxDC = iOptCtxDC;//+drgiExtF[i][j]*3;

			dPredDCAjst[iPredCtxDC] = dPredDC + dAcumDiffAvg[iPredCtxDC];
			if (dPredDCAjst[iPredCtxDC] > iMaxLDC)
				dPredDCAjst[iPredCtxDC] = iMaxLDC;
			if (dPredDCAjst[iPredCtxDC] < 0)
				dPredDCAjst[iPredCtxDC] = 0;
			iCurPredDCAjst = fnRound(dPredDCAjst[iPredCtxDC]);

			for (iIdx = 0; iIdx <= iMaxLDC; iIdx++)
			{
				if (DiffIdxDC[iCurPredDCAjst][iIdx] == drgiQIdxDC[i][j])
					break;
			}

			iCurDiff = TmpIdxDC[iCurPredDCAjst][iIdx];
			//iCurDiff = dAcumDiffAvg[iPredCtxDC] < 0? -iCurDiff:iCurDiff;
			iCurPredDCAjst = fnRound(dPredDCAjst[iPredCtxDC]);
			//if (iCurDiff+iCurPredDCAjst < 0 || iCurDiff+iCurPredDCAjst > iMaxLDC)
			//iCurDiff = -iCurDiff;

			drgiQIdxHDQAjstDCDec[i][j] = iCurDiff + iCurPredDCAjst;
			if (drgiQIdxHDQAjstDCDec[i][j] != drgiQIdxHDQAjstDC[i][j])
				printf("error: DC decoding!\n");

			drgdRecCoeff[i * 8][j * 8] = fnDeQtz(drgiQIdxHDQAjstDC[i][j] - iLLowDC, i, j, 0, 0);

			drgdQIdxHDQDCPred[i][j] = iCurDiff;
			dAcumDiff[iPredCtxDC] += iCurDiff;
			iCtrDiff[iPredCtxDC]++;
			dAcumDiffAvg[iPredCtxDC] = dAcumDiff[iPredCtxDC] / (double)iCtrDiff[iPredCtxDC];

			if (drgiExtFDec[i][j] == 1)
			{
#if MULTIRESOLUTION
				OBF_IMG[i][j] = (unsigned char)255;
#endif
				for (iIdx = 0; iIdx < LyrInfo.ExtNumElemt; iIdx++)
				{
					iFreqRow = LyrInfo.ExtScanRowIdx[iIdx];
					iFreqCol = LyrInfo.ExtScanColIdx[iIdx];

					if (vrgiExtDec[i][j][iFreqRow][iFreqCol][0] == 1)
					{
#if MULTIRESOLUTION
						OTF_IMG[i * 8 + iFreqRow][j * 8 + iFreqCol] = 1;
#endif
						iCurBitLoc = 1;
						while (iCurBitLoc < BIT_LOC_EXT)
						{
							if (vrgiExtDec[i][j][iFreqRow][iFreqCol][iCurBitLoc] == 0)
								break;
							else
								iCurBitLoc++;
						}

						if (iCurBitLoc == BIT_LOC_EXT)
							iCurBitLoc += vrgiHDQIdxDec[i][j][iFreqRow][iFreqCol];

						drgdRecCoeff[i * 8 + iFreqRow][j * 8 + iFreqCol] = fnDeQtzExt(iCurBitLoc * vrgiSign[i][j][iFreqRow][iFreqCol], iFreqRow, iFreqCol);
					}

					else if (vrgiSig[i][j][iFreqRow][iFreqCol][0] == 1)
					{
						iCurBitLoc = 1;
						while (iCurBitLoc < BIT_LOC_SIG)
						{
							if (vrgiSig[i][j][iFreqRow][iFreqCol][iCurBitLoc] == 0)
								break;
							else
								iCurBitLoc++;
						}

						if (iCurBitLoc == BIT_LOC_SIG)
							iCurBitLoc += vrgiHDQIdxDec[i][j][iFreqRow][iFreqCol];

						drgdRecCoeff[i * 8 + iFreqRow][j * 8 + iFreqCol] = fnDeQtzSig(iCurBitLoc * vrgiSign[i][j][iFreqRow][iFreqCol], iFreqRow, iFreqCol);
					}
				}
			}

			else if (drgiCBP[i][j] == 1)
			{
				for (iIdx = 0; iIdx < LyrInfo.SigNumElemt; iIdx++)
				{
					iFreqRow = LyrInfo.SigScanRowIdx[iIdx];
					iFreqCol = LyrInfo.SigScanColIdx[iIdx];

					if (vrgiSig[i][j][iFreqRow][iFreqCol][0] == 1)
					{
						iCurBitLoc = 1;
						while (iCurBitLoc < BIT_LOC_SIG)
						{
							if (vrgiSig[i][j][iFreqRow][iFreqCol][iCurBitLoc] == 0)
								break;
							else
								iCurBitLoc++;
						}

						if (iCurBitLoc == BIT_LOC_SIG)
							iCurBitLoc += vrgiHDQIdxDec[i][j][iFreqRow][iFreqCol];

						drgdRecCoeff[i * 8 + iFreqRow][j * 8 + iFreqCol] = fnDeQtzSig(iCurBitLoc * vrgiSign[i][j][iFreqRow][iFreqCol], iFreqRow, iFreqCol);
					}

				}
			}
		}
	}


	fnEndEncode(GACpEngine, GbspOutputBuffer);

	///printf("\n");
	//printf("********************************************************************************\n");
	////printf("CPU time for DCT: %2.6f seconds \n", time0);
	////printf("CPU time for TCM algorithm: %2.6f seconds \n", time1);
	////printf("CPU time for quantization and entropy coding: %2.6f seconds \n", time2);
	//printf("Total CPU time: %2.6f seconds \n", time0 + time1 + time2);
	//printf("Total bits (excluding overhead) are: %ld\n", (GbspOutputBuffer->m_nBits));

	iRateOvrhd = 63 + 20 + 6;//iLLowDC and iLHighDC: 20-bit + Ttbl: 63*(1-bit) + PSNR index: 6-bit;

	for (u = 0; u < 8; u++)
	{
		for (v = 0; v < 8; v++)
		{
			if ((u + v) > 0)
			{
				if (drgiTTbl[u][v] == 1)
				{
#if ENC_EXTAMP_SEPRT
					if (TAR_PSNR < 3800)
						iRateOvrhd += 20;//Lambda: 8-bit + Gamma: 6-bit + L_prime: 6-bit
					else if (TAR_PSNR < 4400)
						iRateOvrhd += 21;
					else
						iRateOvrhd += 22;
#else
					iRateOvrhd += 14;
#endif
				}
			}
		}
	}
	if (ADPT_DIST_PROFL)
	{
		iRateOvrhd += 14;
	}

	dRate = (float)(GbspOutputBuffer->m_nBits) / (COLS * ROWS);

	//    printf("Average bit test SDQ rate in BPP excluding overheads is: %3.4f\n", dTestRate/512/512/SDQ_LAMBDA);
	//    printf("Average bit test HDQ rate in BPP excluding overheads is: %3.4f\n", dTestRateHDQ/512/512/SDQ_LAMBDA);
	printf("Average bit rate in BPP excluding overheads is: %3.4f\n", dRate);

	bpp_vec.push_back(dRate); 
	std::cout << "bpp_vec "<< bpp_vec[ch -1 ] << std::endl; 
	//dRate += (float)iRateOvrhd / (COLS * ROWS);
	//printf("Average bit rate in BPP including overheads is: %3.4f\n", dRateTmp/ROWS/COLS);
	//printf("Average bit rate in BPP including overheads is: %3.4f\n", dRate);
	//printf("Average mse after SDQ is: %3.4f\n", dTestDist/512/512);

	// fnWtoFile(GbspOutputBuffer, rgcOFName);
	
	float image_psnr; 
	//image_psnr = cal_dist(drgdCoeff, drgdRecCoeff, drgucInput_Y);
	image_psnr = write_buff(drgdCoeff, drgdRecCoeff, drgucInput_Y, ch);
	psnr_vec.push_back(image_psnr); 
	/*FILE* fptxt;
	fptxt = fopen(rgctextOFName, "wb");
	fprintf(fptxt, "PSNR\t%f\n", image_psnr);
	fprintf(fptxt, "BPP\t%f\n", dRate);
	fclose(fptxt);
*/
	/* free space */
	free(GACpEngine);
	GACpEngine = NULL;
	fnFreeBitStream(GbspOutputBuffer);
	free(GbspOutputBuffer);
	GbspOutputBuffer = NULL;




	for (j = 0; j <= iMaxLDC; j++)
	{
		free(DiffIdxDC[j]);
		free(TmpIdxDC[j]);
	}
	free(DiffIdxDC);
	free(RefIdxDC);
	free(TmpIdxDC);
	free(rgiExtQntSeq);
	free(rgiSigQntSeq);
	for (i = 0; i < NUM_QTZ_GAMMA; i++)
	{
		free(Q_cache[i]);
	}
	free(Q_cache);
	fnFreeMemHACACM();
	for (i = 0; i < ROWS; i++)
	{
		delete[] drgdRecCoeff[i];

	}
	delete[] drgdRecCoeff;

	for (i = 0; i < ROWS; i++)
	{
		delete[] drgdCoeff[i];

	}
	delete[] drgdCoeff;
	
}




