#pragma once
#ifndef _HACACM_H_
#define _HACACM_H_

#include "AdapAC.h"
#include <vector>
#define MULTIRESOLUTION 0
#define MULTIRESOLUTION_DC 1//if DC decoding if possible or not
#define DIFFDIST 0//different d for outlier and inlier to test visual quality

#define INTEGER_TCM 0
//#define BAC_EN 0
#define DC_SEP_CTX 0

#define HDQ_EMP_DIST_PROFL 1
#define TEST_RD 0
#define OUTPUT_REC_IMG 1
//#define SDQ 1
//#define ADPT_DIST_PROFL 1
#define NUM_DIST_TRELIS 8
#define JUNE26_SIG_LYR_NEIBER 0

#define SDQ_TEST_SIG

#define ADP_SCN 0

#define NUM_CTX_DC 3
#define SIG_SCH 0//1:use 4 neighbors;0:use 5 neighbors
#define SIG_MERG 1
#define SIG_SEPRT 0//1: encoding of SIG is seperated into two parts, i.e., acumulated 1s and 5 neighbors
#define CTX_INI 0//0:no context initialization scheme
//1:context initialization scheme (only limited context models are initialized)
#define IMMED_NEIBR 0//1: use the immediate neigbour in the scanning layer
//0: use whoever available in the scanning layer that is closest to the current encoding frequency

#define LST2LYR_MERG 0//if merge the contexts for the last 2 context layers or not
#define BIN_SCH 1//if binarilization scheme is performed
#define BIT_LOC_EXT 3//binarilization stops after this number (valid when BIN_SCH = 1)
#define BIT_LOC_SIG 3
#define SIGN_CTX_MODL 0
#define SIGN_AC2AC3_CTX_MODL 1//SIGN_CTX_MODL=1 -> SIGN_AC2AC3_CTX_MODL=1, or use/not use the context modeling based on the neighbouring blocks to encode sign bits of AC2 and AC3
#define EXT_BLK_NEIBR_PREVBIT 0//1:using future block correlation (for BitPos-1) rather than subband correlation (Prof. Yang's suggestion) when coding Ext (BitPos>1)
#define SIG_EXT_HELP 0//1: using already encoded collocated Ext of its subband neighbours info when coding Sig
#define INTER_BLK_SIG_HELP 1//1: using already encoded collocated Sig of its subband neighbours info when coding Sig

#define EXT_WITH_CTX 1

#define COST_FUN 1//1:including the cost for DC into the cost function;
//0: NOT including the cost for DC into the cost function
#define BLK_SIM_AMP_ENC 0//1: use the block similarity coding for the indices inside Yc
#define BLK_SIM_SIGN_ENC 0//1: use the block similarity coding for the indices inside Yc

#define BINARLZ 1//0: arithmetic coding
//1: binary arithmetic coding
#define JNT_EXTSIG 0//1:jointly encode the extension bit and significant bit in EXT blocks
#define ENC_SIG_FRST 0//1:encode AZF first and then ExtF; and encode SIG first and then EXT;0: the other way around
#define JNT_AZFEXTF 0//1:jointly encode the extension flag and all-zero flag for each block
#define ENC_EXTAMP_SEPRT 1//1:encode the extension amplitude for each subband sepreately
#define JNT_SIG_AC2AC3 0//jointly encode the significant bit of AC2 and AC3

#define ENC_EXT 1
#define PRED_EXT 0

#define ENC_METD_EXT 0//0: using arithmetic coding to encode the extension part
//1: using Rice-Golumb coding to encode the extension part

#define NUM_CTX_SIG 13
#define NUM_CTX_AMP 1
#define NUM_CTX_EXTF 16
#define NUM_CTX_EXT 13//maximum number of context models may be used
#define NUM_CTX_EXTAMP 9
#define NUM_CTX_BLKSIMF 5
#define NUM_CTX_BLKSIM_SIGNF 5
#define NUM_CTX_AZF_EXTF 12
#define NUM_CTX_EXTAMP_SEPRT 9
#define NUM_CTX_SIG_AC2AC3 4
#define NUM_CTX_SIGN_AC2AC3 9
#define NUM_CTX_SIGN 9

//define NUM_CTX_AZF based on the target distortion
#if (TAR_PSNR-32)
#define NUM_CTX_AZF 7 //for low rate
#else
#define NUM_CTX_AZF 5 //for high rate
#endif
/* not used */
//extern int drgiAZFSDQ[ROWS / 8][COLS / 8];
//extern int drgiExtFSDQ[ROWS / 8][COLS / 8];
////extern int ctx_ctr[2][16];


//extern int vrgiSign[ROWS / 8][COLS / 8][8][8];
extern std::vector<std::vector<std::vector<std::vector<int>>>> vrgiSign;

//extern int vrgiExt[ROWS / 8][COLS / 8][8][8][BIT_LOC_EXT];
extern int vrgiExt[MAX_ROWS/8][MAX_COLS/8][8][8][BIT_LOC_EXT];


//extern int vrgiExtDec[ROWS / 8][COLS / 8][8][8][BIT_LOC_EXT];
extern int vrgiExtDec[MAX_ROWS / 8][MAX_COLS / 8][8][8][BIT_LOC_EXT];

//extern int vrgiHDQIdxDec[ROWS / 8][COLS / 8][8][8];
extern int vrgiHDQIdxDec[MAX_ROWS / 8][MAX_COLS / 8][8][8];

//extern int vrgiSig[ROWS / 8][COLS / 8][8][8][BIT_LOC_SIG];
extern int vrgiSig[MAX_ROWS / 8][MAX_COLS / 8][8][8][BIT_LOC_SIG];


extern double dRateTmp;
extern int **i2pAmpCtrEmp;

typedef struct {
	int ExtNumElemt;
	int SigNumElemt;
	int SigNumLyr;

	int CtxLyr[63];//7 layers in total
	int CtxLyrSig[63];//the layers for context emerge for Sig and SigExt
	int CtxLyrIdxMaxL[63];
	int CtxLyrIdxMaxLprime[63];

	//int ScanLyrIdxSig[63];

	int LyrIdx_MaxL;
	int LyrIdx_MaxLprime;

	int Acum1Pos;

	int ExtScanRowIdx[63];
	int ExtScanColIdx[63];
	int SigScanRowIdx[63];
	int SigScanColIdx[63];

	int SigLyrNumElemt[16];

	int SwapFlag[14];
	int SigNumCtx[63];
	int SigCtxMdl[63][6];

	int SigAphSize[63];

}Layer;

typedef struct {
	//double min_cost_SDQ;
	double total_cost_HDQ;
	double total_rate_HDQ;
	int opt_amp[8];
	double s_rate[8][4];
	double s_dist[8];
	double s_cost[8][4];
	int s_prvste[8][4];
	int s_SigCtx[8][4];
	int s_opt_amp[8];
	int min_cost_prvste_all;
	int min_cost_prvste_lst3;
	int min_cost_prvste_lst2;
}StateSDQ;

extern Layer LyrInfo;
extern double dTestRate;
extern double dTestRateHDQ;
extern int iErrorNumSDQ;

extern ArithmeticCode* GACpEngine;
extern BitStream* GbspOutputBuffer;

void fnLyrStup();
void fnMemAllocCtr();
void fnEncSign(int iRowIdx, int iColIdx, int iBlkRowIdx, int iBlkColIdx, int iCoef, int iExt);
void fnACItrfc(Context *input_ctx, int input);
void fnACItrfcCoeffDC(Context *input_ctx, int input);
void fnACItrfcVoid(Context *input_ctx, int input);
void fnACmAlpItrfc(Context *input_ctx, int input_sym, int Alp);
void fnACItrfcByps(int input);
Context *fnGetCtxAZF(int iRow, int iCol, int iCurAZF);
Context *fnGetCtxExtF(int iRow, int iCol);
Context *fnGetCtxSigExt(int iAcumSig, int iElemtIdx, int iRow, int iCol, int iBlkRow, int iBlkCol, int iBitPos, int iCurSigExt);
Context *fnGetCtxSig(int iAcumSig, int iElemtIdx, int iRow, int iCol, int iBlkRow, int iBlkCol, int iBitPos, int iCurSig);
Context *fnGetCtxExt(int iElemtIdx, int iRow, int iCol, int iBlkRow, int iBlkCol, int iBitPos);
Context *fnGetCtxSign(int iExt, int iRow, int iCol, int iBlkRowIdx, int iBlkColIdx);
//void fnEncOneBlk(int iCtxDC, int *ipZZIdx, int *ipZZExt, double *dpZZCoef, int iRowIdx, int iColIdx);
void fnEncDC(int iOptCtxDC, int iIdxDC);
int fnHDQ(double dCoef, int iBlkRow, int iBlkCol, int iFreqRow, int iFreqCol);
void fnEncBlkBiLvImg(int *rgdQntSeq, int iBlkRow, int iBlkCol);
void fnEncOneBlk(int *rgdQntSeq, int iBlkRow, int iBlkCol);
void fnColecStatisElemt(int iQntIdx, int iElemtIdx, int iBlkRow, int iBlkCol, int iFreqRow, int iFreqCol);
void fnColecStatisSig(int iCurEXTF, int iBlkRow, int iBlkCol, int iCurBitLoc, int iElemtIdx, int iCurSig);
int fnCalCtxAZF(int iRow, int iCol);
void fnBlkSDQ(int iBlkRow, int iBlkCol, double** drgdOrigCoef);
void fnCalEntrpy();
void fnCtrIni();
int fnTrucCtx(int iBlkRow, int iBlkCol, int iFixCtx, int iBitLoc);
int fnCalFixCtxSig(int iBlkRow, int iBlkCol, int iCurBitLoc, int iElemtIdx);
void fnEndTrellis(int iBlkRow, int iBlkCol, double ** drgdOrigCoef, int iLyrElemtIdxCache[8], int iCurLyrElemtIdx, int *iCBP, StateSDQ* CurLyrTrellis);
void fnInitTrellis(int iBlkRow, int iBlkCol, int iFreqRow, int iFreqCol, int iCurElemtIdx, double dOrigCoef, StateSDQ* CurLyrTrellis);
void fnContTrellis(int iBlkRow, int iBlkCol, int iFreqRow, int iFreqCol, int iLyrElemtIdxCache[8], int iCurLyrElemtIdx, double dOrigCoef, StateSDQ* CurLyrTrellis);
void fnColecStatisBlk(std::vector<std::vector<int>>drgiQtzIdxHDQ, int iBlkRow, int iBlkCol);
void fnFreeMemHACACM();

#endif