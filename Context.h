/*
*  Context.h
*  mjpeg
*
*  Created by Jin Meng on 11-06-08.
*  Copyright 2011 Multimeida Communication Lab, University of Waterloo.
*  All rights reserved.
*/
#pragma once
#ifndef _CONTEXT_H_
#define _CONTEXT_H_

#define NZCTX 5

#define LZCTX 5

/*
*  structure: Context
*  Container storing information of a context
*/
typedef struct{
	int m_iAlpSize;
	long *m_lpCumSum;
} Context;

/*
*  variable: rgnulpCoeffCtx
*  array to store information of each context for coefficeints (abs val) in each layer
*/
extern Context *rgctxpCoeffLayerCtx[15];

/*
*  variable: drgulpNZCount, drgulpLZCount
*  array to store information of each context for nz (none-zero), lz of each layer
*/
extern Context drgulpNZCount[7][NZCTX];
extern Context drgulpLZCount[7][LZCTX];


void fnUpdateCtx(Context *ctxpInst, int iSymbol);
void fnUpdateBACCtx(Context *ctxpInst, int iSymbol);
void fnCreateCtx(Context *ctxpInst, int iAlpSize);
void fnCreateIniCtxEXTF(Context *ctxpInst, int CtxDstrMatrx[2][16], int iAlpSize, int iCtxIdx);
void fnCreateIniCtxExtImg(Context *ctxpInst, int CtxDstrMatrx[6][64], int iAlpSize, int iCtxIdx, int iCoefIdx);
void fnFreeCtx(Context *ctxpInst);

//void fnIntlCount(char *cpCtxFilePath);
void fnIntlCount(char *cpCtxFilePath, int iMaxDCVal);
void fnFreeCount();

int fnCalNZLayer(int ippCoeffblk[8][8]);

extern Context *fnGetCoeffCtx(int ippCoeffBlk[8][8], int iRow, int iCol, int isVal, int nNZLayer);

#endif
