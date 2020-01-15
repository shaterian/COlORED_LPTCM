/*
*  BitStream.c
*  mjpeg
*
*  Created by Jin Meng on 11-06-08.
*  Copyright 2011 __MyCompanyName__. All rights reserved.
*
*/
#pragma warning(disable:4996)
#include <stdio.h>
#include <stdlib.h>
#include "BitStream.h"
#include "global.h"
#include "AdapRecSpace.h"
#include <math.h>
/*
*  function: fnIntlBitStream
*  Initialize an instance of BitStream
*/
void fnIntlBitStream(BitStream *bspInst, long lSize)
{
	if (bspInst->m_cpBuffer != NULL)
	{
		//free(bspInst->m_cpBuffer);
	}

	bspInst->m_cpBuffer = (char *)malloc(lSize * sizeof(char));
	bspInst->m_nCpty = lSize;
	bspInst->m_nBits = 0;
}

/*
*  function: fnWaBit
*  Write a bit into Bitstream buffer
*/
void fnWaBit(BitStream *bspInst, char cBit)
{
	if (bspInst->m_nBits < bspInst->m_nCpty)
	{
		bspInst->m_cpBuffer[bspInst->m_nBits] = cBit;
		bspInst->m_nBits++;
	}
	else
	{
		bspInst->m_cpBuffer = (char *)realloc(bspInst->m_cpBuffer, 2 * bspInst->m_nCpty * sizeof(char));
		bspInst->m_nCpty *= 2;
		bspInst->m_cpBuffer[bspInst->m_nBits] = cBit;
		bspInst->m_nBits++;
	}

}
// reading a bit 
char fnRaBit(BitStream *bspInst)
{
	if (bspInst->m_nBits >= bspInst->m_nCpty)
	{
		return 0;
	}
	else
	{
		bspInst->m_nBits++;
		return bspInst->m_cpBuffer[bspInst->m_nBits - 1];
	}

}

void fnWtoFileInBits(char *cIn, int iNumBits, char *cpOFName)
{
	int i, j, iTmp;
	long iNumBytes = (int)ceil(iNumBits / 8.0);
	int iMul[8] = { 128,64,32,16,8,4,2,1 };
	unsigned char *ucOut;
	FILE *fpOF;

	ucOut = (unsigned char *)malloc((iNumBytes) * sizeof(unsigned char));

	for (j = 0; j<(iNumBytes - 1); j++)
	{
		iTmp = 0;
		for (i = 0; i<8; i++)
			iTmp += iMul[i] * (cIn[j * 8 + i] & 1);

		ucOut[j] = (unsigned char)iTmp;
	}
	iTmp = 0;
	for (i = 0; i<8; i++)
	{
		if ((j * 8 + i) < iNumBits)
			iTmp += iMul[i] * (cIn[j * 8 + i] & 1);
	}
	ucOut[j] = (unsigned char)iTmp;

	fpOF = fopen(cpOFName, "wb");
	fwrite(ucOut, sizeof(unsigned char), iNumBytes, fpOF);
	fclose(fpOF);

	free(ucOut);
}

void fnWtoFile(BitStream *bspInst, char *cpOFName)
{
	int i, j, iTmp, iNumBits, iNumBytes;
	int iMul[8] = { 128,64,32,16,8,4,2,1 };
	unsigned char *ucOut;
	char cTtbl[64];
	FILE* fpOF;

	iNumBits = (int)bspInst->m_nBits;
	iNumBytes = (int)ceil(iNumBits / 8.0);
	fpOF = fopen(cpOFName, "wb");
	fprintf(fpOF, "%d ", iNumBits);
	fprintf(fpOF, "%hu %hu ", (unsigned short)iLLowDC, (unsigned short)iLHighDC);
	fprintf(fpOF, "%hu ", (unsigned short)(TAR_PSNR_ARG));
	//fprintf(fpOF, "%hu ", (unsigned short)fnRound(TAR_PSNR_ARG / 100.0));
	//fprintf(fpOF, "%d %d %d %d ", iNumBits, iLLowDC, iLHighDC, TAR_PSNR_ARG);
	iTmp = 0;
	unsigned char Overhead[192];
	int idx = 0;
	for (i = 0; i<8; i++)
	{
		for (j = 0; j<8; j++)
		{
			cTtbl[iTmp++] = drgiTTbl[i][j];
			if (drgiTTbl[i][j] == 1)
			{
				//Overhead[idx++] = (unsigned char)drgiLTbl[i][j];
				Overhead[idx++] = (unsigned char)drgiLPrimeTbl[i][j];
				Overhead[idx++] = (unsigned char)OptDzoneQtzr.Lambda_idx[i][j];
				Overhead[idx++] = (unsigned char)OptDzoneQtzr.Gamma_idx[i][j];
			}

		}
	}

	ucOut = (unsigned char *)malloc(8 * sizeof(unsigned char));

	for (j = 0; j<8; j++)
	{
		iTmp = 0;
		for (i = 0; i<8; i++)
			iTmp += iMul[i] * (cTtbl[j * 8 + i] & 1);

		ucOut[j] = (unsigned char)iTmp;
	}
	fwrite(ucOut, sizeof(unsigned char), 8, fpOF);
	free(ucOut);
	fwrite(Overhead, sizeof(unsigned char), idx, fpOF);

	ucOut = (unsigned char *)malloc((iNumBytes) * sizeof(unsigned char));

	for (j = 0; j<(iNumBytes - 1); j++)
	{
		iTmp = 0;
		for (i = 0; i<8; i++)
			iTmp += iMul[i] * (bspInst->m_cpBuffer[j * 8 + i] & 1);

		ucOut[j] = (unsigned char)iTmp;
	}
	iTmp = 0;
	for (i = 0; i<8; i++)
	{
		if ((j * 8 + i) < iNumBits)
			iTmp += iMul[i] * (bspInst->m_cpBuffer[j * 8 + i] & 1);
	}
	ucOut[j] = (unsigned char)iTmp;

	fwrite(ucOut, sizeof(unsigned char), iNumBytes, fpOF);
	fclose(fpOF);

	free(ucOut);
}

void fnRfrFile(BitStream *bspInst, char *cpIFName)
{
	FILE* fpIF;
	long lFileSize;
	fpIF = fopen(cpIFName, "rb");

	fseek(fpIF, 0, SEEK_END);
	lFileSize = ftell(fpIF);
	fseek(fpIF, 0, SEEK_SET);

	fnIntlBitStream(bspInst, lFileSize);
	fread(bspInst->m_cpBuffer, sizeof(char), lFileSize, fpIF);

	fclose(fpIF);
}

void fnFreeBitStream(BitStream *bspInst)
{
	free(bspInst->m_cpBuffer);
}

void fnDebugBitStream(BitStream *bspInst, char *cpDFName)
{
	int i;
	FILE* fp;

	fp = fopen(cpDFName, "w");
	for (i = 0; i<bspInst->m_nBits; i++)
	{
		fprintf(fp, "%d", (int)(bspInst->m_cpBuffer[i]));
	}
	fclose(fp);
}
