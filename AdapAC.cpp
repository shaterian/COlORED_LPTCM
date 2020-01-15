/*
*  AdapAC.c
*  mjpeg
*
*  Created by Jin Meng on 11-06-08.
*  Copyright 2011 University of Waterloo. All rights reserved.
*
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "global.h"
#include "AdapAC.h"
//#include "NunfQuanz.h"

#define B_BITS 16

static unsigned short int half = (unsigned short int)(1 << (B_BITS - 1));
static unsigned short int quarter = (unsigned short int)(1 << (B_BITS - 2));
static unsigned short int out_L;   // lower end of the interval in encoder
static unsigned short int out_R;   // range of the interval in encoder
int out_bits_outstanding;

char rgcACDebugrgcIFName[] = "/Users/j4meng/jpeg_ac_debug";
/*
*  function: fnIntlAC
*  Initialize Arithmetic Code Engine
*/
void fnIntlAC(ArithmeticCode *ACpEngine)
{
	ACpEngine->m_llowerbound = 0;
	ACpEngine->m_lupperbound = (1 << PRECISION) - 1;
	ACpEngine->m_k = 0;
}

/*
*  function: fnACEncode
*  Encode a symbol based on a context
*/
void fnACEncode(ArithmeticCode *ACpEngine, Context *ctxpInst, int iSymbol, BitStream *bspInst)
{
 
	long lDeLow{ 0 }, lDeUp{0};
	int i{ 0 };
#if ACDEBUG	
	FILE *fp;
#endif	

     	lDeLow = ACpEngine->m_llowerbound +
		((ACpEngine->m_lupperbound - ACpEngine->m_llowerbound + 1)
			*ctxpInst->m_lpCumSum[iSymbol] / ctxpInst->m_lpCumSum[ctxpInst->m_iAlpSize]);

	lDeUp = ACpEngine->m_llowerbound +
		((ACpEngine->m_lupperbound - ACpEngine->m_llowerbound + 1)
			*ctxpInst->m_lpCumSum[iSymbol + 1] / ctxpInst->m_lpCumSum[ctxpInst->m_iAlpSize]) - 1;

	if (lDeLow > lDeUp)
	{
		printf("error!");
	}

	/* Update context */
	fnUpdateCtx(ctxpInst, iSymbol);

	/* Re-Normalization */
	while ((lDeUp<(1 << (PRECISION - 1)))
		|| (lDeLow >= (1 << (PRECISION - 1)))
		|| ((lDeLow >= (1 << (PRECISION - 2))) && (lDeUp<((1 << (PRECISION - 1)) + (1 << (PRECISION - 2))))))
	{
		if (lDeUp<(1 << (PRECISION - 1)))
		{
			fnWaBit(bspInst, 0);
			for (i = 0; i<ACpEngine->m_k; i++)
			{
				fnWaBit(bspInst, 1);
			}
			lDeUp <<= 1;
			lDeUp++;
			lDeLow <<= 1;
			ACpEngine->m_k = 0;
		}
		else
		{
			if (lDeLow >= (1 << (PRECISION - 1)))
			{
				fnWaBit(bspInst, 1);
				for (i = 0; i<ACpEngine->m_k; i++)
				{
					fnWaBit(bspInst, 0);
				}
				lDeUp = (lDeUp - (1 << (PRECISION - 1))) << 1;
				lDeUp++;
				lDeLow = (lDeLow - (1 << (PRECISION - 1))) << 1;
				ACpEngine->m_k = 0;
			}
			else
			{
				lDeUp = (lDeUp - (1 << (PRECISION - 2))) << 1;
				lDeUp++;
				lDeLow = (lDeLow - (1 << (PRECISION - 2))) << 1;
				ACpEngine->m_k++;
			}
		}
	}

	ACpEngine->m_llowerbound = lDeLow;
	ACpEngine->m_lupperbound = lDeUp;

#if ACDEBUG	
	fp = fopen(rgcACDebugrgcIFName, "a");
	fprintf(fp, "%ld : %ld\n", lDeLow, lDeUp);
	fclose(fp);
#endif	
}

/*
* This routine must be called to initialize the encoding process.
*/
void initialize_arithmetic_encoder()
{
	out_L = 0;
	out_R = half;
	out_bits_outstanding = 0;
}

/*
*  function: fnACEncode
*  Encode a symbol based on a context
*/
void fnBACEncode(ArithmeticCode *ACpEngine, Context *ctxpInst, int iSymbol, BitStream *bspInst)
{
	long c0, c1;        // number of cumulative occurances of bit 0 and 1, respectively
	long cLPS, rLPS;    // count of LPS and corresponding range within previous total range
	int LPS;                    // indicate current less probable symbol, 0 or 1

	c0 = ctxpInst->m_lpCumSum[0];
	c1 = ctxpInst->m_lpCumSum[1];

	if (c0<c1)
	{
		LPS = 0; cLPS = c0;
	}
	else
	{
		LPS = 1; cLPS = c1;
	}

	rLPS = (cLPS*out_R) / (c0 + c1);

	if (iSymbol == LPS)
	{
		out_L += out_R - rLPS;  out_R = rLPS;
	}
	else
		out_R -= rLPS;

	// normalized out_L and out_R
	while (out_R <= quarter)
	{
		if (out_L >= half)
		{
			fnWaBit(bspInst, 1);
			while (out_bits_outstanding > 0)
			{
				fnWaBit(bspInst, 0);
				out_bits_outstanding--;
			}
		}
		else if (out_L + out_R <= half)
		{
			fnWaBit(bspInst, 0);
			while (out_bits_outstanding > 0)
			{
				fnWaBit(bspInst, 1);
				out_bits_outstanding--;
			}
		}
		else
		{
			out_bits_outstanding++;
			out_L -= quarter;
		}
		out_L <<= 1;
		out_R <<= 1;
	}

	/* Update context */
	fnUpdateBACCtx(ctxpInst, iSymbol);
}

/*
*  function: fnACEncodeVoid
*  Encode a symbol based on a context
*/
void fnACEncodeVoid(ArithmeticCode *ACpEngine, Context *ctxpInst, int iSymbol, BitStream *bspInst)
{
	long lDeLow, lDeUp;
	int i;
#if ACDEBUG	
	FILE *fp;
#endif	

	lDeLow = ACpEngine->m_llowerbound +
		((ACpEngine->m_lupperbound - ACpEngine->m_llowerbound + 1)
			*ctxpInst->m_lpCumSum[iSymbol] / ctxpInst->m_lpCumSum[ctxpInst->m_iAlpSize]);

	lDeUp = ACpEngine->m_llowerbound +
		((ACpEngine->m_lupperbound - ACpEngine->m_llowerbound + 1)
			*ctxpInst->m_lpCumSum[iSymbol + 1] / ctxpInst->m_lpCumSum[ctxpInst->m_iAlpSize]) - 1;

	if (lDeLow > lDeUp)
	{
		printf("error!");
	}

	/* Update context */
	fnUpdateCtx(ctxpInst, iSymbol);

	/* Re-Normalization */
	while ((lDeUp<(1 << (PRECISION - 1)))
		|| (lDeLow >= (1 << (PRECISION - 1)))
		|| ((lDeLow >= (1 << (PRECISION - 2))) && (lDeUp<((1 << (PRECISION - 1)) + (1 << (PRECISION - 2))))))
	{
		if (lDeUp<(1 << (PRECISION - 1)))
		{/*
		 fnWaBit(bspInst, 0);
		 for (i=0; i<ACpEngine->m_k; i++)
		 {
		 fnWaBit(bspInst, 1);
		 }*/
			lDeUp <<= 1;
			lDeUp++;
			lDeLow <<= 1;
			ACpEngine->m_k = 0;
		}
		else
		{
			if (lDeLow >= (1 << (PRECISION - 1)))
			{
				/*fnWaBit(bspInst, 1);
				for (i=0; i<ACpEngine->m_k; i++)
				{
				fnWaBit(bspInst, 0);
				}*/
				lDeUp = (lDeUp - (1 << (PRECISION - 1))) << 1;
				lDeUp++;
				lDeLow = (lDeLow - (1 << (PRECISION - 1))) << 1;
				ACpEngine->m_k = 0;
			}
			else
			{
				lDeUp = (lDeUp - (1 << (PRECISION - 2))) << 1;
				lDeUp++;
				lDeLow = (lDeLow - (1 << (PRECISION - 2))) << 1;
				ACpEngine->m_k++;
			}
		}
	}

	ACpEngine->m_llowerbound = lDeLow;
	ACpEngine->m_lupperbound = lDeUp;

#if ACDEBUG	
	fp = fopen(rgcACDebugrgcIFName, "a");
	fprintf(fp, "%ld : %ld\n", lDeLow, lDeUp);
	fclose(fp);
#endif	
}

/*
*  function: fnACEncode
*  Encode a symbol based on a context with modified alphabet size
*/
void fnACEncodeMAlp(ArithmeticCode *ACpEngine, Context *ctxpInst, int iSymbol, BitStream *bspInst, int Alp)
{
	long lDeLow, lDeUp;
	int i;
#if ACDEBUG	
	FILE *fp;
#endif

	assert(iSymbol<Alp);

	lDeLow = ACpEngine->m_llowerbound +
		((ACpEngine->m_lupperbound - ACpEngine->m_llowerbound + 1)
			*ctxpInst->m_lpCumSum[iSymbol] / ctxpInst->m_lpCumSum[Alp]);

	lDeUp = ACpEngine->m_llowerbound +
		((ACpEngine->m_lupperbound - ACpEngine->m_llowerbound + 1)
			*ctxpInst->m_lpCumSum[iSymbol + 1] / ctxpInst->m_lpCumSum[Alp]) - 1;

	if (lDeLow > lDeUp)
	{
		printf("error!");
	}

	/* Update context */
	fnUpdateCtx(ctxpInst, iSymbol);

	/* Re-Normalization */
	while ((lDeUp<(1 << (PRECISION - 1)))
		|| (lDeLow >= (1 << (PRECISION - 1)))
		|| ((lDeLow >= (1 << (PRECISION - 2))) && (lDeUp<((1 << (PRECISION - 1)) + (1 << (PRECISION - 2))))))
	{
		if (lDeUp<(1 << (PRECISION - 1)))
		{
			fnWaBit(bspInst, 0);
			for (i = 0; i<ACpEngine->m_k; i++)
			{
				fnWaBit(bspInst, 1);
			}
			lDeUp <<= 1;
			lDeUp++;
			lDeLow <<= 1;
			ACpEngine->m_k = 0;
		}
		else
		{
			if (lDeLow >= (1 << (PRECISION - 1)))
			{
				fnWaBit(bspInst, 1);
				for (i = 0; i<ACpEngine->m_k; i++)
				{
					fnWaBit(bspInst, 0);
				}
				lDeUp = (lDeUp - (1 << (PRECISION - 1))) << 1;
				lDeUp++;
				lDeLow = (lDeLow - (1 << (PRECISION - 1))) << 1;
				ACpEngine->m_k = 0;
			}
			else
			{
				lDeUp = (lDeUp - (1 << (PRECISION - 2))) << 1;
				lDeUp++;
				lDeLow = (lDeLow - (1 << (PRECISION - 2))) << 1;
				ACpEngine->m_k++;
			}
		}
	}

	ACpEngine->m_llowerbound = lDeLow;
	ACpEngine->m_lupperbound = lDeUp;

#if ACDEBUG	
	fp = fopen(rgcACDebugrgcIFName, "a");
	fprintf(fp, "%ld : %ld\n", lDeLow, lDeUp);
	fclose(fp);
#endif	
}

/*
*  function: fnACEncodebyPass
*  Encode a symbol without context
*/
void fnACEncodebyPass(ArithmeticCode *ACpEngine, int iSymbol, BitStream *bspInst)
{
	long lDeLow, lDeUp;
	int i;
#if ACDEBUG	
	FILE *fp;
#endif

	if (iSymbol == 0)
	{
		lDeLow = ACpEngine->m_llowerbound;
		lDeUp = ACpEngine->m_llowerbound +
			((ACpEngine->m_lupperbound - ACpEngine->m_llowerbound + 1) / 2) - 1;
	}
	else
	{
		lDeLow = ACpEngine->m_llowerbound +
			((ACpEngine->m_lupperbound - ACpEngine->m_llowerbound + 1) / 2);
		lDeUp = ACpEngine->m_lupperbound;
	}

	while ((lDeUp<(1 << (PRECISION - 1)))
		|| (lDeLow >= (1 << (PRECISION - 1)))
		|| ((lDeLow >= (1 << (PRECISION - 2))) && (lDeUp<((1 << (PRECISION - 1)) + (1 << (PRECISION - 2))))))
	{
		if (lDeUp<(1 << (PRECISION - 1)))
		{
			fnWaBit(bspInst, 0);
			for (i = 0; i<ACpEngine->m_k; i++)
			{
				fnWaBit(bspInst, 1);
			}
			lDeUp <<= 1;
			lDeUp++;
			lDeLow <<= 1;
			ACpEngine->m_k = 0;
		}
		else
		{
			if (lDeLow >= (1 << (PRECISION - 1)))
			{
				fnWaBit(bspInst, 1);
				for (i = 0; i<ACpEngine->m_k; i++)
				{
					fnWaBit(bspInst, 0);
				}
				lDeUp = (lDeUp - (1 << (PRECISION - 1))) << 1;
				lDeUp++;
				lDeLow = (lDeLow - (1 << (PRECISION - 1))) << 1;
				ACpEngine->m_k = 0;
			}
			else
			{
				lDeUp = (lDeUp - (1 << (PRECISION - 2))) << 1;
				lDeUp++;
				lDeLow = (lDeLow - (1 << (PRECISION - 2))) << 1;
				ACpEngine->m_k++;
			}
		}
	}

	ACpEngine->m_llowerbound = lDeLow;
	ACpEngine->m_lupperbound = lDeUp;

#if ACDEBUG	
	fp = fopen(rgcACDebugrgcIFName, "a");
	fprintf(fp, "%ld : %ld\n", lDeLow, lDeUp);
	fclose(fp);
#endif
}

/*
*  function: fnEndEncode
*  End Arithmetic Coding Procedure
*/
void fnEndEncode(ArithmeticCode *ACpEngine, BitStream *bspInst)
{
	int i;

	if (ACpEngine->m_llowerbound <= (1 << (PRECISION - 2)))
	{
		fnWaBit(bspInst, 0);
		for (i = 0; i<ACpEngine->m_k + 1; i++)
		{
			fnWaBit(bspInst, 1);
		}
	}
	else
	{
		fnWaBit(bspInst, 1);
		for (i = 0; i<ACpEngine->m_k; i++)
		{
			fnWaBit(bspInst, 0);
		}
	}
}

/*
*  function: fnEncBlk
*  Encode a 8x8 block
*/
void fnEncBlk(ArithmeticCode *ACpEngine, int ippCoeffBlk[8][8], BitStream *bspInst)
{
	int i, j;
	int rgnNZLayer[8];
	int iLastLZLayer;
	int nNZLayer;
	int iExistNZ;

	for (i = 0; i<7; i++)
	{
		for (j = 0; j <= i; j++)
		{
			if ((i + j) == 0)
			{
				fnACEncode(ACpEngine, fnGetCoeffCtx(ippCoeffBlk, j, i - j, 0, 0), abs(ippCoeffBlk[j][i - j]), bspInst);
			}
			else
			{
				fnACEncode(ACpEngine, fnGetCoeffCtx(ippCoeffBlk, j, i - j, 0, 0), MINIMUM(abs(ippCoeffBlk[j][i - j]), ((drgnTrctLv[j][i - j] - 1) / 2)), bspInst);
			}

			if (ippCoeffBlk[j][i - j] != 0)
			{
				fnACEncodebyPass(ACpEngine, (ippCoeffBlk[j][i - j]>0) ? 1 : 0, bspInst);
			}
		}
	}

	for (i = 0; i<8; i++)
	{
		rgnNZLayer[i] = 0;
	}
	iLastLZLayer = 7;

	for (i = 7; i<14; i++)
	{
		for (j = i - 7; j <= 7; j++)
		{
			if (ippCoeffBlk[j][i - j] != 0)
			{
				rgnNZLayer[i - 7] = 1;
				iLastLZLayer = i + 1;
				break;
			}
		}
	}

	if (ippCoeffBlk[7][7] != 0)
	{
		iLastLZLayer = 14;
	}


	nNZLayer = fnCalNZLayer(ippCoeffBlk);
	for (i = 0; i <= (MINIMUM((iLastLZLayer - 7), 6)); i++)
	{
		fnACEncode(ACpEngine, &(drgulpNZCount[i][(int)ceil(log(nNZLayer + 1) / log(2))]), rgnNZLayer[i], bspInst);
		if (rgnNZLayer[i] == 0)
		{
			fnACEncode(ACpEngine, &(drgulpLZCount[i][(int)ceil(log(nNZLayer + 1) / log(2))]), (int)(i == (iLastLZLayer - 7)), bspInst);
		}
	}

	for (i = 7; i<14; i++)
	{
		if (rgnNZLayer[i - 7] == 1)
		{
			iExistNZ = 0;
			for (j = i - 7; j <= 6; j++)
			{
				fnACEncode(ACpEngine, fnGetCoeffCtx(ippCoeffBlk, j, i - j, 0, nNZLayer), (ippCoeffBlk[j][i - j] != 0), bspInst);
				if (ippCoeffBlk[j][i - j] != 0)
				{
					fnACEncode(ACpEngine, fnGetCoeffCtx(ippCoeffBlk, j, i - j, 1, nNZLayer), (abs(ippCoeffBlk[j][i - j])>1), bspInst);
					fnACEncodebyPass(ACpEngine, (ippCoeffBlk[j][i - j]>0) ? 1 : 0, bspInst);
					iExistNZ = 1;
				}
			}
			if (iExistNZ != 0)
			{
				fnACEncode(ACpEngine, fnGetCoeffCtx(ippCoeffBlk, 7, i - 7, 0, nNZLayer), (ippCoeffBlk[7][i - 7] != 0), bspInst);
			}
			if (ippCoeffBlk[7][i - 7] != 0)
			{
				fnACEncode(ACpEngine, fnGetCoeffCtx(ippCoeffBlk, 7, i - 7, 1, nNZLayer), (abs(ippCoeffBlk[7][i - 7])>1), bspInst);
				fnACEncodebyPass(ACpEngine, (ippCoeffBlk[7][i - 7]>0) ? 1 : 0, bspInst);
			}
		}
	}

	if (iLastLZLayer == 14)
	{
		fnACEncode(ACpEngine, fnGetCoeffCtx(ippCoeffBlk, 7, 7, 0, nNZLayer), (ippCoeffBlk[7][7] != 0), bspInst);
		if (ippCoeffBlk[7][7] != 0)
		{
			fnACEncode(ACpEngine, fnGetCoeffCtx(ippCoeffBlk, 7, 7, 1, nNZLayer), (abs(ippCoeffBlk[7][7])>1), bspInst);
			fnACEncodebyPass(ACpEngine, (ippCoeffBlk[7][7]>0) ? 1 : 0, bspInst);
		}
	}
}

/*
*  function: fnRiceGolombCode
*  Encode a symbol using iOrder-th Rice Golomb Code
*/
void fnRiceGolombCode(ArithmeticCode *ACpEngine, int iSymbol, int iOrder, BitStream *bspInst)
{
	int iQuotient, iReminder;
	int i;

	iQuotient = iSymbol >> iOrder;
	iReminder = iSymbol - (iQuotient << iOrder);

	/*
	*  Encode Quotient by unary code
	*/
	for (i = 0; i<iQuotient; i++)
	{
		fnACEncodebyPass(ACpEngine, 1, bspInst);
	}
	fnACEncodebyPass(ACpEngine, 0, bspInst);

	/*
	*  Encode Reminder by iOrder bits binary representation
	*/
	for (i = (iOrder - 1); i >= 0; i--)
	{
		fnACEncodebyPass(ACpEngine, (iReminder >> i) % 2, bspInst);
	}
}

/*
*  function: fnEncExIdx
*  Encode extension indieces within a 8x8 block
*/
void fnEncExIdx(ArithmeticCode *ACpEngine, int ippCoeffBlk[8][8], BitStream *bspInst)
{
	int i, j;
	int index;

	for (i = 0; i<8; i++)
	{
		for (j = 0; j<8; j++)
		{
			if ((abs(ippCoeffBlk[i][j]) >= ((drgnTrctLv[i][j] - 1) / 2)) && ((i + j) != 0))
			{
				index = i + j - 1;
				index = MINIMUM(index, 5);
				fnRiceGolombCode(ACpEngine, (abs(ippCoeffBlk[i][j]) - (drgnTrctLv[i][j] - 1) / 2), rgiRiceGolombOrder[index], bspInst);
			}
		}
	}
}