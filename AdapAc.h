/*
*  AdapAC.h
*  mjpeg
*
*  Created by Jin Meng on 11-06-08.
*  Copyright 2011 Multimeida Communication Lab, University of Waterloo.
*  All rights reserved.
*/
#pragma once
#include "Context.h"
#include "BitStream.h"

/*
*  structure: ArithmeticCode
*  the structure for arithmetic code engine
*/

// variables in context


typedef struct {
	long m_llowerbound;
	long m_lupperbound;
	int m_k;
} ArithmeticCode;

void initialize_arithmetic_encoder();
void fnBACEncode(ArithmeticCode *
	, Context *ctxpInst, int iSymbol, BitStream *bspInst);
void fnIntlAC(ArithmeticCode *ACpEngine);
void fnACEncode(ArithmeticCode *ACpEngine, Context *ctxpInst, int iSymbol, BitStream *bspInst);
void fnACEncodeVoid(ArithmeticCode *ACpEngine, Context *ctxpInst, int iSymbol, BitStream *bspInst);

void fnACEncodeMAlp(ArithmeticCode *ACpEngine, Context *ctxpInst, int iSymbol, BitStream *bspInst, int Alp);

void fnACEncodebyPass(ArithmeticCode *ACpEngine, int iSymbol, BitStream *bspInst);
void fnEndEncode(ArithmeticCode *ACpEngine, BitStream *bspInst);

void fnEncBlk(ArithmeticCode *ACpEngine, int ippCoeffBlk[8][8], BitStream *bspInst);

void fnRiceGolombCode(ArithmeticCode *ACpEngine, int iSymbol, int iOrder, BitStream *bspInst);
void fnEncExIdx(ArithmeticCode *ACpEngine, int ippCoeffBlk[8][8], BitStream *bspInst);




