/*
*  BitStream.h
*  mjpeg
*
*  Created by Jin Meng on 11-06-08.
*  Copyright 2011 Multimeida Communication Lab, University of Waterloo.
*  All rights reserved.
*/

#pragma once
#ifndef  _BITSTREAM_
#define	_BITSTREAM_

#define IntlBfSz 1<<20

/*
*  structure: BitStream
*  the structure for bit buffer
*/
typedef struct {
	long m_nBits;
	long m_nCpty;
	char *m_cpBuffer;
} BitStream;
/*
*  function: fnIntlBitStream
*  Initialize an instance of BitStream
*/
void fnIntlBitStream(BitStream *bspInst, long lSize);
void fnWaBit(BitStream *bspInst, char cBit);
char fnRaBit(BitStream *bspInst);
void fnWtoFile(BitStream *bspInst, char *cpOFName);
void fnRfrFile(BitStream *bspInst, char *cpIFName);
void fnFreeBitStream(BitStream *bspInst);

void fnDebugBitStream(BitStream *bspInst, char *cpDFName);

#endif