//
//  FastDCT.h
//  NewEntropyCoding_Enc
//
//  Created by Chang Sun on 2014-01-07.
//  Copyright (c) 2014 University of Waterloo. All rights reserved.
//

#pragma once
#ifndef _FastDCT_H_
#define _FastDCT_H_

void FDCT(double d[8][8]);
void IDCT(double d[8][8]);

extern double _C[8];
#endif

