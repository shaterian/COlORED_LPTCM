// this class performs the TCM algorithm in order to find the Yk and lamda
//  IntegerTCM.c
//  NewEntropyCoding_Enc
//
//  Created by Chang Sun on 2014-01-08.
//  Copyright (c) 2014 University of Waterloo. All rights reserved.
//

#include "IntegerTCM.h"
#include <iostream>

// it takes in the DCTcoeffs of whole image 
// for each frequency it takes the sequence of that coeffs for that freq for all 
// blocks of the image then it computes the values of a , b , Yc , lambda for each 
// freq as Aval, Bval, drgdYc , drgdLambda and so send this values since this value
// are constant for the whole image so aonly 64 values would be snet 
void fnTCM(double** drgdCoeff)
{
	int* rgiSeq = new int [NUMBLK];
	int peak, i, j, u, v;
	double prob, Yc, lambda;
	FILE *fp;

	for (i = 0; i<8; i++)
	{
		for (j = 0; j<8; j++)
		{
			//            if (i==0 && j==1)
			//                fp = fopen("Lena_AC2.dat","w");
			//            else if (i==1 && j==0)
			//                fp = fopen("Lena_AC3.dat","w");
			//
			int iIdx = 0;
			for (u = 0; u<ROWS / 8; u++)
			{
				for (v = 0; v<COLS / 8; v++)
				{
					rgiSeq[iIdx++] = abs((int)drgdCoeff[8 * u + i][8 * v + j]);

					//                    if (i+j==1)
					//                        fprintf(fp, "%4.4f ",drgdCoeff[8*u+i][8*v+j]);
				}
			}

			//            if (i+j==1)
			//                fclose(fp);
			//printf("%d %d\n", rgiSeq[NUMBLK-1],rgiSeq[NUMBLK-2]);

			TCMprocessOneSequence(rgiSeq, NUMBLK, &peak, &prob, &lambda, &Yc);
			//std::cout << "done tcm" << std::endl; 
			drgdAval[i][j] = (double)peak + 1;
			drgdYc[i][j] = (double)Yc;
			drgdLambda[i][j] = lambda;
			drgdBval[i][j] = prob;
		}
	}
	/*printf("---------------Yc--------------\n");
	for (i = 0; i<8; i++)
	{
		for (j = 0; j<8; j++)
		{
			printf("%6.2f ", drgdYc[i][j]);
		}
		printf("\n");
	}
*/
	/*printf("---------------Lambda--------------\n");

	for (i = 0; i<8; i++)
	{
		for (j = 0; j<8; j++)
		{
			printf("%6.2f ", drgdLambda[i][j]);
		}
		printf("\n");
	}
*/
	//printf( "---------------Lap_Lambda--------------\n");

	//    for (i=0; i<8; i++)
	//    {
	//        for (j=0; j<8; j++)
	//        {
	//            printf( "%6.2f ",drgdLambdaLap[i][j]);
	//        }
	//        printf( "\n");
	//    }

	/*printf("---------------Gamma--------------\n");

	for (i = 0; i<8; i++)
	{
		for (j = 0; j<8; j++)
		{
			printf("%6.2f ", drgdYc[i][j] / drgdLambda[i][j]);
		}
		printf("\n");
	}
*/
	/*printf("---------------A--------------\n");

	for (i = 0; i<8; i++)
	{
		for (j = 0; j<8; j++)
		{
			printf("%6.2f ", drgdAval[i][j]);
		}
		printf("\n");
	}
*/
	//printf("---------------B--------------\n");

	//for (i = 0; i<8; i++)
	//{
	//	for (j = 0; j<8; j++)
	//	{
	//		printf("%6.2f ", drgdBval[i][j]);
	//	}
	//	printf("\n");
	//}

	//printf("\n");
	delete[] rgiSeq;
}

void TCMprocessOneSequence(int *C, int Len, int *peak, double*prob, double *lambda, double*Yc)
{
	int k, MaxPos, StartPoint;
	double MaxLikelyhood, likelyhood;
	tBucket buck[MaxAmp];

	*peak = 0;
	for (k = 0; k<Len; k++)
	{
		int absv = abs(C[k]);
		if (absv>(*peak)) *peak = absv;
	}

	if (*peak == 0 || *peak >= MaxAmp)
	{
		//printf("Input data either all zero or exceed %d (%d)\n", MaxAmp, *peak);

		*lambda = -1;
		*Yc = 0;
		*prob = 1;
		return;
	}

	for (k = 0; k <= (*peak); k++)
	{
		buck[k].count = 0;
		buck[k].AcumAbsAmp = 0;
		buck[k].AcumSampNum = 0;
	}

	for (k = 0; k<Len; k++)
	{
		int cabs = abs(C[k]);
		buck[cabs].count++;
	}
	buck[0].AcumSampNum = buck[0].count;

	for (k = 1; k <= (*peak); k++)
	{
		buck[k].AcumAbsAmp = buck[k - 1].AcumAbsAmp + k * buck[k].count;
		buck[k].AcumSampNum = buck[k - 1].AcumSampNum + buck[k].count;
	}

	// saftey check
	//    if(buck[0].count < (buck[1].count>>1) || buck[1].count < (buck[2].count>>1))
	//    {
	//        *lambda = -1 ;
	//        *Yc = 0;
	//        *prob = 1 ;
	//        return ;
	//    }

	StartPoint = FindStartPoint(buck, *peak, Len);

	//if( StartPoint<=0 )
	//{
	//    *prob = 1 ;
	//    *lambda = -1 ;
	//    *Yc = 0 ;
	//    printf("Data distribution is not like Laplacian1\n");
	//    return ;
	//}

	ComputeLikelyhood(StartPoint, Len, buck, *peak);

	MaxLikelyhood = buck[StartPoint].likelyhood;

	MaxPos = StartPoint;
	for (k = StartPoint + 1; k <= (*peak); k++)
	{
		if (buck[k].count == 0) continue;

		ComputeLikelyhood(k, Len, buck, *peak);

		likelyhood = buck[k].likelyhood;

		if (likelyhood > MaxLikelyhood)
		{
			MaxPos = k;
			MaxLikelyhood = likelyhood;
		}
	}

	if (MaxLikelyhood > MinLikelyhood)
	{
		*prob = buck[MaxPos].prob;
		*lambda = buck[MaxPos].lambda;
		*Yc = MaxPos;
	}
	else  // the search failed
	{
		*prob = 1;
		*lambda = -1;
		*Yc = 0;
		// QQQQ printf("Data distribution is not like Laplacian2\n");
	}
}
// this function finds the value m in the paper which is the first i that is 
// greater than the d 
int FindStartPoint(tBucket *buck, int peak, int N)
{
	int k;
	for (k = peak; k>0; k--)
	{
		if (buck[k].count == 0) continue;

		if (buck[k].AcumSampNum < N* (1.0 - StartPointProb)) break;
	}

	if (buck[0].count>N / 100 && buck[1].count>N / 100 && buck[2].count>N / 100 && buck[3].count>N / 100)
	{
		if (k<3) k = 3;
	}
	else if (buck[0].count>N / 100 && buck[1].count>N / 100 && buck[2].count>N / 100)
	{
		if (k<2) k = 2;
	}
	else
	{
		if (k<1) k = 1;
	}
	return k;
}

// for each yc(which is referd by thePoint here) computes the likelyhood  
// and put lamda and prob and liklyhood in the buk
// of that point (when we are computing the liklyhood for different possile yc)
// buck is the list of each possibel point for yc with their attributs that are 
// the number of values befor and after them and likelyhood and b and lambda of each
void  ComputeLikelyhood(int thePoint, int N, tBucket *buck, int peak)
{

	double N1 = buck[thePoint].AcumSampNum;

	double N2 = N - N1;
	double Yc = thePoint;

	double sumYi = buck[thePoint].AcumAbsAmp;
	double totalNumYi = buck[thePoint].AcumSampNum;

	double lambda = ComputeLambdaGivenYc(Yc, sumYi, totalNumYi);

	double prob = (double)N1 / (double)N; // b 
	// peak is the a in the paper 
	if (lambda>0)
	{
		buck[thePoint].likelyhood = N2 * log(1 - prob) + N1 * log(prob) - N2 * log((peak - Yc)*2.0)
			- N1 * log(1 - exp(-Yc / lambda))
			- N1 * log(2 * lambda) - sumYi / lambda;

		buck[thePoint].lambda = lambda;
		buck[thePoint].prob = prob;
	}
	else
	{
		buck[thePoint].likelyhood = -MinLikelyhood;

		buck[thePoint].lambda = lambda;
		buck[thePoint].prob = 1;
	}
}

double ComputeLambdaGivenYc(double Yc, double sumYi, double totalNumYi)
{
	double c = sumYi / totalNumYi;
	double lambda, lambda_old;

	if (c / Yc >= 0.95) return -1.0; // Too much away from Laplacian

	lambda_old = c;

	lambda = c - Yc * (1.0 - 1.0 / (1.0 - exp(-Yc / lambda_old)));

	{ int k;
	for (k = 0; k<5; k++)
	{
		lambda_old = lambda;

		lambda = c - Yc * (1.0 - 1.0 / (1.0 - exp(-Yc / lambda_old)));
	}
	}
	while (fabs(lambda - lambda_old) > LambdaDelta)
	{
		lambda_old = lambda;

		lambda = c - Yc * (1.0 - 1.0 / (1.0 - exp(-Yc / lambda_old)));
	}

	return lambda;

}
