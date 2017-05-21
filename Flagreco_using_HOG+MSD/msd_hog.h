#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <math.h>
#include "cv.h"
#include "cxcore.h"
#include "highgui.h"
#include "svm.h" 

// RGB - the input original image data is in RGB color space

/*Feature Count*/
int count = 0;
/*the total color quantization  number of HSV color space*/
#define CSA 8  
/*the quantization number of edge orientation*/
#define CSB 6   
/*Defining the size of structure array*/ //40+72
#define MAX_INDEX 76
#define TOT_CLASSES 224
/*No. of bins (defined for use in multiple functions*/
#define BinCount 36     
/*Constant for radian to degree conversion*/ 
#define pi 3.14159265

typedef struct range_model
{
/**
 *	\var double *feature_min
 *	\brief stores minimum feature values
 */
	double *feature_min;
/**
 *	\var double *feature_max
 *	\brief stores maximum feature values
 */	
	double *feature_max;
/**
 *	\var double lower,upper
 *	\brief stores the lowest and highest values specified in range file
 */	
	double lower,upper;

}range_model;



void Mat_color_HSV(double*** HSV, int** img, int colnum1, int colnum2, int colnum3, int wid, int hei);

void mat_color(int*** RGB, double*** HSV, int wid, int hei);

void Mat_ori_HSV(int*** HSV,int** ori, int wid, int hei);

void Map(int** ori, int** img, int** Color, int wid, int hei, int Dx, int Dy);

void microstructure(int** ori, int** ImageX, int** micro, int wid, int hei);

void microdiscriptor(int** ColorX, double* hist, int wid, int hei);

void MSD_feature_extract(int*** HSV,int** ImageX, double* hist, int wid, int hei);

void HSV_color_quant(IplImage*, IplImage*);

void HOG_FeatureEx_FileCreation(IplImage*,char*,FILE*, int);

void AddMat_Bin(CvMat*, CvMat*, CvMat*);

void BlockBinTotal(CvMat*, CvMat*, int,int,int,int, double*);

double BlockBinTotal_Sum(double*);

//Replacement of 'main'
void main_rep(IplImage* IpImage_Color, int* predict_label_seq);

//Scaling function
void scale(double *feature,double *scaledfeature,struct range_model *r,int max_index);