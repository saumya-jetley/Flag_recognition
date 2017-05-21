#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <math.h>
#include "cv.h"
#include "cxcore.h"
#include "highgui.h"
// RGB is the input original image datas in RGB color space

int CSA = 8;  // the total color quantization  number of HSV color space
int CSB = 6;   //the quantization number of edge orientation

//--------------------------------------------
//int Cn1 = 8;   //the quantization number of H 
//int Cn2 = 3;   //the quantization number of S 
//int Cn3 = 3;   //the quantization number of V

/*Images unable to load*/
int failcount = 0;
/*Feature Count*/
int count = 0;
/*No. of bins (defined for use in multiple functions*/

int flag = 0;

#define BinCount 72     
/*Constant for radian to degree conversion*/ 
#define pi 3.14159265

void Mat_color_HSV(double*** HSV, int** img, int colnum1, int colnum2, int colnum3, int wid, int hei);

void mat_color(int*** RGB, double*** HSV, int wid, int hei);

void Mat_ori_HSV(int*** HSV,int** ori, int wid, int hei);

void Map(int** ori, int** img, int** Color, int wid, int hei, int Dx, int Dy);

void microstructure(int** ori, int** ImageX, int** micro, int wid, int hei);

void microdiscriptor(int** ColorX, double* hist, int wid, int hei,int flag);

void MSD_feature_extract(int*** HSV,int** ImageX, double* hist, int wid, int hei,int flag);

void HSV_color_quant(IplImage*, IplImage*);

void HOG_FeatureEx_FileCreation(IplImage*,char*,FILE*, int);

void AddMat_Bin(CvMat*, CvMat*, CvMat*);

void BlockBinTotal(CvMat*, CvMat*, int,int,int,int, int*);

int BlockBinTotal_Sum(int*);
