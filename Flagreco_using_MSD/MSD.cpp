#include "MSD.h"
#include <windows.h>
#include <dirent.h>

//Minimum of two values
int fmin(int a,int b)
{
	if(a < b)
		return a;
	else
		return b;
}

//Maximum of two values
int fmax(int a,int b)
{
	if( a > b)
		return a;
	else
		return b;
}

//*********************************************************************************//
//-----------------------For Color/Texture Histogram(MSD)--------------------------//
//*********************************************************************************//
//RGB to HSV Conversion
//void mat_color(int*** RGB,double*** HSV, int wid, int hei)
/*{
	int i, j;
    //HSV =(double***)calloc(3,sizeof(double**));
	for(i = 0; i < 3;i++)
	{
		HSV[i] = (double**)calloc(wid,sizeof(double*));
		for(j = 0; j < wid;j++)
		{
			HSV[i][j] = (double*)calloc(hei,sizeof(double));
		}
	}
    for (i = 0; i < wid; i++)
    {
        for (j = 0; j < hei; j++)
        {
            double cMax = 255.0;

            int Max, Min, temp;
			
            Max = fmax(RGB[0][i][j], fmax(RGB[1][i][j], RGB[2][i][j]));
            Min = fmin(RGB[0][i][j], fmin(RGB[1][i][j], RGB[2][i][j]));

            temp = Max - Min;

            // V component

            HSV[2][i][j] = Max * 1.0 / cMax;

            //S component
            if (Max > 0)
            {
                HSV[1][i][j] = temp * 1.0 / Max;
            }
            else
            {
                HSV[1][i][j] = 0.0;
            }
            //H component
            if (temp > 0)
            {
                double rr = (Max - RGB[0][i][j]) * 1.0 / temp * 1.0;
                double gg = (Max - RGB[1][i][j]) * 1.0 / temp * 1.0;
                double bb = (Max - RGB[2][i][j]) * 1.0 / temp * 1.0;
                double hh = 0.0;

                if (RGB[0][i][j] == Max)
                {
                    hh = bb - gg;
                }
                else if (RGB[1][i][j] == Max)
                {
                    hh = rr - bb + 2.0;
                }
                else
                {
                    hh = gg - rr + 4.0;
                }
                if (hh < 0)
                {
                    hh = hh + 6;
                }
                HSV[0][i][j] = hh / 6;

            }

            HSV[0][i][j] *= 360.0;

        }
    }
}*/

//Finding the Orientation matrix from the HSV Input
void Mat_ori_HSV(int*** HSV,int** ori, int wid, int hei)
{
	int i,j;
	double gxx = 0.0, gyy = 0.0, gxy = 0.0;

    double rh = 0.0, gh = 0.0, bh = 0.0;
    double rv = 0.0, gv = 0.0, bv = 0.0;

    double theta = 0.0;
    double*** hsv = (double***)calloc(3,sizeof(double**));
	for(i = 0; i < 3;i++)
	{
		hsv[i] = (double**)calloc(wid,sizeof(double*));
		for(j = 0 ;j < wid;j++)
		{
			hsv[i][j] = (double*)calloc(hei,sizeof(double));
		}
	}
	/*ori = (int**)calloc(wid,sizeof(int*));
	for(i = 0;i < wid;i++)
	{
		ori[i] = (int*)calloc(hei,sizeof(int));
	}*/

    for (int i = 0; i < wid; i++)  // HSV based on cylinder 
    {
        for (int j = 0; j < hei; j++)
        {
            hsv[0][i][j] = HSV[1][i][j] * cos((double)HSV[0][i][j]);
            hsv[1][i][j] = HSV[1][i][j] * sin((double)HSV[0][i][j]);
            hsv[2][i][j] = HSV[2][i][j];
        }
    }


    for (int i = 1; i <= wid - 2; i++)
    {
        for (int j = 1; j <= hei - 2; j++)
        {

            //--------------------------------------
            rh = (double)(hsv[0][i - 1][j + 1] + 2 * hsv[0][i][j + 1] + hsv[0][i + 1][j + 1]) - (hsv[0][i - 1][j - 1] + 2 * hsv[0][i][j - 1] + hsv[0][i + 1][j - 1]);
            //gh = (double)(hsv[1, i - 1, j + 1] + 2 * hsv[1, i, j + 1] + hsv[1, i + 1, j + 1]) - (hsv[1, i - 1, j - 1] + 2 * hsv[1, i, j - 1] + hsv[1, i + 1, j - 1]);
            gh = (double)(hsv[1][i - 1][j + 1] + 2 * hsv[1][i][j + 1] + hsv[1][i + 1][j + 1]) - (hsv[1][i - 1][j - 1] + 2 * hsv[1][i][j - 1] + hsv[1][i + 1][j - 1]);
			//bh = (double)(hsv[2, i - 1, j + 1] + 2 * hsv[2, i, j + 1] + hsv[2, i + 1, j + 1]) - (hsv[2, i - 1, j - 1] + 2 * hsv[2, i, j - 1] + hsv[2, i + 1, j - 1]);
			bh = (double)(hsv[0][i - 1][j + 1] + 2 * hsv[2][i][j + 1] + hsv[2][i + 1][j + 1]) - (hsv[2][i - 1][j - 1] + 2 * hsv[2][i][j - 1] + hsv[2][i + 1][j - 1]);
            
			//-----------------------------------------
            rv = (double)(hsv[0][i + 1][j - 1] + 2 * hsv[0][i + 1][j] + hsv[0][i + 1][j + 1]) - (hsv[0][i - 1][j - 1] + 2 * hsv[0][i - 1][j] + hsv[0][i - 1][j + 1]);
            gv = (double)(hsv[1][i + 1][j - 1] + 2 * hsv[1][i + 1][j] + hsv[1][i + 1][j + 1]) - (hsv[1][i - 1][j - 1] + 2 * hsv[1][i - 1][j] + hsv[1][i - 1][j + 1]);
            bv = (double)(hsv[2][i + 1][j - 1] + 2 * hsv[2][i + 1][j] + hsv[2][i + 1][j + 1]) - (hsv[2][i - 1][j - 1] + 2 * hsv[2][i - 1][j] + hsv[2][i - 1][j + 1]);

            //---------------------------------------
            gxx = sqrt(rh * rh + gh * gh + bh * bh);
            gyy = sqrt(rv * rv + gv * gv + bv * bv);
            gxy = rh * rv + gh * gv + bh * bv;

            theta = (acos(gxy / (gxx * gyy + 0.0001)) * 180.0 / 3.14);

            ori[i][j] = theta * CSB / 180.0;

            if (ori[i][j] >= CSB - 1) ori[i][j] = CSB - 1;

        }
    }

	for (i=0; i < 3; ++i) 
	{
        if (hsv[i] != NULL) 
		{
            for (j=0; j < wid; ++j)
                free(hsv[i][j]);
            free(hsv[i]);
        }
    }
    free(hsv);
	hsv = NULL;
        
}

//Obtaining the Quantized HSV Matrix from the HSV input -- NOT IN USE
//void Mat_color_HSV(double*** HSV,int** img, int colnum1, int colnum2, int colnum3, int wid, int hei)
/*{
	int i,j;
	int VI = 0;
	int SI = 0;
	int HI = 0;
	img = (int**)calloc(wid,sizeof(int*));
	for(i = 0; i < wid;i++)
	{
		img[i] = (int*)calloc(hei,sizeof(int));
	}

	for (i = 0; i < wid; i++)
    {
        for (j = 0; j < hei; j++)
        {
			HI = (int)(HSV[0][i][j] * (colnum1 / 360.0));
			if (HI >= colnum1 - 1)
			{
				HI = colnum1 - 1;
			}
            //-------------------------------------

			SI = (int)(HSV[1][i][j] * (colnum2 / 1.0));
            if (SI >= colnum2 - 1)
            {
                SI = colnum2 - 1;
            }
            // -------------------------------------------

			VI = (int)(HSV[2][i][j] * (colnum3 / 1.0));
            if (VI >= colnum3 - 1)
            {
                VI = colnum3 - 1;
            }

           
            //-------------------------------------------
            img[i][j] = (colnum3 * colnum2) * HI + colnum3 * SI + VI;

        }
    }
}*/

//Give the color patch image in "Color" (color map) based on texture similarity from ORI input - pattern analysis from Dx, Dy
void Map(int** ori, int** img, int** Color, int wid, int hei, int Dx, int Dy)
{
	int i,j;
	/*Color = (int**)calloc(wid,sizeof(int*));
	for(i = 0;i < wid;i++)
	{
		Color[i] = (int*)calloc(hei,sizeof(int));
	}*/

    for (int i = 1; i < wid / 3 ; i++)
    {
        for (int j = 1; j < hei / 3 ; j++)
        {
            int WA[9];
            //===========================
            int m = 3 * i + Dx;
            int n = 3 * j + Dy;

            WA[0] = ori[m - 1][n - 1];
            WA[1] = ori[m - 1][n];
            WA[2] = ori[m - 1][n + 1];

            WA[3] = ori[m + 1][n - 1];
            WA[4] = ori[m + 1][n];
            WA[5] = ori[m + 1][n + 1];

            WA[6] = ori[m][n - 1];
            WA[7] = ori[m][n + 1];
            WA[8] = ori[m][n];

            //-------------------------
            if (WA[8] == WA[0])
            {
                Color[m - 1][n - 1] = img[m - 1][n - 1];
            }
            else
            {
                Color[m - 1][n - 1] = -1;
            }
            //--------------------
            if (WA[8] == WA[1])
            {
                Color[m - 1][n] = img[m - 1][n];
            }
            else
            {
                Color[m - 1][n] = -1;
            }
            //----------------------
            if (WA[8] == WA[2])
            {
                Color[m - 1][n + 1] = img[m - 1][n + 1];
            }
            else
            {
                Color[m - 1][n + 1] = -1;
            }
            //----------------------
            if (WA[8] == WA[3])
            {
                Color[m + 1][n - 1] = img[m + 1][n - 1];
            }
            else
            {
                Color[m + 1][n - 1] = -1;

            }
            //-------------------------
            if (WA[8] == WA[4])
            {
                Color[m + 1][n] = img[m + 1][n];
            }
            else
            {
                Color[m + 1][n] = -1;
            }
            //--------------------------
            if (WA[8] == WA[5])
            {
                Color[m + 1][n + 1] = img[m + 1][n + 1];
            }
            else
            {
                Color[m + 1][n + 1] = -1;
            }
            //-----------------------------------------
            if (WA[8] == WA[6])
            {

                Color[m][n - 1] = img[m][n - 1];
            }
            else
            {
                Color[m][n - 1] = -1;
            }
            //----------------------------------------
            if (WA[8] == WA[7])
            {
                Color[m][n + 1] = img[m][n + 1];
            }
            else
            {
                Color[m][n + 1] = -1;
            }
            //------------------------------------------
            if (WA[8] == WA[8]) Color[m][n] = img[m][n];
        }
    }
}

//Consolidate all the color maps together
void microstructure(int** ori, int** ImageX, int** micro, int wid, int hei)
{
	int** ColorA = (int**)calloc(wid,sizeof(int*));
    int** ColorB = (int**)calloc(wid,sizeof(int*));
    int** ColorC = (int**)calloc(wid,sizeof(int*));
    int** ColorD = (int**)calloc(wid,sizeof(int*));
	int i,j;
	
	for(i = 0; i < wid;i++)
	{
		ColorA[i] = (int*)calloc(hei,sizeof(int));		
		ColorB[i] = (int*)calloc(hei,sizeof(int));
		ColorC[i] = (int*)calloc(hei,sizeof(int));
		ColorD[i] = (int*)calloc(hei,sizeof(int));
	}

    Map(ori, ImageX, ColorA, wid, hei, 0, 0);
    Map(ori, ImageX, ColorB, wid, hei, 0, 1);
    Map(ori, ImageX, ColorC, wid, hei, 1, 0);
    Map(ori, ImageX, ColorD, wid, hei, 1, 1);

    //=========the final micro-structure map===============
   /* micro = (int**)calloc(wid,sizeof(int*));
	for(i = 0; i < wid;i++)
	{
		micro[i] = (int*)calloc(hei,sizeof(int));	
	}*/

    for (i = 0; i < wid; i++)
    {
        for (j = 0; j < hei; j++)
        {
            micro[i][j] = fmax(ColorA[i][j], fmax(ColorB[i][j], fmax(ColorC[i][j], ColorD[i][j])));
        }
    }
	
	for (i = 0; i < wid;i++)
	{
		free(ColorA[i]);
		free(ColorB[i]);
		free(ColorC[i]);
		free(ColorD[i]);
	}
	free(ColorA);
	free(ColorB);
	free(ColorC);
	free(ColorD);
	ColorA = NULL;
	ColorB = NULL;
	ColorC = NULL;
	ColorD = NULL;
}

//Build the MSD Histogram
void microdiscriptor(int** ColorX, double* hist, int wid, int hei,int flag)
{
	//int* MS = (int*)calloc(8,sizeof(int));
    //int* HA = (int*)calloc(8,sizeof(int));

	int MS[8] = {0};
	int HA[8] = {0};
	IplImage** CompoBin = (IplImage**)calloc(CSA,sizeof(IplImage*)); 
	int i,j;

	for(i = 0; i < CSA;i++)
	{
		CompoBin[i] = cvCreateImage(cvSize(wid,hei),IPL_DEPTH_8U,1);
		cvSetZero(CompoBin[i]);
	}
	
	//hist = (double*)calloc(CSA,sizeof(double));

   
    //----------------------------------------
    for (i = 0; i < wid - 1; i++)
    {
        for (j = 0; j < hei - 1; j++)
        {
            if (ColorX[i][j] >= 0)
            {
                HA[ColorX[i][j]] += 1;
				if(flag == 0)
					cvSet2D(CompoBin[ColorX[i][j]],j,i,cvScalar(255,0,0,0));
            }
        }
    }
    //----------------------------------------
    for (i = 3; i < 3 * (wid / 3) - 1; i++)
    {
        for (j = 3; j < 3 * (hei / 3) - 1; j++)
        {

            int wa[9];
			int TE1 = 0;
			int m;
            wa[0] = ColorX[i - 1][j - 1];
            wa[1] = ColorX[i - 1][j];
            wa[2] = ColorX[i - 1][j + 1];

            wa[3] = ColorX[i + 1][j - 1];
            wa[4] = ColorX[i + 1][j];
            wa[5] = ColorX[i + 1][j + 1];

            wa[6] = ColorX[i][j - 1];
            wa[7] = ColorX[i][j + 1];
            wa[8] = ColorX[i][j];
            //-------------------------
            

            for ( m = 0; m < 8; m++)
            {
                if ((wa[8] == wa[m]) && (wa[8] >= 0))
                {
                    TE1 = TE1 + 1;
                }

            }
            if (wa[8] >= 0)
            {
                MS[wa[8]] += TE1;
            }


        }
    }
	
	int count = 0;
    // the features vector of MSD 
    for (int i = 0; i < CSA; i++)
    {
        hist[i] = (MS[i] * 1.0) / (8.0 * HA[i] + 0.0001);
		if(flag == 0)
		{
			IplImage* CompoBinDilated = cvCreateImage(cvSize(wid,hei),IPL_DEPTH_8U,1);
			CvMemStorage* storage = cvCreateMemStorage(0);
			CvSeq* contours = 0;
			cvDilate(CompoBin[i],CompoBinDilated);
			//cvShowImage("Image",CompoBin[i]);cvWaitKey(0);
			cvFindContours( CompoBinDilated, storage, &contours, sizeof(CvContour),CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE, cvPoint(0,0) );
			count = 0;
			while(contours != NULL)
			{
				contours = contours->h_next;
				count++;
			}
			hist[CSA + i] = count;
			
			cvReleaseMemStorage(&storage);
			cvReleaseImage(&CompoBinDilated);
		}
		
    }
	
	for(int i = 0; i < CSA;i++)
	{
		cvReleaseImage(&CompoBin[i]);
	}
	/*free(MS);
	free(HA);
	MS = NULL;
	HA = NULL;
	*/
	
}

//Takes HSV and HSV_qnt as input, uses the above funstions, gives Histogram as output 
void MSD_feature_extract(int*** HSV,int** ImageX, double* hist, int wid, int hei,int flag)
{
	int i,j;
	//double*** HSV = (double***)calloc(3,sizeof(double**));
    int** ori = (int**)calloc(wid,sizeof(int*)); // the orientation image
    //int** ImageX = (int**)calloc(wid,sizeof(int*)); // the color index image
	int** micro = (int**)calloc(wid,sizeof(int*));  //micro-structure image

	for(i = 0; i < wid;i++)
	{
		ori[i] = (int*)calloc(hei,sizeof(int));
		//ImageX[i] = (int*)calloc(hei,sizeof(int));
		micro[i] = (int*)calloc(hei,sizeof(int));
	}
   
	/*for(i = 0; i < 3;i++)
	{
		HSV[i] = (double**)calloc(wid,sizeof(double*));
		for(j = 0; j < wid;j++)
		{
			HSV[i][j] = (double*)calloc(hei,sizeof(double));
		}
	}*/

    //mat_color(RGB, HSV, wid, hei); //transform RGB color space to  HSV color space 
	
	 /*cvNamedWindow("myfirstwindow");
	 cvShowImage("myfirstwindow", HSV);
	 cvWaitKey(0);*/

    Mat_ori_HSV(HSV, ori, wid, hei);  // edge orientation detection in HSV color space
                                                                                                                                                                                                                                                                                                                                
    //Mat_color_HSV(HSV,ImageX, Cn1, Cn2, Cn3, wid, hei);  //color quantization in HSV color space     
    
    microstructure(ori, ImageX, micro, wid, hei);   //micro-structure map extraction
	//free the variables
	for(i=0;i<wid;i++)
	{
		free(ori[i]);
	}
	free(ori);
	ori = NULL;
   
    microdiscriptor(micro, hist, wid, hei,flag);    //micro-structure representation

	//free the variables
	for(i=0;i<wid;i++)
	{
		free(micro[i]);
	}
	free(micro);
	micro = NULL;
 }

//*********************************************************************************//
//-----------------------For HOG Feature Extraction--------------------------------//
//*********************************************************************************//
// Quantize the HSV values in the required number of bins (8 for now)
	void HSV_color_quant(IplImage* HSVImg, IplImage* HSV_qnt_mat)
	{
		int h = HSVImg->height;
		int w = HSVImg->width;
		int VI = 0, SI = 0, HI = 0;
		//FILE* fileptr_hsv_chk = fopen("G:\\Flags project\\pick component\\hsv_check.txt","w");
		CvScalar s;

		for(int i = 0;i < h;i++)
		{
			for(int j = 0;j < w;j++)
			{
				s = cvGet2D(HSVImg,i,j);
				HI = (int) (s.val[0]*2);
                SI = (int) ((s.val[1]/255)*100);
                VI = (int) ((s.val[2]/255)*100);
				s.val[0] = 10;
				s.val[1] = 0;
				s.val[2] = 0;
				s.val[3] = 0;
				if (VI < 16)
				{
                    s.val[0] = 0;
                    cvSet2D(HSV_qnt_mat, i, j,s);				
				}
                else
                    {
                        //Check for Saturation - if less than 16, threshold at 31 (Black and white)
                        if (SI < 16)
                        {
                            if (VI < 31)
							{
								s.val[0] = 0;
                                cvSet2D(HSV_qnt_mat,i, j,s);
							}
                            else
							{
								s.val[0] = 1;
                                cvSet2D(HSV_qnt_mat, i, j,s);
							}
                        }
                        //Check for Saturation - if more than 16, threshold at 16 (Black and Color)
                        else
                        {
                            if (VI < 16)
							{
								s.val[0] = 0;	
                                cvSet2D(HSV_qnt_mat, i, j,s);
							}
                            else
                            {
                                if ((HI >= 0 && HI < 30) || (HI >= 330 && HI < 360))
								{
									s.val[0] = 2;
                                    cvSet2D(HSV_qnt_mat, i, j,s);
								}
                                if (HI >= 30 && HI < 90)
								{
									s.val[0] = 3;
                                    cvSet2D(HSV_qnt_mat, i, j,s);
								}
                                if (HI >= 90 && HI < 150)
								{
                                    s.val[0] = 4;
                                    cvSet2D(HSV_qnt_mat, i, j,s);
								}
                                if (HI >= 150 && HI < 210)
								{
                                    s.val[0] = 5;
                                    cvSet2D(HSV_qnt_mat, i, j,s);
								}
                                if (HI >= 210 && HI < 270)
								{
                                    s.val[0] = 6;
                                    cvSet2D(HSV_qnt_mat, i, j,s);
								}
                                if (HI >= 270 && HI < 330)
								{
                                    s.val[0] = 7;
                                    cvSet2D(HSV_qnt_mat, i, j,s);
								}
                            }
                        }
                    }//end-of-main-if
			}
		}
		/*CvScalar hsv;
		for(int i =0; i<h;i++)
		{
			for(int j=0;j<w;j++)
			{
				hsv = cvGet2D(HSV_qnt_mat,i,j);
				fprintf(fileptr_hsv_chk, "%lf ", hsv.val[0]);
			}
			fprintf(fileptr_hsv_chk,"\n");
		}
		fclose(fileptr_hsv_chk);*/
	}



// Extract HOG Features for the Input Image and write them into a file
	void HOG_FeatureEx_FileCreation(IplImage * IpImage_Gray, char * label, FILE* fileptr, int count)
	{
		int BlockSize_h= IpImage_Gray->height/1;	//block height *block size same as image*     
		int BlockSize_w= IpImage_Gray->width/1;		//block width

		int OverlapFactor = 1;						//amount of overlap : half in this case //change1

		/*for display of the gradient images*/
		//IplImage * IpImage_GrayGradient = 0;       //Gradient Image
		//IplImage TempBuffer;                       // temporary buffer

		/*Display the Gray Scale Image*/
		//cvShowImage("Gray Scale Image",IpImage_Gray);cvWaitKey(5000);
		
		/*Creation of the GrayScale and Gradient Matrices(for filter operations)*/
		CvMat* Mat_Gray = cvCreateMat(IpImage_Gray->height,IpImage_Gray->width,CV_8UC1);//to be freed
		cvConvertImage(IpImage_Gray,Mat_Gray);
		
		/*Create a matrix for containing the gray scale gradient (Vertical)*/
		CvMat* Mat_GrayGradient_Vert1 = cvCreateMat(IpImage_Gray->height,IpImage_Gray->width,CV_8UC1);//to be freed
		CvMat* Mat_GrayGradient_Vert2 = cvCreateMat(IpImage_Gray->height,IpImage_Gray->width,CV_8UC1);//to be freed
		CvMat* Mat_GrayGradient_Vert = cvCreateMat(IpImage_Gray->height,IpImage_Gray->width,CV_8UC1);//to be freed
		/*Create a matrix for containing the gray scale gradient (Horizontal)*/
		CvMat* Mat_GrayGradient_Horz1 = cvCreateMat(IpImage_Gray->height,IpImage_Gray->width,CV_8UC1);//to be freed
		CvMat* Mat_GrayGradient_Horz2 = cvCreateMat(IpImage_Gray->height,IpImage_Gray->width,CV_8UC1);//to be freed
		CvMat* Mat_GrayGradient_Horz = cvCreateMat(IpImage_Gray->height,IpImage_Gray->width,CV_8UC1);//to be freed
		

		/* Create the kernel matrix -- for vertical edges*/
		float DataMV[] = {1,0,-1,2,0,-2,1,0,-1};
		CvMat Mat_KernelV = cvMat(3,3,CV_32FC1,DataMV);
		/*Perform 2D Gradient Filtering -- vertical edges*/
		cvFilter2D(Mat_Gray,Mat_GrayGradient_Vert1,&Mat_KernelV);

		DataMV[0] = -1;
		DataMV[1] = 0;
		DataMV[2] = 1;
		DataMV[3] = -2;
		DataMV[4] = 0;
		DataMV[5] = 2;
		DataMV[6] = -1;
		DataMV[7] = 0;
		DataMV[8] = 1;
		Mat_KernelV = cvMat(3,3,CV_32FC1,DataMV);
		/*Perform 2D Gradient Filtering -- vertical edges*/
		cvFilter2D(Mat_Gray,Mat_GrayGradient_Vert2,&Mat_KernelV);

		AddMat_Bin(Mat_GrayGradient_Vert1, Mat_GrayGradient_Vert2, Mat_GrayGradient_Vert);
		

		/*Create the Kernel matrix -- for horizontal edges*/
		float DataMH[] = {1,2,1,0,0,0,-1,-2,-1};
		CvMat Mat_KernelH = cvMat(3,3,CV_32FC1,DataMH);
		/*Perform 2D Gradient Filtering -- horizontal filtering*/
		cvFilter2D(Mat_Gray,Mat_GrayGradient_Horz1,&Mat_KernelH);

		DataMH[0] = -1;
		DataMH[1] = -2;
		DataMH[2] = -1;
		DataMH[3] = 0;
		DataMH[4] = 0;
		DataMH[5] = 0;
		DataMH[6] = 1;
		DataMH[7] = 2;
		DataMH[8] = 1;
		Mat_KernelH = cvMat(3,3,CV_32FC1,DataMH);
		/*Perform 2D Gradient Filtering -- horizontal filtering*/
		cvFilter2D(Mat_Gray,Mat_GrayGradient_Horz2,&Mat_KernelH);
		
		AddMat_Bin(Mat_GrayGradient_Horz1, Mat_GrayGradient_Horz2, Mat_GrayGradient_Horz);


		/*Dont want display 
		//-----------------------------strictly for display----------------------------------------//
		//Convert Mat to IplImage for display	
		IpImage_GrayGradient = cvGetImage(Mat_GrayGradient_Vert,&TempBuffer);
		// Display the gradient image
		cvShowImage("Gray Gradient Image (Vertical edges)",IpImage_GrayGradient);cvWaitKey(5000);
		
		//Convert Mat to IplImage for display	
		IpImage_GrayGradient = cvGetImage(Mat_GrayGradient_Horz,&TempBuffer);
		// Display the gradient image
		cvShowImage("Gray Gradient Image (Horizontal Edges)",IpImage_GrayGradient);cvWaitKey(5000);
		//-------------------------------------------------------------------------------------------//
		*/

		/*Magnitude and Angle Matrix*/
		CvMat * Mat_GradientMag = cvCreateMat(IpImage_Gray->height,IpImage_Gray->width,CV_8UC1);//to be freed
		CvMat * Mat_GradientTheta = cvCreateMat(IpImage_Gray->height,IpImage_Gray->width,CV_8UC1);//to be freed
		
		for(int i=0;i<IpImage_Gray->height;i++)
			for(int j=0;j<IpImage_Gray->width;j++)
			{
				double gradx = cvGetReal2D(Mat_GrayGradient_Horz,i,j);
				double grady = cvGetReal2D(Mat_GrayGradient_Vert,i,j);
				cvSetReal2D(Mat_GradientMag,i,j,sqrt(pow(gradx,2) + pow(grady,2)));\
				//Keep the gradient values in 0-180 degree range
				if(grady<0)
				{
					grady = grady*(-1);
					gradx = gradx*(-1);
				}
				cvSetReal2D(Mat_GradientTheta,i,j,(atan2(grady, gradx)+pi)*180/pi);
			}
		
		/*Declaration of the 2D array of pointers - each to the histogram of a block*/
		int ySteps = (OverlapFactor*((IpImage_Gray->height/BlockSize_h)-1))+1;
		int xSteps = (OverlapFactor*((IpImage_Gray->width/BlockSize_w)-1))+1;
		int *** Arr_BlockBinTotal = (int***) calloc(ySteps,sizeof(int**));//to be freed
		for(int i=0;i<ySteps;i++)
		{	
			//*(Arr_BlockBinTotal+i) = (int**) calloc(xSteps, sizeof(int*));
			Arr_BlockBinTotal[i] = (int**) calloc(xSteps, sizeof(int*));
			for(int j=0;j<xSteps;j++)
			{
				//*(*(Arr_BlockBinTotal+i)+j) = (int*) calloc(BinCount, sizeof(int));
				Arr_BlockBinTotal[i][j] = (int*) calloc(BinCount, sizeof(int));
			}
		}
		
		/*(I/P)Gradient Magnitude, Gradient Theta, BlockSize, x, y -->(O/P)Histogram for the cell*/
		// Assumption:: square window and the dimensions of the image are whole multiple of window size//
		for(int indy=0; indy<ySteps; indy++)
			for(int indx=0; indx<xSteps; indx++)
			{
				int XStartingPoint = indx*(BlockSize_w/OverlapFactor);
				int YStartingPoint = indy*(BlockSize_h/OverlapFactor);
				BlockBinTotal(Mat_GradientMag, Mat_GradientTheta, BlockSize_w, BlockSize_h, XStartingPoint, YStartingPoint, Arr_BlockBinTotal[indy][indx]);
			}

		/*Write the HOG Features of the image with the label*/
		/*concatenate the FV referenced in each point in 2-D Grid and send as whole FV*/
		
		for(int i=0; i<ySteps; i++)
			for(int j=0; j<xSteps; j++)
			{
				int sum = BlockBinTotal_Sum(Arr_BlockBinTotal[i][j]);
				for(int k=0; k<BinCount; k++)
				{	
					float FeatureValue_norm = (float)(*((Arr_BlockBinTotal[i][j])+k))/((float)sum+0.0001);
					if(FeatureValue_norm!=0.0)
					{
						fprintf(fileptr,"%d:%f ",++count,FeatureValue_norm);
					}
					else
					++count;
				}
			}
		fprintf(fileptr,"\n");

		/*free all the pointers here*/
		for(int i=0;i<ySteps;i++)
		{
			for(int j=0;j<xSteps;j++)
			{
				free(Arr_BlockBinTotal[i][j]);
			}
			free(Arr_BlockBinTotal[i]);
		}
		free(Arr_BlockBinTotal);
		Arr_BlockBinTotal = NULL;
		/*free all the arrays*/
		cvReleaseMat(&Mat_Gray);
		cvReleaseMat(&Mat_GrayGradient_Vert);
		cvReleaseMat(&Mat_GrayGradient_Vert1);
		cvReleaseMat(&Mat_GrayGradient_Vert2);
		cvReleaseMat(&Mat_GrayGradient_Horz);
		cvReleaseMat(&Mat_GrayGradient_Horz1);
		cvReleaseMat(&Mat_GrayGradient_Horz2);
		cvReleaseMat(&Mat_GradientMag);
		cvReleaseMat(&Mat_GradientTheta);
		Mat_Gray=NULL;
		Mat_GrayGradient_Vert=NULL;
	    Mat_GrayGradient_Vert1=NULL;
		Mat_GrayGradient_Vert2=NULL;
		Mat_GrayGradient_Horz=NULL;
		Mat_GrayGradient_Horz1=NULL;
		Mat_GrayGradient_Horz2=NULL;
		Mat_GradientMag=NULL;
		Mat_GradientTheta=NULL;
	}



// Average the matrices
	void AddMat_Bin(CvMat* Mat1, CvMat* Mat2, CvMat* Mat_tot)
	{
		int h = Mat1->height;
		int w = Mat1->width;

		CvScalar s1,s2,s;

		for(int i = 0; i < h ; i++)
		{
			for(int j = 0 ; j < w ; j++)
			{
				s1 = cvGet2D(Mat1,i,j);		
				s2 = cvGet2D(Mat2,i,j);

				
				if(s1.val[0]>s2.val[0])
					s.val[0] = s1.val[0];
				else 
					s.val[0] = s2.val[0];

				cvSet2D(Mat_tot,i,j,s);
			}
		}
	}






// Calculate Histogram of each specified block
	void BlockBinTotal(CvMat* M_Magnitude, CvMat* M_Theta, int BlockSize_w, int BlockSize_h, int x, int y, int* EachBinTotal)
	{ 
			
		for(int i=0; i<BlockSize_w; i++)
			for(int j=0; j<BlockSize_h; j++)
			{
				double theta = cvGetReal2D(M_Theta,y+j,x+i);
				int groupno = floor(theta/5);
				int temp = cvGetReal2D(M_Magnitude,y+j,x+i);
				*(EachBinTotal+groupno) += cvGetReal2D(M_Magnitude,y+j,x+i); 
			}
	}




// Normalise the Histogram of each Block---------------------------------//
	int BlockBinTotal_Sum(int * inputptr)
	{
		float sum=0;
		for(int i=0;i<BinCount;i++)
		{
			sum += *(inputptr+i);
		}
		return (sum);
	}
//*********************************************************************************//
//---------------------------------Main Function-----------------------------------//
//*********************************************************************************//
int main(void)
{
	int count1 = 0; //for component count
	FILE* fileptr = fopen("G:\\flags_all\\mag_only\\MSD+HOG180_test_features.txt","w");
	FILE* fileptr_comp = fopen("G:\\flags_all\\mag_only\\MSD+HOG_test_features_others.txt","w");
	
	DIR *dir_m, *dir_sub;
	struct dirent *ent_m, *ent_sub;
	char dirpath_m[] = "G:\\flags_all\\mag_only\\flags_test";
	char fullpath[500],filename[200];
    dir_m = opendir (dirpath_m);
	if (dir_m != NULL) 
	{
		//Read the folders in the main directory
		while ((ent_m = readdir (dir_m)) != NULL) 
		{
			if(strcmp(".",ent_m->d_name) && strcmp("..",ent_m->d_name)  && strcmp("Thumbs.db",ent_m->d_name))
			{
				strcpy(fullpath,dirpath_m);
				strcat(fullpath,"\\");
				strcat(fullpath,ent_m->d_name);
				
				dir_sub = opendir(fullpath);
				if(dir_sub != NULL)
				{
					//Read the image files within each folder
					while((ent_sub = readdir (dir_sub)) != NULL)
					{
						if(strcmp(".",ent_sub->d_name) && strcmp("..",ent_sub->d_name)  && strcmp("Thumbs.db",ent_sub->d_name))
						{
							strcpy(filename,fullpath);
							strcat(filename,"\\");
							strcat(filename,ent_sub->d_name);
							
							IplImage* IpImage_Color = cvLoadImage(filename,CV_LOAD_IMAGE_COLOR);
							
							IplImage* IpImage_Resized = 0;
							IplImage* IpImage_HSV = 0;
							IplImage* HSV_qnt  = 0;
							flag = 0;
							//Unable to load
							if(!IpImage_Color)
							{
								fprintf(stderr,"\nFailed to load input image : %s->%s",ent_m->d_name,ent_sub->d_name);
								failcount++;
								continue;
								//return -1;
								
							}

							int w = IpImage_Color->width;
							int h = IpImage_Color->height;
							int i , j;

							//Resizing
							//Fit the image into a max size (1000 : the lengthiest side)
							int maxside;
							if(w>h)
								maxside = w;
							else
								maxside = h;
                            if (maxside > 1000)
                            {
                                float as_ratio = (((float)maxside) / 1000);

                                h = (int)floor(h / as_ratio);
                                w = (int)floor(w / as_ratio);
							}
							IpImage_Resized = cvCreateImage(cvSize(w,h),IPL_DEPTH_8U,3);
							cvResize(IpImage_Color, IpImage_Resized,CV_INTER_CUBIC);
							
							//IpImage_Resized = cvCreateImage(cvSize(w,h),IPL_DEPTH_8U,3);
							//cvCopy(IpImage_Color, IpImage_Resized);


							//Convert from RGB to HSV
							IpImage_HSV = cvCreateImage(cvSize(w,h), IPL_DEPTH_8U, 3);
							cvCvtColor(IpImage_Resized,IpImage_HSV,CV_BGR2HSV);
							
							//Quantize the HSV Values
							HSV_qnt = cvCreateImage(cvSize(w,h), IPL_DEPTH_8U, 1);
							HSV_color_quant(IpImage_HSV, HSV_qnt);
							//cvShowImage("quantized",HSV_qnt);cvWaitKey(0);
							
							//Fetch the Class ID - label
							char label[10]= "";
							strcpy(label,ent_m->d_name);

							int part_w = w / 3;
                            int part_h = h / 3;

							int*** image = (int***)calloc(3,sizeof(int**));
							int*** image1 = (int***)calloc(3,sizeof(int**));
							int*** image2 = (int***)calloc(3,sizeof(int**));
							int*** image3 = (int***)calloc(3,sizeof(int**));
							int*** image4 = (int***)calloc(3,sizeof(int**));
							int** ImageX = (int**)calloc(w,sizeof(int*));
							int** ImageX1 = (int**)calloc(2*part_w,sizeof(int*));
							int** ImageX2 = (int**)calloc(w-part_w,sizeof(int*));
							int** ImageX3 = (int**)calloc(2*part_w,sizeof(int*));
							int** ImageX4 = (int**)calloc(w-part_w,sizeof(int*));
							
							double* hist;
							if(flag == 1)
								hist = (double*)calloc(CSA,sizeof(double)); // the features vector of MSD 
							else
								hist = (double*)calloc(2 * CSA,sizeof(double));// // the features vector of MSD + component count
							
							double* hist1 = (double*)calloc(CSA,sizeof(double));
							double* hist2 = (double*)calloc(CSA,sizeof(double));
							double* hist3 = (double*)calloc(CSA,sizeof(double));
							double* hist4 = (double*)calloc(CSA,sizeof(double));

							//call the function-- MSD Feature Extraction and creation of SVM sample file
							
							//Allocate memory to the variables
							
							//For HSV 3D array
							for(i = 0; i < 3;i++)
							{
								image[i] = (int**)calloc(w,sizeof(int*));
								image1[i] = (int**)calloc(2*part_w,sizeof(int*));
								image2[i] = (int**)calloc(w-part_w,sizeof(int*));
								image3[i] = (int**)calloc(2*part_w,sizeof(int*));
								image4[i] = (int**)calloc(w-part_w,sizeof(int*));
								for(j = 0; j < w;j++)
								{
									image[i][j] = (int*)calloc(h,sizeof(int));
									if(j<2*part_w)
									{
										image1[i][j] = (int*)calloc(2*part_h,sizeof(int));
										image3[i][j] = (int*)calloc(h-part_h,sizeof(int));
									}
									if(j<(w-part_w))
									{
										image2[i][j] = (int*)calloc(2*part_h,sizeof(int));
										image4[i][j] = (int*)calloc(h-part_h,sizeof(int));
									}
								}
							}
							//For HSV Quant 1D array
							for(j = 0; j < w;j++)
								{
									ImageX[j] = (int*)calloc(h,sizeof(int));
									if(j<2*part_w)
									{
										ImageX1[j] = (int*)calloc(2*part_h,sizeof(int));
										ImageX3[j] = (int*)calloc(h-part_h,sizeof(int));
									}
									if(j<(w-part_w))
									{
										ImageX2[j] = (int*)calloc(2*part_h,sizeof(int));
										ImageX4[j] = (int*)calloc(h-part_h,sizeof(int));
									}
								}
							//Put values in the array
							CvScalar color;
							CvScalar color_qnt;
							for(i = 0; i < w;i++)
							{
								for(j = 0;j < h;j++)
								{
									color = cvGet2D(IpImage_HSV,j,i);
									color_qnt = cvGet2D(HSV_qnt,j,i);

									image[0][i][j] = color.val[0];
									image[1][i][j] = color.val[1];
									image[2][i][j] = color.val[2];
									ImageX[i][j] = color_qnt.val[0];
									if (i < 2 * part_w && j<2*part_h)
                                    {
                                     image1[0][i][j] = color.val[0];
                                     image1[1][i][j] = color.val[1];
                                     image1[2][i][j] = color.val[2];
									 ImageX1[i][j] = color_qnt.val[0];
                                    }
                                    if ((i >= part_w && i < w) && (j < 2 * part_h))
                                    {
                                     image2[0][i - part_w][j] = color.val[0];
                                     image2[1][i - part_w][j] = color.val[1];
                                     image2[2][i - part_w][j] = color.val[2];
									 ImageX2[i - part_w][j] = color_qnt.val[0];
                                    }
                                    if ((i < 2 * part_w) && (j >= part_h && j < h))
                                    {
                                     image3[0][i][j - part_h] = color.val[0];
                                     image3[1][i][j - part_h] = color.val[1];
                                     image3[2][i][j - part_h] = color.val[2];
									 ImageX3[i][j - part_h] = color_qnt.val[0];
                                    }
                                    if ((i >= part_w && i < w) && (j >= part_h && j < h))
                                    {
                                     image4[0][i - part_w][j - part_h] = color.val[0];
                                     image4[1][i - part_w][j - part_h] = color.val[1];
                                     image4[2][i - part_w][j - part_h] = color.val[2];
									 ImageX4[i - part_w][j - part_h] = color_qnt.val[0];
                                    }
								}
							}

							// Call the function for MSD based Histogram evaluation
							// And Write MSD features to the file

							
							//0.
							
							MSD_feature_extract(image,ImageX,hist,w,h,flag);
							flag = 1;
							fprintf(fileptr,"%s ",label);
							fprintf(fileptr_comp,"%s ",label);
							count = 0;
							for (int iHIstIndex = 0; iHIstIndex < CSA; iHIstIndex++)
							{	
								if(hist[iHIstIndex]!=0.0)
								{
									fprintf(fileptr,"%d:%f ",++count,hist[iHIstIndex]);
								}
								else
								++count;
							}
							count1 = 0;
							for (int iHIstIndex = CSA; iHIstIndex < 2 * CSA; iHIstIndex++)
							{	
								if(hist[iHIstIndex]!=0.0)
								{
									fprintf(fileptr_comp,"%d:%f ",++count1,hist[iHIstIndex]);
								}
								else
								++count1;
							}
							fprintf(fileptr_comp,"\n");

							free(hist);							
							for(i=0;i<w;i++)
							{
								free(ImageX[i]);
							}
							free(ImageX);
							for(i=0;i<3;i++)
							{
								for(j=0;j<w;j++)
								{
									free(image[i][j]);
								}
								free(image[i]);
							}
							free(image);
							

							//1.
							MSD_feature_extract(image1,ImageX1,hist1,2*part_w,2*part_h,flag);
							for (int iHIstIndex = 0; iHIstIndex < CSA; iHIstIndex++)
							{	
								if(hist1[iHIstIndex]!=0.0)
								{
									fprintf(fileptr,"%d:%f ",++count,hist1[iHIstIndex]);
								}
								else
								++count;
							}
							free(hist1);
							for(i=0;i<2*part_w;i++)
							{
								free(ImageX1[i]);
							}
							free(ImageX1);
							for(i=0;i<3;i++)
							{
								for(j=0;j<2*part_w;j++)
								{
									free(image1[i][j]);
								}
								free(image1[i]);
							}
							free(image1);


							//2.
							MSD_feature_extract(image2,ImageX2,hist2,w-part_w,2*part_h,flag);
							for (int iHIstIndex = 0; iHIstIndex < CSA; iHIstIndex++)
							{	
								if(hist2[iHIstIndex]!=0.0)
								{
									fprintf(fileptr,"%d:%f ",++count,hist2[iHIstIndex]);
								}
								else
								++count;
							}
							free(hist2);
							for(i=0;i<w-part_w;i++)
							{
								free(ImageX2[i]);
							}
							free(ImageX2);
							for(i=0;i<3;i++)
							{
								for(j=0;j<w-part_w;j++)
								{
									free(image2[i][j]);
								}
								free(image2[i]);
							}
							free(image2);

							//3.
							MSD_feature_extract(image3, ImageX3, hist3, 2*part_w, h-part_h,flag);
							for (int iHIstIndex = 0; iHIstIndex < CSA; iHIstIndex++)
							{	
								if(hist3[iHIstIndex]!=0.0)
								{
									fprintf(fileptr,"%d:%f ",++count,hist3[iHIstIndex]);
								}
								else
								++count;
							}
							free(hist3);
							for(i=0;i<2*part_w;i++)
							{
								free(ImageX3[i]);
							}
							free(ImageX3);
							for(i=0;i<3;i++)
							{
								for(j=0;j<2*part_w;j++)
								{
									free(image3[i][j]);
								}
								free(image3[i]);
							}
							free(image3);

							//4.
							MSD_feature_extract(image4, ImageX4, hist4, w-part_w, h-part_h,flag);
							for (int iHIstIndex = 0; iHIstIndex < CSA; iHIstIndex++)
							{	
								if(hist4[iHIstIndex]!=0.0)
								{
									fprintf(fileptr,"%d:%f ",++count,hist4[iHIstIndex]);
								}
								else
								++count;
							}

							free(hist4);
							for(i=0;i<w-part_w;i++)
							{
								free(ImageX4[i]);
							}
							free(ImageX4);
							for(i=0;i<3;i++)
							{
								for(j=0;j<w-part_w;j++)
								{
									free(image4[i][j]);
								}
								free(image4[i]);
							}
							free(image4);
							
							/*//clean all the memory
							free(hist);
							free(hist1);
							free(hist2);
							free(hist3);
							free(hist4);
							for(i=0;i<w;i++)
							{
								free(ImageX[i]);
								if(i<2*part_w)
								{
									free(ImageX1[i]);
									free(ImageX3[i]);
								}
								if(i<(w-part_w))
								{
									free(ImageX2[i]);
									free(ImageX4[i]);
								}
							}
							free(ImageX);
							free(ImageX1);
							free(ImageX2);
							free(ImageX3);
							free(ImageX4);

							for(i=0;i<3;i++)
							{
								for(j=0;j<w;j++)
								{
									free(image[i][j]);
									if(j<2*part_w)
									{
										free(image1[i][j]);
										free(image3[i][j]);
									}
									if(j<(w-part_w))
									{
										free(image2[i][j]);
										free(image4[i][j]);
									}
								}
								free(image[i]);
								free(image1[i]);
								free(image2[i]);
								free(image3[i]);
								free(image4[i]);
							}
							free(image);
							free(image1);
							free(image2);
							free(image3);
							free(image4);
							*/
							
							// Call the function for HOG descriptor generation
							HOG_FeatureEx_FileCreation(HSV_qnt,label,fileptr, count);
							
							//Release Image Data
							cvReleaseImage(&IpImage_Resized);
							cvReleaseImage(&IpImage_Color);
							cvReleaseImage(&IpImage_HSV);
							cvReleaseImage(&HSV_qnt);
							
							IpImage_Resized = NULL;
							IpImage_Color = NULL;
							IpImage_HSV = NULL;
							HSV_qnt = NULL;

							image = NULL;
							image1 = NULL;
							image2 = NULL;
							image3 = NULL;
							image4 = NULL;

							ImageX = NULL;
							ImageX1 = NULL;
							ImageX2 = NULL;
							ImageX3 = NULL;
							ImageX4 = NULL;

							hist = NULL;
							hist1 = NULL;
							hist2 = NULL;
							hist3 = NULL;
							hist4 = NULL;
							
						}
					}
				}
			}
		}
	}
	fclose(fileptr_comp);
	fclose(fileptr);
	return 0;
}