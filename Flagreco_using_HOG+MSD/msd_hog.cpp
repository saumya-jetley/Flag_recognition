#include "msd_hog.h"
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
           
            gh = (double)(hsv[1][i - 1][j + 1] + 2 * hsv[1][i][j + 1] + hsv[1][i + 1][j + 1]) - (hsv[1][i - 1][j - 1] + 2 * hsv[1][i][j - 1] + hsv[1][i + 1][j - 1]);
			
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


//Give the color patch image in "Color" (color map) based on texture similarity from ORI input - pattern analysis from Dx, Dy
void Map(int** ori, int** img, int** Color, int wid, int hei, int Dx, int Dy)
{
	int i,j;
	
    for (int i = 1; i < wid / 3 ; i++)
    {
        for (int j = 1; j < hei / 3 ; j++)
        {
            int WA[9];
            //===========================
            int m = 3 * i + Dx;
            int n = 3 * j + Dy;

            WA[0] = ori[m - 1][n - 1];
            WA[1] = ori[m - 1][ n];
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

    //=========the final micro-structure map===============//
   
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

//Build the MSD Histogram//
void microdiscriptor(int** ColorX, double* hist, int wid, int hei)
{
	int MS[8] = {0};
	int HA[8] = {0};

	int i,j;
	   
    //----------------------------------------
    for (i = 0; i < wid - 1; i++)
    {
        for (j = 0; j < hei - 1; j++)
        {
            if (ColorX[i][j] >= 0)
            {
                HA[ColorX[i][j]] += 1;
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

    // the features vector of MSD //
    for (int i = 0; i < CSA; i++)
    {
        hist[i] = (MS[i] * 1.0) / (8.0 * HA[i] + 0.0001);
    }

}

//Takes HSV and HSV_qnt as input, uses the above functions, gives Histogram as output// 
void MSD_feature_extract(int*** HSV,int** ImageX, double* hist, int wid, int hei)
{
	int i,j;
	
    int** ori = (int**)calloc(wid,sizeof(int*)); // the orientation image
    
	int** micro = (int**)calloc(wid,sizeof(int*));  //micro-structure image

	for(i = 0; i < wid;i++)
	{
		ori[i] = (int*)calloc(hei,sizeof(int));
		
		micro[i] = (int*)calloc(hei,sizeof(int));
	}
   
	Mat_ori_HSV(HSV, ori, wid, hei);  // edge orientation detection in HSV color space
                                                                                                                                                                                                                                                                                                                                
    microstructure(ori, ImageX, micro, wid, hei);   //micro-structure map extraction
	
	//free the variables
	for(i=0;i<wid;i++)
	{
		free(ori[i]);
	}
	free(ori);
	ori = NULL;
   
    microdiscriptor(micro, hist, wid, hei);    //micro-structure representation

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

		CvScalar s;

		for(int i = 0;i < h;i++)
		{
			for(int j = 0;j < w;j++)
			{
				s = cvGet2D(HSVImg,i,j);
				HI = (int) (s.val[0]*2);
                		SI = (int) ((s.val[1]/255)*100);
		                VI = (int) ((s.val[2]/255)*100);
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
	                        //Check for Value - if more than 16, threshold at 16 (Black and Color)
        	                else
                	        {
                        	    if (VI < 16)
					{
					s.val[0] = 0;	
                                	cvSet2D(HSV_qnt_mat, i, j,s);
					}
	                            else
        	                    {
                	                if ((HI >= 0 && HI < 30) || (HI > 330 && HI <= 359))
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
	}

// Extract HOG Features for the Input Image and write them into a file
	void HOG_FeatureEx_FileCreation(IplImage * IpImage_Gray, double* hog_feat)
	{
		int BlockSize_h= IpImage_Gray->height/1;	//block height *block size same as image*     
		int BlockSize_w= IpImage_Gray->width/1;		//block width

		int OverlapFactor = 1;						//amount of overlap : half in this case //change1

		/*Creation of the GrayScale and Gradient Matrices(for filter operations)*/
		CvMat* Mat_Gray = cvCreateMat(IpImage_Gray->height,IpImage_Gray->width,CV_8UC1);//to be freed
//1.LINUX Changes
		cvConvertImage(IpImage_Gray,Mat_Gray);
		
		/*Create a matrix for containing the gray scale gradient (Vertical)*/
		CvMat* Mat_GrayGradient_Vert = cvCreateMat(IpImage_Gray->height,IpImage_Gray->width,CV_16SC1);//to be freed
		/*Create a matrix for containing the gray scale gradient (Horizontal)*/
		CvMat* Mat_GrayGradient_Horz = cvCreateMat(IpImage_Gray->height,IpImage_Gray->width,CV_16SC1);//to be freed
		
//2. LINUX Changes
		/* Create the kernel matrix -- for vertical edges*/
		float DataMV[] = {1,0,-1,2,0,-2,1,0,-1};
		CvMat Mat_KernelV = cvMat(3,3,CV_32FC1,DataMV);
		/*Perform 2D Gradient Filtering -- vertical edges*/
		cvFilter2D(Mat_Gray,Mat_GrayGradient_Vert,&Mat_KernelV);

		/*Create the Kernel matrix -- for horizontal edges*/
		float DataMH[] = {1,2,1,0,0,0,-1,-2,-1};
		CvMat Mat_KernelH = cvMat(3,3,CV_32FC1,DataMH);
		/*Perform 2D Gradient Filtering -- horizontal filtering*/
		cvFilter2D(Mat_Gray,Mat_GrayGradient_Horz,&Mat_KernelH);

		/*Magnitude and Angle Matrix*/
		CvMat * Mat_GradientMag = cvCreateMat(IpImage_Gray->height,IpImage_Gray->width,CV_32FC1);//to be freed
		CvMat * Mat_GradientTheta = cvCreateMat(IpImage_Gray->height,IpImage_Gray->width,CV_32FC1);//to be freed
		
		for(int i=0;i<IpImage_Gray->height;i++)
			for(int j=0;j<IpImage_Gray->width;j++)
			{
				double gradx = cvGetReal2D(Mat_GrayGradient_Horz,i,j);
				double grady = cvGetReal2D(Mat_GrayGradient_Vert,i,j);
				cvSetReal2D(Mat_GradientMag,i,j,sqrt(pow(gradx,2) + pow(grady,2)));
				//Keep the gradient values in 0-180 degree range
				if(grady<0)
				{
					grady = grady*(-1);
					gradx = gradx*(-1);
				}
				cvSetReal2D(Mat_GradientTheta,i,j,(atan2(grady, gradx))*180/pi);
			}
		
		/*Declaration of the 2D array of pointers - each to the histogram of a block*/
		int ySteps = (OverlapFactor*((IpImage_Gray->height/BlockSize_h)-1))+1;
		int xSteps = (OverlapFactor*((IpImage_Gray->width/BlockSize_w)-1))+1;
		double *** Arr_BlockBinTotal = (double***) calloc(ySteps,sizeof(double**));//to be freed
		for(int i=0;i<ySteps;i++)
		{	
			//*(Arr_BlockBinTotal+i) = (int**) calloc(xSteps, sizeof(int*));
			Arr_BlockBinTotal[i] = (double**) calloc(xSteps, sizeof(double*));
			for(int j=0;j<xSteps;j++)
			{
				//*(*(Arr_BlockBinTotal+i)+j) = (int*) calloc(BinCount, sizeof(int));
				Arr_BlockBinTotal[i][j] = (double*) calloc(BinCount, sizeof(double));
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

		/*Write the HOG Features in the array*/
		/*concatenate the FV referenced in each point in 2-D Grid and send as whole FV*/
		
		for(int i=0; i<ySteps; i++)
			for(int j=0; j<xSteps; j++)
			{
				double sum = BlockBinTotal_Sum(Arr_BlockBinTotal[i][j]);
				double FeatureValue_norm;
				for(int k=0; k<BinCount; k++)
				{	
					FeatureValue_norm = (*((Arr_BlockBinTotal[i][j])+k))/(sum+0.0001);
					hog_feat[k] = FeatureValue_norm;
				}
			}

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
		cvReleaseMat(&Mat_GrayGradient_Horz);
		cvReleaseMat(&Mat_GradientMag);
		cvReleaseMat(&Mat_GradientTheta);
		Mat_Gray=NULL;
		Mat_GrayGradient_Vert=NULL;
		Mat_GrayGradient_Horz=NULL;
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
	void BlockBinTotal(CvMat* M_Magnitude, CvMat* M_Theta, int BlockSize_w, int BlockSize_h, int x, int y, double* EachBinTotal)
	{ 
		double theta;
		int groupno;
		int groupdiv = 180/BinCount;	
		for(int i=0; i<BlockSize_w; i++)
			for(int j=0; j<BlockSize_h; j++)
			{
				theta = cvGetReal2D(M_Theta,y+j,x+i);
				if(theta == 180)
					theta = 0;
				groupno = floor(theta/groupdiv);
				//double cee = cvGetReal2D(M_Magnitude,y+j,x+i);
				*(EachBinTotal+groupno) += cvGetReal2D(M_Magnitude,y+j,x+i); 
			}
	}


// Normalise the Histogram of each Block---------------------------------//
	double BlockBinTotal_Sum(double * inputptr)
	{
		double sum=0;
		for(int i=0;i<BinCount;i++)
		{
			sum += *(inputptr+i);
		}
		return (sum);
	}

	


//*******************************Scaling function**********************************//

void scale(double *feature,double *scaledfeature,struct range_model *r,int max_index)
{ 
	int i;
	for(i=0;i<max_index;i++)
	{
		if(r->feature_max[i+1]==r->feature_min[i])
		{
			scaledfeature[i]=-100;
			continue;
		}
		if(feature[i] == r->feature_min[i+1])
		{
			scaledfeature[i] = r->lower;
		}
		else if(feature[i] == r->feature_max[i+1])
		{
			scaledfeature[i] = r->upper;
		}
		else
		{
			scaledfeature[i] = r->lower + (r->upper-r->lower) * (feature[i]-r->feature_min[i+1])/(r->feature_max[i+1]-r->feature_min[i+1]);
		}
		if(scaledfeature[i]==0)
			scaledfeature[i]=-100;
	}
}

//**********************Replacement of main****************************************//

	void main_rep(IplImage* IpImage_Color, int* predict_label_seq)
	{
		/*------------------------------------------------*/
		/* Load the range file*/
		char range_name[200] = {0};
		strcpy(range_name,"D:\\SVM_training_tool\\SVM_Tools\\msd+hog_180\\fortesting\\range.txt");			//change		
		FILE* fp = fopen(range_name,"r");
		if(fp==NULL)
		{
			fprintf(stderr,"can't open range file %s\n", range_name);
			getch();
			exit(1);
		}
		int idx;
		double fmin, fmax;
		struct range_model* range;
		range = new struct range_model;
		range->feature_max =new double[MAX_INDEX+1];
		range->feature_min =new double[MAX_INDEX+1];
		for(int i=0;i<=MAX_INDEX;i++)
		{
			range->feature_max[i] = 0.00;
			range->feature_min[i] = 0.00;
		}
		if(range->feature_max == NULL || range->feature_min == NULL)
		{
			fprintf(stderr,"can't allocate enough memory\n");
			getch();
			exit(1);
		}
		if (fgetc(fp) == 'x') 
		{
			fscanf(fp, "\n%lf %lf\n", &range->lower, &range->upper);
			while(fscanf(fp,"%d %lf %lf\n",&idx,&fmin,&fmax)==3)
			{
				range->feature_min[idx] = fmin;
				range->feature_max[idx] = fmax;
			}
		}
		fclose(fp);
	/*---------------------------------------------------------------*/

	/*----------------------------------------------------------------*/
		//Load the model file
		char model_name[200];
		struct svm_model* model;

		strcpy(model_name,"D:\\SVM_training_tool\\SVM_Tools\\msd+hog_180\\fortesting\\train.model");
		if((model=svm_load_model_BIN(model_name))==0)//only for sirs model file with a flag for probability
		{
			fprintf(stderr,"can't open model file %s.\n",model_name);
			getch();
			exit(1);
		}
	/*---------------------------------------------------------------*/

		IplImage* IpImage_Resized = 0;
		IplImage* IpImage_HSV = 0;
		IplImage* HSV_qnt  = 0;

		int w = IpImage_Color->width;
		int h = IpImage_Color->height;
		int i , j;
		int part_w, part_h;
		int iHIstIndex;
		int Index;

		//feature vector
		double featvec[MAX_INDEX] = {0};
		double hog_feat[BinCount] = {0};
		double scaledfeatvec[MAX_INDEX] = {0};
		struct svm_node* featseq = (struct svm_node*)malloc(sizeof(struct svm_node)*(MAX_INDEX+1));	
		double prob_estimates[TOT_CLASSES];
		int predict_label;

		//Resizing
		//Fit the image into a max size (1000 : the lengthiest side)
		int maxside;
		if(w>h)
			maxside = w;
		else
			maxside = h;
        if (maxside > 600)
        {
            float as_ratio = (((float)maxside) / 600);

            h = (int)floor(h / as_ratio);
            w = (int)floor(w / as_ratio);
		}
		
		IpImage_Resized = cvCreateImage(cvSize(w,h),IPL_DEPTH_8U,3);
		cvResize(IpImage_Color, IpImage_Resized,CV_INTER_CUBIC);
		
		//Convert from RGB to HSV
		IpImage_HSV = cvCreateImage(cvSize(w,h), IPL_DEPTH_8U, 3);
		cvCvtColor(IpImage_Resized,IpImage_HSV,CV_BGR2HSV);
		
		//Quantize the HSV Values
		HSV_qnt = cvCreateImage(cvSize(w,h), IPL_DEPTH_8U, 1);
		HSV_color_quant(IpImage_HSV, HSV_qnt);
		//cvShowImage("quantized",HSV_qnt);cvWaitKey(0);

		part_w = w / 3;
        	part_h = h / 3;

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
		
		double* hist = (double*)calloc(CSA,sizeof(double)); // the features vector of MSD 
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

				image[0][i][j] = color.val[0]; %H-channel
				image[1][i][j] = color.val[1]; %S-channel
				image[2][i][j] = color.val[2]; %V-channel
				ImageX[i][j] = color_qnt.val[0];
				if (i < 2 * part_w && j<2*part_h) %Top-left
                		{
		                 image1[0][i][j] = color.val[0];
                		 image1[1][i][j] = color.val[1];
		                 image1[2][i][j] = color.val[2];
				 ImageX1[i][j] = color_qnt.val[0];
                		}
		                if ((i >= part_w && i < w) && (j < 2 * part_h)) %Top-right
                		{
		                 image2[0][i - part_w][j] = color.val[0];
                		 image2[1][i - part_w][j] = color.val[1];
		                 image2[2][i - part_w][j] = color.val[2];
				 ImageX2[i - part_w][j] = color_qnt.val[0];
                		}
		                if ((i < 2 * part_w) && (j >= part_h && j < h)) %Bottom-left
                		{
		                 image3[0][i][j - part_h] = color.val[0];
                		 image3[1][i][j - part_h] = color.val[1];
		                 image3[2][i][j - part_h] = color.val[2];
				 ImageX3[i][j - part_h] = color_qnt.val[0];
		                }
                		if ((i >= part_w && i < w) && (j >= part_h && j < h))	%Bottom-right
		                {
		                 image4[0][i - part_w][j - part_h] = color.val[0];
		                 image4[1][i - part_w][j - part_h] = color.val[1];
		                 image4[2][i - part_w][j - part_h] = color.val[2];
				 ImageX4[i - part_w][j - part_h] = color_qnt.val[0];
                		}
			}
		}

		// Call the function for MSD based Histogram evaluation
		// And Write MSD features to the feature vector, scale them, and then write to the structure array

		//0.
		MSD_feature_extract(image,ImageX,hist,w,h);
		count = 0;
		for (iHIstIndex = 0; iHIstIndex < CSA; iHIstIndex++)
		{	
			featvec[count++] = hist[iHIstIndex];
		}
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
		MSD_feature_extract(image1,ImageX1,hist1,2*part_w,2*part_h);
		for (iHIstIndex = 0; iHIstIndex < CSA; iHIstIndex++)
		{	
			featvec[count++] = hist1[iHIstIndex];
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
		MSD_feature_extract(image2,ImageX2,hist2,w-part_w,2*part_h);
		for (iHIstIndex = 0; iHIstIndex < CSA; iHIstIndex++)
		{	
			featvec[count++] = hist2[iHIstIndex];
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
		MSD_feature_extract(image3, ImageX3, hist3, 2*part_w, h-part_h);
		for (iHIstIndex = 0; iHIstIndex < CSA; iHIstIndex++)
		{	
			featvec[count++] = hist3[iHIstIndex];
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
		MSD_feature_extract(image4, ImageX4, hist4, w-part_w, h-part_h);
		for (iHIstIndex = 0; iHIstIndex < CSA; iHIstIndex++)
		{	
			featvec[count++] = hist4[iHIstIndex];
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
		
		//5.
		// Call the function for HOG descriptor generation
		HOG_FeatureEx_FileCreation(HSV_qnt, hog_feat);
		for (iHIstIndex = 0; iHIstIndex < BinCount; iHIstIndex++)
		{	
			featvec[count++] = hog_feat[iHIstIndex];
		}

		
		//Release Image Data
		cvReleaseImage(&IpImage_Resized);
		cvReleaseImage(&IpImage_HSV);
		cvReleaseImage(&HSV_qnt);
		
		IpImage_Resized = NULL;
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

		//Scale the feature vector	
		scale(featvec,scaledfeatvec,range,MAX_INDEX);
			
		//Put the feature values in the structure array
		Index = 0;
		for(j=0;j<MAX_INDEX;j++)
		{
			if(scaledfeatvec[j]!=-100)
			{
				featseq[Index].index=j+1;
				featseq[Index].value = scaledfeatvec[j];
				Index++;
			}
		}
		featseq[Index].index=-1;

		//Calling the svm_predict function to get the recognized class index
		//predict_label = (int)svm_predict(model, featseq);
		
		predict_label = (int)svm_predict_probability(model,featseq,prob_estimates);
		
		//sort the label_seq in the decreasing order of prob_estimates
		int label_seq[] = {1,10,100,101,102,103,104,105,106,107,108,109,11,110,111,112,113,114,115,116,117,118,119,12,120,121,122,123,124,125,126,127,128,129,13,130,131,132,133,134,135,136,137,138,139,14,140,141,142,143,144,145,146,147,148,149,15,150,151,152,153,154,155,156,157,158,159,16,160,161,162,163,164,165,166,167,168,169,17,170,171,172,173,174,175,176,177,178,179,18,180,181,182,183,184,185,186,187,188,189,19,190,191,192,193,194,195,196,197,198,199,2,20,200,201,202,203,204,205,206,207,208,209,21,210,211,212,213,214,215,216,217,218,219,22,220,221,222,223,224,23,24,25,26,27,28,29,3,30,31,32,33,34,35,36,37,38,39,4,40,41,42,43,44,45,46,47,48,49,5,50,51,52,53,54,55,56,57,58,59,6,60,61,62,63,64,65,66,67,68,69,7,70,71,72,73,74,75,76,77,78,79,8,80,81,82,83,84,85,86,87,88,89,9,90,91,92,93,94,95,96,97,98,99};
		double temp = 0;
		for (int i = 0; i < TOT_CLASSES - 1; i++)
        {
            for (int j = i + 1; j < TOT_CLASSES; j++)
            {
               if (prob_estimates[j] > prob_estimates[i])
               {
                  temp = prob_estimates[i];
                  prob_estimates[i] = prob_estimates[j];
                  prob_estimates[j] = temp;

				  temp = label_seq[i];
                  label_seq[i] = label_seq[j];
                  label_seq[j] = (int)temp;
                }
             }
         }
			
		//clear memory
		free(featseq);
		svm_destroy_model_bet(model);

		//Put the label sequence in the pass-by-ref variable
		for(int i = 0; i<10; i++)
			predict_label_seq[i] = label_seq[i];
		
	}
//*********************************************************************************//
//---------------------------------Main Function-----------------------------------//
//*********************************************************************************//
int main(void)
{
	int predict_label_seq[10];

	DIR *dir_m, *dir_sub;
	struct dirent *ent_m, *ent_sub;
	char dirpath_m[] = "G:\\flags_all\\mag_only\\flags_train";
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
							//Unable to load
							if(!IpImage_Color)
							{
								fprintf(stderr,"\nFailed to load input image!");
								getch();
								return -1;
							}
							
							//Call the function for predicting the class id of the image
							main_rep(IpImage_Color, predict_label_seq);
							cvReleaseImage(&IpImage_Color);
							IpImage_Color = NULL;

						}
					}
					closedir(dir_sub);
				}
			}
		}closedir(dir_m);
	}
	return 0;
}