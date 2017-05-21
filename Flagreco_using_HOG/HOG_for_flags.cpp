/*This code is used to generate the HOG Feature Vector for the Input Image (based on the paper by Dalal-2005)and
is followed by creation of SVM sample file as per the LIBSVM standard format.
*/

#include"HOG_for_flags.h"

/*No. of bins (defined for use in multiple functions*/
#define BinCount 72     
/*Constant for radian to degree conversion*/ 
#define pi 3.14159265
/*Function Declarations*/
void BlockBinTotal(CvMat*, CvMat*, int,int,int,int, int*);
int BlockBinTotal_Sum(int*);
void HOG_FeatureEx_FileCreation(IplImage*,char*,FILE*);
void HSV_color_quant(IplImage*, IplImage*);
void AddMat_Bin(CvMat*, CvMat*, CvMat*);

int main()
{
	/*Filename variables */
	char filename[261] = {}; 
	char fullpath[100];
	int failcount=0;
	/* LInk the file pointer to the sample file */
	FILE * fileptr = fopen("G:\\flags_all\\mag_only\\MSD_features_test_HOGreal.txt","w");
	if(fileptr == NULL)
		return(-1);
	/*Read all the files in the folder*/
	DIR *dir_m, *dir_sub;
	struct dirent *ent_m, *ent_sub;
	char dirpath_m[] = "G:\\flags_all\\mag_only\\flags_train";
    dir_m = opendir (dirpath_m);
	if (dir_m != NULL) 
	{
		/*Read the folders in the main directory*/
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
					while((ent_sub = readdir (dir_sub)) != NULL)
					{
						if(strcmp(".",ent_sub->d_name) && strcmp("..",ent_sub->d_name)  && strcmp("Thumbs.db",ent_sub->d_name))
						{
							/*Images being used*/
							IplImage * IpImage_Color = 0; //to be freed 
							IplImage * IpImage_HSV = 0; //to be freed
							IplImage* HSV_qnt  = 0;//to be freed
							IplImage* IpImage_Resized = 0;//to be freed

							strcpy(filename,fullpath);
							strcat(filename,"\\");
							strcat(filename,ent_sub->d_name);
							IpImage_Color = cvLoadImage(filename,CV_LOAD_IMAGE_COLOR);
							if(!IpImage_Color)
							{
								//fprintf(stderr,"failed to load input image\n");
								//return -1;
								failcount++;
								continue;
								
							}
							int w = IpImage_Color->width;
							int h = IpImage_Color->height;

							//Fit the image into a max size (1000 : the lengthiest side)
							int maxside;
							if(w>h)
								maxside = w;
							else
								maxside = h;
                            if (maxside > 1000)
                            {
                                int as_ratio = ceil(((float)maxside) / 1000);

                                h = h / as_ratio;
                                w = w / as_ratio;
							}

							IpImage_Resized = cvCreateImage(cvSize(w,h),IPL_DEPTH_8U,3);
							cvResize(IpImage_Color, IpImage_Resized,CV_INTER_CUBIC);

							/*Color quantize the input image - considering values of H, S and V*/
							IpImage_HSV = cvCreateImage(cvSize(w,h), IPL_DEPTH_8U, 3);
							HSV_qnt = cvCreateImage(cvSize(w,h), IPL_DEPTH_8U, 1);
							
							//Convert from RGB to HSV
							cvCvtColor(IpImage_Resized,IpImage_HSV,CV_BGR2HSV);
							
							//Quantize the HSV Values
							HSV_color_quant(IpImage_HSV, HSV_qnt);
							//cvShowImage("quantized",HSV_qnt);cvWaitKey(0);
							
							/*Extract the 'label' of the file*/
							//strcpy(filename,ent_sub->d_name);
							//char delimiter[] = "_"; int pos; char label[10]= "";
							char label[10]= "";
							//pos = strcspn(filename,delimiter);
							//strncpy(label,filename,pos);
							//label[pos]= '\0';
							strcpy(label,ent_m->d_name);
							
							/*call the function-- HOG Feature Extraction and creation of SVM sample file*/
							HOG_FeatureEx_FileCreation(HSV_qnt,label,fileptr);
							
							//Release Data
							cvReleaseImage(&IpImage_Resized);
							cvReleaseImage(&IpImage_Color);
							cvReleaseImage(&IpImage_HSV);
							cvReleaseImage(&HSV_qnt);
							IpImage_Resized = NULL;
							IpImage_Color = NULL;
							IpImage_HSV = NULL;
							HSV_qnt = NULL;

						}						
					}
					closedir(dir_sub);
				}
			}
		}
	closedir (dir_m);
	} 
	else 
	{
	/*could not open directory*/
	perror ("");
	return EXIT_FAILURE;	
	}
	printf("SUCCESSFUL!!");
	printf("\n%dFail Count", failcount);
	getch();
}

	

/*     ================================================================================                 */
//******************************     Functions' Definitions     ****************************************//
/*     ================================================================================                 */ 

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
                    s.val[0] = 0*30;
                    cvSet2D(HSV_qnt_mat, i, j,s);				
				}
                else
                    {
                        //Check for Saturation - if less than 16, threshold at 31 (Black and white)
                        if (SI < 16)
                        {
                            if (VI < 31)
							{
								s.val[0] = 0*30;
                                cvSet2D(HSV_qnt_mat,i, j,s);
							}
                            else
							{
								s.val[0] = 1*30;
                                cvSet2D(HSV_qnt_mat, i, j,s);
							}
                        }
                        //Check for Saturation - if more than 16, threshold at 16 (Black and Color)
                        else
                        {
                            if (VI < 16)
							{
								s.val[0] = 0*30;	
                                cvSet2D(HSV_qnt_mat, i, j,s);
							}
                            else
                            {
                                if ((HI >= 0 && HI < 30) || (HI > 330 && HI <= 359))
								{
									s.val[0] = 2*30;
                                    cvSet2D(HSV_qnt_mat, i, j,s);
								}
                                if (HI >= 30 && HI < 90)
								{
									s.val[0] = 3*30;
                                    cvSet2D(HSV_qnt_mat, i, j,s);
								}
                                if (HI >= 90 && HI < 150)
								{
                                    s.val[0] = 4*30;
                                    cvSet2D(HSV_qnt_mat, i, j,s);
								}
                                if (HI >= 150 && HI < 210)
								{
                                    s.val[0] = 5*30;
                                    cvSet2D(HSV_qnt_mat, i, j,s);
								}
                                if (HI >= 210 && HI < 270)
								{
                                    s.val[0] = 6*30;
                                    cvSet2D(HSV_qnt_mat, i, j,s);
								}
                                if (HI >= 270 && HI < 330)
								{
                                    s.val[0] = 7*30;
                                    cvSet2D(HSV_qnt_mat, i, j,s);
								}
                            }
                        }
                    }//end-of-main-if
			}
		}		
	}
	
//-------------Function to extract HOG Features for the Input Image and write it in a file---------------//
	

	void HOG_FeatureEx_FileCreation(IplImage * IpImage_Gray,char * label, FILE* fileptr)
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
				cvSetReal2D(Mat_GradientMag,i,j,sqrt(pow(gradx,2) + pow(grady,2)));
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
		int count = 0;
		fprintf(fileptr,"%s ",label);
		for(int i=0; i<ySteps; i++)
			for(int j=0; j<xSteps; j++)
			{
				int sum = BlockBinTotal_Sum(Arr_BlockBinTotal[i][j]);
				for(int k=0; k<BinCount; k++)
				{	
					float FeatureValue_norm = (float)(*((Arr_BlockBinTotal[i][j])+k))/((float)sum+0.0001);
					if(FeatureValue_norm!=0.0)
					{
						fprintf(fileptr,"%d:%f ",count++,FeatureValue_norm);
					}
					else
						count++;
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

//---------------------Function to Add the matrices and threshold the same--------------------------------//
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

				s.val[0] = (s1.val[0]+s2.val[0])/2;
				cvSet2D(Mat_tot,i,j,s);
			}
		}
	}


//---------------------Function to calculate Histogram of each specified block---------------------------//

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



//-------------------Function to Normalise the Histogram of each Block---------------------------------//

	int BlockBinTotal_Sum(int * inputptr)
	{
		float sum=0;
		for(int i=0;i<BinCount;i++)
		{
			sum += *(inputptr+i);
		}
		return (sum);
	}
