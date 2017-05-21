#include"flags_v1_header.h"

using namespace std;

void MakeBlank(int, IplImage*);
void AddEltWoDup(char[], char, float);
void LabelCount(CvMat*, float*);
char* ColorSetAlongLine(CvArr* , CvPoint, CvPoint);
void DuplicateFeatureSet(char[223][41]);
int minimum(int a, int b, int c);
//int levenshtein_distance(char s[31], char t[41]);
int levenshtein_distance(char *s, char *t);
int getRecognizedFlag(char featureSet[41]);
void readFeatures();
char allFlagFeatures[223][41];
int main()
{	
	//Feature Set
	char* ColorSet;
	char FeatureSet[41]={};
	/*File for Writing the Color Sets*/
	FILE* filesave = fopen("--path_to_file--\\flags\\Test.txt","w");
	if(filesave==NULL)
	{
		return(9);
	}
	char WholeFeatureSet[223][41];
	int count = 0;
	/*Filename variables */
	char filename[261] = {}; 
	char fullpath[100];
	/*Image being read (Colored)*/
	IplImage * IpImage_flag = 0; 
	IplImage * IpImage_flag_r = 0;
	IplImage * IpImage_flag_g = 0;
	IplImage * IpImage_flag_b = 0;
	/*Read all the files in the folder*/
	DIR *dir;
	struct dirent *ent;
	char dirpath[] = "D:\\Test data";
    dir = opendir (dirpath);
	readFeatures();
	if (dir != NULL) 
	{
		/*Fetch the image in the folder*/
		while ((ent = readdir (dir)) != NULL) 
		{
			if(strcmp(".",ent->d_name) && strcmp("..",ent->d_name) && strcmp("Thumbs.db",ent->d_name))
			{
				strcpy(fullpath,"");
				//strcat(fullpath,"--path_to_data--\\flags_large\\");
				strcat(fullpath,"D:\\Test data\\");
				IpImage_flag = cvLoadImage(strcat(fullpath,ent->d_name),CV_LOAD_IMAGE_COLOR);
				/*strcat(fullpath,"--path_to_data--\\flags_style1_large\\Czech_Republic.png");
				IpImage_flag = cvLoadImage(fullpath,CV_LOAD_IMAGE_COLOR);*/
				if(!IpImage_flag)
				{
					fprintf(stderr,"failed to load input image\n");
					return -1;
				}
			
			//Convert from RGB to HSV
			IplImage* IpImage_flag_HSV = cvCreateImage(cvGetSize(IpImage_flag),IPL_DEPTH_8U,3);
			cvCvtColor(IpImage_flag, IpImage_flag_HSV,CV_BGR2HSV);

			// IMage Dimensions
			float w = IpImage_flag->width;
			float h = IpImage_flag->height;

			////Along the diagonal
			//CvPoint pt1 = cvPoint(0,0);
			//CvPoint pt2 = cvPoint(w-1, h-1);
 		//	ColorSet = ColorSetAlongLine(IpImage_flag_HSV, pt1, pt2);
			//strcpy(FeatureSet, ColorSet);
			
			////Along the 3/5th above-diagonal
			//CvPoint pt1 = cvPoint(0,0);
			//CvPoint pt2 = cvPoint(w-1, (int)(2*(h-1))/5);
 		//	ColorSet = ColorSetAlongLine(IpImage_flag_HSV, pt1, pt2);
			//strcpy(FeatureSet, ColorSet);

			////Along the 3/5th below-diagonal
			//pt1 = cvPoint(0,0);
			//pt2 = cvPoint((int)(2*(w-1))/5, h-1);
 		//	ColorSet = ColorSetAlongLine(IpImage_flag_HSV, pt1, pt2);
			//strcat(FeatureSet, ColorSet);

			////Along the 2/5th below-diagonal
			//pt1 = cvPoint(0,0);
			//pt2 = cvPoint((int)(3*(w-1))/5, h-1);
 		//	ColorSet = ColorSetAlongLine(IpImage_flag_HSV, pt1, pt2);
			//strcat(FeatureSet, ColorSet);

			////Along the 2/5th below-diagonal
			//pt1 = cvPoint(0,0);
			//pt2 = cvPoint((int)(3*(w-1))/5, h-1);
 		//	ColorSet = ColorSetAlongLine(IpImage_flag_HSV, pt1, pt2);
			//strcat(FeatureSet, ColorSet);
			
			//Along the Horizontal Line (a little below the upper boundary)
			CvPoint pt1 = cvPoint(0,2);
			CvPoint pt2 = cvPoint(w-1,2);	
			ColorSet = ColorSetAlongLine(IpImage_flag_HSV, pt1, pt2);
			strcpy(FeatureSet, ColorSet);
			strcat(FeatureSet, " ");
			//fprintf(filesave,"%d \t %s\n", flagNo,FeatureSet);
			
			//Along the Horizontal Line (a little below the middle-1 boundary)
			pt1 = cvPoint(0,(int)((h-1)/4)+2);
			pt2 = cvPoint(w-1,(int)((h-1)/4)+2);	
			ColorSet = ColorSetAlongLine(IpImage_flag_HSV, pt1, pt2);
			strcat(FeatureSet, ColorSet);
			strcat(FeatureSet, " ");
			//fprintf(filesave,"%d \t %s\n", flagNo,FeatureSet);

			//Along the Horizontal Line (a little below the middle-2 boundary)
			pt1 = cvPoint(0,(int)((h-1)/2)+2);
			pt2 = cvPoint(w-1,(int)((h-1)/2)+2);	
			ColorSet = ColorSetAlongLine(IpImage_flag_HSV, pt1, pt2);
			strcat(FeatureSet, ColorSet);
			strcat(FeatureSet, " ");
			//fprintf(filesave,"%d \t %s\n", flagNo,FeatureSet);

			//Along the Horizontal Line (a little below the middle-3 boundary)
			pt1 = cvPoint(0,(int)(3*(h-1)/4)+2);
			pt2 = cvPoint(w-1,(int)(3*(h-1)/4)+2);	
			ColorSet = ColorSetAlongLine(IpImage_flag_HSV, pt1, pt2);
			strcat(FeatureSet, ColorSet);
			strcat(FeatureSet, " ");
			//fprintf(filesave,"%d \t %s\n", flagNo,FeatureSet);

			//Along the Horizontal Line (a little above the lower boundary)
			pt1 = cvPoint(0,h-3);
			pt2 = cvPoint(w-1,h-3);	
			ColorSet = ColorSetAlongLine(IpImage_flag_HSV, pt1, pt2);
			strcat(FeatureSet, ColorSet);
			strcat(FeatureSet, " ");
			//fprintf(filesave,"%d \t %s\n", flagNo,FeatureSet);

			//Along the vertical
			pt1 = cvPoint((int)((w-1)/2),0);
			pt2 = cvPoint((int)((w-1)/2),h-1);
 			ColorSet = ColorSetAlongLine(IpImage_flag_HSV, pt1, pt2);
			strcat(FeatureSet, ColorSet);
			strcat(FeatureSet, " ");
			//fprintf(filesave,"%d \t %s\n", flagNo,FeatureSet);

			//Along the vertical
			pt1 = cvPoint((int)((w-1)/10),0);
			pt2 = cvPoint((int)((w-1)/10),h-1);
 			ColorSet = ColorSetAlongLine(IpImage_flag_HSV, pt1, pt2);
			strcat(FeatureSet, ColorSet);

			int flagNo = getRecognizedFlag(FeatureSet);
			//fprintf(filesave,"\n%s",ent->d_name);
			
			//fprintf(filesave,"%d\n",flagNo);
			printf("\t%d\n",flagNo);
			//strcpy(WholeFeatureSet[count++],FeatureSet);

			}
		}
		fclose(filesave);
		getch();
		//DuplicateFeatureSet(WholeFeatureSet);
		closedir (dir);
	} 
	else 
	{
	/*could not open directory*/
	perror ("");
	return EXIT_FAILURE;	
	}
}


/*--------------- To find the color set along a given line---------------- */
char* ColorSetAlongLine(CvArr* IpImage_flag_HSV, CvPoint pt1, CvPoint pt2)
{
	char ColorSet[11] = {};
	/* Iteratively fetch the pixel values along the line and fill the "ColorSet" array */
	float B=0,M=0,W=0,R=0,Y=0,G=0,I=0,C=0,P=0;
	//To iteratively fetch the pixel values along a line
	float max_buffer;
	CvLineIterator iterator;
	max_buffer = cvInitLineIterator(IpImage_flag_HSV,pt1,pt2,&iterator,8,0);
	// for cvKMeans2 KC
	int CountValidElts = 0;
	CvMat* BufferBigger = cvCreateMat(max_buffer,3,CV_32FC1);
				
	for(int j=0; j<max_buffer; j++)
	{
		// HSV Values
	/*	printf("%d,", iterator.ptr[0]); 
		printf("%d,", iterator.ptr[1]); 
		printf("%d\n", iterator.ptr[2]);*/
			
		int h = iterator.ptr[0];	
		int s = iterator.ptr[1];	
		int v = iterator.ptr[2];

		if(s<13)
		{ // grey scale image 
			if(v<120){
				//add black (0)- if not already present
				B++;
				AddEltWoDup(ColorSet, 'B', B/max_buffer);}
			//else if(v<154){
			//	//add grey	(1)- if not already present
			//	AddEltWoDup(ColorSet, 'M');
			//	M++;}
			else{
				// add white(2)- if not already present
				W++;
				AddEltWoDup(ColorSet, 'W', W/max_buffer);}
		}

		else if(s>13)
		{
			// color image
			if(v<26){
				//add black (0)- if not already present
				B++;
				AddEltWoDup(ColorSet, 'B', B/max_buffer);}
			else
			{
				//valid for hue based separation
				CvScalar HsvVals;
				
				HsvVals.val[0] = h;
				cvSet2D(BufferBigger,CountValidElts,0,HsvVals);
				
				HsvVals.val[0] = s;
				cvSet2D(BufferBigger,CountValidElts,1,HsvVals);
				
				HsvVals.val[0] = v;
				cvSet2D(BufferBigger,CountValidElts,2,HsvVals);
				
				CountValidElts += 1;
			}
		}

		CV_NEXT_LINE_POINT(iterator); //Step to the next pixel
	}// end for
	
	// Hue based separation
	CvScalar PixelVal;int h;
	int dimset=1;int nClusters = 10;
	CvMat * Buffer = cvCreateMat(CountValidElts,dimset,CV_32FC1);
	CvMat* Centres = cvCreateMat(nClusters,dimset,CV_32FC1);
	CvMat * Labels = cvCreateMat(CountValidElts,1,CV_32S);
	for(int i=0;i<CountValidElts;i++)
		for(int j=0;j<dimset;j++)
		{
			PixelVal = cvGet2D(BufferBigger,i,j);
			cvSet2D(Buffer,i,j,PixelVal);
		}
	if(CountValidElts>10)
	{
		cvKMeans2(Buffer, nClusters, Labels,cvTermCriteria(CV_TERMCRIT_EPS,20,(double)0.5),1,0,0, Centres);		
		float * lblcnt = (float*)malloc(nClusters*sizeof(float));
		for(int i=0; i<nClusters; i++)
			lblcnt[i]=0;
		LabelCount(Labels, lblcnt);
		for(int i=0; i<nClusters; i++)
		{
			PixelVal = cvGet2D(Centres,i,0);
			h = PixelVal.val[0];
							
			if(h>=0 && h<15 || h>=165 && h<180){
				// add Red (3)- if not already present 
				R += (int)lblcnt[i];
				AddEltWoDup(ColorSet, 'R', R/max_buffer);}
			else if(h>=15 && h<45){
				// add Yellow (4)- if not already present
				Y += (int)lblcnt[i];
				AddEltWoDup(ColorSet, 'Y', Y/max_buffer);}
			else if(h>=45 && h<75){
				// add Green (5)- if not already present
				G += (int)lblcnt[i];
				AddEltWoDup(ColorSet, 'G', G/max_buffer);}
			else if(h>=75 && h<105){
				// add Cyan (6)- if not already present
				C += (int)lblcnt[i];
				AddEltWoDup(ColorSet, 'C', C/max_buffer);}
			else if(h>=105 && h<135){
				// add Blue (7)- if not already present
				I += (int)lblcnt[i];
				AddEltWoDup(ColorSet, 'I', I/max_buffer);}
			else if(h>=135 && h<165){
				// add Magenta (8)- if not already present
				P += (int)lblcnt[i];
				AddEltWoDup(ColorSet, 'P', P/max_buffer);}
		}
	}	
	//Order the color set//

	string clrseq = ColorSet;
	sort(clrseq.begin(), clrseq.end());
	strcpy(ColorSet, clrseq.c_str());
	return(ColorSet);

	/*cvShowImage("flag image",IpImage_flag);
	cvWaitKey(0);*/
}

void MakeBlank(int lmt, IplImage* Image)
{
	CvScalar Blank;
	Blank.val[0]=0;

	for(int i=Image->width-lmt;i<Image->width;i++)
		for(int j=0;j<Image->height;j++)
				cvSet2D(Image,j,i,Blank);

	for(int i=0;i<lmt;i++)
		for(int j=0;j<Image->height;j++)
				cvSet2D(Image,j,i,Blank);

	for(int i=0;i<Image->width;i++)
		for(int j=0;j<lmt;j++)
				cvSet2D(Image,j,i,Blank);

	for(int i=0;i<Image->width;i++)
		for(int j=Image->height-lmt;j<Image->height;j++)
				cvSet2D(Image,j,i,Blank);
}

void AddEltWoDup(char CharArray[], char elt, float count)
{
	int len = strlen(CharArray);
	int flag=0;
	if(!(count<0.02))
	{		
		for(int i=0; i<=len;i++)
		{
			if(CharArray[i]==elt)
				flag=1;
		}
		if(!flag)
			CharArray[len] = elt;
	}
}

void LabelCount(CvMat* Label, float * lblcnt)
{
	CvSize dim = cvGetSize(Label);
	CvScalar PixelVal;
	for(int i=0; i<dim.height; i++)
	{
		for(int j=0; j<dim.width; j++)
		{
			PixelVal = cvGet2D(Label, i, j);
			lblcnt[(int)PixelVal.val[0]]++;

		}
	}
}

//void DuplicateFeatureSet(char features[223][41])
//{
//	int ret=-1;
//	int check[223] = {0};
//	filesave = fopen("--path_to_file--\\ColorSetData.txt","a");
//	for(int i=0;i<223;i++)
//	{
//		if(check[i]==1)
//			continue;
//		fprintf(filesave,"\n%d",i+1);
//		for(int j=i+1;j<223;j++)
//		{	
//			ret = strcmp(features[i],features[j]);
//			if(!ret)
//			{
//				fprintf(filesave , "----%d",j+1);
//				check[j] = 1;
//			}
//		}
//
//	}
//	fclose(filesave);
//}

int minimum(int a, int b, int c)
{
    int min = a;
    if (b < min)
        min = b;
    if (c < min)
        min = c;
    return min;
}
//int levenshtein_distance(char s[41], char t[41])
//{
//    //Step 1
//    int k, i, j, n, m, cost, distance;
//    int *d;
//    n = strlen(s);
//    m = strlen(t);
//    if (n != 0 && m != 0)
//    {
//        d=(int*)malloc((sizeof(int))*(m+1)*(n+1));
//        //d = new int[(m + 1) * (n + 1)];
//        m++;
//        n++;
//        //Step 2	
//        for (k = 0; k < n; k++)
//            d[k] = k;
//        for (k = 0; k < m; k++)
//            d[k * n] = k;
//        //Step 3 and 4	
//        for (i = 1; i < n; i++)
//            for (j = 1; j < m; j++)
//            {
//                //Step 5
//                if (s[i - 1] == t[j - 1])
//                    cost = 0;
//                else
//                    cost = 2;
//                //Step 6			 
//                d[j * n + i] = minimum(d[(j - 1) * n + i] + 1, d[j * n + i - 1] + 1, d[(j - 1) * n + i - 1] + cost);
//            }
//        distance = d[n * m - 1];
//        //free(d);
//        return distance;
//    }
//    else
//        return -1; //a negative return value means that one or both strings are empty.
//}

int levenshtein_distance(char *s, char *t)
{
    //Step 1
    int k, i, j, n, m, cost, distance;
    int *d;
    n = strlen(s);
    m = strlen(t);
    if (n != 0 && m != 0)
    {
        d=(int*)malloc((sizeof(int))*(m+1)*(n+1));
        //d = new int[(m + 1) * (n + 1)];
        m++;
        n++;
        //Step 2	
        for (k = 0; k < n; k++)
            d[k] = k;
        for (k = 0; k < m; k++)
            d[k * n] = k;
        //Step 3 and 4	
        for (i = 1; i < n; i++)
            for (j = 1; j < m; j++)
            {
                //Step 5
                if (s[i - 1] == t[j - 1])
                    cost = 0;
                else
                    cost = 2;
                //Step 6			 
                d[j * n + i] = minimum(d[(j - 1) * n + i] + 1, d[j * n + i - 1] + 1, d[(j - 1) * n + i - 1] + cost);
            }
        distance = d[n * m - 1];
        //free(d);
        return distance;
    }
    else
        return -1; //a negative return value means that one or both strings are empty.
}


void readFeatures()
{
	int count = 0;
	int countFeatures = 0;
	char colorset[31];
	//char *directionFeatureSet
	FILE *fr = fopen("--path_to_file--\\ColorSetEachDirection.txt","rb");
	int noOfFeatures = 0;
	while(fgets(colorset,sizeof(colorset), fr) != NULL)
	{
		if (count == 0)
		{
			noOfFeatures = atoi(colorset);
			count = 1;
		}
		else
		{
			int lengthOfFeatureSet = strlen(colorset);
			for(int i = 0; i < lengthOfFeatureSet;i++)
			{
				allFlagFeatures[countFeatures][i] = colorset[i];
			}
			
			countFeatures++;
		}
	}
	fclose(fr);
}

int getRecognizedFlag(char testFeatureSet[41])
{
	int arrEditDistance[223];
	int minDistance = 100;
	int minIndex = 0;
	for(int i = 0; i < 223;i++)
	{
		char *directionTestFeatureSet[7];
		char *temp;
		char *directionReferenceFeatureSet;
		int count = 0;
		temp = strtok(testFeatureSet," \n");
		while(temp != NULL)
		{
			directionTestFeatureSet[count] = temp;
			temp = strtok(NULL," \n");
			count++;
		}
		directionReferenceFeatureSet = strtok(allFlagFeatures[i]," \n");
		count = 0;
		while(directionReferenceFeatureSet != NULL)
		{
			int distance = levenshtein_distance(allFlagFeatures[i],directionTestFeatureSet[count]);
			arrEditDistance[i] += distance;
			/*if(arrEditDistance[i] < minDistance)
			{
				minDistance = arrEditDistance[i];
				minIndex = i;
			}*/
			//directionTestFeatureSet = strtok(NULL," \n");
			directionReferenceFeatureSet = strtok(NULL," \n");
		}
		if(arrEditDistance[i] < minDistance)
		{
			minDistance = arrEditDistance[i];
			minIndex = i;
		}
	}
	return minIndex;
}
