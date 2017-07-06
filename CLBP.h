/*
声明类以及变量
*/

#pragma once

#include "lib.h"

#define DIM_OF_LBP 256
using namespace std;

/*
输入： 1，图片名字
       2，图片指针 IplImage*
	   
输出： 深度是DIM_OF_LBP的double数组 outPut 

*/


class CLBP
{
public:
	CLBP(void);
	~CLBP(void);

private:
int bin2value( int a[], int n )//2进制转换为10进制
{
    int value = 0;
    for( int i=0; i<n; ++i )
        value |= a[i]<<i;
    return value;
}	
	
void ComputingPro(int r,int c,int **a)//求得数组LBP中的每个元素概率
{

	int count[DIM_OF_LBP];    //分配一个数组，求频数
	for(int t=0;t<DIM_OF_LBP;t++)
			count[t]=0;
	
	for(int i=0;i<DIM_OF_LBP;i++)
	{
		for(int k=0;k<r-2;k++)
		{
			for(int j=0;j<c-2;j++)
			{
				if(i==a[k][j])
					count[i]++;
			}
		}
		outPut[i]=count[i]/((double)((r-2)*(c-2)));//计算频率
	}
}	
void clear(){
	if (m_cvImageFlag==1) {cvReleaseImage(&m_cvImage); m_cvImageFlag=0;}
	if (m_cvImageRBGFlag==1) {cvReleaseImage(&m_cvImageRBG); m_cvImageRBGFlag=0;}

	if(LBPFlag==1){
	  for(int i=0;i<width;i++)
	  {
		  delete [] LBPData[i];
		  LBPData[i]=NULL;
	  }
	  delete [] LBPData;
	  LBPData=NULL;	
	  LBPFlag=0;
	}

	height=0;
	width=0;	
}

public:
void setup(string src_image_fl){
	clear();
	m_cvImageRBG = cvLoadImage(src_image_fl.c_str(), 1);// 读取图片
	m_cvImageRBGFlag=1;

	height=m_cvImageRBG->height;
	width=m_cvImageRBG->width;	

	m_cvImage = cvCreateImage(cvSize(m_cvImageRBG->width, m_cvImageRBG->height), 8, 1);
    m_cvImageFlag=1;

	LBPData=new int*[width];    //分配一个指针数组，将其首地址保存在LBP中  
	for(int i=0;i<width;i++)             //为指针数组的每个元素分配一个数组
	{LBPData[i]=new int[height];}
	LBPFlag=1;

	cvCvtColor(m_cvImageRBG, m_cvImage, CV_BGR2GRAY);

	
}
void setup(IplImage* temp_image){
	m_cvImageRBG=temp_image;	
	m_cvImageRBGFlag=0;
	
	if(height==m_cvImageRBG->height&&width==m_cvImageRBG->width){
	}else{
	clear();	
	m_cvImage = cvCreateImage(cvSize(m_cvImageRBG->width, m_cvImageRBG->height), 8, 1);
    m_cvImageFlag=1;	
	height=m_cvImageRBG->height;
	width=m_cvImageRBG->width;	

	LBPData=new int*[width];    //分配一个指针数组，将其首地址保存在LBP中  
	for(int i=0;i<width;i++)             //为指针数组的每个元素分配一个数组
	{LBPData[i]=new int[height];}
	LBPFlag=1;

	}
	cvCvtColor(m_cvImageRBG, m_cvImage, CV_BGR2GRAY);
	
}

void transform(){
	
	if(m_cvImageFlag==0) {
		cout<<"CLBP::transform m_cvImageFlag==0"<<endl;
		exit(0);
	}




	/*创建动态二维指针表示LBP特征*/


	int center=0;
	int center_lbp=0;
	for (int row=1; row<m_cvImage->width-1; row++)
	{
		for (int col=1; col<m_cvImage->height-1; col++)
		{
			int lbp[8]={0,0,0,0,0,0,0,0};//保存局部lbp特征值
			center = (int)cvGetReal2D(m_cvImage, row, col);
			center_lbp=0; 
			if (center >= cvGetReal2D(m_cvImage, row-1, col-1))
			{
				center_lbp= 1;
				lbp[0]=center_lbp;
			}
			else 
			{
				center_lbp =0;
				lbp[0]=center_lbp;
			}

			if (center >= cvGetReal2D(m_cvImage, row-1, col))
			{
				center_lbp= 1;
				lbp[1]=center_lbp;
			}
			else 
			{
				center_lbp =0; 
				lbp[1]=center_lbp;
			}

			if (center >= cvGetReal2D(m_cvImage, row-1, col+1))
			{
				center_lbp= 1;
				lbp[2]=center_lbp;	
			}
			else 
			{
				center_lbp =0; 
				lbp[2]=center_lbp;
			}

			if (center >= cvGetReal2D(m_cvImage, row, col-1))
			{
				center_lbp= 1;
				lbp[3]=center_lbp;
			}
			else 
			{
				center_lbp = 0; 
				lbp[3]=center_lbp;	
			}

			if (center >= cvGetReal2D(m_cvImage, row, col+1))
			{
				center_lbp= 1;
				lbp[4]=center_lbp;
			}
			else 
			{
				center_lbp = 0; 
				lbp[4]=center_lbp;
			}

			if (center >= cvGetReal2D(m_cvImage, row+1, col-1))
			{
				center_lbp= 1;
				lbp[5]=center_lbp;	
			}
			else 
			{
				center_lbp = 0; 
				lbp[5]=center_lbp;
			}

			if (center >= cvGetReal2D(m_cvImage, row+1, col))
			{
				center_lbp=1;
				lbp[6]=center_lbp;
			}
			else 
			{
				center_lbp =0;
				lbp[6]=center_lbp;
			}

			if (center >= cvGetReal2D(m_cvImage, row+1, col+1))
			{
				center_lbp= 1;
				lbp[7]=center_lbp;	
			}
			else 
			{
				center_lbp = 0;
				lbp[7]=center_lbp;
			}
			if(row-1<0||row-1>=m_cvImageRBG->width)cout<<row-1<<"  x "<<m_cvImageRBG->width<<endl;
			if(col-1<0||col-1>=m_cvImageRBG->height)cout<<col-1<<" y  "<<m_cvImageRBG->height<<endl;
			LBPData[row-1][col-1]=bin2value(lbp,sizeof(lbp)/sizeof(lbp[0]));
		}
	}

	ComputingPro(m_cvImage->width,m_cvImage->height,LBPData);

	

	
}

void print(){
	cout<<"输出LBP特征值：（概率）"<<endl;
	for(int t=0;t<DIM_OF_LBP;t++)
	{
		cout<<outPut[t]<<" ";
		if((t+1)%15==0)
			cout<<endl;
	}	
}

	void	LBP(string src_image_fl);
	void   LBP(IplImage* temp_image);

private:
	IplImage*			m_cvImageRBG;
	IplImage*			m_cvImage;
	int m_cvImageRBGFlag;	
	int m_cvImageFlag;	
	int height;
	int width;
	int **LBPData;
	int LBPFlag;
public:
	double outPut[DIM_OF_LBP];

};




/*
对LBP特征向量进行提取的步骤
（1）首先将检测窗口划分为N×N的小区域（cell）；（N 的大小主要是看图片的像素，也就是说一个cell就是一个像素）
（2）对于cell中的像素，将相邻的8个像素的灰度值与其进行比较，若周围像素值大于中心像素值，则该像素点的位置被标记为1，否则为0。
		这样，3*3邻域内的8个点经比较可产生8位二进制数（可转换为10进制数），即得到该窗口中心像素点的LBP值；
（3）然后计算每个cell的直方图，即每个数字（假定是十进制数LBP值）出现的频率；然后对该直方图进行归一化处理（此处就是求概率）。
（4）最后将得到的每个cell的统计直方图进行连接成为一个特征向量，也就是整幅图的LBP纹理特征向量；
*/


CLBP::CLBP(void)
{
	m_cvImageFlag=0;
	m_cvImageRBGFlag=0;
	height=0;
	width=0;
}
CLBP::~CLBP(void)
{
clear();
}





void CLBP::LBP(string src_image_fl){
	IplImage* temp_image = cvLoadImage(src_image_fl.c_str(), 1);// 读取图片
	
	LBP(temp_image);
	cvReleaseImage(&temp_image);
}


void CLBP::LBP(IplImage* temp_image)
{
	
	m_cvImage = cvCreateImage(cvSize(temp_image->width, temp_image->height), 8, 1);

	if (temp_image->nChannels == 3)
	{
		cvCvtColor(temp_image, m_cvImage, CV_BGR2GRAY);
	}

	
	int r=m_cvImage->width,c=m_cvImage->height;


	/*创建动态二维指针表示LBP特征*/
	int **LBP=new int*[r-2];    //分配一个指针数组，将其首地址保存在LBP中  
	for(int i=0;i<r-2;i++)             //为指针数组的每个元素分配一个数组
		LBP[i]=new int[c-2];

	int center=0;
	int center_lbp=0;
	for (int row=1; row<c-1; row++)
	{
		for (int col=1; col<r-1; col++)
		{
			int lbp[8]={0,0,0,0,0,0,0,0};//保存局部lbp特征值
			center = (int)cvGetReal2D(m_cvImage, row, col);
			center_lbp=0; 
			if (center >= cvGetReal2D(m_cvImage, row-1, col-1))
			{
				center_lbp= 1;
				lbp[0]=center_lbp;
			}
			else 
			{
				center_lbp =0;
				lbp[0]=center_lbp;
			}

			if (center >= cvGetReal2D(m_cvImage, row-1, col))
			{
				center_lbp= 1;
				lbp[1]=center_lbp;
			}
			else 
			{
				center_lbp =0; 
				lbp[1]=center_lbp;
			}

			if (center >= cvGetReal2D(m_cvImage, row-1, col+1))
			{
				center_lbp= 1;
				lbp[2]=center_lbp;	
			}
			else 
			{
				center_lbp =0; 
				lbp[2]=center_lbp;
			}

			if (center >= cvGetReal2D(m_cvImage, row, col-1))
			{
				center_lbp= 1;
				lbp[3]=center_lbp;
			}
			else 
			{
				center_lbp = 0; 
				lbp[3]=center_lbp;	
			}

			if (center >= cvGetReal2D(m_cvImage, row, col+1))
			{
				center_lbp= 1;
				lbp[4]=center_lbp;
			}
			else 
			{
				center_lbp = 0; 
				lbp[4]=center_lbp;
			}

			if (center >= cvGetReal2D(m_cvImage, row+1, col-1))
			{
				center_lbp= 1;
				lbp[5]=center_lbp;	
			}
			else 
			{
				center_lbp = 0; 
				lbp[5]=center_lbp;
			}

			if (center >= cvGetReal2D(m_cvImage, row+1, col))
			{
				center_lbp=1;
				lbp[6]=center_lbp;
			}
			else 
			{
				center_lbp =0;
				lbp[6]=center_lbp;
			}

			if (center >= cvGetReal2D(m_cvImage, row+1, col+1))
			{
				center_lbp= 1;
				lbp[7]=center_lbp;	
			}
			else 
			{
				center_lbp = 0;
				lbp[7]=center_lbp;
			}
			LBP[row-1][col-1]=bin2value(lbp,sizeof(lbp)/sizeof(lbp[0]));
		}
	}

	ComputingPro(r,c,LBP);
	cout<<"输出LBP特征值：（概率）"<<endl;
	for(int t=0;t<DIM_OF_LBP;t++)
	{
		cout<<outPut[t]<<" ";
		if((t+1)%15==0)
			cout<<endl;
	}
	
	/*释放二维动态数组LBP*/
	  for(int i=0;i<r-2;i++)
	  {
		  delete [c-2]LBP[i];
		  LBP[i]=NULL;
	  }
	  delete [r-2]LBP;
	  LBP=NULL;

	 system("pause");

	cvWaitKey(0);

} 
