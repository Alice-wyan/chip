#include"lib.h"
//!  a CommomFunction class
class  CommomFunction
{
public:
double getrmse(vector<double> xx)
{
        double sum=0;
		double ave=0;
	for(int a=0;a<xx.size();a++)//先算出平均值
	{

		ave+=xx[a];
	}

	for(int i=0;i<xx.size();i++)
	{
		
			//cout<<"ave== "<<ave/xx.size()<<endl;
			//cout<<"sum "<<ave/n<<endl;
			sum+=(xx[i]-(ave/xx.size()))*(xx[i]-(ave/xx.size()));	
	}
	
  double deta=sqrt((double)(sum/xx.size()));
  //cout<<"deta== "<<deta<<endl;
  return deta;
}

float * FitLine(IplImage* dst,vector<double> detectpoint)//通过点像素点坐标拟合直线
{


			float *liner = new float[4];///定义数组

			CvMemStorage* storage1 = cvCreateMemStorage(0);///创建存储空间

			CvSeq* point_seqr = cvCreateSeq( CV_32FC2,sizeof(CvSeq), sizeof(CvPoint2D32f), storage1 );///创建序列

			for(int s=0;s<detectpoint.size()-1;s=s+2)
			{
			cvSeqPush(point_seqr, &cvPoint2D32f(detectpoint[s],detectpoint[s+1]));///把像素点精确坐标push给point_seqr
			}


			cvFitLine(point_seqr,CV_DIST_HUBER ,0,0.001,0.001,liner);///直线拟合函数   x 方向
			//cvFitLine(point_seql,CV_DIST_L2,0,0.001,0.001,linel);///把ydetectpoint的像素点进行直线拟合   y 方向
		

			
			double x1=liner[2]; double y1=liner[3];/// line[0--3] 分别为 (vx, vy,x0, y0)


			//cout<<"liner[0]="<<liner[0]<<"  "<<liner[1]<<endl;
			//cout<<"linel[0]="<<linel[0]<<"  "<<linel[1]<<endl;
			
			//cout<<"x1="<<x1<<"  "<<y1<<endl;
			//cout<<"x2="<<x2<<"  "<<y2<<endl;
			//cout<<"angle1 "<<angle1<<"  x= "<<x1<<"   "<<y1<<endl;
			//cout<<"angle2 "<<angle2<<"  x=  "<<x2<<"   "<<y2<<endl;
			double k1=liner[1]/liner[0];double b11=liner[3]-k1*liner[2];


//.....................将拟合直线画出来.............start.............................//
	double cos_theta1 = liner[0];
    double sin_theta1 = liner[1];
    double x01 = liner[2], y01 = liner[3];




	double k11 = sin_theta1 / cos_theta1;
	double b1 = y01 - k11 * x01;
	double x11 = 0;
    double y11 = k11 * x11 + b1;
	Mat dd=dst;

	line(dd, Point(x01,y01), Point(x11,y11), cv::Scalar(255), 1);
	// DrawLine(dst, phi, rho, cv::Scalar(0));
//.....................将拟合直线画出来.............end.............................//


			cvReleaseMemStorage (&storage1);
			return liner;


}

vector<double> getcovergepoint(float *line1,float* line2)//求取拟合后两直线的交点
{
	        vector<double> xy;
			xy.clear();
	        //line[0--3] 分别为 (vx, vy,x0, y0)
	        double x1=line1[2]; double y1=line1[3];/// line[0--3] 分别为 (vx, vy,x0, y0)
			double x2=line2[2]; double y2=line2[3];

			//cout<<"liner[0]="<<liner[0]<<"  "<<liner[1]<<endl;
			//cout<<"linel[0]="<<linel[0]<<"  "<<linel[1]<<endl;
			
			//cout<<"x1="<<x1<<"  "<<y1<<endl;
			//cout<<"x2="<<x2<<"  "<<y2<<endl;
			//cout<<"angle1 "<<angle1<<"  x= "<<x1<<"   "<<y1<<endl;
			//cout<<"angle2 "<<angle2<<"  x=  "<<x2<<"   "<<y2<<endl;
			double k1=line1[1]/line1[0];double b11=line1[3]-k1*line1[2];
			double k2=line2[1]/line2[0];double b22=line2[3]-k2*line2[2];///确定直线斜率及直线方程
			xy.push_back((y2-y1+k1*x1-k2*x2)/(k1-k2));///利用两直线求角点			
			xy.push_back(k1*(xy[0]-x1)+y1);
   
            return xy;

       

}



};