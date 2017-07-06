//
//linetocorner函数旨在利用houghlinesp函数检测图片直线后得出直线两端端点，
//利用端点坐标在两端点之间画探测线，通过探测线获取的 边缘像素坐标直线拟合，通过拟合后的直线求取出两直线的交点
//通过交点和利用直线的角度再设置探测线进行边缘像素点的探测，后再进行直线拟合后求取出精确坐标的交点。
//


#include"lib.h"

class LineToDetect{

public:


void linetocorner(IplImage *src,IplImage *dst,vector<Vec4i> lines)//linetocorner函数输入原图，输入划线后便于观察的图，和检测的直线参数
	{

		TestDetectPoint  tdp;//调用设置探测线位置的类
		
		CommomFunction cf;
		double degree1=0;double degree2=0;
		
		
		vector<double> b;//存储直线截距
		b.clear();
	
		vector<float *> line;//存储拟合直线参数
		line.clear();

		vector<float *> linee;//存储拟合直线参数
		linee.clear();
		vector<double> degree;
		degree.clear();
		int distance=0;	//拟合直线的长度
		for(int i=0;i<lines.size();i++)//循环直线的数量
		{
			double k=0; double degree=0; 
			Vec4i l= lines[i]; 
						
			CvPoint ll;//x坐标小的直线端点
			CvPoint lr; //x坐标大的直线端点   当x坐标相等时，则比较y坐标
			if(l[0]<l[2]){ ll.x=l[0];ll.y=l[1];lr.x=l[2];lr.y=l[3];} //端点ll 必须是直线上x坐标小的端点  这样会使角度范围在(-90~90)
			else if(l[0]==l[2]){if(l[1]<l[3]){ll.x=l[0];ll.y=l[1];lr.x=l[2];lr.y=l[3];}else{ll.x=l[2];ll.y=l[3];lr.x=l[0];lr.y=l[1];}}
			else{ ll.x=l[2];ll.y=l[3];lr.x=l[0];lr.y=l[1];}

			 if(ll.x==lr.x){ degree=-90;}			
			 else
			 { k=(l[1]-l[3])/(l[0]-l[2]);
			   if(k==0){degree=0;}
			   else{degree=atan((double)(k))/CV_PI*180.0;}}//degree 为直线角度
			
			
			 distance=sqrt((double)((l[0]-l[2])*(l[0]-l[2])+(l[1]-l[3])*(l[1]-l[3])));//两端点之间的长度

			// line.push_back( tdp.getdetcet(src,dst,ll,150,degree));//调用探测线利用探测线求取边缘像素点并拟合直线
	
		}

		//if(line.size()==2)
		//{
		//  vector<double> xy=cf.getcovergepoint(line[0],line[1]);//获取两拟合直线的交点


		//     fstream ph;

	 //        ph.clear();

		//     ph.open("cx.txt",ios::app);
		//	 fstream phy;

	 //        phy.clear();

		//     phy.open("cy.txt",ios::app);

		//	 CvPoint p;
		//	 p.x=xy[0];p.y=xy[1];
		//	// cvCircle(dst,p,10, CV_RGB(0,255,0), 1,8,0);//在dst上画点

		//	 ph<<xy[0]<<endl;
		//	 phy<<xy[1]<<endl;

		//
		//for(int d=0;d<line.size();d++)
		//{
		//	double k=line[d][1]/line[d][0];
		//	double deg=atan((double)(k))/CV_PI*180.0;
		//	cout<<"degree=="<<deg<<endl;
		//	degree.push_back(deg);
		//}

		//CvPoint point;
		//point.x=xy[0];
		//point.y=xy[1];
		////point.x=205;
		////point.y=400;
		//cout<<"pp=="<<xy[0]<<"  "<<point.x<<"  "<<point.y<<endl;
		////cvCircle(dst,point,10, CV_RGB(0,255,0), 1,8,0);//在dst上画交点--绿色
		//double tmp=0;
		//if(abs(degree[0])<abs(degree[1])){tmp=degree[0];degree[0]=degree[1];degree[1]=tmp;}//degree[0]是角度比较大的角度
  //    
		////cout<<"degree[0]="<<degree[0]<<"  "<<degree[1]<<endl;

		//tdp.getdetectedge(src,dst,point,distance-40,degree[0],degree[1]);
		//}
		







		//for(int j=0;j<sitaa.size();j++)
		//{
		//	double b1=lines[j][1]-sitaa[j]*lines[j][0];

		//	b.push_back(b1);
		//}

		//double x=(b[1]-b[0])/(sitaa[0]-sitaa[1]);//交代x
		//double y=sitaa[0]*x+b[0];//交点y
		//

		//cout<<"sita="<<sitaa[0]<<"  "<<sitaa[1]<<endl;
		//
		//cout<<"x=="<<x<<"  "<<"y=="<<y<<endl;
		//
		//CvPoint xy;
		//xy.x=x;
		//xy.y=y;

		////cvCircle(dst,xy,10, CV_RGB(0,255,0), 1,8,0);//在dst上画点
		//
		//degree1=atan((double)(1/sitaa[0]))/CV_PI*180;
		//degree2=atan((double)(1/sitaa[1]))/CV_PI*180;

		//double tmp=0;
		//if(degree1<degree2){tmp=degree1;degree1=degree2;degree2=tmp;}
		//
		//
		//cout<<"degree="<<degree1<<"  "<<degree2<<endl; 
		//
		//tdp.getdetectedge(src,dst,xy,92,(4.08562-6));
		

		cvShowImage( "dst", dst );///  show  the  IplImage *dst
	}


















};