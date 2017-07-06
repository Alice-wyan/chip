#include"lib.h"
//#include"getcorner.h"

//用于通过确定的角点和直线的角度来确定 探测线的角度。通过探测线求取出边缘像素点坐标，
//然后将这些 边缘像素点拟合，拟合后求取出拟合的两直线的交点---即为角点的精确坐标




class TestDetecteCorner
{
  public:

	  void getcorner(IplImage *src,IplImage *dst)
	  {
		  DetectLine dll;
		  Getcorner gc;

		  vector<double> l1;
		  l1.clear();
		  vector<double> l2;
		  l2.clear();
		  vector<double> lx;
		  lx.clear();

		  vector<double> l3;
		  l3.clear();
		  vector<double> l4;
		  l4.clear();


		  vector<double> line;
		  line.clear();

		  vector<double> XY;
		  XY.clear();
		  CvPoint p;
		  p.x=281;p.y=265;
				 int i=10; int j=90;

				 CvPoint po,p1;
				 po.x=p.x+i*cos((60)*CV_PI/180.0); po.y=p.y+i*sin((60)*CV_PI/180.0);//通过拟合像素点角度来获取需要被探测的像素点
				 p1.x=p.x+j*cos((60)*CV_PI/180.0); p1.y=p.y+j*sin((60)*CV_PI/180.0);//通过拟合像素点角度来获取需要被探测的像素点

				 l1=dll.detectlinedegree(po,(30),5,src,dst);
				 lx.push_back(l1[0]);//存储y方向的x坐标
				 lx.push_back(l1[1]);//存储y方向的y坐标
				 l1x.push_back(l1[0]);
				 l1y.push_back(l1[1]);

				 l2=dll.detectlinedegree(p1,(30),5,src,dst);
				 lx.push_back(l2[0]);//存储y方向的x坐标
				 lx.push_back(l2[1]);//存储y方向的y坐标
				 l2x.push_back(l2[0]);
				 l2y.push_back(l2[1]);


				 CvPoint p2,p3;
				 p2.x=p.x+i*cos((30)*CV_PI/180.0); p2.y=p.y-i*sin((30)*CV_PI/180.0);//通过拟合像素点角度来获取需要被探测的像素点
				 p3.x=p.x+j*cos((30)*CV_PI/180.0); p3.y=p.y-j*sin((30)*CV_PI/180.0);//通过拟合像素点角度来获取需要被探测的像素点

				 l3=dll.detectlinedegree(p2,(-60),5,src,dst);
				 line.push_back(l3[0]);//存储y方向的x坐标
				 line.push_back(l3[1]);//存储y方向的y坐标
				 l3x.push_back(l3[0]);
				 l3y.push_back(l3[1]);

				 l4=dll.detectlinedegree(p3,(-60),5,src,dst);
				 line.push_back(l4[0]);//存储y方向的x坐标
				 line.push_back(l4[1]);//存储y方向的y坐标
				 l4x.push_back(l4[0]);
				 l4y.push_back(l4[1]);



			  
				 //line.push_back(325);
				 //line.push_back(236);

				 //line.push_back(343);
				 //line.push_back(224);

             double k1=(l2[1]-l1[1])/(l2[0]-l1[0]);
			 double b1=l2[1]-k1*l2[0];

			 double k2=(line[3]-line[1])/(line[2]-line[0]);
			 double b2=line[3]-k2*line[2];

			 XY.push_back((b2-b1)/(k1-k2));
			 XY.push_back(k1*XY[0]+b1);

			 //cout<<"pp="<<XY[0]<<"  "<<XY[1]<<endl;


			// if((l1.size()!=0)&&(l2.size()!=0))///确定有角点才执行
			//{
			//	XY=gc.FitLine(dst,lx, line);
			//	//gc.lineFit(xdetectpoint);

			//}
          	CvPoint xy;
			xy.x=XY[0];
			xy.y=XY[1];
			cvCircle(dst,xy,10, CV_RGB(255,0,0), 1,8,0);//在dst上画点
			X.push_back(XY[0]);
			Y.push_back(XY[1]);


	  }






public:
///
///公共变量
///
vector<double> l1x;//存储y方向像素点
vector<double> l1y;//存储x方向像素点
vector<double> l2x;///定义常量与变量并初始化。
vector<double> l2y;///存储边缘像素坐标

vector<double> l3x;//存储y方向像素点
vector<double> l3y;//存储x方向像素点
vector<double> l4x;///定义常量与变量并初始化。
vector<double> l4y;///存储边缘像素坐标



vector<double> X;
vector<double> Y;



};
