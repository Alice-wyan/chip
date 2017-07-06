#include"lib.h"

//!  a DetectLine class
class DetectLine
{

public:
	//构造函数
	DetectLine(){};

	
vector<CvScalar> getrgb(double x,double y,IplImage* src)
{
		// CvScalar s1;
		// s1.val[0]=100;
		// s1.val[1]=200;
		// s1.val[2]=100;
		//  CvScalar s2;
		// s2.val[0]=200;
		// s2.val[1]=100;
		// s2.val[2]=100;
		//  CvScalar s3;
		// s3.val[0]=200;
		// s3.val[1]=100;
		// s3.val[2]=100;
		// CvScalar s4;
		// s4.val[0]=200;
		// s4.val[1]=100;
		// s4.val[2]=100;
	 //cvSet2D(src,(int)(y),int(x),s1);
	 // cvSet2D(src,(int)(y+1),int(x),s2);
	 //  cvSet2D(src,(int)y,int(x+1),s3);
	 //   cvSet2D(src,(int)(y+1),int(x+1),s4);

	 vector<CvScalar> ps;
	 ps.clear();
	 ps.push_back(cvGet2D(src,(int)(y),int(x)));
	 ps.push_back(cvGet2D(src,(int)(y+1),int(x)));
	 ps.push_back(cvGet2D(src,(int)y,int(x+1)));
	 ps.push_back(cvGet2D(src,(int)(y+1),int(x+1)));

	 return ps;

}

//vector<double>pointsx,vector<double>pointsy  探测线所有像素点xy坐标集合
void getpixelrgb(vector<double>pointsx,vector<double>pointsy,IplImage* src)//获取探测线所有像素点的rgb值
{
	   bb.clear();///vector clear
	   gg.clear();
	   rr.clear();
	 CvPoint* pp=new CvPoint[4];

	  int pmin ;//距离次小值
	  int min;//距离最小值

	                	  	

	 //求取出像素点的rgb值
		for(int i=0;i<pointsx.size();i++)
		{
			vector<CvScalar>prgb=getrgb(pointsx[i],pointsy[i],src);//获取像素点周围的4个像素点的rgb；
			if((pointsx[i]==int(pointsx[i]))&&(pointsy[i]==int(pointsy[i]))){bb.push_back(prgb[0].val[0]);gg.push_back(prgb[0].val[1]);rr.push_back(prgb[0].val[2]);}
			//if((pointsx[i]==int(pointsx[i]))&&(pointsy[i]==int(pointsy[i]))){bb.push_back(prgb[0].val[0]);gg.push_back(0);rr.push_back(0);}
			else
			{
				pp[0].x=int(pointsx[i]);pp[0].y=int(pointsy[i]);
				pp[1].x=int(pointsx[i]);pp[1].y=int(pointsy[i]+1);
				pp[2].x=int(pointsx[i]+1);pp[2].y=int(pointsy[i]);
				pp[3].x=int(pointsx[i]+1);pp[3].y=int(pointsy[i]+1);

			 vector<double> L;//存储像素点距离周围四个像素点所有的距离
			 L.clear();

				for(int j=0;j<4;j++)//计算出像素点离周围4个像素点的距离
				{	  
				   L.push_back((pointsx[i]-pp[j].x)*(pointsx[i]-pp[j].x)+(pointsy[i]-pp[j].y)*(pointsy[i]-pp[j].y));
				   //cout<<"L="<<L[j]<<endl;
				}

				    double sum=0;
					double b=0;
					double g=0;
					double r=0;

             for(int k=0;k<L.size();k++)
			 {
				 sum+=L[k];

			 }
			// cout<<"sum="<<sum<<endl;
			 for(int k=0;k<prgb.size();k++)
			 {
				b+=prgb[k].val[0]*(L[k]/sum);				
				g+=prgb[k].val[1]*(L[k]/sum);
				//g+=0;
				//r+=0;
				r+=prgb[k].val[2]*(L[k]/sum);
			 }
			 bb.push_back(b);
			 gg.push_back(g);
			 rr.push_back(r);
			}
		}

			//fstream ph;

	  //      ph.clear();

		 //ph.open("ppd.txt",ios::app);
	//for(int i=0;i<bb.size();i++){ph<<"bb== "<<bb[i]<<endl;}

}

///
///函数getedgepix 用于获取像素点的边缘像素点
///输入 被探测的像素点坐标
///
vector<double> getedgepix(CvPoint point)
{

	  vector<double> XY;/// vector<double> XY is to store the xy coordinate of accurate point
	  XY.clear();///vector clear


 	 vector<double> bgr;///用于存储bgr的加权平均值
	 bgr.clear();///vector clear

	 CvScalar ss;///define the object of CvScalar

	  vector<double> dd;///vector<double> dd 用来存储像素rgb的差值
	  dd.clear();//清空

	  int aa=0;//用于判断情况
	  double m=0;//用于判断情况

	 double x=0;double y=0;//边缘像素坐标
	 int p1=0;int p2=0; //rbg差值最大值和次大值对应的像素点
	 int a=0;//区分两种差值情况设的变量



		 for(int i=0;i<bb.size();i++)///改变bgr值是为了增强信号
	 {
		 double aveage=(bb[i]+gg[i]+rr[i]);///探测线上每个像素点的rgb的加权平均值
		 if(aveage!=0)
		 {bgr.push_back((bb[i]*(bb[i]/aveage)+gg[i]*(gg[i]/aveage)+rr[i]*(rr[i]/aveage)));}///bgr存储着bgr的加权平均值
		 else{bgr.push_back(0);}
   		 ss.val[0]=ss.val[1]=ss.val[2]=bgr[i];///把探测线上每个像素点的的bgr分别设置同一个值即为bgr的加权平均值
		 //cout<<"bgr=="<<i<<"=="<<bgr[i]<<endl;
		 //cvSet2D(detect,points[i].y,points[i].x,s1);///把原图探测线上的像素点设置为改变了rgb值的像素点
	 }
     



	for(int i=0;i<bb.size()-1;i++)///for 用于存储探测线上相邻像素的差值
	 {
		 
		 dd.push_back(abs(bgr[i]-bgr[i+1]));///dd为探测线上相邻像素的rgb的差值
		// cout<<"dd=="<<dd[i]<<endl;
	 }

	 

	for(int i=0;i<dd.size();i++)///取rgb差值的最大值
	{

		if(dd[i]>m){ m=dd[i];}///选取出像素之间的rgb值的差值的最大值
		/**aa  用来判断探测线上的点的rgb值是否有变化如果没有变化则说明探测线上并没有边缘点
		*如果aa 的size==dd.size()-1  则说明并没有边缘点在探测线上
		*只有当aa<dd.size()-1 时才说明探测线上的相邻两像素之间的rgb差值大于2  即有边缘点存在
		*/
		if(dd[i]<2){aa++;}
	}
	
	
	if(aa<=(dd.size()-1))///进入此if循环内的说明有边缘点
	{
		
        /**
		方法一：
		双重for循环用来判断像素的rgb差值中有几个最大值
		*如果有两个最大值 则需要两个最大值对应的坐标的中值
		*若只有一个最大值，则提取出接近最大值周围的次大值，即坐标为出两者的加权平均数。
	    *若次大值远小于最大值，则直接利用最大值的像素点为边缘像素点
		*/
		/**
		方法二：
		双重for循环用来判断像素的rgb差值中有几个最大值
		若只有一个最大值，则提取出差值最大的两个像素点坐标，即坐标为出两者的加权平均数。
		*/

		for(int i=0;i<dd.size();i++)
		{
			
			for(int j=(dd.size()-1);j>i;j--)
			{
			  if((dd[i]==m)&&(dd[j]==m))///说明rgb差值中有两个都等于最大值
				  {
					 a=1;
					{
						 x=(double)(pointsx[i]+pointsx[j+1])/2;///提取出两个最大值的中值即为边缘坐标点的精确坐标
						 y=(double)(pointsy[i]+pointsy[j+1])/2;
					 }

				   XY.push_back(x);///把边缘点x坐标传递给vectorXY[0]
				   XY.push_back(y);///把边缘点y坐标传递给vectorXY[1]
				   break;
				  }
			    }
		}
	

			if(a==0)///则说明只有一个最大值
			{
		
				for(int i=0;i<dd.size();i++)
				{
					if(dd[i]==m)///选取出哪个坐标点为最大值
					{
						p1=i+1;///p1即为拥有rgb差值最大的坐标点
						p2=i;//采用方法二  差值最大的两个像素点坐标
					}
				 }

					double mm=0;
					for(int c=0;c<dd.size();c++)
					{
						if(c!=p1-1)
						{
							
							if(dd[c]>mm){mm=dd[c];}///选取出差值中出最大值之外的次大值
						}
						
					}
				//	
					{
					/*for(int h=0;h<dd.size();h++)
					{
							if(dd[h]==mm)
						{
							p2=h+1;///p2即为rgb差值为次大值的坐标点
						
						}
					}*/

					double ave=(bgr[p1]+bgr[p2]);
				   //  cout<<"p1=="<<p1<<"  "<<p2<<endl;
					x=(pointsx[p1]*(bgr[p1]/ave))+(pointsx[p2]*(bgr[p2]/ave));///求取出差值最大值与次大值的加权平均值即为边缘点坐标
					y=(pointsy[p1]*(bgr[p1]/ave))+(pointsy[p2]*(bgr[p2]/ave));

					}
					/*else
					{
						x=pointsx[p1];y=pointsy[p1];
					}*/
				   XY.push_back(x);///把边缘点x坐标传递给vectorXY[0]
				   XY.push_back(y);///把边缘点y坐标传递给vectorXY[1]
				   			
	            }
	}
	else{
		   cout<<"ERROR=="<<point.x<<"  "<<point.y<<"输入坐标中间没有边缘点"<<endl;				
	       XY.push_back(point.x);///把边缘点x坐标传递给vectorXY[0]
		   XY.push_back(point.y);///把边缘点y坐标传递给vectorXY[1]
           
	    }
	return  XY;


}



///
///CvPoint point int型像素点坐标   vector<double>pointsx,vector<double>pointsy存储探测线像素点xy坐标 int distance探测线长度
///
vector<double> detectline(CvPoint point,vector<double>pointsx,vector<double>pointsy,int distance,IplImage* src,IplImage* dst)//src 是原图
{


	//.............................变量定义区.............................................//
	  CvPoint X,Y;	//定义探测线左右端点以便于在dst上画线
	  vector<double> XY;/// vector<double> XY is to store the xy coordinate of accurate point
	  XY.clear();///vector clear

//....................................................................................................//

	 X.x=pointsx[0];X.y=pointsy[0];//XY分别为两端端点，转换为int 数据类型仅为了画直线好看 
	Y.x=pointsx[2*distance];Y.y=pointsy[2*distance];//Y为探测线右端断点

	cvLine(dst,X,Y,CV_RGB(0,255,0),1,8,0);///将XY坐标连线 	
	getpixelrgb(pointsx,pointsy,src);//获取探测线像素点rgb值


	  XY=getedgepix(point);//获取边缘坐标	 
	  return XY;///return   边缘点xy坐标

}
///
///CvPoint x即探测线的中间像素点    double degree即探测线的角度（负值为向下，正值向上） IplImage *src--输入原图  IplImage *dst--为画上探测线的图片
///
 vector<double> detectlinedegree(CvPoint x,double degree,int distance,IplImage *src,IplImage *dst)
{
	
	vector<double>XY;//定义集合用于存储探测到的边缘像素的坐标
	XY.clear();//清空集合
	///
	///调用getdetectpoints来获取探测线上所有像素点坐标
	///
	getdetectpoints(x,degree,distance);	 
	///
	///调用detectline函数来获取边缘像素点坐标并存储在集合XY中
	///
	XY=detectline(x,pointsx,pointsy,distance,src,dst);
	
	return XY;


}

 ///
 ///getdetectpoints函数用于获取探测线上所有的像素点坐标
 ///CvPoint x 输入探测线的中点坐标  degree为探测线角度  distance为探测线长度的一半
 ///
void getdetectpoints(CvPoint x,double degree,int distance)
{
	  ///
	  ///探测线必须是直线，则由于探测线还可能是斜的，则探测线上的点必须是double型数据
	  ///
		double* x1=new double[distance]; double* y1=new double[distance];//为探测线的左端像素点
        double* x2=new double[distance]; double* y2=new double[distance];//为探测线的右端像素点

		pointsx.clear();//存储探测线左侧端点
		pointsy.clear();//存储探测线右侧端点

    if(degree<=0)//当探测线角度小于或等于0时
	{
		for(int i=0;i<distance;i++)//把x1设为x像素的上半部分，x2为x的下半部分
		{
			x1[i]=x.x-((distance-i)*cos(degree*CV_PI/180.0));//利用探测线角度与探测线长度来获取探测线上所有像素点坐标
			y1[i]=x.y+((distance-i)*sin(degree*CV_PI/180.0));//利用探测线角度与探测线长度来获取探测线上所有像素点坐标

			x2[i]=x.x+((i+1)*cos(degree*CV_PI/180.0));//利用探测线角度与探测线长度来获取探测线上所有像素点坐标
			y2[i]=x.y-((i+1)*sin(degree*CV_PI/180.0));//利用探测线角度与探测线长度来获取探测线上所有像素点坐标
		
		}
	}
	else//当探测线角度为正时
	{    for(int i=0;i<distance;i++)//把x1设为x像素的上半部分，x2为x的下半部分
		{

			x1[i]=x.x+((distance-i)*cos(degree*CV_PI/180.0));//利用探测线角度与探测线长度来获取探测线上所有像素点坐标
			y1[i]=x.y-((distance-i)*sin(degree*CV_PI/180.0));//利用探测线角度与探测线长度来获取探测线上所有像素点坐标

			x2[i]=x.x-((i+1)*cos(degree*CV_PI/180.0));//利用探测线角度与探测线长度来获取探测线上所有像素点坐标
			y2[i]=x.y+((i+1)*sin(degree*CV_PI/180.0));//利用探测线角度与探测线长度来获取探测线上所有像素点坐标

		}

	}

		for(int i=0;i<distance;i++)//把探测线上像素点都存储在集合pointsx、pointsy
		{
			pointsx.push_back(x1[i]);
			pointsy.push_back(y1[i]);
		}
		  pointsx.push_back(x.x);///把探测线中点像素点存储在集合pointsx、pointsy
		  pointsy.push_back(x.y);
		for(int k=0;k<distance;k++)//把探测线上像素点都存储在集合pointsx、pointsy
		{
			pointsx.push_back(x2[k]);
			pointsy.push_back(y2[k]);
		}


		delete x1;delete x2;///删除动态数组空间
		delete y1; delete y2;

}

private:
///
///私有变量  
///
	vector<double> pointsx;///存储探测线上所有像素点坐标
	vector<double> pointsy;///存储探测线上所有像素点坐标
	 vector<double> bb;///存储着原图的探测线上的每个像素点的b值 
     vector<double> gg;///存储着原图的探测线上的每个像素点的g值 
	 vector<double> rr;///存储着原图的探测线上的每个像素点的r值 

};