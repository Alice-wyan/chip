#include"lib.h"
#include"findcorner.h"
//#include"detectline.h"
#include"getdistanceofcorner.h"
#define R2D    57.295779513082323
#define Deg2Rad( Deg ) (Deg)/ R2D  // 角度到弧度
#define Rad2Deg( Rad ) (Rad)* R2D  // 弧度到角度

class FindCorner;

//!  a  Getcorner class
class Getcorner
{
//!  a constructor of Getcorner
public:Getcorner(){}

///@brief   the setup function
/*! the setup function is to  transform IplImage* frame to    CSixSword  class 
* to get the contour of CEdgeHitChain chain
*/
///@param   IplImage* frame input IplImage* 
void  setup(IplImage* frame)///输入原图
{  
	  if(!frame){cout<<"no frame !";}/// if  no frame   will cout  no frame
      ss.setup(frame);///  setup frame to   CSixSword  class
	  fr.set(frame);///transform frame to CFrame class  for read and draw contour in picture
	  ss.transform();/// invoke transform function of CSixSword  class  to get the contour of CEdgeHitChain chain
	 
 }

///@brief   the getcornerpoint function
/*! the getcornerpoint function's  object to  get coernerpoint  */
///@param   IplImage* frame input IplImage* 
void getcornerpoint(IplImage* frame, IplImage* dst,int T1,int T2,int T3,int T4)//T1与T2是判断哪三个点来构成三角形计算出余弦值  T1与T2即为欧式距离
{

	  CvSeq* temp;///define the  CvSeq* temp (contour)
	  FindCorner fc;///define object  fc of FindCorner class
	  DetectLine dl;///define object dl  of DetectLine class
	  
	  points.clear();///存储着过滤后的图片中的所有角点坐标
	  accos.clear();
	  vector<double> allcorner;///存储着所有轮廓的所有候选角点坐标
	  allcorner.clear();/// clear the vector

	  vector<double> allacos;///存储着左右轮廓的所有候选角点的反余弦值
	  allacos.clear();/// clear the vector

	 vector<CvPoint>  intpoints;///中间变量  存储着过滤后的所有角点坐标  转换为CvPoint 用来在图上标记出角点
	 intpoints.clear();///clear the vector


	 IplImage *src=cvCloneImage(frame);///  frame clone to src  for  CSixSword ss to setup
	
	 setup(src);///CSixSword ss setup  IplImage *src  to get CEdgeHitChain chain(contour)


	  cout<<"size=="<<ss.v.size()<<endl;///cout the size of contour 

          /// for 循环循环轮廓数   分别对每个轮廓求取角点
	 	  	for(unsigned int k=0;k<ss.v.size();k++)
			{
               	 vector<CvPoint> cornerpoints;///存储着单个轮廓的所有候选角点
				 vector<double> corneracos;////存储着单个轮廓的所有候选角点的反余弦值
				 cornerpoints.clear();///clear the vector
				 corneracos.clear();///clear the vector
				  vector<double> XY ;
			      XY.clear();

			  /* vector<CvPoint> point;///存储着单个轮廓上的所有像素点
	           point.clear();/// clear the vector*/
			
			  temp=transform(ss.v[k]);/// transform the contour of CEdgeHitChain chain  to contour of CvSeq* 
			 
			  ss.v[k].print(fr);///draw the contour in IplImage *src 
			  
			  /**  invoke the findcorner function of FindCorner class
			  * input  the contour of CvSeq* 
			  * input  int T1
			  *   T1 是 设置double L 的最小范围阈值  即  L>T1
			  *input  int T2
			  *   T2  是设置double L的最大范围阈值   即L<=T2
			  *input   int T3
			  *   T3 是设置个角度阈值，abs(90-acos)<T3  是比较于90度的角度值阈值
			  */
			 fc.findcorner(dst,temp,T1,T2,T3);//获取轮廓的角点
			 
			  for(int i=0;i<fc.corner.size();i++)///  通过for循环把每个轮廓的候选角点均存储在集合中
			  {
				  cornerpoints.push_back(fc.corner[i]);///  using the circulation to put every corner in the vector of cornerpoints
				  corneracos.push_back(fc.accos[i]);/// using the circulation to put every corner's the value of acos in the vector of cornerpoints
			  }

			  for(int i=0;i<cornerpoints.size();i++)
			  {
				 // cout<<"cornerpoints"<<corneracos[i]<<"   "<<cornerpoints[i].x<<"    "<<cornerpoints[i].y<<endl;
				  cvCircle(dst,cornerpoints[i], 5, CV_RGB(255,0,0), 1,8,0);
			  }

			  fc.setup(temp); ///调用findcontour的类的setup函数来获得当前轮廓上所有的点
			  //cout<<"point.size== "<<fc.points.size()<<endl;

		
			if(cornerpoints.size()!=0){XY=detectp(fc.points,cornerpoints,corneracos,frame,dst);}
					
			
			 if(XY.size()!=0)///XY有值
			 {
			 for(int i=0;i<XY.size()-2;i=i+3)
			 { 
			    allcorner.push_back(XY[i]);///把x坐标存储给allcorner   为什么不用CvPoint  是因为CvPoint只支持int 
				allcorner.push_back(XY[i+1]);///把x坐标存储给allcorner
				allacos.push_back(XY[i+2]);///把角点反余弦值存储给allacos
				
			 }
			 }

			 for(int i=0;i<fc.points.size()-1;i++){cvLine( src,fc.points[i], fc.points[i+1], CV_RGB(255,0,0),1, 8, 0 );}/// 通过轮廓上的点来画出轮廓		
 
		      cvReleaseMemStorage (&tempStorage);///释放transform函数上创建的内存。
			}

			//! invoke the filtercorner function
			/** input the threshold value T3、T4 
			*T3是阈值  即在最终得到每个轮廓的角点坐标后需要选出abs(90-acos)<T3的角度阈值
			*T4是阈值，即在T4单位像素之内  筛选出acos值大的作为角点保留且设置为角点为90度必须保留
			*/
			if(allcorner.size()!=0)
			{
			    filter(allcorner,allacos,T4);///过滤角点
			}


          if(points.size()!=0)
		  {
			for(int i=0;i<points.size()-1;i=i+2)
			{
				if((points[i]>0)&&(points[i+1]>0))///把points存储给intpoints是方便在图上画出最终角点。
				{
				CvPoint pp;
				pp.x=points[i];
				pp.y=points[i+1];
				intpoints.push_back(pp);
				}
            }	

//为了判断稳定性则需要通过下面的比较选择想要留取的点坐标
			//此处想 保留y坐标最大的角点坐标
			//double pointxx=0; double pointyy=0;
			//for(int i=1;i<points.size();i=i+2)
			//{
			//	if(abs(points[i]-400)<5&&(abs(points[i-1]-292)<5)){pointyy=points[i];pointxx=points[i-1];}

			//}

			//fstream pppx;

	  //      pppx.clear();

		 //pppx.open("pppx.txt",ios::app);
	
		 //pppx<<pointxx<<endl;

	  //    fstream pppy;

	  //      pppy.clear();

		 //pppy.open("pppy.txt",ios::app);
	
		 //pppy<<pointyy<<endl;



			if(intpoints.size()!=0)
			{
		   	  for(int i=0;i<intpoints.size();i++)///for  circulation to draw the corenerpoint in picture 
				{

					cout<<"intpoints==  "<<accos[i]<<"   "<<intpoints[i].x<<"   "<<intpoints[i].y<<endl;
					cvCircle(dst,intpoints[i],10, CV_RGB(0,0,255), 1,8,0);/// 在图片上画出角点，以角点为中心，5个像素为半径画圆。
				} 
		    }
		 }
		
      
       
		cvNamedWindow( "dst", CV_WINDOW_AUTOSIZE );	 ///define a window
        cvShowImage( "dst", dst );///  show  the  IplImage *dst
		cvNamedWindow( "srcc", CV_WINDOW_AUTOSIZE );	 ///define a window
        cvShowImage( "srcc", src );///  show  the  IplImage *dst

		cvReleaseImage(&src);

}


vector<double>  FitLine(  IplImage* dst,vector<double> xdetectpoint, vector<double> ydetectpoint)//输入图片与探测后的角点两边的边缘像素点坐标
{
	       vector<double> xy;//存储直线拟合后的两直线的交点
		   xy.clear();

			float *liner = new float[4];///定义数组
			float *linel = new float[4];
			CvMemStorage* storage1 = cvCreateMemStorage(0);///创建存储空间
			CvMemStorage* storage2 = cvCreateMemStorage(0);
			CvSeq* point_seqr = cvCreateSeq( CV_32FC2,sizeof(CvSeq), sizeof(CvPoint2D32f), storage1 );///创建序列
			CvSeq* point_seql = cvCreateSeq( CV_32FC2,sizeof(CvSeq), sizeof(CvPoint2D32f), storage2 );
			for(int s=0;s<xdetectpoint.size()-1;s=s+2)
			{
			cvSeqPush(point_seqr, &cvPoint2D32f(xdetectpoint[s],xdetectpoint[s+1]));///把像素点精确坐标push给point_seqr
			}
			for(int s=0;s<ydetectpoint.size()-1;s=s+2)
			{
			cvSeqPush(point_seql, &cvPoint2D32f(ydetectpoint[s],ydetectpoint[s+1]));///把像素点精确坐标push给point_seql
			}

			cvFitLine(point_seqr,CV_DIST_HUBER ,0,0.001,0.001,liner);///直线拟合函数   x 方向
			//cvFitLine(point_seql,CV_DIST_L2,0,0.001,0.001,linel);///把ydetectpoint的像素点进行直线拟合   y 方向
			cvFitLine(point_seql,CV_DIST_HUBER ,0,0.001,0.001,linel);
			//double kr=atan(liner[1]/liner[0]);
			//double kl=atan(linel[1]/linel[0]);
			//double angle1= Rad2Deg(kr);
			//double angle2=Rad2Deg(kl);

			
			
			double x1=liner[2]; double y1=liner[3];/// line[0--3] 分别为 (vx, vy,x0, y0)
			double x2=linel[2]; double y2=linel[3];

			//cout<<"liner[0]="<<liner[0]<<"  "<<liner[1]<<endl;
			//cout<<"linel[0]="<<linel[0]<<"  "<<linel[1]<<endl;
			
			//cout<<"x1="<<x1<<"  "<<y1<<endl;
			//cout<<"x2="<<x2<<"  "<<y2<<endl;
			//cout<<"angle1 "<<angle1<<"  x= "<<x1<<"   "<<y1<<endl;
			//cout<<"angle2 "<<angle2<<"  x=  "<<x2<<"   "<<y2<<endl;
			double k1=liner[1]/liner[0];double b11=liner[3]-k1*liner[2];
			double k2=linel[1]/linel[0];double b22=linel[3]-k2*linel[2];///确定直线斜率及直线方程


			//cout<<"k1=="<<k1<<"  "<<b11<<"   "<<k2<<"  "<<b22<<"  "<<endl;

			xy.push_back((y2-y1+k1*x1-k2*x2)/(k1-k2));///利用两直线求角点			
			xy.push_back(k1*(xy[0]-x1)+y1);

			
			double cos_theta1 = liner[0];
			double sin_theta1 = liner[1];
			double x01 = liner[2], y01 = liner[3];

			double cos_theta = linel[0];
			double sin_theta = linel[1];
			double x0 = linel[2], y0 = linel[3];


				 double k = sin_theta / cos_theta;
			   double b = y0 - k * x0;
			 double x = 0;
			double y = k * x + b;


				 double k11 = sin_theta1 / cos_theta1;
			   double b1 = y01 - k11 * x01;
			 double x11 = 0;
			double y11 = k11 * x11 + b1;
			Mat dd=dst;

	//line(dd, Point(x0,y0), Point(x,y), cv::Scalar(0,255,255), 1);
	//line(dd, Point(x01,y01), Point(x11,y11), cv::Scalar(0,255,255), 1);
	// DrawLine(dst, phi, rho, cv::Scalar(0));



			cvReleaseMemStorage (&storage1);
			cvReleaseMemStorage (&storage2);


            return xy;

}



int lineFit(vector<double> xdetectpoint)
{
     int size = points.size();
	 double a=0;double b=0;double c=0;
	
     double x_mean = 0;
     double y_mean = 0;
     for(int i = 0; i < xdetectpoint.size()-1; i=i+2)
     {
         x_mean += xdetectpoint[i];
         y_mean +=xdetectpoint[i+1];
     }
     x_mean /= size;
     y_mean /= size; //至此，计算出了 x y 的均值

     double Dxx = 0, Dxy = 0, Dyy = 0;

     for(int i = 0; i < size; i++)
     {
         Dxx += (xdetectpoint[i] - x_mean) * (xdetectpoint[i] - x_mean);
         Dxy += (xdetectpoint[i] - x_mean) * (xdetectpoint[i+1] - y_mean);
         Dyy += (xdetectpoint[i+1] - y_mean) * (xdetectpoint[i+1] - y_mean);
     }
     double lambda = ( (Dxx + Dyy) - sqrt( (Dxx - Dyy) * (Dxx - Dyy) + 4 * Dxy * Dxy) ) / 2.0;
     double den = sqrt( Dxy * Dxy + (lambda - Dxx) * (lambda - Dxx) );
      if(fabs(den) < 1e-5)
     {
         if( fabs(Dxx / Dyy - 1) < 1e-5) //这时没有一个特殊的直线方向，无法拟合
         {
             return 0;
         }
         else
         {
             a = 1;
             b = 0;
             c = - x_mean;
         }
     }
     else
     {
         a = Dxy / den;
         b = (lambda - Dxx) / den;
         c = - a * x_mean - b * y_mean;
     }
   // cout<<"a="<<a<<"  "<<b<<"  "<<c<<endl;
} 







double getdtectdegree(vector<CvPoint> accurate)
{
	double degree=0;
	for(int g=0;g<accurate.size();g++)//确定y轴角度
	{
		if(g<4)///如果是顶点则选择g与g+1之差来确定斜率
		{

			if(accurate[g].x==accurate[g+4].x){degree=0;}//如果两端点的x坐标相等则说明探测线角度为0度

			if(accurate[g].x!=accurate[g+4].x)///如果两像素点x坐标不相等则通过斜率公式计算后求出探测线角度即可。
			{
							
				double slopey=(-1)*(accurate[g].y-accurate[g+4].y)/(accurate[g].x-accurate[g+4].x);
				degree=atan((double)(-1/slopey))/CV_PI*180;

			}
									 
		}
		else
		{

			if(accurate[g].x==accurate[g-4].x){degree=0;}//如果两端点的x坐标相等则说明探测线角度为0度

			if(accurate[g].x!=accurate[g-4].x)///如果两像素点x坐标不相等则通过斜率公式计算后求出探测线角度即可。
			{							
				double slopey=(-1)*(accurate[g].y-accurate[g-4].y)/(accurate[g].x-accurate[g-4].x);			                     
				degree=atan((double)(-1/slopey))/CV_PI*180;
			}

		}
	}
		return degree;

}

vector<double> comparedis(CvPoint  coener,vector<CvPoint> cornerpoints)
{
	vector<double> cd;//角点之间的距离
	cd.clear();
	vector<double>  dist;
	dist.clear();
	
	  for(int d=0;d<cornerpoints.size();d++)//计算角点之间的相互距离
	{
		if((abs(coener.x-cornerpoints[d].x)>5)||(abs(coener.y-cornerpoints[d].y)>5))
		{ double dist=sqrt((double)((coener.x-cornerpoints[d].x)*(coener.x-cornerpoints[d].x)+(coener.y-cornerpoints[d].y)*(coener.y-cornerpoints[d].y)));
		cd.push_back(dist);}
	}
	  double min=cd[0];//距当前角点最小距离
	  double smin=cd[1];//距当前角点距离次小值
	  for(int i=0;i<cd.size();i++)
	  {
		  if(min>cd[i]){min=cd[i];}
	  }

	   for(int i=0;i<cd.size();i++)
	  {
		  if(abs(cd[i]-min)>5)
		  {
		  if(smin>cd[i]){smin=cd[i];}
		  }
	  }
	   dist.push_back(min);
	   dist.push_back(smin);
       
	   return dist;	   

}


vector<double > detectp(vector<CvPoint> point, vector<CvPoint> cornerpoints,vector<double> corneracos,IplImage *src, IplImage *dst)
{
	//..........变量定义区.......................//
	     //int L1=5;int L2=25;
	      int L1,L2x,L2y;//确定探测线的总长
		 DetectLine dl;///创建DetectLine 类的对象 dl
		// CvPoint  point1,point2,point3,point4;///定义CvPoint 对象  即两条直线上的四个点
		 vector<double > corner;
		 corner.clear();///定义存储两直线交点的xy坐标的集合
		 

		  double x=0;double y=0;///
		  int k=1; int errord=5;
	//.........................................................................................................................................//

    //通过角点之间的长度来确定探测线总共的长度L1-L2
		  vector<double> ds;
		  ds.clear();

   //......................................................................................................................................................//

			for(int i=0;i<cornerpoints.size();i++)///双重循环来确定两直线的交点坐标，即角点的精确坐标
		   {
			   
			   			  vector<double> xdetectpoint;///定义常量与变量并初始化。
						  vector<double> ydetectpoint;//存储y方向探测线坐标
						  xdetectpoint.clear();//清空集合
						  ydetectpoint.clear();//清空集合
						  vector<double> dist;
		                  dist.clear();
						  
                       dist=comparedis(cornerpoints[i],cornerpoints);
					   L1=5;L2x=dist[1];L2y=dist[0];  
						//cout<<"L2x=="<<L2x<<"  "<<L2y<<endl;  
             //if(k==1)
			{
			   for(int j=0;j<point.size();j++)///利用角点和轮廓上的点来进行循环旨在分别找出角点左右两端的像素点 point--轮廓上所有点坐标
			   {

				  				   
				  // ...........................................................................................................//
				   		vector<CvPoint> dispoint;///   vector<CvPoint> dispoint 用来储存进行余弦定理求夹角的三个像素点坐标		
						vector<double> disL;///   vector<double> disL 用来储存进行余弦定理求夹角的三个像素点之间的距离		 
						dispoint.clear();///清空vector<CvPoint>  dispoint 集合		
						disL.clear(); ///清空vector<double>  dispoint 集合
						
                      
				          vector<CvPoint> xaccurate;//存储x方向被探测像素点
						  vector<CvPoint> yaccurate;//存储y方向被探测像素点
						  xaccurate.clear();
						  yaccurate.clear();

						  double degreey;//y方向探测线角度
						  double degreex;//x方向探测线角度
						
						   vector<double>xyx;
						   vector<double>xyy;
						   xyx.clear();
						   xyy.clear();


						   vector<double> XY;//直线交点坐标
						   XY.clear();

						   
                //..................................................................................//
						 
		     // if((abs(cornerpoints[i].x-202)<5)&&(abs(cornerpoints[i].y-339)<5))
				{
				  // k=k+1;
				if((point[j].x==cornerpoints[i].x)&&(point[j].y==cornerpoints[i].y))///如果轮廓上的点坐标等于角点坐标则进入循环
				  {
					  
		  				fstream pppx;

						pppx.clear();

					 pppx.open("pppx.txt",ios::app);
	

					  fstream pppy;

						pppy.clear();

					 pppy.open("pppy.txt",ios::app);
					  //.....................................................................................................................//
						for (int d=point.size()-1; d>=0; d--)///for循环 像素点的末端来进行与正在遍历的像素点的两者之间的欧氏距离
						{		 			
							double diss=((cornerpoints[i].x-point[d].x)*(cornerpoints[i].x-point[d].x)+(cornerpoints[i].y-point[d].y)*(cornerpoints[i].y-point[d].y));///double diss 为两像素点之间的距离  x*x+y*y
			
							double L=sqrt(diss);///double L 为两像素之间的欧式距离 sqrt(x*x+y*y)

								if((L>=L1)&&(L<=L2x))///if 条件来进行筛选  选出两像素之间的欧式距离在 T1<L<=T2的两像素点 
								{
									disL.push_back(L);///符合条件的欧氏距离存储在vector<double> disL
									dispoint.push_back(point[d]);///符合条件的像素点坐标存储在 vector<CvPoint> dispoint	
									//pp<<L<<"  "<<point[d].x<<"  "<<point[d].y<<endl;
								}
						 }  
							//pppx<<"point "<<point[j].x<<"  "<<point[j].y<<endl;
							//pppy<<"point "<<point[j].x<<"  "<<point[j].y<<endl;
							for(int g=0;g<dispoint.size();g++)
							{
								//if((abs(dispoint[g].x-point[j].x)<10)&&(abs(dispoint[g].y-point[j].y)>=L1)&&(dispoint[g].y<point[j].y))
								//if((abs(dispoint[g].x-point[j].x)<5)&&(abs(dispoint[g].y-point[j].y)>=L1))
								//{yaccurate.push_back(dispoint[g]);pppy<<dispoint[g].x<<"  "<<dispoint[g].y<<endl;
								//}//说明是y轴方向像素点
								//if((abs(dispoint[g].x-point[j].x)>=L1)&&(abs(dispoint[g].y-point[j].y)<5))
								//{xaccurate.push_back(dispoint[g]);pppx<<dispoint[g].x<<"  "<<dispoint[g].y<<endl;
								//}//说明是x轴方向像素点
                                 //.......................试试xy完全理想化...........................................//
								
								if((abs(dispoint[g].x-point[j].x)<=L1)&&(abs(dispoint[g].y-point[j].y)>L1)&&(abs(dispoint[g].y-point[j].y)<(L2y-5)))
								{yaccurate.push_back(dispoint[g]);//pppy<<dispoint[g].x<<"  "<<dispoint[g].y<<endl;
								}//说明是y轴方向像素点
								if((abs(dispoint[g].x-point[j].x)>L1)&&((abs(dispoint[g].x-point[j].x)<=(L2x-5))&&(abs(dispoint[g].y-point[j].y)<=L1)))
								{xaccurate.push_back(dispoint[g]);//pppx<<dispoint[g].x<<"  "<<dispoint[g].y<<endl;
								}//说明是x轴方向像素点

							}

							if((yaccurate.size()>=4)&&(xaccurate.size()>=4))
							{
								int dky=0; int dkx=0;
								if(((yaccurate.size()>=8)&&(xaccurate.size()>=8))){ dkx=4;}//两个像素点之间的距离来确定探测线角度
								if((yaccurate.size()<8)&&(xaccurate.size()>=8)){ dky=yaccurate.size()/2;}
								if((xaccurate.size()<8)&&(yaccurate.size()>=8)){dkx=xaccurate.size()/2;}
							for(int g=0;g<yaccurate.size();g++)//确定y轴角度
							{
								
								if(g<dky)///如果是顶点则选择g与g+1之差来确定斜率
								{   
																
									if(yaccurate[g].x==yaccurate[g+dky].x){degreey=0;}//如果两端点的x坐标相等则说明探测线角度为0度

									if(yaccurate[g].x!=yaccurate[g+dky].x)///如果两像素点x坐标不相等则通过斜率公式计算后求出探测线角度即可。
									{							          
									  double slopey=(-1)*(yaccurate[g].y-yaccurate[g+dky].y)/(yaccurate[g].x-yaccurate[g+dky].x);
									  degreey=atan((double)(-1/slopey))/CV_PI*180;//确定探测线角度
									}
									
								}
								else
								{   
									if(yaccurate[g].x==yaccurate[g-dky].x){degreey=0;}//如果两端点的x坐标相等则说明探测线角度为0度

									if(yaccurate[g].x!=yaccurate[g-dky].x)///如果两像素点x坐标不相等则通过斜率公式计算后求出探测线角度即可。
									{							
									   double slopey=(-1)*(yaccurate[g].y-yaccurate[g-dky].y)/(yaccurate[g].x-yaccurate[g-dky].x);			                     
									    degreey=atan((double)(-1/slopey))/CV_PI*180;//确定探测线角度
									}

								}
								
								
								//pppy<<"degreey=="<<degreey<<endl;
							    xyy=dl.detectlinedegree(yaccurate[g],degreey,5,src,dst);///调用detectlinedegree函数求取出像素点精确xy坐标
								ydetectpoint.push_back(xyy[0]);///把求取的像素点精确坐标的xy坐标存储在ydetectpoint
								ydetectpoint.push_back(xyy[1]);
								p.push_back(yaccurate[g]);//存储被探测的像素点坐标
								pd.push_back(xyy[0]);//存储探测后的像素点坐标
								pd.push_back(xyy[1]);
							
								//pppy<<xyy[0]<<"   "<<xyy[1]<<endl;
							}
							
							for(int k=0;k<xaccurate.size();k++)///确定x轴方向的像素点的角度
							{
							
								if(k<dkx)
								{   
									 if(xaccurate[k].x==xaccurate[k+dkx].x){degreex=-90;}/////如果两端点的x坐标相等则说明探测线角度为90度
									 else ///如果两像素点x坐标不相等则通过斜率公式计算后求出探测线角度即可
									 {
										
									double slopex=(-1)*(xaccurate[k].y-xaccurate[k+dkx].y)/(xaccurate[k].x-xaccurate[k+dkx].x);//求取直线斜率
									 degreex=atan((double)(-1/slopex))/CV_PI*180;//换算成角度
									 }								
								}
								else///非顶点像素点
								{ 
									if(xaccurate[k].x==xaccurate[k-dkx].x){degreex=-90;}//如果两端点的x坐标相等则说明探测线角度为90度
									else//如果两像素点x坐标不相等则通过斜率公式计算后求出探测线角度即可
									{
										
							         double slopex=(-1)*(xaccurate[k].y-xaccurate[k-dkx].y)/(xaccurate[k].x-xaccurate[k-dkx].x);
									 degreex=atan((double)(-1/slopex))/CV_PI*180;
									}
									 
								}
								//pppx<<"degreex=="<<degreex<<endl;
								xyx=dl.detectlinedegree(xaccurate[k],degreex,5,src,dst);//调用detectlinedegree函数求取像素点精确xy坐标
								xdetectpoint.push_back(xyx[0]);//把求取的像素点精确坐标的xy坐标存储在xdetectpoint
							    xdetectpoint.push_back(xyx[1]);
								p.push_back(xaccurate[k]);//存储被探测的像素点坐标
								pd.push_back(xyx[0]);//存储探测后的边缘像素点坐标
								pd.push_back(xyx[1]);
								pppx<<xyx[0]<<"  "<<xyx[1]<<endl;
						}
//
//
					}//for
					    if((xdetectpoint.size()!=0)&&(ydetectpoint.size()!=0))///确定有角点才执行
						{
							XY=FitLine(dst,xdetectpoint, ydetectpoint);
						}
						//cout<<"XY=="<<XY[0]<<"  "<<XY[1]<<"   cornerpoints="<<cornerpoints[i].x<<" "<<cornerpoints[i].y<<endl;
///.........................................................................................................................//
						if(XY.size()!=0)
						{
					 if((abs(XY[0]-cornerpoints[i].x)<errord)&&(abs(XY[1]-cornerpoints[i].y)<errord))//如果计算的精确角点坐标与原角点坐标相差大于5则返回原角点坐标
						{
							corner.push_back(XY[0]);
							corner.push_back(XY[1]);
							corner.push_back(corneracos[i]);
						}

					else
						{
		
						    cout<<"探测数据ERROR"<<endl;
							corner.push_back(point[j].x);
							corner.push_back(point[j].y);
							corner.push_back(corneracos[i]);
						}
						}else{
							corner.push_back(point[j].x);
							corner.push_back(point[j].y);
							corner.push_back(corneracos[i]);


						}
			    }///if结束
				}
			  }	
		     }///j for循环结束
				  
	      }///i for循环结束
		  
		  return corner;

}




///@brief   the filtercorner function
/*! the filtercorner function  is  to  filter coernerpoint 
*可能有很多个角点在一起 则需要进行过滤
* 设置阈值-- 即在一定区域内选取最大的acos值为角点并保存
*/
///@param   int T  a threshold value
void filter(vector<double> allcorner,vector<double> allcos,int T4)
{
	vector<double>  pp;/// vector pp is to store  the cornerpoint of mismatch condition 
	pp.clear();/// vector clear
	points.clear();
	accos.clear();

	bool* flag=new bool[allcorner.size()];/// define a  bool array to mark the cornerpoint which  match condition
	////!  double circulation is to select  which cornerpoint match condition  in  threshold value area
	///** 双重for循环来进行角点筛选  如果相邻两像素在T4阈值内 则保留acos值大的角点
	//*/
	for(int k=0;k<allcorner.size()-1;k=k+2)///解决相邻角点的问题 若两角点之间的距离小于T4 则保留acos最大的
	{ 

		for(int j=allcorner.size()-2;j>k;j=j-2)///若夹角值为90度则必须保留则本函数值正对直角
		{
			int max=0;
				if((abs(allcorner[j]-allcorner[k])<T4)&&(abs(allcorner[j+1]-allcorner[k+1])<T4))/// 判断相邻两角点是否在阈值区域内
				{

					if((abs(allcos[k/2]-90))>=(abs(allcos[j/2]-90)))
					{max=k;}
					else{max=j;}

					flag[max]=true;///把acos值小的角点在对数组做标记后删除
					
			     }
		}			
	}

	
		///! for circulation 
		/** the for circulation is to put the  points that bigger value of acos  in vector points
		* 把acos 值大的角点存储进vector<CvPoint> points*/
		for(int i=0;i<allcorner.size()-1;i=i+2)
		{
			if(flag[i]==true)/// acos 值小的点
			{
             pp.push_back(allcorner[i]);
			 pp.push_back(allcorner[i+1]);
			 /// the vector pp is to store the the cornerpoint of mismatch condition
			}
			else
			{
				points.push_back(allcorner[i]);
				points.push_back(allcorner[i+1]);
				accos.push_back(allcos[i/2]);
			}///put the the cornerpoint of match condition in vector points
			
		}


}



	///@brief   the transform function
	/*! the transform function  is  to  transform CEdgeHitChain chain to CvSeq* 
	*把数据类型CEdgeHitChain 的轮廓转换为CvSeq* 的轮廓
	*/
	///@param   input CEdgeHitChain chain
    ///@return  the contour of CvSeq* 
	CvSeq* transform(CEdgeHitChain chain)//自定义的transform函数
	{
		
		CEdgeHit *tempHit;///define the  object tempHit  of CEdgeHit class
		CvPoint tempPoint;/// define the object of CvPoint 
	    //CvMemStorage* tempStorage;//CvMemStorage 动态内存存储
		tempStorage = cvCreateMemStorage(0);///初始化内存
		CvSeq* tempContour = cvCreateSeq(CV_SEQ_POLYGON, sizeof(CvSeq), sizeof(CvPoint), tempStorage);//创建CvSeq*序列
		
		tempHit = chain.first;///把CEdgeHitChain 的首地址给CEdgeHit *tempHit
		do ///循环
		{
			tempPoint.x = tempHit->x;///把CEdgeHitChain上的点全部赋值给CvPoint tempPoint
			tempPoint.y = tempHit->y;
			//cout<<"tempHit->x"<<tempHit->x<<"  "<<tempHit->y<<endl;
			//cout<<"tempPoint"<<tempPoint.x<<"  "<<tempPoint.y<<endl;
			cvSeqPush(tempContour ,&tempPoint);///把tempPoint 一个一个写入CvSeq* tempContour
			tempHit = tempHit->next;///继续下一个点
		}while (tempHit!=chain.last->next);//直到CEdgeHitChain尾端跳出循环

		//cvReleaseMemStorage (&tempStorage);
		return tempContour;///返回CvSeq*序列轮廓
	}



public:///  the object of public
    CSixSword ss;
	CFrame fr;
	 CvMemStorage* tempStorage;//CvMemStorage 动态内存存储
	vector<double>  points;//最后输出的角点坐标
	vector<double> accos;//最后输出的角点的反余弦值


	vector<CvPoint> p;
	vector<CvPoint> ppx;
	vector<double> pd;
};