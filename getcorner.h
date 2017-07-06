#include"lib.h"
#include"findcorner.h"
//#include"detectline.h"
#include"getdistanceofcorner.h"
#define R2D    57.295779513082323
#define Deg2Rad( Deg ) (Deg)/ R2D  // �Ƕȵ�����
#define Rad2Deg( Rad ) (Rad)* R2D  // ���ȵ��Ƕ�

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
void  setup(IplImage* frame)///����ԭͼ
{  
	  if(!frame){cout<<"no frame !";}/// if  no frame   will cout  no frame
      ss.setup(frame);///  setup frame to   CSixSword  class
	  fr.set(frame);///transform frame to CFrame class  for read and draw contour in picture
	  ss.transform();/// invoke transform function of CSixSword  class  to get the contour of CEdgeHitChain chain
	 
 }

///@brief   the getcornerpoint function
/*! the getcornerpoint function's  object to  get coernerpoint  */
///@param   IplImage* frame input IplImage* 
void getcornerpoint(IplImage* frame, IplImage* dst,int T1,int T2,int T3,int T4)//T1��T2���ж��������������������μ��������ֵ  T1��T2��Ϊŷʽ����
{

	  CvSeq* temp;///define the  CvSeq* temp (contour)
	  FindCorner fc;///define object  fc of FindCorner class
	  DetectLine dl;///define object dl  of DetectLine class
	  
	  points.clear();///�洢�Ź��˺��ͼƬ�е����нǵ�����
	  accos.clear();
	  vector<double> allcorner;///�洢���������������к�ѡ�ǵ�����
	  allcorner.clear();/// clear the vector

	  vector<double> allacos;///�洢���������������к�ѡ�ǵ�ķ�����ֵ
	  allacos.clear();/// clear the vector

	 vector<CvPoint>  intpoints;///�м����  �洢�Ź��˺�����нǵ�����  ת��ΪCvPoint ������ͼ�ϱ�ǳ��ǵ�
	 intpoints.clear();///clear the vector


	 IplImage *src=cvCloneImage(frame);///  frame clone to src  for  CSixSword ss to setup
	
	 setup(src);///CSixSword ss setup  IplImage *src  to get CEdgeHitChain chain(contour)


	  cout<<"size=="<<ss.v.size()<<endl;///cout the size of contour 

          /// for ѭ��ѭ��������   �ֱ��ÿ��������ȡ�ǵ�
	 	  	for(unsigned int k=0;k<ss.v.size();k++)
			{
               	 vector<CvPoint> cornerpoints;///�洢�ŵ������������к�ѡ�ǵ�
				 vector<double> corneracos;////�洢�ŵ������������к�ѡ�ǵ�ķ�����ֵ
				 cornerpoints.clear();///clear the vector
				 corneracos.clear();///clear the vector
				  vector<double> XY ;
			      XY.clear();

			  /* vector<CvPoint> point;///�洢�ŵ��������ϵ��������ص�
	           point.clear();/// clear the vector*/
			
			  temp=transform(ss.v[k]);/// transform the contour of CEdgeHitChain chain  to contour of CvSeq* 
			 
			  ss.v[k].print(fr);///draw the contour in IplImage *src 
			  
			  /**  invoke the findcorner function of FindCorner class
			  * input  the contour of CvSeq* 
			  * input  int T1
			  *   T1 �� ����double L ����С��Χ��ֵ  ��  L>T1
			  *input  int T2
			  *   T2  ������double L�����Χ��ֵ   ��L<=T2
			  *input   int T3
			  *   T3 �����ø��Ƕ���ֵ��abs(90-acos)<T3  �ǱȽ���90�ȵĽǶ�ֵ��ֵ
			  */
			 fc.findcorner(dst,temp,T1,T2,T3);//��ȡ�����Ľǵ�
			 
			  for(int i=0;i<fc.corner.size();i++)///  ͨ��forѭ����ÿ�������ĺ�ѡ�ǵ���洢�ڼ�����
			  {
				  cornerpoints.push_back(fc.corner[i]);///  using the circulation to put every corner in the vector of cornerpoints
				  corneracos.push_back(fc.accos[i]);/// using the circulation to put every corner's the value of acos in the vector of cornerpoints
			  }

			  for(int i=0;i<cornerpoints.size();i++)
			  {
				 // cout<<"cornerpoints"<<corneracos[i]<<"   "<<cornerpoints[i].x<<"    "<<cornerpoints[i].y<<endl;
				  cvCircle(dst,cornerpoints[i], 5, CV_RGB(255,0,0), 1,8,0);
			  }

			  fc.setup(temp); ///����findcontour�����setup��������õ�ǰ���������еĵ�
			  //cout<<"point.size== "<<fc.points.size()<<endl;

		
			if(cornerpoints.size()!=0){XY=detectp(fc.points,cornerpoints,corneracos,frame,dst);}
					
			
			 if(XY.size()!=0)///XY��ֵ
			 {
			 for(int i=0;i<XY.size()-2;i=i+3)
			 { 
			    allcorner.push_back(XY[i]);///��x����洢��allcorner   Ϊʲô����CvPoint  ����ΪCvPointֻ֧��int 
				allcorner.push_back(XY[i+1]);///��x����洢��allcorner
				allacos.push_back(XY[i+2]);///�ѽǵ㷴����ֵ�洢��allacos
				
			 }
			 }

			 for(int i=0;i<fc.points.size()-1;i++){cvLine( src,fc.points[i], fc.points[i+1], CV_RGB(255,0,0),1, 8, 0 );}/// ͨ�������ϵĵ�����������		
 
		      cvReleaseMemStorage (&tempStorage);///�ͷ�transform�����ϴ������ڴ档
			}

			//! invoke the filtercorner function
			/** input the threshold value T3��T4 
			*T3����ֵ  �������յõ�ÿ�������Ľǵ��������Ҫѡ��abs(90-acos)<T3�ĽǶ���ֵ
			*T4����ֵ������T4��λ����֮��  ɸѡ��acosֵ�����Ϊ�ǵ㱣��������Ϊ�ǵ�Ϊ90�ȱ��뱣��
			*/
			if(allcorner.size()!=0)
			{
			    filter(allcorner,allacos,T4);///���˽ǵ�
			}


          if(points.size()!=0)
		  {
			for(int i=0;i<points.size()-1;i=i+2)
			{
				if((points[i]>0)&&(points[i+1]>0))///��points�洢��intpoints�Ƿ�����ͼ�ϻ������սǵ㡣
				{
				CvPoint pp;
				pp.x=points[i];
				pp.y=points[i+1];
				intpoints.push_back(pp);
				}
            }	

//Ϊ���ж��ȶ�������Ҫͨ������ıȽ�ѡ����Ҫ��ȡ�ĵ�����
			//�˴��� ����y�������Ľǵ�����
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
					cvCircle(dst,intpoints[i],10, CV_RGB(0,0,255), 1,8,0);/// ��ͼƬ�ϻ����ǵ㣬�Խǵ�Ϊ���ģ�5������Ϊ�뾶��Բ��
				} 
		    }
		 }
		
      
       
		cvNamedWindow( "dst", CV_WINDOW_AUTOSIZE );	 ///define a window
        cvShowImage( "dst", dst );///  show  the  IplImage *dst
		cvNamedWindow( "srcc", CV_WINDOW_AUTOSIZE );	 ///define a window
        cvShowImage( "srcc", src );///  show  the  IplImage *dst

		cvReleaseImage(&src);

}


vector<double>  FitLine(  IplImage* dst,vector<double> xdetectpoint, vector<double> ydetectpoint)//����ͼƬ��̽���Ľǵ����ߵı�Ե���ص�����
{
	       vector<double> xy;//�洢ֱ����Ϻ����ֱ�ߵĽ���
		   xy.clear();

			float *liner = new float[4];///��������
			float *linel = new float[4];
			CvMemStorage* storage1 = cvCreateMemStorage(0);///�����洢�ռ�
			CvMemStorage* storage2 = cvCreateMemStorage(0);
			CvSeq* point_seqr = cvCreateSeq( CV_32FC2,sizeof(CvSeq), sizeof(CvPoint2D32f), storage1 );///��������
			CvSeq* point_seql = cvCreateSeq( CV_32FC2,sizeof(CvSeq), sizeof(CvPoint2D32f), storage2 );
			for(int s=0;s<xdetectpoint.size()-1;s=s+2)
			{
			cvSeqPush(point_seqr, &cvPoint2D32f(xdetectpoint[s],xdetectpoint[s+1]));///�����ص㾫ȷ����push��point_seqr
			}
			for(int s=0;s<ydetectpoint.size()-1;s=s+2)
			{
			cvSeqPush(point_seql, &cvPoint2D32f(ydetectpoint[s],ydetectpoint[s+1]));///�����ص㾫ȷ����push��point_seql
			}

			cvFitLine(point_seqr,CV_DIST_HUBER ,0,0.001,0.001,liner);///ֱ����Ϻ���   x ����
			//cvFitLine(point_seql,CV_DIST_L2,0,0.001,0.001,linel);///��ydetectpoint�����ص����ֱ�����   y ����
			cvFitLine(point_seql,CV_DIST_HUBER ,0,0.001,0.001,linel);
			//double kr=atan(liner[1]/liner[0]);
			//double kl=atan(linel[1]/linel[0]);
			//double angle1= Rad2Deg(kr);
			//double angle2=Rad2Deg(kl);

			
			
			double x1=liner[2]; double y1=liner[3];/// line[0--3] �ֱ�Ϊ (vx, vy,x0, y0)
			double x2=linel[2]; double y2=linel[3];

			//cout<<"liner[0]="<<liner[0]<<"  "<<liner[1]<<endl;
			//cout<<"linel[0]="<<linel[0]<<"  "<<linel[1]<<endl;
			
			//cout<<"x1="<<x1<<"  "<<y1<<endl;
			//cout<<"x2="<<x2<<"  "<<y2<<endl;
			//cout<<"angle1 "<<angle1<<"  x= "<<x1<<"   "<<y1<<endl;
			//cout<<"angle2 "<<angle2<<"  x=  "<<x2<<"   "<<y2<<endl;
			double k1=liner[1]/liner[0];double b11=liner[3]-k1*liner[2];
			double k2=linel[1]/linel[0];double b22=linel[3]-k2*linel[2];///ȷ��ֱ��б�ʼ�ֱ�߷���


			//cout<<"k1=="<<k1<<"  "<<b11<<"   "<<k2<<"  "<<b22<<"  "<<endl;

			xy.push_back((y2-y1+k1*x1-k2*x2)/(k1-k2));///������ֱ����ǵ�			
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
     y_mean /= size; //���ˣ�������� x y �ľ�ֵ

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
         if( fabs(Dxx / Dyy - 1) < 1e-5) //��ʱû��һ�������ֱ�߷����޷����
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
	for(int g=0;g<accurate.size();g++)//ȷ��y��Ƕ�
	{
		if(g<4)///����Ƕ�����ѡ��g��g+1֮����ȷ��б��
		{

			if(accurate[g].x==accurate[g+4].x){degree=0;}//������˵��x���������˵��̽���߽Ƕ�Ϊ0��

			if(accurate[g].x!=accurate[g+4].x)///��������ص�x���겻�����ͨ��б�ʹ�ʽ��������̽���߽Ƕȼ��ɡ�
			{
							
				double slopey=(-1)*(accurate[g].y-accurate[g+4].y)/(accurate[g].x-accurate[g+4].x);
				degree=atan((double)(-1/slopey))/CV_PI*180;

			}
									 
		}
		else
		{

			if(accurate[g].x==accurate[g-4].x){degree=0;}//������˵��x���������˵��̽���߽Ƕ�Ϊ0��

			if(accurate[g].x!=accurate[g-4].x)///��������ص�x���겻�����ͨ��б�ʹ�ʽ��������̽���߽Ƕȼ��ɡ�
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
	vector<double> cd;//�ǵ�֮��ľ���
	cd.clear();
	vector<double>  dist;
	dist.clear();
	
	  for(int d=0;d<cornerpoints.size();d++)//����ǵ�֮����໥����
	{
		if((abs(coener.x-cornerpoints[d].x)>5)||(abs(coener.y-cornerpoints[d].y)>5))
		{ double dist=sqrt((double)((coener.x-cornerpoints[d].x)*(coener.x-cornerpoints[d].x)+(coener.y-cornerpoints[d].y)*(coener.y-cornerpoints[d].y)));
		cd.push_back(dist);}
	}
	  double min=cd[0];//�൱ǰ�ǵ���С����
	  double smin=cd[1];//�൱ǰ�ǵ�����Сֵ
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
	//..........����������.......................//
	     //int L1=5;int L2=25;
	      int L1,L2x,L2y;//ȷ��̽���ߵ��ܳ�
		 DetectLine dl;///����DetectLine ��Ķ��� dl
		// CvPoint  point1,point2,point3,point4;///����CvPoint ����  ������ֱ���ϵ��ĸ���
		 vector<double > corner;
		 corner.clear();///����洢��ֱ�߽����xy����ļ���
		 

		  double x=0;double y=0;///
		  int k=1; int errord=5;
	//.........................................................................................................................................//

    //ͨ���ǵ�֮��ĳ�����ȷ��̽�����ܹ��ĳ���L1-L2
		  vector<double> ds;
		  ds.clear();

   //......................................................................................................................................................//

			for(int i=0;i<cornerpoints.size();i++)///˫��ѭ����ȷ����ֱ�ߵĽ������꣬���ǵ�ľ�ȷ����
		   {
			   
			   			  vector<double> xdetectpoint;///���峣�����������ʼ����
						  vector<double> ydetectpoint;//�洢y����̽��������
						  xdetectpoint.clear();//��ռ���
						  ydetectpoint.clear();//��ռ���
						  vector<double> dist;
		                  dist.clear();
						  
                       dist=comparedis(cornerpoints[i],cornerpoints);
					   L1=5;L2x=dist[1];L2y=dist[0];  
						//cout<<"L2x=="<<L2x<<"  "<<L2y<<endl;  
             //if(k==1)
			{
			   for(int j=0;j<point.size();j++)///���ýǵ�������ϵĵ�������ѭ��ּ�ڷֱ��ҳ��ǵ��������˵����ص� point--���������е�����
			   {

				  				   
				  // ...........................................................................................................//
				   		vector<CvPoint> dispoint;///   vector<CvPoint> dispoint ��������������Ҷ�����нǵ��������ص�����		
						vector<double> disL;///   vector<double> disL ��������������Ҷ�����нǵ��������ص�֮��ľ���		 
						dispoint.clear();///���vector<CvPoint>  dispoint ����		
						disL.clear(); ///���vector<double>  dispoint ����
						
                      
				          vector<CvPoint> xaccurate;//�洢x����̽�����ص�
						  vector<CvPoint> yaccurate;//�洢y����̽�����ص�
						  xaccurate.clear();
						  yaccurate.clear();

						  double degreey;//y����̽���߽Ƕ�
						  double degreex;//x����̽���߽Ƕ�
						
						   vector<double>xyx;
						   vector<double>xyy;
						   xyx.clear();
						   xyy.clear();


						   vector<double> XY;//ֱ�߽�������
						   XY.clear();

						   
                //..................................................................................//
						 
		     // if((abs(cornerpoints[i].x-202)<5)&&(abs(cornerpoints[i].y-339)<5))
				{
				  // k=k+1;
				if((point[j].x==cornerpoints[i].x)&&(point[j].y==cornerpoints[i].y))///��������ϵĵ�������ڽǵ����������ѭ��
				  {
					  
		  				fstream pppx;

						pppx.clear();

					 pppx.open("pppx.txt",ios::app);
	

					  fstream pppy;

						pppy.clear();

					 pppy.open("pppy.txt",ios::app);
					  //.....................................................................................................................//
						for (int d=point.size()-1; d>=0; d--)///forѭ�� ���ص��ĩ�������������ڱ��������ص������֮���ŷ�Ͼ���
						{		 			
							double diss=((cornerpoints[i].x-point[d].x)*(cornerpoints[i].x-point[d].x)+(cornerpoints[i].y-point[d].y)*(cornerpoints[i].y-point[d].y));///double diss Ϊ�����ص�֮��ľ���  x*x+y*y
			
							double L=sqrt(diss);///double L Ϊ������֮���ŷʽ���� sqrt(x*x+y*y)

								if((L>=L1)&&(L<=L2x))///if ����������ɸѡ  ѡ��������֮���ŷʽ������ T1<L<=T2�������ص� 
								{
									disL.push_back(L);///����������ŷ�Ͼ���洢��vector<double> disL
									dispoint.push_back(point[d]);///�������������ص�����洢�� vector<CvPoint> dispoint	
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
								//}//˵����y�᷽�����ص�
								//if((abs(dispoint[g].x-point[j].x)>=L1)&&(abs(dispoint[g].y-point[j].y)<5))
								//{xaccurate.push_back(dispoint[g]);pppx<<dispoint[g].x<<"  "<<dispoint[g].y<<endl;
								//}//˵����x�᷽�����ص�
                                 //.......................����xy��ȫ���뻯...........................................//
								
								if((abs(dispoint[g].x-point[j].x)<=L1)&&(abs(dispoint[g].y-point[j].y)>L1)&&(abs(dispoint[g].y-point[j].y)<(L2y-5)))
								{yaccurate.push_back(dispoint[g]);//pppy<<dispoint[g].x<<"  "<<dispoint[g].y<<endl;
								}//˵����y�᷽�����ص�
								if((abs(dispoint[g].x-point[j].x)>L1)&&((abs(dispoint[g].x-point[j].x)<=(L2x-5))&&(abs(dispoint[g].y-point[j].y)<=L1)))
								{xaccurate.push_back(dispoint[g]);//pppx<<dispoint[g].x<<"  "<<dispoint[g].y<<endl;
								}//˵����x�᷽�����ص�

							}

							if((yaccurate.size()>=4)&&(xaccurate.size()>=4))
							{
								int dky=0; int dkx=0;
								if(((yaccurate.size()>=8)&&(xaccurate.size()>=8))){ dkx=4;}//�������ص�֮��ľ�����ȷ��̽���߽Ƕ�
								if((yaccurate.size()<8)&&(xaccurate.size()>=8)){ dky=yaccurate.size()/2;}
								if((xaccurate.size()<8)&&(yaccurate.size()>=8)){dkx=xaccurate.size()/2;}
							for(int g=0;g<yaccurate.size();g++)//ȷ��y��Ƕ�
							{
								
								if(g<dky)///����Ƕ�����ѡ��g��g+1֮����ȷ��б��
								{   
																
									if(yaccurate[g].x==yaccurate[g+dky].x){degreey=0;}//������˵��x���������˵��̽���߽Ƕ�Ϊ0��

									if(yaccurate[g].x!=yaccurate[g+dky].x)///��������ص�x���겻�����ͨ��б�ʹ�ʽ��������̽���߽Ƕȼ��ɡ�
									{							          
									  double slopey=(-1)*(yaccurate[g].y-yaccurate[g+dky].y)/(yaccurate[g].x-yaccurate[g+dky].x);
									  degreey=atan((double)(-1/slopey))/CV_PI*180;//ȷ��̽���߽Ƕ�
									}
									
								}
								else
								{   
									if(yaccurate[g].x==yaccurate[g-dky].x){degreey=0;}//������˵��x���������˵��̽���߽Ƕ�Ϊ0��

									if(yaccurate[g].x!=yaccurate[g-dky].x)///��������ص�x���겻�����ͨ��б�ʹ�ʽ��������̽���߽Ƕȼ��ɡ�
									{							
									   double slopey=(-1)*(yaccurate[g].y-yaccurate[g-dky].y)/(yaccurate[g].x-yaccurate[g-dky].x);			                     
									    degreey=atan((double)(-1/slopey))/CV_PI*180;//ȷ��̽���߽Ƕ�
									}

								}
								
								
								//pppy<<"degreey=="<<degreey<<endl;
							    xyy=dl.detectlinedegree(yaccurate[g],degreey,5,src,dst);///����detectlinedegree������ȡ�����ص㾫ȷxy����
								ydetectpoint.push_back(xyy[0]);///����ȡ�����ص㾫ȷ�����xy����洢��ydetectpoint
								ydetectpoint.push_back(xyy[1]);
								p.push_back(yaccurate[g]);//�洢��̽������ص�����
								pd.push_back(xyy[0]);//�洢̽�������ص�����
								pd.push_back(xyy[1]);
							
								//pppy<<xyy[0]<<"   "<<xyy[1]<<endl;
							}
							
							for(int k=0;k<xaccurate.size();k++)///ȷ��x�᷽������ص�ĽǶ�
							{
							
								if(k<dkx)
								{   
									 if(xaccurate[k].x==xaccurate[k+dkx].x){degreex=-90;}/////������˵��x���������˵��̽���߽Ƕ�Ϊ90��
									 else ///��������ص�x���겻�����ͨ��б�ʹ�ʽ��������̽���߽Ƕȼ���
									 {
										
									double slopex=(-1)*(xaccurate[k].y-xaccurate[k+dkx].y)/(xaccurate[k].x-xaccurate[k+dkx].x);//��ȡֱ��б��
									 degreex=atan((double)(-1/slopex))/CV_PI*180;//����ɽǶ�
									 }								
								}
								else///�Ƕ������ص�
								{ 
									if(xaccurate[k].x==xaccurate[k-dkx].x){degreex=-90;}//������˵��x���������˵��̽���߽Ƕ�Ϊ90��
									else//��������ص�x���겻�����ͨ��б�ʹ�ʽ��������̽���߽Ƕȼ���
									{
										
							         double slopex=(-1)*(xaccurate[k].y-xaccurate[k-dkx].y)/(xaccurate[k].x-xaccurate[k-dkx].x);
									 degreex=atan((double)(-1/slopex))/CV_PI*180;
									}
									 
								}
								//pppx<<"degreex=="<<degreex<<endl;
								xyx=dl.detectlinedegree(xaccurate[k],degreex,5,src,dst);//����detectlinedegree������ȡ���ص㾫ȷxy����
								xdetectpoint.push_back(xyx[0]);//����ȡ�����ص㾫ȷ�����xy����洢��xdetectpoint
							    xdetectpoint.push_back(xyx[1]);
								p.push_back(xaccurate[k]);//�洢��̽������ص�����
								pd.push_back(xyx[0]);//�洢̽���ı�Ե���ص�����
								pd.push_back(xyx[1]);
								pppx<<xyx[0]<<"  "<<xyx[1]<<endl;
						}
//
//
					}//for
					    if((xdetectpoint.size()!=0)&&(ydetectpoint.size()!=0))///ȷ���нǵ��ִ��
						{
							XY=FitLine(dst,xdetectpoint, ydetectpoint);
						}
						//cout<<"XY=="<<XY[0]<<"  "<<XY[1]<<"   cornerpoints="<<cornerpoints[i].x<<" "<<cornerpoints[i].y<<endl;
///.........................................................................................................................//
						if(XY.size()!=0)
						{
					 if((abs(XY[0]-cornerpoints[i].x)<errord)&&(abs(XY[1]-cornerpoints[i].y)<errord))//�������ľ�ȷ�ǵ�������ԭ�ǵ�����������5�򷵻�ԭ�ǵ�����
						{
							corner.push_back(XY[0]);
							corner.push_back(XY[1]);
							corner.push_back(corneracos[i]);
						}

					else
						{
		
						    cout<<"̽������ERROR"<<endl;
							corner.push_back(point[j].x);
							corner.push_back(point[j].y);
							corner.push_back(corneracos[i]);
						}
						}else{
							corner.push_back(point[j].x);
							corner.push_back(point[j].y);
							corner.push_back(corneracos[i]);


						}
			    }///if����
				}
			  }	
		     }///j forѭ������
				  
	      }///i forѭ������
		  
		  return corner;

}




///@brief   the filtercorner function
/*! the filtercorner function  is  to  filter coernerpoint 
*�����кܶ���ǵ���һ�� ����Ҫ���й���
* ������ֵ-- ����һ��������ѡȡ����acosֵΪ�ǵ㲢����
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
	///** ˫��forѭ�������нǵ�ɸѡ  ���������������T4��ֵ�� ����acosֵ��Ľǵ�
	//*/
	for(int k=0;k<allcorner.size()-1;k=k+2)///������ڽǵ������ �����ǵ�֮��ľ���С��T4 ����acos����
	{ 

		for(int j=allcorner.size()-2;j>k;j=j-2)///���н�ֵΪ90������뱣���򱾺���ֵ����ֱ��
		{
			int max=0;
				if((abs(allcorner[j]-allcorner[k])<T4)&&(abs(allcorner[j+1]-allcorner[k+1])<T4))/// �ж��������ǵ��Ƿ�����ֵ������
				{

					if((abs(allcos[k/2]-90))>=(abs(allcos[j/2]-90)))
					{max=k;}
					else{max=j;}

					flag[max]=true;///��acosֵС�Ľǵ��ڶ���������Ǻ�ɾ��
					
			     }
		}			
	}

	
		///! for circulation 
		/** the for circulation is to put the  points that bigger value of acos  in vector points
		* ��acos ֵ��Ľǵ�洢��vector<CvPoint> points*/
		for(int i=0;i<allcorner.size()-1;i=i+2)
		{
			if(flag[i]==true)/// acos ֵС�ĵ�
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
	*����������CEdgeHitChain ������ת��ΪCvSeq* ������
	*/
	///@param   input CEdgeHitChain chain
    ///@return  the contour of CvSeq* 
	CvSeq* transform(CEdgeHitChain chain)//�Զ����transform����
	{
		
		CEdgeHit *tempHit;///define the  object tempHit  of CEdgeHit class
		CvPoint tempPoint;/// define the object of CvPoint 
	    //CvMemStorage* tempStorage;//CvMemStorage ��̬�ڴ�洢
		tempStorage = cvCreateMemStorage(0);///��ʼ���ڴ�
		CvSeq* tempContour = cvCreateSeq(CV_SEQ_POLYGON, sizeof(CvSeq), sizeof(CvPoint), tempStorage);//����CvSeq*����
		
		tempHit = chain.first;///��CEdgeHitChain ���׵�ַ��CEdgeHit *tempHit
		do ///ѭ��
		{
			tempPoint.x = tempHit->x;///��CEdgeHitChain�ϵĵ�ȫ����ֵ��CvPoint tempPoint
			tempPoint.y = tempHit->y;
			//cout<<"tempHit->x"<<tempHit->x<<"  "<<tempHit->y<<endl;
			//cout<<"tempPoint"<<tempPoint.x<<"  "<<tempPoint.y<<endl;
			cvSeqPush(tempContour ,&tempPoint);///��tempPoint һ��һ��д��CvSeq* tempContour
			tempHit = tempHit->next;///������һ����
		}while (tempHit!=chain.last->next);//ֱ��CEdgeHitChainβ������ѭ��

		//cvReleaseMemStorage (&tempStorage);
		return tempContour;///����CvSeq*��������
	}



public:///  the object of public
    CSixSword ss;
	CFrame fr;
	 CvMemStorage* tempStorage;//CvMemStorage ��̬�ڴ�洢
	vector<double>  points;//�������Ľǵ�����
	vector<double> accos;//�������Ľǵ�ķ�����ֵ


	vector<CvPoint> p;
	vector<CvPoint> ppx;
	vector<double> pd;
};