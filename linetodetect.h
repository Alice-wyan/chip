//
//linetocorner����ּ������houghlinesp�������ͼƬֱ�ߺ�ó�ֱ�����˶˵㣬
//���ö˵����������˵�֮�仭̽���ߣ�ͨ��̽���߻�ȡ�� ��Ե��������ֱ����ϣ�ͨ����Ϻ��ֱ����ȡ����ֱ�ߵĽ���
//ͨ�����������ֱ�ߵĽǶ�������̽���߽��б�Ե���ص��̽�⣬���ٽ���ֱ����Ϻ���ȡ����ȷ����Ľ��㡣
//


#include"lib.h"

class LineToDetect{

public:


void linetocorner(IplImage *src,IplImage *dst,vector<Vec4i> lines)//linetocorner��������ԭͼ�����뻮�ߺ���ڹ۲��ͼ���ͼ���ֱ�߲���
	{

		TestDetectPoint  tdp;//��������̽����λ�õ���
		
		CommomFunction cf;
		double degree1=0;double degree2=0;
		
		
		vector<double> b;//�洢ֱ�߽ؾ�
		b.clear();
	
		vector<float *> line;//�洢���ֱ�߲���
		line.clear();

		vector<float *> linee;//�洢���ֱ�߲���
		linee.clear();
		vector<double> degree;
		degree.clear();
		int distance=0;	//���ֱ�ߵĳ���
		for(int i=0;i<lines.size();i++)//ѭ��ֱ�ߵ�����
		{
			double k=0; double degree=0; 
			Vec4i l= lines[i]; 
						
			CvPoint ll;//x����С��ֱ�߶˵�
			CvPoint lr; //x������ֱ�߶˵�   ��x�������ʱ����Ƚ�y����
			if(l[0]<l[2]){ ll.x=l[0];ll.y=l[1];lr.x=l[2];lr.y=l[3];} //�˵�ll ������ֱ����x����С�Ķ˵�  ������ʹ�Ƕȷ�Χ��(-90~90)
			else if(l[0]==l[2]){if(l[1]<l[3]){ll.x=l[0];ll.y=l[1];lr.x=l[2];lr.y=l[3];}else{ll.x=l[2];ll.y=l[3];lr.x=l[0];lr.y=l[1];}}
			else{ ll.x=l[2];ll.y=l[3];lr.x=l[0];lr.y=l[1];}

			 if(ll.x==lr.x){ degree=-90;}			
			 else
			 { k=(l[1]-l[3])/(l[0]-l[2]);
			   if(k==0){degree=0;}
			   else{degree=atan((double)(k))/CV_PI*180.0;}}//degree Ϊֱ�߽Ƕ�
			
			
			 distance=sqrt((double)((l[0]-l[2])*(l[0]-l[2])+(l[1]-l[3])*(l[1]-l[3])));//���˵�֮��ĳ���

			// line.push_back( tdp.getdetcet(src,dst,ll,150,degree));//����̽��������̽������ȡ��Ե���ص㲢���ֱ��
	
		}

		//if(line.size()==2)
		//{
		//  vector<double> xy=cf.getcovergepoint(line[0],line[1]);//��ȡ�����ֱ�ߵĽ���


		//     fstream ph;

	 //        ph.clear();

		//     ph.open("cx.txt",ios::app);
		//	 fstream phy;

	 //        phy.clear();

		//     phy.open("cy.txt",ios::app);

		//	 CvPoint p;
		//	 p.x=xy[0];p.y=xy[1];
		//	// cvCircle(dst,p,10, CV_RGB(0,255,0), 1,8,0);//��dst�ϻ���

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
		////cvCircle(dst,point,10, CV_RGB(0,255,0), 1,8,0);//��dst�ϻ�����--��ɫ
		//double tmp=0;
		//if(abs(degree[0])<abs(degree[1])){tmp=degree[0];degree[0]=degree[1];degree[1]=tmp;}//degree[0]�ǽǶȱȽϴ�ĽǶ�
  //    
		////cout<<"degree[0]="<<degree[0]<<"  "<<degree[1]<<endl;

		//tdp.getdetectedge(src,dst,point,distance-40,degree[0],degree[1]);
		//}
		







		//for(int j=0;j<sitaa.size();j++)
		//{
		//	double b1=lines[j][1]-sitaa[j]*lines[j][0];

		//	b.push_back(b1);
		//}

		//double x=(b[1]-b[0])/(sitaa[0]-sitaa[1]);//����x
		//double y=sitaa[0]*x+b[0];//����y
		//

		//cout<<"sita="<<sitaa[0]<<"  "<<sitaa[1]<<endl;
		//
		//cout<<"x=="<<x<<"  "<<"y=="<<y<<endl;
		//
		//CvPoint xy;
		//xy.x=x;
		//xy.y=y;

		////cvCircle(dst,xy,10, CV_RGB(0,255,0), 1,8,0);//��dst�ϻ���
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