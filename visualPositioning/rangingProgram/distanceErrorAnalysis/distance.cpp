#include<iostream>
#include<math.h>
#include<fstream>
#include<cstdlib>
double PI=3.1415926;
using namespace std;
int main()
{
	double H=17;//������߶�					
	double y1=13;//�������ֱ�ӽ�ͶӰ�ڵ�����������					
	double y12=159;//�������ֱ�ӽ�ͶӰ�ڵ������Զ����					
	double a=52.6/180*PI;//�������ֱ�ӽ����ƽ��y��ļнǵ����ֵ					
	double b=6.1/180*PI;//�������ֱ�ӽ����ƽ��y��ļнǵ���Сֵ
	double c=47.5/180*PI;//�����ˮƽ�ӽ��ڵ����ϵ�ͶӰ��y��ļн�					
	double Sx=320;//����ƽ������������������xֵ
	double Sy=240;//����ƽ������������������yֵ

	double L;//����Ŀ�����������ڵ���ͶӰ֮��ľ���					
	double angle;//p����С��ǰ���������ɼн�					
	double u,v;//���ص���ͼƬ�е�����
	double ce_x,ce_y;//�����������������ϵ�е�����
	double real_x,real_y;//�����������������ϵ�е�ʵ������
	double ce_a,ce_b;//���������������ϵ�е�����					
	double real_a,real_b;//���������������ϵ�е�ʵ������
	double L1,num=0,jun=0;

	ifstream fin("����1.txt",ios::in);
	  if(!fin){cerr<<"�޷���"<<endl;exit(1);}
	   while(fin>>real_a>>real_b>>real_x>>real_y>>u>>v)
	   {
		//cout<<"real_a  real_b  real_x  real_y  u  v:"<<endl;
		//cin>>real_a>>real_b>>real_x>>real_y>>u>>v;
		ce_y=H*tan(PI/2-a/2-b/2-atan((Sy-2*v)/Sy*tan(a/2-b/2)));
		double ce_x_1=ce_y*tan((Sx/2-u)/(Sx/2)*c);
		double ce_x_2=ce_y*tan((Sx/2-u)/v);
		double ce_x_3=u*sqrt(y1*y1+H*H)/sqrt(H*H+y1*y1+u*u-cos(PI/2-a/2+b/2)*2*u*sqrt(H*H+y1*y1));
		ce_x=(ce_x_1+ce_x_2+ce_x_3)/3;
		ce_a=real_x-ce_x;
		ce_b=real_y-ce_y;
		//cout<<"ʵ��("<<real_a<<","<<real_b<<")  "<<"����("<<ce_a<<","<<ce_b<<")"<<endl;
		L=sqrt(ce_x*ce_x+ce_y*ce_y);
		L1=sqrt((real_a-real_x)*(real_a-real_x)+(real_b-real_y)*(real_b-real_y));
		cout<<"L="<<L<<"    ,ʵ��"<<L1<<"     ��"<<L1-L<<"     ����ʣ�"<<(L1-L)/L1<<endl;
		if(((L1-L)/L1)<0.3&&((L1-L)/L1)>0)
		{	jun+=(L1-L)/L1;
			num++;
		}
			if(((L1-L)/L1)<0)
		{	jun-=(L1-L)/L1;
			num++;
		}
	   }
	   cout<<"ƽ������ʣ�"<<jun/num<<endl;
	return 0;
}
