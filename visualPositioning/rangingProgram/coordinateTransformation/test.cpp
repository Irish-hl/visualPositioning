#include<iostream>
#include<cstdio>
#include<iomanip>
#include<windows.h>
#include<math.h>
#include<fstream>
#include<cstdlib>
double PI=3.1415926;
double highRate=0.20;
double lowRate=0.05;
using namespace std;
int main()
{
	double H=17;//������߶�					
	double y1=13;//�������ֱ�ӽ�ͶӰ�ڵ�����������					
	double y12=159;//�������ֱ�ӽ�ͶӰ�ڵ������Զ����					
	double a=52.6/180*PI;//�������ֱ�ӽ����ƽ��y��ļнǵ����ֵ					
	double b=6.1/180*PI;//�������ֱ�ӽ����ƽ��y��ļнǵ���Сֵ
	double c=85.0/180*PI;//�����ˮƽ�ӽ��ڵ����ϵ�ͶӰ��y��ļн�					
	double Sx=320;//����ƽ������������������xֵ
	double Sy=240;//����ƽ������������������yֵ

	double L;//����Ŀ�����������ڵ���ͶӰ֮��ľ���					
	//double angle;//p����С��ǰ���������ɼн�					
	double u,v;//���ص���ͼƬ�е�����
	double ce_x,ce_y;//�����������������ϵ�е�����
	double real_x,real_y;//�����������������ϵ�е�ʵ������
	//double ce_a,ce_b;//���������������ϵ�е�����					
	double real_a,real_b;//���������������ϵ�е�ʵ������
	double L1,num=0,jun=0;

	ifstream fin("����1.txt",ios::in);
	if(!fin){cerr<<"�޷���"<<endl;exit(1);}
		while(fin>>real_a>>real_b>>real_x>>real_y>>u>>v)
		{
			//ce_y=H*tan((PI/2-a+(a-b)/2+(v-Sy/2)*(a-b)/2/(Sy/2)));
			//��������ģ�͵�ce_y
			ce_y=H*tan((PI/2-a+(a-b)/2+atan((v-Sy/2)/(Sy/2)*tan((a-b)/2))));
			double ce_x_1=ce_y*tan((Sx/2-u)/(Sx/2)*c);
			double ce_x_2=ce_y*tan((Sx/2-u)/v);
			double ce_x_3=u*sqrt(y1*y1+H*H)/sqrt(H*H+y1*y1+u*u-cos(PI/2-a/2+b/2)*2*u*sqrt(H*H+y1*y1));
			double ce_x_4=ce_y*tan((u-Sx/2)*c/(Sx/2));
			//��������ģ�͵�ce_x
			double ce_x_5=ce_y*tan(atan((u-Sx/2)/(Sx/2)*tan(c/2)));
			ce_x=(ce_x_1+ce_x_2+ce_x_3+ce_x_4+ce_x_5)/3;
			//���������ϵ��Ŀ�������㵽ԭ��ľ���
			L=sqrt(ce_x_5*ce_x_5+ce_y*ce_y);
			//��������ϵ��Ŀ�������㵽���������ϵԭ��ľ���
			L1=sqrt((real_a-real_x)*(real_a-real_x)+(real_b-real_y)*(real_b-real_y));
			//������ʣ�-rate��rate�������ݺ�ɫ��ʾ
			if(fabs((L1-L)/L1)>=highRate)
			{
				SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_RED);
				cout<<"�����("<<setprecision(2)<<real_a<<","<<setprecision(2)<<real_b<<"),  �����㣨"<<setprecision(2)<<real_x<<","<<setprecision(2)<<real_y<<"),  ���ص�("<<setprecision(3)<<u<<","<<setprecision(3)<<v<<"),  ""��þ��룺"<<setprecision(5)<<L<<",  ʵ�ʾ��룺"<<setprecision(5)<<L1<<",  ��"<<setprecision(5)<<L1-L<<",  ����ʣ�"<<setprecision(5)<<(L1-L)/L1<<endl;	
				SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),7);
			}
			else if(fabs((L1-L)/L1)<=lowRate)
			{
				SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_GREEN);
				cout<<"�����("<<setprecision(2)<<real_a<<","<<setprecision(2)<<real_b<<"),  �����㣨"<<setprecision(2)<<real_x<<","<<setprecision(2)<<real_y<<"),  ���ص�("<<setprecision(3)<<u<<","<<setprecision(3)<<v<<"),  ""��þ��룺"<<setprecision(5)<<L<<",  ʵ�ʾ��룺"<<setprecision(5)<<L1<<",  ��"<<setprecision(5)<<L1-L<<",  ����ʣ�"<<setprecision(5)<<(L1-L)/L1<<endl;
				SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),7);
			}
			else
				cout<<"�����("<<setprecision(2)<<real_a<<","<<setprecision(2)<<real_b<<"),  �����㣨"<<setprecision(2)<<real_x<<","<<setprecision(2)<<real_y<<"),  ���ص�("<<setprecision(3)<<u<<","<<setprecision(3)<<v<<"),  ""��þ��룺"<<setprecision(5)<<L<<",  ʵ�ʾ��룺"<<setprecision(5)<<L1<<",  ��"<<setprecision(5)<<L1-L<<",  ����ʣ�"<<setprecision(5)<<(L1-L)/L1<<endl;
			//ͳ�������λ�ڣ�-rate��rate��֮�������ʾ���ֵƽ��ֵ
			if(fabs((L1-L)/L1)<highRate)
			{	
				jun+=fabs(L1-L)/L1;
				num++;
			}
		}
	    cout<<"ƽ������ʣ�"<<jun/num<<endl;
		fin.close();
	return 0;
}