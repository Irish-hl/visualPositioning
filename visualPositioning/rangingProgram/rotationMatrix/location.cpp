#include <iostream>
#include <iomanip>
#include<math.h>
using namespace std;

int main()
{
	double p_cx,p_cy;
	double p_wx,p_wy;
	double q_cx,q_cy;
	double q_wx,q_wy;
	double 
	double cos,sin,tx,ty;

	cout<<"�������һ��������������ϵ�µ�����"<<endl;
	cin>>p_cx>>p_cy;
	cout<<"�������һ�������������ϵ�µ�����"<<endl;
	cin>>p_wx>>p_wy;
	cout<<"������ڶ���������������ϵ�µ�����"<<endl;
	cin>>q_cx>>q_cy;
	cout<<"������ڶ��������������ϵ�µ�����"<<endl;
	cin>>q_wx>>q_wy;
	//(p_wx-q_wx)*cos-(p_wy-q_wy)*sin==(p_wx-q_wx)*sin+(p_wy-q_wy)*cos;
	sin=sqrt(1-cos*cos);
	cos=(p_wx-q_wx+p_wy-q_wy)/(p_wx-q_wx-p_wy+q_wy)*sin;
	cout<<"��ת����Ϊ:"<<endl<<cos<<" "<<-sin<<endl<<sin<<" "<<cos<<endl;
}





