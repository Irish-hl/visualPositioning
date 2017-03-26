#ifndef _matrix_h
#define _matrix_h
#include <iostream>
#include <iomanip>
#include <cmath>
#define eps (1e-06)
#define KMAX 60  //����������
using namespace std;
int cifang(int i,int j); //ȫ�ֺ�������-1��(i+j)�η�
int chartoint(const char* p);//�ַ�תΪ������
int getRank(double *Array,int row,int col);//ȫ�ֺ���������
int IsZero(double x){ return fabs(x) <= 0.00001?1:0;} //ȫ�ֺ������ж��Ƿ�Ϊ��
int jacobieigen(double *a,double *u,int jie);//a ��u ����jie*jie�ķ���
void getEigen(double *Array,int dim);//ȫ�ֺ������������ֵ������������Array��dim*dim����
void show(double array[],int jie); //ȫ�ֺ�������ʾ����
void transpose(double array[],int jie); //ȫ�ֺ�����ת��
bool IsSymmetricMatrix(double *Array, int dim);//�ж��ǶԳƾ���
bool notfindchar(const char* p,int index);
bool Is0Array( double *Array,int row,int col); //ȫ�ֺ������ж��Ƿ���0����
bool chardouble(const char *p);//ȫ�ֺ������ж��ַ��Ƿ�����ת��Ҫ��
double det(double array[],int Jie); //ȫ�ֺ�����������ʽ
double chartodouble(const char* p);//ȫ�ֺ������ַ���תdouble��
double* doInverse(double Array[],int row);//ȫ�ֺ���������

int cifang(int i,int j) //ȫ�ֺ�������-1��(i+j)�η�
{ 
	if( i < 0 || j < 0 )
	{
		cerr<<"i,j�Ƿ�!" << endl;
		exit(-1);
	}
	int temp = i+j;
	if( temp%2 == 0)
		return 1;
	else
		return -1;
	}
int chartoint(const char* p)//�ַ�תΪ������
{
	if( p == NULL)
		return 0;
	int length = strlen(p);
	int INT = 0,temp =0;
	if( p[0] != '-' )
	{
		for( int i = 0; i < length; i++)
		{
			temp = p[i] - '0';
			INT = INT*10 + temp;
		}
		return INT;
	}
	else
	{
		for( int j = 1; j < length; j++)
		{
			temp = p[j] - '0';
			INT = INT*10 + temp;
		}
		return (-INT);
	}
}
int getRank(double *Array,int row,int col) //ȫ�ֺ���������
{
	int i,j,k,l,i1,j1,main_row,main_col,rank;
	double main_element,temp;
	for( i = 0; i < row; i++ ) //��ѭ��
	{ 
		main_element = *(Array + i*col + i );//�Խ�Ԫ��
		main_row = i;//��Ԫ���꣺���Խ�Ԫ(i,i)
		main_col=i;
		for( j = i; j < row; j++ ) //ѭ����Ѱ�ҵ�ǰ(i,i)Ԫ֮������Ԫ������Ϊ����(main_row,mai_col)
			for( k = i; k < col; k++ )
			{ 
				if( fabs(*( Array + j*col + k )) >= fabs(main_element) && (i != j) )
				{ 
					main_element = *(Array + j*col + k);
					main_row=j;
					main_col=k;
				}
			}
		for(l = 0; l < col; l++ )
		{ 
			temp = *(Array + main_row*col + l);
			*(Array + main_row*col + l) = *(Array +i*col + l);
			*(Array +i*col + l) = temp;
		}
		for(l = 0; l < row; l++ )
		{ 
			temp = *(Array + l*col + main_col);
			*(Array + l*col + main_col) = *(Array + l*col + i);
			*(Array + l*col + i) = temp;
		}
		if( IsZero(main_element) == 1)
			break;
		for( i1 = i+1; i1 < row; i1++ )
		{ 
			temp = *(Array + i1*col + i);
			for( j1 = 0; j1 < col; j1++ )
				*(Array + i1*col + j1) = *(Array + i1*col +j1) - *(Array + i*col + j1)*temp/main_element;
		}
	}
	rank = 0;
	for( i = 0; i < row; i++ )
	{ 
		for( j = 0; j < col; j++ )
		{
			if( IsZero( *(Array + i*col + j)) == 0 )
			{
				rank++;
				break;
			}
			else
				continue;
		}
	}
	return rank;
}

int jacobieigen(double *a,double *u,int jie)//a ��u ����jie*jie�ķ���
{ 
	int i,j,k,p,q;
	double d,m,x,y,sn,cn,w;
	for( i = 0; i < jie; i++ ) //������λ����
	{ 
		(*(u+i*jie+i)) = 1;
		for( j = 0; j < jie; j++ )
		{ 
			if( i!=j )
				(*(u+i*jie+j)) = 0;
		}
	}
	k = 1;
	while(1)
	{
		m = 0;
		for( i = 0; i <= jie-1; i++ ) //ѡȡ����ֵ���ĶԽ���Ԫ��
		{ 
			for( j = 0; j <= i-1; j++ )
			{ 
				d = fabs((*(a+i*jie+j)));
				if( (i!=j) && (d>m) )
				{ 
					m = d;
					p = i;
					q = j;
				}
			}
		}
		if( m < eps )  //���㾫��Ҫ����������
			return(1);
		if( k > KMAX )  //������������������
			return(-1);
		k = k + 1;
		x = -(*(a+p*jie+q));
		y = ( (*(a+q*jie+q)) - (*(a+p*jie+p)) )/2.0;
		w = x/sqrt( x*x + y*y );
		if( y < 0 )
			w = -w;
		sn = 1 + sqrt( 1 - w*w );
		sn = w/sqrt( 2*sn );
		cn = sqrt( 1 - sn*sn );
		m = (*(a+p*jie+p));  //�������A����Ԫ��
		(*(a+p*jie+p)) = m*cn*cn + (*(a+q*jie+q))*sn*sn + (*(a+p*jie+q))*w;
		(*(a+q*jie+q)) = m*sn*sn + (*(a+q*jie+q))*cn*cn - (*(a+p*jie+q))*w;
		(*(a+p*jie+q)) = 0;
		(*(a+q*jie+p)) = 0;
		for( j = 0; j < jie; j++ )
		{ 
			if( (j!=p) && (j!=q) )
			{ 
				m = (*(a+p*jie+j));
				(*(a+p*jie+j)) = m*cn + (*(a+q*jie+j))*sn;
				(*(a+q*jie+j)) = -m*sn + (*(a+q*jie+j))*cn;
			}
		}
		for( i = 0; i < jie; i++ )
		{
			if( (i!=p)&&(i!=q) )
			{ 
				m = (*(a+i*jie+p));
				(*(a+i*jie+p)) = m*cn + (*(a+i*jie+q))*sn;
				(*(a+i*jie+q)) = -m*sn + (*(a+i*jie+q))*cn;
			}
		}
		for( i = 0; i < jie; i++ )
		{ 
			m = (*(u+i*jie+p));
			(*(u+i*jie+p)) = m*cn + (*(u+i*jie+q))*sn;
			(*(u+i*jie+q)) = -m*sn + (*(u+i*jie+q))*cn;
		}
	}
}

void getEigen(double *Array,int dim)//ȫ�ֺ������������ֵ������������Array��dim*dim����
{ 
	double *a = new double[dim*dim],*b = new double[dim*dim],*v = new double[dim*dim];
	int i,j,k;
	for( i = 0; i < dim; i++ )
		for( j = 0; j < dim; j++ )
			*(a+i*dim+j) = *(Array+i*dim+j);
	if( !IsSymmetricMatrix(a,dim) )
	{
		cerr << "��������ֵ����������\n";
		return; 
	}
	k = jacobieigen(a,v,dim);
	if( k == 1 )
	{ 
		cout << "\n����ֵ������\n";
		for( i = 0; i < dim; i++ )
		{ 
			cout << i+1 << ":";
			for( j = 0; j < dim; j++ )
				if( i == j )
					printf("%11f",(*(a+i*dim+j)));
			cout << endl;
		}
		cout << endl;
		for( i = 0; i < dim; i++ )
			for( j = 0; j < dim; j++ )
				*(b+i*dim+j) = *(v+i*dim+j);
		cout << "��Ӧ����������\n";
		for( i = 0; i < dim; i++ )
		{ 
			printf(" %d:",i+1);
			for( j = 0; j < dim; j++ )
				printf(" %11f",*(b+i*dim+j));
			cout << endl;
		}
	}
	else
		cout << "dimoEigenValue\n";
}

void transpose(double array[],int jie) //ȫ�ֺ�����ת��
{ 
	if( jie <= 0){cerr<< "����Ƿ�!" << endl;exit(-1);}
	double temp;
	for( int i = 0; i < jie; i++ )
		for( int j = 0; j < i; j++ )
		{ 
			temp = *(array + i*jie + j);
			*(array + i*jie + j) = *(array + j*jie + i );
			*(array + j*jie + i) = temp;
		}
}

void show(double array[],int jie) //ȫ�ֺ������������
{ 
	if(jie <= 0)
	{ 
		cerr << " �Ƿ�" << endl;
		exit(-1);
	}
	cout << "�������Ϊ: " << endl;
	int i, j;
	for( i = 0; i < jie; i++ )
	{ 
		for( j = 0; j < jie; j++ )
			cout << *(array + i*jie + j) << '\t';
		cout << endl;
	}
}
bool IsSymmetricMatrix(double *Array, int dim)//�ж��ǶԳƾ���
{ 
	if( dim <= 0 )
	{
		cerr << "�����Ƿ�,�˳�!" << endl;
		exit(-1);
	}
	int i,j;
	for( i = 0; i < dim; i++ )
		for( j = 0; j <= i; j++ )
		{ 
			if( *(Array+i*dim+j) != *(Array+j*dim+i) )
			{ 
				cout << "���ǶԳƾ���" << endl;
				return false;
			}
		}
	cout << "�ǶԳƾ���!" << endl;
	return true;
}

bool Is0Array( double *Array,int row,int col)//ȫ�ֺ������ж��Ƿ���0����
{ 
	if( row <= 0 || col<= 0 )
	{
		cerr <<"����ֵ�Ƿ�,�˳�!\n";
		exit(-1);
	}
	int i,j;
	double temp1 = 0,temp2 = 0;
	for( i = 0; i < row; i++ )
		for( j = 0; j < col; j++ )
		{
			temp2 = fabs(*(Array + i*col +j));
			temp1 += temp2;
		}
	if(temp1 < 0.00001)
	{
		cout << "��" << row << '*' << col << "��0����\n";
		return true;
	}
	return false;
}

bool notfindchar(const char* p,int index)
{
	if( p == NULL )
		return true;
	int length = strlen(p);
	for( int i = index; i <= length-1; i++ )
		if(p[i]<'0' || p[i] >'9')
			return false;
	return true;
}

double det(double array[],int Jie) //ȫ�ֺ�����������ʽ
{ 
	if( Jie <= 0 )
	{
		cerr << "��С��0�����0!" << endl;
		return 0;
	}
	else if( Jie == 1)
		return array[0];
	else
	{ 
		int i,j,k,tag;
		double *subArray[500];  //����̶�ֵ?
		for( i = 0; i < Jie; i++ )
			subArray[i] = new double[(Jie-1)*(Jie-1)];
		for( i = 0; i < Jie; i++ )
			for( j = 0; j < Jie-1; j++ )
				for( k = 0; k < Jie-1; k++ )
					*(subArray[i] + j*(Jie-1) + k) = 0;
		for( i = 0; i < Jie; i++ )
			for( j = 0; j < Jie-1; j++ )
				for( k = 0; k < Jie-1; k++ )
				{ 
					if( k < i )
						*(subArray[i] + j*(Jie-1) + k) = *(array + (j+1)*Jie + k );
					else
						*(subArray[i] + j*(Jie-1) + k) = *(array + (j+1)*Jie + k+1 );     
				}
		double temp= 0;
		tag = 1;
		for( i = 0 ; i < Jie; i++)
		{ 
			temp += tag * det(subArray[i],Jie-1) * array[i];
			tag *= -1;
		}
		return temp;
	}
}

double *doInverse(double Array[],int row)//ȫ�ֺ���������  
{ 
	double d_det = det(Array,row);
	if(d_det == 0)
	{
		cerr << "����ʽΪ0,�����������!" << endl;
		exit(-1);
	} 
	int i,j,k,h,subRow = row-1;
	double *subArray[1000];      //����̶�ֵ?
	for( i = 0; i < 500; i++ )
	{ 
		if( i < row*row )
		{
			subArray[i] = new double[subRow*subRow];
			for( j = 0; j < subRow; j++ )
				for( k = 0; k < subRow; k++ )
					*(subArray[i] + j*subRow + k) = 0;
		}
		else
			subArray[i] = NULL;
	}
	for( i= 0; i < row; i++ )
		for( h = 0; h < row; h++ )
			for( j = 0; j < subRow; j++ )
				for( k = 0; k < subRow; k++ )
				{
					if( j < i && k < h )
						*(subArray[i*row+h] + j*subRow + k) = *(Array + j*row + k );
					if( j < i && k >= h )
						*(subArray[i*row+h] + j*subRow + k) = *(Array + j*row + k+1 );
					if( j >= i && k < h )
						*(subArray[i*row+h] + j*subRow + k) = *(Array + (j+1)*row + k );
					if( j >= i && k >= h)
						*(subArray[i*row+h] + j*subRow + k) = *(Array + (j+1)*row + k+1 );
				}
	double *tempArray = new double[row*row]; //������ʱ���飬��ֵ�󲢷���
	for( i = 0; i < row*row; i++ )    //��ʼ��
		tempArray[i] = 0; 
	for( i = 0; i < row; i++ )
		for( j = 0; j < row; j++ )
			*(tempArray + i*row + j) = cifang(i,j)*det(subArray[i*row+j],subRow)/d_det;//��������(i,j)Ԫ
	transpose(tempArray,row); //����ת�ú�İ���
	return tempArray;
}
bool chardouble(const char *p)
{
	if( p == NULL )
		return true;//����?
	int length=strlen(p),count1=0,count2=0,i;
	for( i = 0; i < length; i++ )
	{
		if( p[i] =='.')
			count1++;
		if( p[i] < '0' || p[i] > '9' )
			count2++;
	}
	if( count2 == 0)
		return true;
	if( count2 == 1  && p[0] == '-' && p[1]!='\0' )
		return true;
	if( count2 == 1 && count1 == 1 && p[0]!='.' )
		return true;
	if( count2 == 2 && count1 == 1 && p[0] == '-' && p[1] !='.' )
		return true;
	return false;
}

double chartodouble(const char* p)
{
	if( p == NULL )
		return 0;
	double DOUBLE = 0,DOUBLE1 = 0;
	int length=strlen(p),count1=0,count2=0,temp = 0,i,j=-1;
	for( i = 0; i < length; i++ )
	{
		if( p[i] =='.')
		{
			count1++;
			j=i;
		}
		if( p[i] < '0' || p[i] > '9' )
			count2++;
	}
	if( count2 == 0)
	{
		for( i = 0;i< length;i++)
		{
			temp = p[i]-'0';
			DOUBLE = DOUBLE*10 + temp;
		}
		return DOUBLE;
	}
	if( count2 == 1  && p[0] == '-'&& p[1]!='\0' )
	{
		for( i = 1;i < length;i++)
		{
			temp = p[i] - '0';
			DOUBLE = DOUBLE*10 + temp;
		}
		return -DOUBLE;
	}
	if( count2 == 1 && count1 == 1 && p[0]!='.' )
	{
		for( i = 0;i < j;i++)
		{
			temp = p[i] - '0';
			DOUBLE = DOUBLE*10 + temp;
		}
		for( i = length-1; i > j;i--)
		{
			temp = p[i] - '0';
			DOUBLE1 = DOUBLE1*0.1+temp;
		}
		return (DOUBLE+0.1*DOUBLE1);
	}
	if( count2 == 2 && count1 == 1 && p[0] == '-' && p[1] != '.' )
	{
		for( i = 1;i < j;i++)
		{
			temp = p[i] - '0';
			DOUBLE = DOUBLE*10 + temp;
		}
		for( i = length-1; i > j;i--)
		{
			temp = p[i] - '0';
			DOUBLE1 = DOUBLE1*0.1+temp;
		}
		return -(DOUBLE+0.1*DOUBLE1);
	}
	return 0;
}
#endif