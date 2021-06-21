#include "Constants.h"
#include <math.h>
#include <fstream>
#include <thread>
#include "Functions.h"
#include <iostream>
#include <string>
using namespace std;

extern double Arr_Noise[10000][2];
extern double Arr_StateValue[Nx][Ny];
extern double Arr_StateValueNew[Nx][Ny];
extern State Arr_State[Nx][Ny];

void Func_Init()//��ʼ������
{
	//�������д������"noise.dat"�ļ������ļ�����10000����[0,0]Ϊ��������[[0.01,0],[0,0.01]]ΪЭ�������Ķ�ά��̬�ֲ��������
	ifstream ifs("noise.dat", ios::binary | ios::in);
	ifs.read((char*)Arr_Noise, sizeof(double) * 10000 * 2);
	ifs.close();

	//��ʼ��Arr_State����ʾ�����������״̬
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			Arr_State[i][j].x = xd + gridx * i;
			Arr_State[i][j].y = yd + gridy * j;
		}
	}
}

//״̬ת�Ʒ���
State Func_Transition(State s0, double u, int noiseindex)
{
	double x0 = s0.x, y0 = s0.y;
	double x1 = x0 + 0.1 * y0;
	double y1 = y0 + 0.1 * (-x0 + x0 * x0 * x0 + u);
	//����Arr_Noise�����е�Ԫ��
	x1 = x1 + Arr_Noise[noiseindex][0];
	y1 = y1 + Arr_Noise[noiseindex][1];

	//��ֹ״̬�뿪������
	if (x1>=xu)
	{
		x1 = xu - 0.0001;
	}

	if (x1 <= xd)
	{
		x1 = xd + 0.0001;
	}

	if (y1 >= yu)
	{
		y1 = yu - 0.0001;
	}

	if (y1 <= yd)
	{
		y1 = yd + 0.0001;
	}

	return State{ x1,y1 };
}

bool Func_IsTagetSet(State s0, int time)//�ж�״̬�Ƿ�����Ŀ�꼯
{
	double cx = 0.5 * cos(0.1 * PI * time);
	double cy = 0.5 * sin(0.1 * PI * time);

	if (s0.x<=cx+0.25 && s0.x >= cx - 0.25 && s0.y <= cy + 0.25 && s0.y >= cy - 0.25)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool Func_IsBarrierSet(State s0, int time)//�ж�״̬�Ƿ����ڽ�ֹ��
{
	double cx = 0.5 * cos(0.2 * PI * time+PI);
	double cy = 0.5 * sin(0.2 * PI * time+PI);

	if (s0.x <= cx + 0.25 && s0.x >= cx - 0.25 && s0.y <= cy + 0.25 && s0.y >= cy - 0.25)
	{
		return true;
	}
	else
	{
		return false;
	}
}

double Func_IndicatorA(State s0, int time)//Ŀ�꼯��ָʾ����
{
	if (Func_IsTagetSet(s0,time))
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

double Func_IndicatorB(State s0, int time)//��ֹ����ָʾ����
{
	if (Func_IsBarrierSet(s0, time))
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

double Func_InterPolation(State s0)//�������Բ�ֵ
{
	double ialpha_double = ((s0.x - xd) / gridx);
	double itheta_double = ((s0.y - yd) / gridy);

	int ialpha = int(ialpha_double);
	int itheta = int(itheta_double);

	double v00 = Arr_StateValue[ialpha][itheta];
	double v01 = Arr_StateValue[ialpha][itheta + 1];
	double v10 = Arr_StateValue[ialpha + 1][itheta];
	double v11 = Arr_StateValue[ialpha + 1][itheta + 1];


	//if (v000==0 && v001 == 0 && v010==0 && v011 == 0 && v100 == 0 && v101 == 0 && v110 == 0 && v111 == 0)
	//{
	//	return 0;
	//}

	double dalpha0 = ialpha_double - ialpha;
	double dalpha1 = 1 - dalpha0;
	double dtheta0 = itheta_double - itheta;
	double dtheta1 = 1 - dtheta0;


	double V00 = dalpha0 * dtheta0;
	double V01 = dalpha0 * dtheta1;
	double V10 = dalpha1 * dtheta0;
	double V11 = dalpha1 * dtheta1;


	//double VTotal = gridx * gridy * gridz;
	return (v00 * V11 + v01 * V10 + v10 * V01 + v11 * V00);
}

double Func_Expectation(State s0, double u)//���ؿ��޷���������
{
	double SUM = 0;
	for (int i = 0; i < 10000; i++)
	{
		State newstate = Func_Transition(s0, u, i);
		double value = Func_InterPolation(newstate);
		SUM += value;
	}
	return SUM / 10000;
}

double Func_OptExpectation(State s0)//�������е�u������������
{
	double maxvalue = 0;
	for (int i = 0; i < Nu; i++)
	{
		double u = ud + gridu * i;
		double value = Func_Expectation(s0, u);
		if (maxvalue<value)
		{
			maxvalue = value;
		}
	}
	return maxvalue;
}

double Func_OptValue(State s0, int time)//����min{O_B(s),max[I_A(s),max E(V_k(F))]}
{
	if (time==Ntimestep)
	{
		if (Func_IsBarrierSet(s0, time))
		{
			return 0;
		}

		if (Func_IsTagetSet(s0, time))
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}

	if (Func_IsBarrierSet(s0,time))
	{
		return 0;
	}

	if (Func_IsTagetSet(s0,time))
	{
		return 1;
	}
	return Func_OptExpectation(s0);
}

void Func_SavaData(string filename0, string filename1)//��������
{
	//���������Զ�������ʽ��������
	ofstream ofs(filename0, ios::binary | ios::out);
	ofs.write((const char*)Arr_StateValue, sizeof(double) * Nx * Ny);
	ofs.close();

	//�����ݱ���Ϊcsv�ļ�����ȡcsv�ļ�������reshape֮���������Խ�����ת��Ϊ201X201�Ķ�ά����
	ofstream ofs1;
	ofs1.open(filename1);
	for (int i = 0; i < Nx; i += 1)
	{
		for (int j = 0; j < Ny; j += 1)
		{
			ofs1 << Arr_StateValue[i][j] << endl;
		}
	}
	ofs1.close();
}

void Func_RecursionST(int threadid, int time)//�����̵߳ĵݹ飬���ܹ�10�߳�Ϊ����ÿ���̳߳�ʼindexΪ0,1,...,9������ÿ������10����index>=201*201ʱѭ������
{
	int index = threadid;
	while (true)
	{
		if (index >= Nx * Ny)
		{
			break;
		}

		if (index%1000==0)
		{
			cout << index << endl;
		}

		int ix = index / Ny;
		int iy = (index - ix * Ny);
		Arr_StateValueNew[ix][iy] = Func_OptValue(Arr_State[ix][iy], time);
		index += threadnum;
	}
}

void Func_RecursionMT(int time)//10���̵߳ĵݹ飬ÿ�εݹ��������Arr_StateValueNew������Arr_StateValue
{
	thread t0(Func_RecursionST, 0, time);
	thread t1(Func_RecursionST, 1, time);
	thread t2(Func_RecursionST, 2, time);
	thread t3(Func_RecursionST, 3, time);
	thread t4(Func_RecursionST, 4, time);
	thread t5(Func_RecursionST, 5, time);
	thread t6(Func_RecursionST, 6, time);
	thread t7(Func_RecursionST, 7, time);
	thread t8(Func_RecursionST, 8, time);
	thread t9(Func_RecursionST, 9, time);

	t0.join();
	t1.join();
	t2.join();
	t3.join();
	t4.join();
	t5.join();
	t6.join();
	t7.join();
	t8.join();
	t9.join();

	memcpy(Arr_StateValue, Arr_StateValueNew, sizeof(double) * Nx * Ny);
}