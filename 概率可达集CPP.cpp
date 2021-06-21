// 概率可达集CPP.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <string>
#include "Constants.h"
#include "Functions.h"
using namespace std;


string filenames0[22] = { "statevalue0.dat","statevalue1.dat","statevalue2.dat","statevalue3.dat","statevalue4.dat","statevalue5.dat","statevalue6.dat","statevalue7.dat","statevalue8.dat","statevalue9.dat","statevalue10.dat",
"statevalue11.dat","statevalue12.dat","statevalue13.dat","statevalue14.dat","statevalue15.dat","statevalue16.dat","statevalue17.dat","statevalue18.dat","statevalue19.dat","statevalue20.dat","statevalue21.dat" };


string filenames1[22] = { "statevalue0.csv","statevalue1.csv","statevalue2.csv","statevalue3.csv","statevalue4.csv","statevalue5.csv","statevalue6.csv","statevalue7.csv","statevalue8.csv","statevalue9.csv","statevalue10.csv",
"statevalue11.csv","statevalue12.csv","statevalue13.csv","statevalue14.csv","statevalue15.csv","statevalue16.csv","statevalue17.csv","statevalue18.csv","statevalue19.csv","statevalue20.csv","statevalue21.csv" };


double Arr_Noise[10000][2];
double Arr_StateValue[Nx][Ny];
double Arr_StateValueNew[Nx][Ny];
State Arr_State[Nx][Ny];


int main()
{
	Func_Init();
	//从第21个时刻开始递归至0时刻，个时间步保存一次数据
	for (int i = Ntimestep; i >= 0; i--)
	{
		Func_RecursionMT(i);
		Func_SavaData(filenames0[i], filenames1[i]);
		cout << i << endl;
	}
}


