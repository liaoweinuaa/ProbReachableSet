#pragma once

//以下四行代码定义计算域
const double xd = -1;
const double xu = 1;
const double yu = 1;
const double yd = -1;
const double PI = 3.141592653;

//以下两行代码定义网格数量
const int Nx = 201;
const int Ny = 201;

//时间步的数量
const int Ntimestep = 21;

const double gridx = (xu - xd) / (Nx - 1);
const double gridy = (yu - yd) / (Ny - 1);

//控制输入的上下限
const double ud = -0.5;
const double uu = 0.5;
const int Nu = 21;
const double gridu = (uu - ud) / (Nu - 1);


//线程数
const int threadnum = 10;


//状态结构体s=[x,y]^T
struct State
{
	double x;
	double y;
};




























