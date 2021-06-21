#pragma once

//�������д��붨�������
const double xd = -1;
const double xu = 1;
const double yu = 1;
const double yd = -1;
const double PI = 3.141592653;

//�������д��붨����������
const int Nx = 201;
const int Ny = 201;

//ʱ�䲽������
const int Ntimestep = 21;

const double gridx = (xu - xd) / (Nx - 1);
const double gridy = (yu - yd) / (Ny - 1);

//���������������
const double ud = -0.5;
const double uu = 0.5;
const int Nu = 21;
const double gridu = (uu - ud) / (Nu - 1);


//�߳���
const int threadnum = 10;


//״̬�ṹ��s=[x,y]^T
struct State
{
	double x;
	double y;
};




























