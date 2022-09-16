#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <chrono>
#include <string>
#include "newmat.h"
#define  varName(x) #x 
using namespace std;
double **Stage1, **Stage2, **Stage3, **Stage4;
double deltaS, deltaT, T;
double Srange;
int Sstep;
int Tstep;
double Strike, S0, Trigger, sigma, rd, rf;
double Muti;

double max(double a, double b) {
	return (b < a) ? a : b;
}
void FiniteDifference(double**StageN)
{
	Matrix K(Sstep + 1, Sstep + 1);
	ColumnVector TempVector(Sstep + 1);
	Matrix K_1(Sstep + 1, Sstep + 1);
	ColumnVector KVector(Sstep + 1);
	ColumnVector K_1Vector(Sstep + 1);
	for (int i = 1; i <= Sstep + 1; i++)
	{
		for (int j = 1; j <= Sstep + 1; j++)
		{
			K(i, j) = 0;
			K_1(i, j) = 0;
		}
	}
	for (int i = 1; i <= Sstep + 1; i++)
	{
		KVector(i) = 0;
		K_1Vector(i) = 0;
		TempVector(i) = 0;
	}

	K_1(1, 1) = 1;
	K(1, 1) = 1;
	K_1(Sstep + 1, Sstep + 1) = 1;
	K(Sstep + 1, Sstep + 1) = 1;

	for (int i = 2; i <= Sstep; i++)
	{
		K(i, i - 1) = (-0.25*(rd - rf)* (i - 1) + 0.25* sigma*sigma * (i - 1) *(i - 1))*deltaT;
		K(i, i) = 1 - 0.5*(rd + sigma * sigma*(i - 1)*(i - 1))*deltaT;
		K(i, i + 1) = (0.25*(rd - rf) *(i - 1) + 0.25* sigma*sigma *(i - 1) *(i - 1))*deltaT;
		K_1(i, i - 1) = -(-0.25*(rd - rf)* (i - 1) + 0.25* sigma*sigma * (i - 1) *(i - 1))*deltaT;
		K_1(i, i) = 1 + 0.5*(rd + sigma * sigma*(i - 1)*(i - 1))*deltaT;
		K_1(i, i + 1) = -(0.25*(rd - rf)*(i - 1) + 0.25* sigma*sigma * (i - 1) *(i - 1))*deltaT;
	}

	for (int j = Tstep - 1; j >= 0; j--)
	{
		for (int i = 1; i <= Sstep + 1; i++)
		{
			KVector(i) = StageN[i - 1][j + 1];
			K_1Vector(i) = 0;
			TempVector(i) = 0;
		}
		TempVector = K * KVector;
		K_1Vector = K_1.i() * TempVector;
		for (int i = 1; i <= Sstep + 1; i++)
		{
			StageN[i - 1][j] = K_1Vector(i);
		}
	}
}
void Final(double**Final)
{
	for (int i = 0; i <= Sstep; i++)
	{
		for (int j = 0; j <= Tstep; j++)
		{
			Final[i][j] = -9999;
		}
	}//Set initial condition
	for (int j = 0; j <= Tstep; j++)
	{
		Final[0][j] = 0;
		Final[Sstep][j] = 0;
	}
	for (int i = 0; i <= Sstep; i++)
	{
		if (i > (Trigger * Muti))
		{
			Final[i][Tstep] = 0;
		}
		else
		{
			Final[i][Tstep] = 31.375;
		}
	}
}
void Midway1(double**Former, double**Later)
{
	for (int i = 0; i <= Sstep; i++)
	{
		for (int j = 0; j <= Tstep; j++)
		{
			Later[i][j] = -9999;
		}
	}
	//Set initial condition
	for (int j = 0; j <= Tstep; j++)
	{
		Later[0][j] = 0;
		Later[Sstep][j] = 0;
	}
	for (int i = 0; i <= Sstep; i++)
	{
		if (i > (Trigger * Muti))
		{
			Later[i][Tstep] = Former[i][0];
		}
		else
		{
			Later[i][Tstep] = 0;
		}
	}
}
void Midway2(double**Former, double**Later)
{
	for (int i = 0; i <= Sstep; i++)
	{
		for (int j = 0; j <= Tstep; j++)
		{
			Later[i][j] = -9999;
		}
	}
	//Set initial condition
	for (int j = 0; j <= Tstep; j++)
	{
		Later[0][j] = 0;
		Later[Sstep][j] = 0;
	}
	for (int i = 0; i <= Sstep; i++)
	{
		if (i > (Trigger * Muti))
		{
			Later[i][Tstep] = 0;
		}
		else if (i > (Strike*Muti))
		{
			Later[i][Tstep] = 31.375+Former[i][0];
		}
		else
		{
			Later[i][Tstep] = 0;
		}
	}
}
void Output(double**Matrix, string name)
{
	ofstream oFile;
	oFile.open(name, ios::out | ios::ate);
	oFile << ",";
	for (int i = 0; i <= Sstep; i++)
	{
		oFile << i << ",";
	}
	oFile << endl;
	for (int j = 0; j <= Tstep; j++)
	{
		oFile << j << ",";
		for (int i = 0; i <= Sstep; i++)
		{
			oFile << Matrix[i][j] << ",";
		}
		oFile << endl;
	}
	oFile << endl;
	oFile.close();
}

int main(int argc, char* argv[])
{
	Muti = 10;
	Sstep = 20 * Muti;
	Tstep = 100;
	Strike = 5.7458;
	Trigger = 1.55*Strike;
	T = 0.25;
	Srange = 20;
	deltaS = Srange / Sstep;
	deltaT = T / Tstep;
	rd = 0.02242997;
	rf = 0.19590;
	sigma = 0.2763623;
	S0 = 5.7627;
	Stage4 = new double*[Sstep + 1];
	Stage3 = new double*[Sstep + 1];
	Stage2 = new double*[Sstep + 1];
	Stage1 = new double*[Sstep + 1];
	for (int i = 0; i <= Sstep; i++)
	{
		Stage4[i] = new double[Tstep + 1];
		Stage3[i] = new double[Tstep + 1];
		Stage2[i] = new double[Tstep + 1];
		Stage1[i] = new double[Tstep + 1];
	}

	Final(Stage4);
	FiniteDifference(Stage4);
	Midway1(Stage4, Stage3);
	FiniteDifference(Stage3);
	Midway2(Stage3, Stage2);
	FiniteDifference(Stage2);
	Midway2(Stage2, Stage1);
	FiniteDifference(Stage1);

	string csv = ".csv";
	string Filename4 = varName(Stage4);
	string name4 = Filename4 + csv;
	string Filename3 = varName(Stage3);
	string name3 = Filename3 + csv;
	string Filename2 = varName(Stage2);
	string name2 = Filename2 + csv;
	string Filename1 = varName(Stage1);
	string name1 = Filename1 + csv;
	Output(Stage4, name4);
	Output(Stage3, name3);
	Output(Stage2, name2);
	Output(Stage1, name1);
	cout << Stage1[(int)(S0*Muti)][0] << endl;
	system("pause");
}
