#include<iostream>
#include<vector>
#include<fstream>

using namespace std;
#include "convection_solver.h"

int main()
{
	initialize_parameter();

	flow_initialization();

	for (iter = 0; iter < numberOfTimeSteps; ++iter)
	{
		load_qField();

		time_marching();

		//boundary_condition();
		boundary_condition_periodic();

		compute_residual();

		output_residual();
		
		physicalTime += dt;
	}

	output_results(outFile);
	
	compute_exact_solution();

	return 0;
}

void compute_exact_solution()
{
	vector<double> exactSolution(numberOfGridPoints);
}

void output_results( string fileName )
{
	cout << "dumping results..." << endl;
	fstream file;
	file.open(fileName, ios_base::out );

	file << "TITLE     = \"results\"" << endl;
	file << "VARIABLES = \"x\", \"qField\"" << endl;

	for (int iNode = 0; iNode < numberOfGridPoints; ++iNode)
	{
		file << xCoordinates[iNode] << "\t" << qField_N1[iNode] << endl;
	}	

	file.close();

	cout << "done!" << endl;
}

void output_residual()
{
	int residualOutPut = 200;
	if (iter % (numberOfTimeSteps/residualOutPut) == 0)
	{
		cout << "\titer " << "\tresidual" << endl;
	}

	if (iter % residualOutPut == 0)
	{
		cout << "\t" << iter << "\t" << residual << endl;
	}
}

void load_qField()
{
	qField_M1 = qField;
	qField = qField_N1;
}

void compute_residual()
{
	for (int iNode = 0; iNode < numberOfGridPoints; ++iNode)
	{
		residual += (qField_N1[iNode] - qField[iNode]) * (qField_N1[iNode] - qField[iNode]);
	}

	residual = sqrt(residual / numberOfGridPoints);  //L2
}

void time_marching_CTCS()
{
	//BC
	qField_N1[0] = qField_M1[0] - sigma * (qField[1] - qField_LB[1]);

	//interior
	for (int iNode = 1; iNode < numberOfGridPoints - 1; ++iNode)
	{
		qField_N1[iNode] = qField_M1[iNode] - sigma * (qField[iNode + 1] - qField[iNode - 1]);
	}
}

void time_marching_1st_upwind()
{	
	qField_N1[0] = qField[0] - sigma * (qField[0] - qField_LB[1]);

	for (int iNode = 1; iNode < numberOfGridPoints - 1; ++iNode)
	{
		qField_N1[iNode] = qField[iNode] - sigma * (qField[iNode] - qField[iNode-1]);
	}
}

void time_marching_2nd_upwind()
{
	qField_N1[0] = qField[0] - 0.5 * sigma * (3.0 * qField[0] - 4.0 * qField_LB[1] + qField_LB[0]);
	qField_N1[1] = qField[1] - 0.5 * sigma * (3.0 * qField[1] - 4.0 * qField[0]    + qField_LB[1]);

	for (int iNode = 2; iNode < numberOfGridPoints - 1; ++iNode)
	{
		qField_N1[iNode] = qField[iNode] - 0.5 * sigma * (3.0*qField[iNode] - 4.0*qField[iNode - 1] + qField[iNode-2]);
	}
}

void time_marching_lax_wendroff()
{
	qField_N1[0] = qField[0] - 0.5 * sigma * (qField[1] - qField_LB[1])
		+ 0.5 * sigma * sigma * (qField[1] - 2.0 * qField[0] + qField_LB[1]);

	for (int iNode = 1; iNode < numberOfGridPoints - 1; ++iNode)
	{
		qField_N1[iNode] = qField[iNode] - 0.5 * sigma * (qField[iNode+1] - qField[iNode-1]) 
			             + 0.5 * sigma * sigma *(qField[iNode+1] -2.0 * qField[iNode]+ qField[iNode-1]);
	}
}

void time_marching_beam_warming()
{
	qField_N1[0] = qField[0] - sigma * (qField[0] - qField_LB[1])
		+ 0.5 * sigma * (1.0 - sigma) * (qField[0] - 2.0 * qField_LB[1] + qField_LB[0]);

	qField_N1[1] = qField[1] - sigma * (qField[1] - qField[0])
		+ 0.5 * sigma * (1.0 - sigma) * (qField[1] - 2.0 * qField[0] + qField_LB[1]);

	for (int iNode = 2; iNode < numberOfGridPoints - 1; ++iNode)
	{
		qField_N1[iNode] = qField[iNode] - sigma * (qField[iNode] - qField[iNode-1]) 
                         + 0.5*sigma*(1.0-sigma)*(qField[iNode]-2.0*qField[iNode-1]+qField[iNode-2]);
	}
}

void boundary_condition()
{
	qField_N1[0] = 0;		//�߽�ֵ
	qField_N1[numberOfGridPoints - 1] = 2.0 * qField_N1[numberOfGridPoints - 2] - qField_N1[numberOfGridPoints - 3];  //�ұ߽�ֵ��ͨ�����õ���

	qField_LB[0] = 2.0 * qField_N1[0] - qField_N1[1];	//���õ������ֵ
	qField_LB[1] = 2.0 * qField_N1[0] - qField_N1[1];
	qField_RB[0] = 2.0 * qField_N1[numberOfGridPoints - 1] - qField_N1[numberOfGridPoints - 2];
	qField_RB[1] = 2.0 * qField_N1[numberOfGridPoints - 1] - qField_N1[numberOfGridPoints - 2];

	//qField_N1[0] = 0;
	//qField_N1[numberOfGridPoints - 1] = 2.0 * qField_N1[numberOfGridPoints - 2] - qField_N1[numberOfGridPoints - 3];
}

void boundary_condition_periodic()
{
	qField_N1[numberOfGridPoints - 1] = 2.0 * qField[numberOfGridPoints - 2] - qField[numberOfGridPoints - 3];  //�ұ߽�ֵ��ͨ�����õ���
	qField_RB[0] = 2.0 * qField_N1[numberOfGridPoints - 1] - qField_N1[numberOfGridPoints - 2];//���õ������ֵ
	qField_RB[1] = 2.0 * qField_N1[numberOfGridPoints - 1] - qField_N1[numberOfGridPoints - 2];

	qField_N1[0] = qField_N1[numberOfGridPoints - 1];
	qField_LB[0] = qField_N1[numberOfGridPoints - 2];
	qField_LB[1] = qField_N1[numberOfGridPoints - 1];

	//qField_RB[0] = 2.0 * qField[numberOfGridPoints - 2] - qField[numberOfGridPoints - 3];
	//qField_RB[1] = 2.0 * qField[numberOfGridPoints - 2] - qField[numberOfGridPoints - 3];

	//qField_LB[0] = qField[numberOfGridPoints - 2];
	//qField_LB[1] = qField[numberOfGridPoints - 1];

	//qField[numberOfGridPoints - 1] = 2.0 * qField[numberOfGridPoints - 2] - qField[numberOfGridPoints - 3];
	//qField[0] = qField[numberOfGridPoints - 1];

	//qField_M1[numberOfGridPoints - 1] = 2.0 * qField_M1[numberOfGridPoints - 2] - qField_M1[numberOfGridPoints - 3];
	//qField_M1[0] = qField_M1[numberOfGridPoints - 1];

	//qField_N1[numberOfGridPoints - 1] = 2.0 * qField_N1[numberOfGridPoints - 2] - qField_N1[numberOfGridPoints - 3];
	//qField_N1[0] = qField_N1[numberOfGridPoints - 1];
}

void flow_initialization()
{
	qField.resize	(numberOfGridPoints);
	qField_N1.resize(numberOfGridPoints);
	qField_M1.resize(numberOfGridPoints);

	qField_LB.resize(2);  //���׸�ʽ��Ҫ2�������
	qField_RB.resize(2);

	cout << "please choose inflow type..." << endl;
	cout << "0--��Ծ������\t1--sine������\t2-�ֶκ�����" << endl;
	cin >> inflowType;

	if (inflowType==0)
	{
		flow_initialization_inflow0();
	}
	else if (inflowType == 1)
	{
		flow_initialization_inflow1();
	}
	else if (inflowType == 2)
	{
		flow_initialization_inflow2();
	}
}

void flow_initialization_inflow2()
{
	if (iter == 0)
	{
		for (int iNode = 0; iNode < numberOfGridPoints; ++iNode)
		{
			double xNode = xCoordinates[iNode];
			if (xNode >= 0.0 && xNode <= 0.2)
			{
				qField[iNode] = 0.0;
			}
			else if(xNode > 0.2 && xNode <= 0.5)
			{
				qField[iNode] = sin((xNode - 0.2) * 10.0 * PI);
			}
			else if (xNode > 0.5 && xNode <= 0.7)
			{
				qField[iNode] = 7.5 * (xNode - 0.5);
			}
			else if (xNode > 0.7 && xNode <= 1.0)
			{
				qField[iNode] = -1.0;
			}
		}
	}

	qField_M1 = qField;
	qField_N1 = qField;

	output_results("results-accurate.dat");

	boundary_condition();
}

void flow_initialization_inflow0()
{
	if (iter == 0)
	{
		for (int iNode = 0; iNode < numberOfGridPoints; ++iNode)
		{
			double xNode = xCoordinates[iNode];
			if (xNode >= 0.25 && xNode <= 0.75)
			{
				qField[iNode] = 1.0;
			}
			else
			{
				qField[iNode] = 0.0;
			}
		}
	}

	qField_M1 = qField;
	qField_N1 = qField;

	output_results("results-accurate.dat");

	boundary_condition_periodic();
}

void flow_initialization_inflow1()
{
	if (iter == 0)
	{
		for (int iNode = 0; iNode < numberOfGridPoints; ++iNode)
		{
			double xNode = xCoordinates[iNode];
			if (xNode >= 0 && xNode <= 2.0*PI)
			{
				qField[iNode] = sin(xNode);
			}
			else
			{
				qField[iNode] = 0.0;
			}
		}
	}

	qField_M1 = qField;
	qField_N1 = qField;

	output_results("results-accurate.dat");

	boundary_condition_periodic();
}

void initialize_parameter()
{		
	//cout << "Enter number of grid points..." << endl;
	//cin >> numberOfGridPoints;
	cout << "numberOfGridPoints = " << numberOfGridPoints << endl;

	//cout << "Enter totalTime..." << endl;
	//cin >> totalTime;
	cout << "totalTime = " << totalTime << endl;

	set_time_march_method();

	generate_grid_1D( numberOfGridPoints );

	int iter_min = int( totalTime * coeff_a / ds );

	cout << "Enter number of time steps..." << "iter > " << iter_min << endl;
	//cin >> numberOfTimeSteps;
	cout << "numberOfTimeSteps = " << numberOfTimeSteps << endl;
	
	dt = totalTime / numberOfTimeSteps;
	sigma = coeff_a * dt / ds;
	cout << "sigma = " << sigma << endl;
}

void generate_grid_1D( int numberOfGridPoints )
{
	double startCoord = 0.0; 
	//double endCoord = 2.0 * PI;
	ds = ( endCoord - startCoord ) / ( numberOfGridPoints - 1 );

	xCoordinates.resize(numberOfGridPoints);
	for (int iNode = 0; iNode < numberOfGridPoints; ++iNode)
	{
		xCoordinates[iNode] = startCoord + ds * iNode;
	}
	//return xCoordinates;
}

void set_time_march_method()
{
	cout << "1--CTCS;\t2--1st_upwind;\t3--2nd_upwind;\t4--Lax_Wendroff;\t5--Beam_Warming, please choose!" << endl;
	int time_march_method;
	cin >> time_march_method;
	if (time_march_method == 1)
	{
		time_marching = &time_marching_CTCS;
		cout << "time marching method is CTCS!" << endl;
		outFile = "results-CTCS.dat";		
	}
	else if (time_march_method == 2)
	{
		time_marching = &time_marching_1st_upwind;
		cout << "time marching method is 1st_upwind!" << endl;
		outFile = "results-1st.dat";
	}
	else if (time_march_method == 3)
	{
		time_marching = &time_marching_2nd_upwind;
		cout << "time marching method is 2nd_upwind!" << endl;
		outFile = "results-2nd.dat";
	}
	else if (time_march_method == 4)
	{
		time_marching = &time_marching_lax_wendroff;
		cout << "time marching method is lax_wendroff!" << endl;
		outFile = "results-LW.dat";
	}
	else if (time_march_method == 5)
	{
		time_marching = &time_marching_beam_warming;
		cout << "time marching method is beam_warming!" << endl;
		outFile = "results-BW.dat";
	}
	else
	{
		cout << "invalid time marching method, program ends!" << endl;
		exit(1);
	}
}

//����������
double minmod_limiter(double a, double b)
{
	if (a * b <= 0)
		return 0;
	else
	{
		if ((fabs(a) - fabs(b)) > 0)
			return b;
		else
			return a;
	}
}

