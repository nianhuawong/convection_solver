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

		boundary_condition();

		compute_residual();

		output_residual();
	}

	output_results(outFile);

	return 0;
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
	if (iter % 10 == 0)
	{
		cout << "\titer " << "\tresidual" << endl;
	}

	cout << "\t" << iter << "\t" << residual << endl;
}

void load_qField()
{
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

void time_marching_1st_upwind()
{	
	for (int iNode = 1; iNode < numberOfGridPoints - 1; ++iNode)
	{
		qField_N1[iNode] = qField[iNode] - sigma * (qField[iNode] - qField[iNode-1]);
	}
}

void time_marching_lax_wendroff()
{
	for (int iNode = 1; iNode < numberOfGridPoints - 1; ++iNode)
	{
		qField_N1[iNode] = qField[iNode] - 0.5 * sigma * (qField[iNode+1] - qField[iNode-1]) 
			             + 0.5 * sigma * sigma *(qField[iNode+1] -2.0 * qField[iNode]+ qField[iNode-1]);
	}
}

void time_marching_beam_warming()
{
	for (int iNode = 2; iNode < numberOfGridPoints - 1; ++iNode)
	{
		qField_N1[iNode] = qField[iNode] - sigma * (qField[iNode] - qField[iNode-1]) 
                         + 0.5*sigma*(1.0-sigma)*(qField[iNode]-2.0*qField[iNode-1]+qField[iNode-2]);
	}
}

void boundary_condition()
{
	qField[0] = 0;
	qField[numberOfGridPoints - 1] = 2.0 * qField[numberOfGridPoints - 2] - qField[numberOfGridPoints - 3];

	qField_N1[0] = 0;
	qField_N1[numberOfGridPoints - 1] = 2.0 * qField[numberOfGridPoints - 2] - qField[numberOfGridPoints - 3];
}

void flow_initialization()
{
	//初值
	qField.resize(numberOfGridPoints);
	qField_N1.resize(numberOfGridPoints);

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

	qField_N1 = qField;

	boundary_condition();
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
	double endCoord = 1.0;
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
	cout << "1--1st_upwind;\t2--Lax_Wendroff;\t3--Beam_Warming, please choose!" << endl;
	int time_march_method;
	cin >> time_march_method;

	if (time_march_method == 1)
	{
		time_marching = &time_marching_1st_upwind;
		cout << "time marching method is 1st_upwind!" << endl;
		outFile = "results-1st.dat";
	}
	else if (time_march_method == 2)
	{
		time_marching = &time_marching_lax_wendroff;
		cout << "time marching method is lax_wendroff!" << endl;
		outFile = "results-LW.dat";
	}
	else if (time_march_method == 3)
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

//限制器函数
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

