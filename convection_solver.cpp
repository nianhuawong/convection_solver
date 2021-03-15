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

		//time_marching_1st_upwind();
		//time_marching_lax_wendroff();
		time_marching_beam_warming();   //格式还有bug

		boundary_condition();

		compute_residual();

		output_residual();
	}

	output_results();

	return 0;
}

void output_results()
{
	cout << "dumping results..." << endl;
	fstream file;
	file.open("results.dat", ios_base::out );

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
	//qField[0] = 0;
	//qField[numberOfGridPoints - 1] = -1.0;

	//qField_N1[0] = 0;
	//qField_N1[numberOfGridPoints - 1] = -1.0;
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
	//cout << "numberOfGridPoints = " << numberOfGridPoints << endl;

	//cout << "Enter totalTime..." << endl;
	//cin >> totalTime;
	//cout << "totalTime = " << totalTime << endl;
	numberOfGridPoints = 101;
	totalTime = 1.0;

	generate_grid_1D( numberOfGridPoints );

	int iter_min = int( totalTime * coeff_a / ds );

	cout << "Enter number of time steps..." << "iter > " << iter_min << endl;
	//cin >> numberOfTimeSteps;
	numberOfTimeSteps = 40;
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