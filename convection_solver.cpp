#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
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

		if (inflowType == 3)
		{
			boundary_condition();
		}
		else
		{
			boundary_condition_periodic();
		}	

		compute_residual();

		output_residual();
		
		physicalTime += dt;
	}

	output_results(outFile, qField_N1);
	
	compute_exact_solution();

	return 0;
}

void compute_exact_solution()
{
	vector<double> exactSolution(numberOfTotalPoints);
	for (int iNode = 0; iNode < numberOfGridPoints; ++iNode)
	{
		double ds = coeff_a * physicalTime;
		double xNode = xCoordinates[iNode];
		
		int nodeIndex = iNode + numberOfGhostPoints;

		if (xNode >= 0.0 + ds && xNode <= 0.2 + ds)
		{
			exactSolution[nodeIndex] = 0.0;
		}
		else if (xNode > 0.2 + ds && xNode <= 0.5 + ds)
		{
			xNode -= ds;
			exactSolution[nodeIndex] = sin((xNode - 0.2) * 10.0 * PI);
		}
		else if (xNode > 0.5 + ds && xNode <= 0.7 + ds)
		{
			xNode -= ds;
			exactSolution[nodeIndex] = 7.5 * (xNode - 0.5);
		}
		else if (xNode > 0.7 + ds && xNode <= 1.0 + ds)
		{
			exactSolution[nodeIndex] = -1.0;
		}
	}

	output_results("results-accurate.dat", exactSolution);
}

void output_results( string fileName, vector< double >& qField_out)
{
	cout << "dumping results..." << endl;
	fstream file;
	file.open(fileName, ios_base::out );

	file << "TITLE     = \"results\"" << endl;
	file << "VARIABLES = \"x\", \"u\"" << endl;

	file << setiosflags(ios::right);
	file << setiosflags(ios::scientific);
	file << setprecision(15);

	for (int iNode = 0; iNode < numberOfGridPoints; ++iNode)
	{
		int nodeIndex = iNode + numberOfGhostPoints;
		file << xCoordinates[iNode] << "\t" << qField_out[nodeIndex] << endl;
	}	

	file.close();

	cout << "done!" << endl;
}

void output_residual()
{
	int residualOutPut = 40;
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
	qField_M1	= qField;
	qField		= qField_N1;
}

void compute_residual()
{
	for (int iNode = 0; iNode < numberOfGridPoints; ++iNode)
	{
		int nodeIndex = iNode + numberOfGhostPoints;
		residual += (qField_N1[nodeIndex] - qField[nodeIndex]) * (qField_N1[nodeIndex] - qField[nodeIndex]);
	}

	residual = sqrt(residual / numberOfGridPoints);  //L2
}

void time_marching_CTCS()
{
	//BC
	//iNode = numberOfGhostPoints��iNode = boundaryIndex

	//�����߽�
	for (int iNode = numberOfGhostPoints; iNode <= boundaryIndex; ++iNode)
	{
		qField_N1[iNode] = qField_M1[iNode] - sigma * (qField[iNode + 1] - qField[iNode - 1]);
	}
}

void time_marching_1st_upwind()
{	
	for (int iNode = numberOfGhostPoints; iNode <= boundaryIndex; ++iNode)
	{
		qField_N1[iNode] = qField[iNode] - sigma * (qField[iNode] - qField[iNode-1]);
	}
}

void time_marching_2nd_upwind()
{
	for (int iNode = numberOfGhostPoints; iNode <= boundaryIndex; ++iNode)
	{
		qField_N1[iNode] = qField[iNode] - 0.5 * sigma * (3.0 * qField[iNode] - 4.0 * qField[iNode - 1] + qField[iNode-2]);
	}
}

void time_marching_lax_wendroff()
{
	for (int iNode = numberOfGhostPoints; iNode <= boundaryIndex; ++iNode)
	{
		qField_N1[iNode] = qField[iNode] - 0.5 * sigma *		(qField[iNode+1] -		qField[iNode-1])           
										 + 0.5 * sigma * sigma *(qField[iNode+1] -2.0 * qField[iNode]+ qField[iNode-1]);
	}
}

void time_marching_lax_wendroff_TVD()
{
	double sigma1 = (sigma + abs(sigma)) / 2.0;
	double sigma2 = (sigma - abs(sigma)) / 2.0;
	double sigma3 = 0.5 * abs(sigma) * (1.0 - abs(sigma));

	for (int iNode = numberOfGhostPoints; iNode <= boundaryIndex; ++iNode)
	{
		double dum1 = qField[iNode	  ] - qField[iNode - 1] + SMALL;
		double dum2 = qField[iNode - 1] - qField[iNode - 2] + SMALL;
		double dup1 = qField[iNode + 1] - qField[iNode	  ] + SMALL;
		double dup2 = qField[iNode + 2] - qField[iNode + 1] + SMALL;

		double ita_p1, ita_m1;
		if ( coeff_a > 0 )
		{//r-
			ita_p1 = dum1 / dup1;
			ita_m1 = dum2 / dum1;
		}
		else
		{//r+
			ita_p1 = dup2 / dup1;
			ita_m1 = dup1 / dum1;
		}

		double fai_p1 = limiter_fun(ita_p1, 1.0);
		double fai_m1 = limiter_fun(ita_m1, 1.0);
		
		qField_N1[iNode] = qField[iNode] - sigma1 * dum1  - sigma2 * dup1 - sigma3 * (fai_p1 * dup1 - fai_m1 * dum1);
	}
}

void compute_rhs_weno(vector< double >& qField, vector<double>& rhs)
{
	//j�㴦��ͨ��
	vector<double> fluxVector(numberOfTotalPoints);
	for (int iNode = numberOfGhostPoints; iNode <= boundaryIndex; ++iNode)
	{
		fluxVector[iNode] = coeff_a * qField[iNode];
	}

	vector<double> q1(numberOfTotalPoints);
	vector<double> q2(numberOfTotalPoints);
	vector<double> q3(numberOfTotalPoints);
	for (int iNode = numberOfGhostPoints; iNode <= boundaryIndex; ++iNode)
	{
		q1[iNode] =  1.0 / 3.0 * fluxVector[iNode - 2] - 7.0 / 6.0 * fluxVector[iNode - 1] + 11.0 / 6.0 * fluxVector[iNode    ];
		q2[iNode] = -1.0 / 6.0 * fluxVector[iNode - 1] + 5.0 / 6.0 * fluxVector[iNode    ] + 1.0 / 3.0  * fluxVector[iNode + 1];
		q3[iNode] =  1.0 / 3.0 * fluxVector[iNode    ] + 5.0 / 6.0 * fluxVector[iNode + 1] - 1.0 / 6.0  * fluxVector[iNode + 2];
	}

	vector<double> IS1(numberOfTotalPoints);
	vector<double> IS2(numberOfTotalPoints);
	vector<double> IS3(numberOfTotalPoints);
	for (int iNode = numberOfGhostPoints; iNode <= boundaryIndex; ++iNode)
	{
		IS1[iNode] = 13.0 / 12.0 * pow( (fluxVector[iNode - 2] - 2.0 * fluxVector[iNode - 1] +		 fluxVector[iNode])		, 2) 
				    + 1.0 / 4.0  * pow( (fluxVector[iNode - 2] - 4.0 * fluxVector[iNode - 1] + 3.0 * fluxVector[iNode])		, 2);

		IS2[iNode] = 13.0 / 12.0 * pow( (fluxVector[iNode - 1] - 2.0 * fluxVector[iNode    ] +		 fluxVector[iNode + 1])	, 2) 
				    + 1.0 / 4.0  * pow( (fluxVector[iNode - 1] -       fluxVector[iNode + 1])								, 2);

		IS3[iNode] = 13.0 / 12.0 * pow( (fluxVector[iNode    ] - 2.0 * fluxVector[iNode + 1] +		 fluxVector[iNode + 2])	, 2) 
				    + 1.0 / 4.0  * pow( (fluxVector[iNode + 2] - 4.0 * fluxVector[iNode + 1] + 3.0 * fluxVector[iNode])		, 2);
	}

	double C1 = 1.0 / 10.0, C2 = 3.0 / 5.0, C3 = 3.0 / 10.0;
	double eps = 1e-6;

	//j+1/2����ͨ��
	for (int iNode = numberOfGhostPoints; iNode <= boundaryIndex; ++iNode)
	{
		double a1 = C1 / pow((eps + IS1[iNode]), 2);
		double a2 = C2 / pow((eps + IS2[iNode]), 2);
		double a3 = C3 / pow((eps + IS3[iNode]), 2);

		double w1 = a1 / (a1 + a2 + a3);
		double w2 = a2 / (a1 + a2 + a3);
		double w3 = a3 / (a1 + a2 + a3);

		fluxVector[iNode] = w1 * q1[iNode] + w2 * q2[iNode] + w3 * q3[iNode];
	}

	for (int iNode = numberOfGhostPoints; iNode <= boundaryIndex; ++iNode)
	{
		rhs[iNode] = - ( fluxVector[iNode] - fluxVector[iNode - 1]) / ds;
	}
}

void time_marching_weno_RK3()
{
	vector<double> u0(numberOfTotalPoints);
	u0 = qField;

	vector<double> rhs0(numberOfTotalPoints);
	compute_rhs_weno(u0,rhs0);

	vector<double> u1(numberOfTotalPoints);
	for (int iNode = numberOfGhostPoints; iNode <= boundaryIndex; ++iNode)
	{
		u1[iNode] = u0[iNode] + dt * rhs0[iNode];
	}

	vector<double> rhs1(numberOfTotalPoints);
	compute_rhs_weno(u1, rhs1);

	vector<double> u2(numberOfTotalPoints);
	for (int iNode = numberOfGhostPoints; iNode <= boundaryIndex; ++iNode)
	{
		u2[iNode] = 3.0 / 4.0 * u0[iNode] + 1.0 / 4.0 * u1[iNode] + 1.0 / 4.0 * dt * rhs1[iNode];
	}

	vector<double> rhs2(numberOfTotalPoints);
	compute_rhs_weno(u2, rhs2);

	vector<double> u3(numberOfTotalPoints);
	for (int iNode = numberOfGhostPoints; iNode <= boundaryIndex; ++iNode)
	{
		u3[iNode] = 1.0 / 3.0 * u0[iNode] + 2.0 / 3.0 * u2[iNode] + 2.0 / 3.0 * dt * rhs2[iNode];
	}

	qField_N1 = u3;
}

void compute_rhs_wcns(vector< double >& qField, vector<double>& rhs)
{
	vector<double> g1(numberOfTotalPoints);
	vector<double> g2(numberOfTotalPoints);
	vector<double> g3(numberOfTotalPoints);

	vector<double> s1(numberOfTotalPoints);
	vector<double> s2(numberOfTotalPoints);
	vector<double> s3(numberOfTotalPoints);
	for (int iNode = numberOfGhostPoints; iNode <= boundaryIndex; ++iNode)
	{
		g1[iNode] = (       qField[iNode - 2] - 4.0 * qField[iNode - 1] + 3.0 * qField[iNode    ]) / 2.0 / ds;
		g2[iNode] = (       qField[iNode + 1] -       qField[iNode - 1]                          ) / 2.0 / ds;
		g3[iNode] = (-3.0 * qField[iNode    ] + 4.0 * qField[iNode + 1] -       qField[iNode + 2]) / 2.0 / ds;

		s1[iNode] = (       qField[iNode - 2] - 2.0 * qField[iNode - 1] +		qField[iNode    ]) / ds / ds;
		s2[iNode] = (       qField[iNode - 1] - 2.0 * qField[iNode    ] +		qField[iNode + 1]) / ds / ds;
		s3[iNode] = (		qField[iNode    ] - 2.0 * qField[iNode + 1] +       qField[iNode + 2]) / ds / ds;
	}

	vector<double> IS1(numberOfTotalPoints);
	vector<double> IS2(numberOfTotalPoints);
	vector<double> IS3(numberOfTotalPoints);
	for (int iNode = numberOfGhostPoints; iNode <= boundaryIndex; ++iNode)
	{
		IS1[iNode] = pow((ds * g1[iNode]), 2) + pow((ds * ds * s1[iNode]), 2);
		IS2[iNode] = pow((ds * g2[iNode]), 2) + pow((ds * ds * s2[iNode]), 2);
		IS3[iNode] = pow((ds * g3[iNode]), 2) + pow((ds * ds * s3[iNode]), 2);
	}

	double C1 = 1.0 / 16.0, C2 = 10.0 / 16.0, C3 = 5.0 / 16.0;
	double eps = 1e-6;

	//j+1/2���ı�����Roeͨ��
	vector<double> fluxVector(numberOfTotalPoints);
	for (int iNode = numberOfGhostPoints; iNode <= boundaryIndex; ++iNode)
	{
		double a1 = C1 / pow((eps + IS1[iNode]), 2);
		double a2 = C2 / pow((eps + IS2[iNode]), 2);
		double a3 = C3 / pow((eps + IS3[iNode]), 2);

		double w1 = a1 / (a1 + a2 + a3);
		double w2 = a2 / (a1 + a2 + a3);
		double w3 = a3 / (a1 + a2 + a3);

		double qField_Left1  = qField[iNode] + ds * g1[iNode] / 2.0 + ds * ds * s1[iNode] / 8.0;
		double qField_Left2  = qField[iNode] + ds * g2[iNode] / 2.0 + ds * ds * s2[iNode] / 8.0;
		double qField_Left3  = qField[iNode] + ds * g3[iNode] / 2.0 + ds * ds * s3[iNode] / 8.0;

		double qField_Right1 = qField[iNode + 1] - ds * g1[iNode + 1] / 2.0 + ds * ds * s1[iNode + 1] / 8.0;
		double qField_Right2 = qField[iNode + 1] - ds * g2[iNode + 1] / 2.0 + ds * ds * s2[iNode + 1] / 8.0;
		double qField_Right3 = qField[iNode + 1] - ds * g3[iNode + 1] / 2.0 + ds * ds * s3[iNode + 1] / 8.0;

		//j+1/2��������ֵ
		double qField_Left  = w1 * qField_Left1  + w2 * qField_Left2  + w3 * qField_Left3;
		double qField_Right = w1 * qField_Right1 + w2 * qField_Right2 + w3 * qField_Right3;

		double fL = coeff_a * qField_Left ;
		double fR = coeff_a * qField_Right;
		double du = qField_Right - qField_Left;

		fluxVector[iNode] = 0.5 * (fL + fR - abs(coeff_a) * du);   //Roe��ʽ
	}

	for (int iNode = numberOfGhostPoints; iNode <= boundaryIndex; ++iNode)
	{
		rhs[iNode] = -(fluxVector[iNode] - fluxVector[iNode - 1]) / ds;
	}
}

void time_marching_wcns_RK3()
{
	vector<double> u0(numberOfTotalPoints);
	u0 = qField;

	vector<double> rhs0(numberOfTotalPoints);
	compute_rhs_wcns(u0, rhs0);

	vector<double> u1(numberOfTotalPoints);
	for (int iNode = numberOfGhostPoints; iNode <= boundaryIndex; ++iNode)
	{
		u1[iNode] = u0[iNode] + dt * rhs0[iNode];
	}

	vector<double> rhs1(numberOfTotalPoints);
	compute_rhs_wcns(u1, rhs1);

	vector<double> u2(numberOfTotalPoints);
	for (int iNode = numberOfGhostPoints; iNode <= boundaryIndex; ++iNode)
	{
		u2[iNode] = 3.0 / 4.0 * u0[iNode] + 1.0 / 4.0 * u1[iNode] + 1.0 / 4.0 * dt * rhs1[iNode];
	}

	vector<double> rhs2(numberOfTotalPoints);
	compute_rhs_wcns(u2, rhs2);

	vector<double> u3(numberOfTotalPoints);
	for (int iNode = numberOfGhostPoints; iNode <= boundaryIndex; ++iNode)
	{
		u3[iNode] = 1.0 / 3.0 * u0[iNode] + 2.0 / 3.0 * u2[iNode] + 2.0 / 3.0 * dt * rhs2[iNode];
	}

	qField_N1 = u3;
}

void time_marching_lax_wendroff_TVD_RK3()
{
	vector<double> u0(numberOfTotalPoints);
	u0 = qField;

	vector<double> rhs0(numberOfTotalPoints);
	for (int iNode = numberOfGhostPoints; iNode <= boundaryIndex; ++iNode)
	{
		double dum1 = u0[iNode	  ] - u0[iNode - 1] + SMALL;
		double dum2 = u0[iNode - 1] - u0[iNode - 2] + SMALL;
		double dup1 = u0[iNode + 1] - u0[iNode	  ] + SMALL;

		double ita_p1 = dum1 / dup1;
		double ita_m1 = dum2 / dum1;

		double fai_p1 = limiter_fun(ita_p1, 1.0);
		double fai_m1 = limiter_fun(ita_m1, 1.0);

		rhs0[iNode] = - sigma * dum1 - 0.5 * sigma * (1.0 - sigma) * (fai_p1 * dup1 - fai_m1 * dum1);
	}

	vector<double> u1(numberOfTotalPoints);
	for (int iNode = numberOfGhostPoints; iNode <= boundaryIndex; ++iNode)
	{
		u1[iNode] = u0[iNode] + rhs0[iNode];
	}

	vector<double> rhs1(numberOfTotalPoints);
	for (int iNode = numberOfGhostPoints; iNode <= boundaryIndex; ++iNode)
	{
		double dum1 = u1[iNode	  ] - u1[iNode - 1] + SMALL;
		double dum2 = u1[iNode - 1] - u1[iNode - 2] + SMALL;
		double dup1 = u1[iNode + 1] - u1[iNode	  ] + SMALL;

		double ita_p1 = dum1 / dup1;
		double ita_m1 = dum2 / dum1;

		double fai_p1 = limiter_fun(ita_p1, 1.0);
		double fai_m1 = limiter_fun(ita_m1, 1.0);

		rhs1[iNode] = - sigma * dum1 - 0.5 * sigma * (1.0 - sigma) * (fai_p1 * dup1 - fai_m1 * dum1);
	}

	vector<double> u2(numberOfTotalPoints);
	for (int iNode = numberOfGhostPoints; iNode <= boundaryIndex; ++iNode)
	{
		u2[iNode] = 3.0 / 4.0 * u0[iNode] + 1.0 / 4.0 * u1[iNode] + 1.0 / 4.0 * rhs1[iNode];
	}
	
	vector<double> rhs2(numberOfTotalPoints);
	for (int iNode = numberOfGhostPoints; iNode <= boundaryIndex; ++iNode)
	{
		double dum1 = u2[iNode	  ] - u2[iNode - 1] + SMALL;
		double dum2 = u2[iNode - 1] - u2[iNode - 2] + SMALL;
		double dup1 = u2[iNode + 1] - u2[iNode	  ] + SMALL;

		double ita_p1 = dum1 / dup1;
		double ita_m1 = dum2 / dum1;

		double fai_p1 = limiter_fun(ita_p1, 1.0);
		double fai_m1 = limiter_fun(ita_m1, 1.0);

		rhs2[iNode] = - sigma * dum1 - 0.5 * sigma * (1.0 - sigma) * (fai_p1 * dup1 - fai_m1 * dum1);
	}

	vector<double> u3(numberOfTotalPoints);
	for (int iNode = numberOfGhostPoints; iNode <= boundaryIndex; ++iNode)
	{
		u3[iNode] = 1.0 / 3.0 * u0[iNode] + 2.0 / 3.0 * u2[iNode] + 2.0 / 3.0  * rhs2[iNode];
	}

	qField_N1 = u3;
}

void time_marching_beam_warming()
{
	for (int iNode = numberOfGhostPoints; iNode <= boundaryIndex; ++iNode)
	{
		qField_N1[iNode] = qField[iNode] - sigma *  (qField[iNode] -	   qField[iNode-1]) 
					- 0.5 * sigma * (1.0 - sigma) * (qField[iNode] - 2.0 * qField[iNode-1] + qField[iNode-2]);
	}
}

void boundary_condition()
{	
	qField_N1[0] = -qField_N1[4];
	qField_N1[1] = -qField_N1[3];
	qField_N1[2] = 0.0;

	qField_N1[ghostIndex  ]  = 2.0 * qField_N1[boundaryIndex] - qField_N1[boundaryIndex - 1];
	qField_N1[ghostIndex+1]  = 2.0 * qField_N1[boundaryIndex] - qField_N1[boundaryIndex - 2];

	qField[0] = - qField[4];			
	qField[1] = - qField[3];	
	qField[2] = 0.0;

	qField[ghostIndex    ]	= 2.0 * qField[boundaryIndex] - qField[boundaryIndex - 1];
	qField[ghostIndex + 1]	= 2.0 * qField[boundaryIndex] - qField[boundaryIndex - 2];

	//qField[0] = 0.0;
	//qField[1] = 0.0;
	//qField[2] = 0.0;
	//qField[boundaryIndex] = -1.0;
	//qField[ghostIndex]	  = -1.0;
	//qField[ghostIndex+1]  = -1.0;

	//qField_N1[0] = 0.0;
	//qField_N1[1] = 0.0;
	//qField_N1[2] = 0.0;
	//qField_N1[boundaryIndex]	= -1.0;
	//qField_N1[ghostIndex]		= -1.0;
	//qField_N1[ghostIndex + 1]	= -1.0;
}

void boundary_condition_periodic()
{
	qField_N1[ghostIndex]   = 2.0 * qField_N1[boundaryIndex] - qField_N1[boundaryIndex - 1];  
	qField_N1[ghostIndex+1] = 2.0 * qField_N1[boundaryIndex] - qField_N1[boundaryIndex - 2];

	qField_N1[0] = qField_N1[ghostIndex + 1];
	qField_N1[1] = qField_N1[ghostIndex    ];
}

void flow_initialization()
{
	numberOfTotalPoints = numberOfGridPoints + 2 * numberOfGhostPoints;//���׸�ʽ��Ҫ2�������
	ghostIndex			= numberOfTotalPoints - numberOfGhostPoints;
	boundaryIndex		= ghostIndex - 1;

	qField.resize	(numberOfTotalPoints);
	qField_N1.resize(numberOfTotalPoints);   //n+1ʱ�䲽��Qֵ
	qField_M1.resize(numberOfTotalPoints);	 //n-1ʱ�䲽��Qֵ�����Ĳ�ָ�ʽ���õ�n-1ʱ�䲽

	cout << "========================================="	 << endl;
	cout << "please choose inflow type..."				 << endl;
	cout << "1--��Ծ������\t2--sine������\t3-�ֶκ�����" << endl;
	cin >> inflowType;

	if (inflowType==1)
	{
		flow_initialization_inflow1();
	}
	else if (inflowType == 2)
	{
		flow_initialization_inflow2();
	}
	else if (inflowType == 3)
	{
		flow_initialization_inflow3();
	}
}

void flow_initialization_inflow3()
{
	if (iter == 0)
	{
		for (int iNode = 0; iNode < numberOfGridPoints; ++iNode)
		{
			double xNode = xCoordinates[iNode];
			int nodeIndex = iNode + numberOfGhostPoints;

			if (xNode >= 0.0 && xNode <= 0.2)
			{
				qField[nodeIndex] = 0.0;
			}
			else if(xNode > 0.2 && xNode <= 0.5)
			{
				qField[nodeIndex] = sin((xNode - 0.2) * 10.0 * PI);
			}
			else if (xNode > 0.5 && xNode <= 0.7)
			{
				qField[nodeIndex] = 7.5 * (xNode - 0.5);
			}
			else if (xNode > 0.7 && xNode <= 1.0)
			{
				qField[nodeIndex] = -1.0;
			}
		}
	}

	qField_M1 = qField;
	qField_N1 = qField;

	//output_results("results-accurate.dat");

	boundary_condition();
}

void flow_initialization_inflow1()
{
	if (iter == 0)
	{
		for (int iNode = 0; iNode < numberOfGridPoints; ++iNode)
		{
			double xNode = xCoordinates[iNode];

			int nodeIndex = iNode + numberOfGhostPoints;

			if (xNode >= 0.25 && xNode <= 0.75)
			{
				qField[nodeIndex] = 1.0;
			}
			else
			{
				qField[nodeIndex] = 0.0;
			}
		}
	}

	qField_M1 = qField;
	qField_N1 = qField;

	//output_results("results-accurate.dat");

	boundary_condition_periodic();
}

void flow_initialization_inflow2()
{
	if (iter == 0)
	{
		for (int iNode = 0; iNode < numberOfGridPoints; ++iNode)
		{
			double xNode = xCoordinates[iNode];
			int nodeIndex = iNode + numberOfGhostPoints;
			if (xNode >= 0 && xNode <= 2.0*PI)
			{
				qField[nodeIndex] = sin(xNode);
			}
			else
			{
				qField[nodeIndex] = 0.0;
			}
		}
	}

	qField_M1 = qField;
	qField_N1 = qField;

	//output_results("results-accurate.dat");

	boundary_condition_periodic();
}

void initialize_parameter()
{		
	//cout << "Enter number of grid points..." << endl;
	//cin >> numberOfGridPoints;
	cout << "numberOfGridPoints = \t" << numberOfGridPoints << endl;

	//cout << "Enter totalTime..." << endl;
	//cin >> totalTime;
	cout << "totalTime = \t\t" << totalTime << endl;
	cout << endl;

	set_time_march_method();
	set_limiter();

	generate_grid_1D( numberOfGridPoints );

	int iter_min = int( totalTime * coeff_a / ds );

	cout << "Enter number of time steps..." << "iter > " << iter_min << endl;
	//cin >> numberOfTimeSteps;
	cout << "numberOfTimeSteps = \t" << numberOfTimeSteps << endl;
	
	dt = totalTime / numberOfTimeSteps;
	sigma = coeff_a * dt / ds;
	cout << "sigma = \t\t" << sigma << endl;
	cout << endl;
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
	cout << "========================================="			<< endl;
	cout << "1--CTCS;                  2--1st_upwind;"			<< endl;
	cout << "3--2nd_upwind;            4--Lax_Wendroff;"		<< endl;
	cout << "5--Beam_Warming;          6--lax_wendroff_TVD;  "	<< endl;
	cout << "7--lax_wendroff_TVD_RK3;  8--weno-RK3;"			<< endl;
	cout << "9-WCNS_RK3"     << endl;
	cout << "please choose!" << endl;
	
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
	else if (time_march_method == 6)
	{
		time_marching = &time_marching_lax_wendroff_TVD;
		cout << "time marching method is lax_wendroff_TVD!" << endl;
		outFile = "results-LWTVD.dat";
	}
	else if (time_march_method == 7)
	{
		time_marching = &time_marching_lax_wendroff_TVD_RK3;
		cout << "time marching method is lax_wendroff_TVD_RK3!" << endl;
		outFile = "results-LWTVDRK3.dat";
	}
	else if (time_march_method == 8)
	{
		time_marching = &time_marching_weno_RK3;
		cout << "time marching method is weno_RK3!" << endl;
		outFile = "results-wenoRK3.dat";
	}
	else if (time_march_method == 9)
	{
		time_marching = &time_marching_wcns_RK3;
		cout << "time marching method is wcns_RK3!" << endl;
		outFile = "results-wcnsRK3.dat";
	}
	else
	{
		cout << "invalid time marching method, program ends!" << endl;
		exit(1);
	}
	cout << endl;
}

void set_limiter()
{
	if ((time_march_method != 6) && (time_march_method != 7))
	{
		return;
	}
	cout << "========================================="				<< endl;
	cout << "1--minmod;\t2--vanleer;\t3--superbee; please choose!"	<< endl;
	int limiter_type;
	cin >> limiter_type;
	if (limiter_type == 1)
	{
		limiter_fun = &minmod_limiter;
		cout << "limiter is minmod!" << endl;
	}
	else if (limiter_type == 2)
	{
		limiter_fun = &vanleer_limiter;
		cout << "limiter is vanleer!" << endl;
	}
	else if (limiter_type == 3)
	{
		limiter_fun = &superbee_limiter;
		cout << "limiter is superbee!" << endl;
	}
	else
	{
		cout << "invalid limiter, program ends!" << endl;
		exit(1);
	}
	cout << endl;
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

double vanleer_limiter(double a, double)
{
	return (a + abs(a)) / (1.0 + abs(a));
}

double superbee_limiter(double a, double)
{
	double tmp1 = min(2.0 * a, 1.0);
	double tmp2 = min(a, 2.0);

	return max(0.0, max(tmp1,tmp2));
}
