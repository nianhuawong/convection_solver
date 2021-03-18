#pragma once

const double coeff_a = 1.0;
const double PI = 3.14159;

int iter;
int numberOfGridPoints = 201;
int numberOfTimeSteps  = 10000;
int inflowType = 0;

double totalTime = 4.0*PI;
double dt;
double ds;
double sigma;

double residual;
string outFile;

vector< double > qField;
vector< double > qField_N1;
vector< double > qField_M1;
vector< double > xCoordinates;

void initialize_parameter();
void flow_initialization();
void flow_initialization_inflow0();
void flow_initialization_inflow1();
void flow_initialization_inflow2();
void load_qField();

typedef void (*Time_Marching_Pointer)(); //自定义void类型的指针函数
Time_Marching_Pointer time_marching;
void set_time_march_method();

void time_marching_CTCS();
void time_marching_1st_upwind();
void time_marching_2nd_upwind();
void time_marching_lax_wendroff();
void time_marching_beam_warming();

void boundary_condition();
void boundary_condition_periodic();
void compute_residual();
void output_residual();
void output_results(string fileName);

void generate_grid_1D(int numberOfGridPoints);

double minmod_limiter(double a, double b);