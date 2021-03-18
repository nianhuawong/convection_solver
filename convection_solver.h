#pragma once

const double coeff_a = 0.1;
const double PI = 3.14159;

int iter;
int numberOfGridPoints = 101;
int numberOfTimeSteps  = 40;

double totalTime = 1.0;
double dt;
double ds;
double sigma;

double residual;
string outFile;

vector< double > qField;
vector< double > qField_N1;
vector< double > xCoordinates;

void initialize_parameter();
void flow_initialization();
void load_qField();

typedef void (*Time_Marching_Pointer)(); //自定义void类型的指针函数
Time_Marching_Pointer time_marching;

void time_marching_1st_upwind();
void time_marching_lax_wendroff();
void time_marching_beam_warming();

void boundary_condition();
void boundary_condition_periodic();
void compute_residual();
void output_residual();
void output_results(string fileName);

void generate_grid_1D(int numberOfGridPoints);
void set_time_march_method();
double minmod_limiter(double a, double b);