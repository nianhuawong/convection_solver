#pragma once
const double PI = 3.1415926;

int iter;
int numberOfTimeSteps  = 20000;
int inflowType = 0;

//int numberOfGridPoints = 801;
//const double coeff_a = 1.0;//方波参数
//double endCoord  = 1.0;
//double totalTime = 1.0;

//int numberOfGridPoints = 81; 
//const double coeff_a = 1;		//sine函数的参数
//double endCoord  = 2.0 * PI;
//double totalTime = 4.0 * PI;  //CTCS格式和L-W格式色散特性相同,为负色散，取时间t=32pi，其他均为4pi

int numberOfGridPoints = 801; 
const double coeff_a = 0.2;//分段函数参数
double endCoord = 1.0;
double totalTime = 1.0;

double dt;
double ds;
double sigma;
double physicalTime = 0.0;
double residual;
string outFile;

vector< double > qField;
vector< double > qField_N1;
vector< double > qField_M1;
vector< double > qField_LB;
vector< double > qField_RB;
vector< double > xCoordinates;

void initialize_parameter();
void flow_initialization();
void flow_initialization_inflow1();
void flow_initialization_inflow2();
void flow_initialization_inflow3();
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
void compute_exact_solution();