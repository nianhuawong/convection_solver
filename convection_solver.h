#pragma once

const double coeff_a = 0.2;
const double PI = 3.14159;

int iter;
int numberOfGridPoints;
int numberOfTimeSteps;

double totalTime;
double dt;
double ds;
double sigma;

double residual;

vector< double > qField;
vector< double > qField_N1;
vector< double > xCoordinates;

void initialize_parameter();
void flow_initialization();
void load_qField();

void time_marching_1st_upwind();
void time_marching_lax_wendroff();
void time_marching_beam_warming();

void boundary_condition();
void compute_residual();
void output_residual();
void output_results();

void generate_grid_1D(int numberOfGridPoints);
