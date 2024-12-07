#pragma once

#include"ode_solver.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);
matrix df1(matrix, matrix = NAN, matrix = NAN);
matrix ff2T(matrix, matrix = NAN, matrix = NAN);
double M(double alpha, double omega, double k1, double k2);
matrix df2(double t, matrix Y, matrix ud1, matrix ud2);
matrix ff2R(matrix x, matrix ud1, matrix ud2);
matrix ff3T(matrix, matrix = NAN, matrix = NAN);
matrix df3(double t, matrix Y, matrix = NAN, matrix = NAN);
matrix ff3R(matrix x, matrix = NAN, matrix = NAN);
matrix f1(double, matrix, matrix = NAN, matrix = NAN);
bool check_constraints_2D(matrix, matrix);
matrix f2(matrix, matrix = NAN, matrix = NAN);
matrix ff4T(matrix x, matrix, matrix);
matrix gf4T(matrix x, matrix, matrix);
matrix hf4T(matrix = NAN, matrix = NAN, matrix = NAN);
double sigmoid(matrix X, matrix theta);
matrix logisticCostFunction(matrix x, matrix ud1, matrix ud2);
matrix computeGradient(matrix X, matrix Y, matrix theta);
matrix computeHessian(matrix X, matrix Y, matrix theta);
double computeAccuracy(matrix X, matrix Y, matrix theta);