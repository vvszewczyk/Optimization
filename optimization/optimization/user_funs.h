#pragma once

#include"ode_solver.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);
matrix df1(matrix, matrix = NAN, matrix = NAN);
matrix df2(matrix, matrix = NAN, matrix = NAN);
double M(double alpha, double omega, double k1, double k2);
matrix Q(double t, matrix Y, matrix ud1, matrix ud2);
matrix f1(double, matrix, matrix = NAN, matrix = NAN);
matrix f2(matrix, matrix = NAN, matrix = NAN);