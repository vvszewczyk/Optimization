#include"user_funs.h"
#define PI 3.14

matrix ff0T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);
	int n = get_len(Y[0]);
	double teta_max = Y[1](0, 0);
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));
	Y[0].~matrix();
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);
	double m = 1, l = 0.5, b = 0.5, g = 9.81;
	double I = m*pow(l, 2);
	dY(0) = Y(1);
	dY(1) = ((t <= ud2(1))*ud2(0) - m*g*l*sin(Y(0)) - b*Y(1)) / I;
	return dY;
}

matrix df1(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = -cos(0.1 * x(0)) * exp(-pow(0.1 * x(0) - (2 * PI), 2)) + 0.002 * pow(0.1 * x(0), 2);
	return y;
}

matrix f1(double t, matrix Y, matrix ud1, matrix ud2)
{
	const double PA = 0.5;  // Pole podstawy zbiornika A [m^2]
	const double TA = 90;   // Temperatura w zbiorniku A [°C]
	

	const double a = 0.98;  // Wspó³czynnik lepkoœci cieczy
	const double b = 0.63;  // Wspó³czynnik zwê¿enia strumienia
	const double g = 9.81;  // Przyspieszenie ziemskie [m/s^2]


	const double PB = 1.0;  // Pole podstawy zbiornika B [m^2]
	const double DB = 0.00365665;  // Pole przekroju otworu w zbiorniku B [m^2]
	const double Fin = 0.01;  // Szybkoœæ nap³ywu wody do zbiornika B [m^3/s]
	const double Tin = 20;    // Temperatura nap³ywaj¹cej wody [°C]


	matrix dY(3, 1);

	const double FAout = (Y(0) > 0) ? (a * b * m2d(ud2) * sqrt(2 * g * Y(0) / PA)) : 0;
	const double FBout = (Y(1) > 1) ? (a * b * DB * sqrt(2 * g * Y(1) / PB)) : 0;

	dY(0) = -FAout;
	dY(1) = FAout + Fin - FBout;
	dY(2) = (Fin / Y(1)) * (Tin - Y(2)) + (FAout / Y(1)) * (TA - Y(2));

	return dY;
}

matrix f2(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(3, new double[3] {5, 1, 20});
	matrix* Y = solve_ode(f1, 0, 1, 2000, Y0, ud1, x);
	int n = get_len(Y[0]);
	double max = Y[1](0, 2);
	for (size_t i = 0; i < n; i++)
	{
		if (max < Y[1](i, 2))
		{
			max = Y[1](i, 2);
		}
		y = abs(max - 50);
	}
	
	return y;
}

matrix df2(matrix x, matrix ud2, matrix ud1)
{
	double x1 = x(0, 0);
	double x2 = x(1, 0);
	matrix result(1, 1);
	result(0, 0) = pow(x1, 2) + pow(x2, 2) - cos(2.5 * PI * x1) - cos(2.5 * PI * x2) + 2;
	return result;
}



