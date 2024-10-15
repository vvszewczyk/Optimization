#include"opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
			Xopt.fit_fun(ff, ud1, ud2);
			if (Xopt.y < epsilon)
			{
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = 0;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2]{ 0,0 };
		//Tu wpisz kod funkcji

		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2) 
{
    try 
	{
        int k = 1;
        double fik_prev = 1, fik = 1, bufor;

		// Znalezienie najmniejszej liczby spełniającej poniższy warunek
        while (fik <= (b - a) / epsilon)
		{
            k++;
			bufor = fik;
			fik = fik + fik_prev;
			fik_prev = bufor;
        }

        // a(0) = a
		double ai = a;

		// b(0) = b
		double	bi = b;

        // c(0) = b(0) - (fi(k-1) / fi(k)) * (b(0) * a(0))
        double ci = bi - (fik_prev / fik) * (bi - ai);

		// d(0) = a(0) + b(0) - c(0)
        double di = ai + bi - ci;

        solution Xopt;

        for (int i = 0; i <= k - 3; i++) 
		{
            matrix x_c(1, 1, ci);
            matrix x_d(1, 1, di);
            double fc = m2d(ff(x_c, ud1, ud2));
            double fd = m2d(ff(x_d, ud1, ud2));

            if (fc < fd) 
			{
				// a(i+1) = a(i)
                ai = ai;
				// b(i+1) = d(i)
                bi = di;
            } 
			else 
			{
				// a(i+1) = c(i)
                ai = ci;
				// b(i+1) = b(i)
                bi = bi;
            }

            // Aktualizacja dla następnej iteracji
			bufor = fik_prev;
			fik_prev = fik - fik_prev;    // fi(k-1) -> fi(k-2)
			fik = bufor;                   // fi(k) -> fi(k-1)

			// c(i+1) = b(i+1) - (fi(k-i-2) / fi(k-i-1)) * (b(i+1) - a(i+1))
            ci = bi - (fik_prev / fik) * (bi - ai);

            // d(i+1) = a(i+1) + b(i+1) - c(i+1)
            di = ai + bi - ci;
        }

        // x* = c(i+1)
        Xopt.x = matrix(1, 1, ci);

        Xopt.fit_fun(ff, ud1, ud2);

        return Xopt;
    } 
	catch (string ex_info) 
	{
        throw ("solution fib(...):\n" + ex_info);
    }
}




solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//int i = 0;
//double ai = a, bi = b, ci = c;
//double di = 0.0, diMinus1 = 0.0;

//do {
//	double l = (ff(ai) * (pow(bi, 2) - pow(ci, 2))) + (f(bi) * (pow(ci, 2) - pow(ai, 2))) + (f(ci) * (pow(ai, 2) - pow(bi, 2));
//	double m = (ff(ai) * (bi - ci)) + (f(bi) * (ci - ai)) + (f(ci) * (ai - bi));

//	if (m <= 0) {
//		return -1; // error
//	}

//	di = 0.5 * l / m;

//	if (ai < di && di < ci) {
//		if (f(di) < f(ci)) {
//			ai = ai;
//			ci = di;
//			bi = ci;
//		}
//		else {
//			ai = di;
//			ci = ci;
//			bi = bi;
//		}
//	}
//	else if (ci < di && di < bi) {
//		if (f(di) < f(ci)) {
//			ai = ci;
//			ci = di;
//			bi = bi;
//		}
//		else {
//			ai = ai;
//			ci = ci;
//			bi = di;
//		}
//	}
//	else {
//		return -1; // error
//	}

//	i++;

//	if (fcalls > Nmax) {
//		return -1; // error
//	}
//} while ((bi - ai) < epsilon || std::abs(di - diMinus1) < gamma);

//return di;

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
