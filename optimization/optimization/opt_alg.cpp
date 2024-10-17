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
		auto* p = new double[2] { 0, 0 };
		auto i = 0;
		auto x1 = x0 + d;
		matrix fx0 = ff(matrix(1, &x0), ud1, ud2);
		matrix fx1 = ff(matrix(1, &x1), ud1, ud2);
		if (m2d(fx1) == m2d(fx0)) {
			p[0] = x0;
			p[1] = x1;
			return p;
		}
		if (m2d(fx1) > m2d(fx0)) {
			d = -d;
			x1 = x0 + d;
			fx1 = ff(matrix(1, &x1), ud1, ud2);

			if (m2d(fx1) >= m2d(fx0)) {
				p[0] = x1;
				p[1] = x0 - d;
				return p;
			}
		}
		while (solution::f_calls < Nmax) {
			++i;
			double xi_next = x0 + pow(alpha, i) * d;
			matrix fxi_next = ff(matrix(1, &xi_next), ud1, ud2);
			if (m2d(fxi_next) >= m2d(fx1)) {
				if (d > 0) {
					p[0] = x0 + pow(alpha, i - 1) * d;
					p[1] = xi_next;
				}
				else {
					p[0] = xi_next;
					p[1] = x0 + pow(alpha, i - 1) * d;
				}
				return p;
			}
			x1 = xi_next;
			fx1 = fxi_next;
		}
		throw string("Metoda ekspansji osiagnela maksymalna liczbe iteracji");
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
		solution Xopt;  // Zmienna przechowująca rozwiązanie
		Xopt.ud = b - a; // Dodanie długosci przedzialu
		
		solution A(a), B(b), C; // Krok2: a(0), b(0)
		C.x = (a + b) / 2;  // Krok2: c(0) = (a + b) / 2

		solution D, D_old(a); // inicjalizacja d i kopia pomocnicza

		A.fit_fun(ff, ud1, ud2);  // sledzenie wywolan A
		B.fit_fun(ff, ud1, ud2);  // sledzenie wywolan B
		C.fit_fun(ff, ud1, ud2);  // sledzenie wywolan C

		double l, m;

		for (int i = 0; i < Nmax; i++)  // Krok1/krok3: repeat
		{
			l = m2d(A.y * (pow(B.x, 2) - pow(C.x, 2)) + B.y * (pow(C.x, 2) - pow(A.x, 2)) + C.y * (pow(A.x, 2) - pow(B.x, 2))); // Krok 4: oblicz l
			m = m2d(A.y * (B.x - C.x) + B.y * (C.x - A.x) + C.y * (A.x - B.x)); // Krok 5: oblicz m

			if (m <= 0) // Krok 6: Sprawdzanie, czy m <= 0
			{
				Xopt = D_old;  // Zwrócenie poprzedniego rozwiązania, jeśli m jest niepoprawne
				Xopt.flag = -1; // krok 7: Flaga bledu (error)
				return Xopt;   // Zakończenie algorytmu z błędem
			} // krok 8

			D.x = 0.5 * l / m; // Krok 9: Obliczanie nowego punktu d(i)
			D.fit_fun(ff, ud1, ud2);  // Osledzenie wywolan D

			if (A.x < D.x && D.x < C.x)  //krok 10: a(i) < d(i) < c(i)
			{
				if (D.y < C.y)  // krok 11: f(d(i)) < f(c(i))
				{
					B = C;      // Aktualizacja b(i+1) = c(i)
					C = D;      // Aktualizacja c(i+1) = d(i)
				}
				else
				{
					A = D;      // Aktualizacja a(i+1) = d(i)
				}
			}
			else if (C.x <= D.x && D.x <= B.x)  // Sprawdzanie, czy d(i) jest pomiędzy c(i) a b(i)
			{
				if (D.y < C.y)  // Jeśli f(d(i)) < f(c(i))
				{
					A = C;      // Aktualizacja a(i+1) = c(i)
					C = D;      // Aktualizacja c(i+1) = d(i)
				}
				else
				{
					B = D;      // Aktualizacja b(i+1) = d(i)
				}
			}
			else  // Jeśli d(i) nie mieści się w przedziale [a, b]
			{
				Xopt = D_old;  // Zwrócenie poprzedniego rozwiązania
				Xopt.flag = 2; // Flaga błędu
				return Xopt;   // Zakończenie algorytmu z błędem
			}

			// Zapisanie długości przedziału do wyniku
			Xopt.ud.add_row((B.x - A.x)());

			// Krok 39: Warunki zakończenia
			if (B.x - A.x < epsilon || abs(D.x() - D_old.x()) < gamma)
			{
				Xopt = D;      // Zwrócenie optymalnego punktu d(i)
				Xopt.flag = 0; // Flaga sukcesu
				break;         // Zakończenie algorytmu
			}

			// Krok 36: Sprawdzanie maksymalnej liczby iteracji
			if (solution::f_calls > Nmax)
			{
				Xopt = D;      // Zwrócenie optymalnego punktu d(i)
				Xopt.flag = 1; // Flaga przekroczenia liczby iteracji
				break;         // Zakończenie algorytmu
			}

			D_old = D;  // Zapisanie poprzedniego punktu d(i)
		}

		return Xopt;  // Zwrócenie optymalnego rozwiązania
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);  // Obsługa błędów
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
