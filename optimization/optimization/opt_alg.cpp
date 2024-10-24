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
		solution X0(x0), X1(x0 + d);
		X0.fit_fun(ff, ud1, ud2);
		X1.fit_fun(ff, ud1, ud2);
		if (X1.y == X0.y) {
			p[0] = m2d(X0.x);
			p[1] = m2d(X1.x);
			return p;
		}
		else if (X1.y > X0.y) {
			d = -d;
			X1.x = x0 + d;
			X1.fit_fun(ff, ud1, ud2);
			if (X1.y >= X0.y) {
				p[0] = m2d(X0.x - d);
				p[1] = m2d(X1.x);
				return p;
			}
		}
		solution X2;
		while (true) {
			++i;
			X2.x = x0 + pow(alpha, i) * d;
			X2.fit_fun(ff, ud1, ud2);
			if(X2.y >= X1.y || solution::f_calls > Nmax)
				break;
			X0 = X1;
			X1 = X2;
		}
		if (d > 0) {
			p[0] = m2d(X0.x);
			p[1] = m2d(X2.x);
			return p;
		}
		p[0] = m2d(X2.x);
		p[1] = m2d(X0.x);
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
		double bi = b;

		// c(0) = b(0) - (fi(k-1) / fi(k)) * (b(0) - a(0))
		double ci = bi - (fik_prev / fik) * (bi - ai);

		// d(0) = a(0) + b(0) - c(0)
		double di = ai + bi - ci;

		solution Xopt;

		std::cout << "i = 0: Range = " << bi - ai << std::endl;

		for (int i = 0; i <= k - 3; i++)
		{
			matrix x_c(1, 1, ci);
			matrix x_d(1, 1, di);

			double fc = m2d(ff(x_c, ud1, ud2));
			double fd = m2d(ff(x_d, ud1, ud2));

			solution::f_calls += 2;

			if (fc < fd)
			{
				ai = ai;
				bi = di;
			}
			else
			{
				ai = ci;
				bi = bi;
			}

			// Aktualizacja dla następnej iteracji
			bufor = fik_prev;
			fik_prev = fik - fik_prev;
			fik = bufor;

			ci = bi - (fik_prev / fik) * (bi - ai);
			di = ai + bi - ci;

			std::cout << "i = " << i + 1 << ": Range = " << bi - ai << std::endl;
		}

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

        // Krok 2: inicjalizacja a(0), b(0)
        solution A(a), B(b), C; 
        C.x = (a + b) / 2;  // c(0) = (a + b) / 2

        solution D, D_old(a); // d(i) i kopia pomocnicza

        A.fit_fun(ff, ud1, ud2);  // Wywołanie funkcji dla A
        B.fit_fun(ff, ud1, ud2);  // Wywołanie funkcji dla B
        C.fit_fun(ff, ud1, ud2);  // Wywołanie funkcji dla C

        double l, m;
        int i = 0;  // Inicjalizacja licznika iteracji
		std::cout << "i = " << i << ": Range = " << m2d(B.x) - m2d(A.x) << std::endl;

        // Pętla główna optymalizacji Lagrange'a
        while (i < Nmax) // repeat
        {
            // Krok 4: Obliczanie licznika l
            l = m2d(A.y * (pow(m2d(B.x), 2) - pow(m2d(C.x), 2)) + B.y * (pow(m2d(C.x), 2) - pow(m2d(A.x), 2)) + C.y * (pow(m2d(A.x), 2) - pow(m2d(B.x), 2)));

            // Krok 5: Obliczanie mianownika m
            m = m2d(A.y * (m2d(B.x) - m2d(C.x)) + B.y * (m2d(C.x) - m2d(A.x)) + C.y * (m2d(A.x) - m2d(B.x)));

            // Sprawdzenie, czy m <= 0
			if (m <= 0)
            {
                Xopt = D_old;  // Zwrócenie poprzedniego rozwiązania
                Xopt.flag = -1; // Flaga błędu (error)
                return Xopt;    // Zakończenie algorytmu
            }

            // Krok 9: Obliczanie nowego punktu d(i)
            D.x = 0.5 * l / m;
            D.fit_fun(ff, ud1, ud2);  // Wywołanie funkcji dla D

			// Krok 10-19: Aktualizacja przedziału w zależności od pozycji d(i)
			if (A.x < D.x && D.x < C.x)  // Krok 10: a(i) < d(i) < c(i)
			{
				if (D.y < C.y)  // Krok 11: f(d(i)) < f(c(i))
				{
					// A = A    // krok 12: a(i+1) = a(i)
					B = C;      // krok 14: b(i+1) = c(i)
					C = D;      // krok 13: c(i+1) = d(i) odwrotnie ponieważ nadpisujemy C
				}
				else  // Krok 15
				{
					A = D;      // krok 16: a(i+1) = d(i)
					// C = C    // krok 17: c(i+1) = c(i)
					// D = D    // krok 18: d(i+1) = d(i)
				}
			}// krok 19 & 20
			else if (C.x < D.x && D.x < B.x)  // Krok 21: c(i) < d(i) < b(i)
			{
				if (D.y < C.y)  // Krok 22: f(d(i)) < f(c(i))
				{
					A = C;      // krok 23: a(i+1) = c(i)
					C = D;      // krok 24: c(i+1) = d(i)
					//B = B     // krok 25: b(i+1) = b(i)
				}
				else  // Krok 26
				{
					// A = A;   // krok 27: a(i+1) = a(i)
					// C = C;   // krok 28: c(i+1) = c(i)
					B = D;      // krok 29: b(i+1) = d(i)
				}
			} // krok 30
			else  // Krok 31: Jeśli d(i) nie mieści się w przedziale [a, b]
			{
				Xopt = D_old;  // Zwrócenie poprzedniego rozwiązania
				Xopt.flag = -1; // Krok 32: Flaga błędu (error)
				if (m2d(B.x) - m2d(A.x) < 1) Xopt.flag = 0;
				return Xopt;    // Zakończenie algorytmu
			} // Krok 33 & 24

            // Wypisywanie długości przedziału
            std::cout << "i = " << i + 1 << ": Range = " << m2d(B.x) - m2d(A.x) << std::endl;

            // Warunek zakończenia
            if (m2d(B.x) - m2d(A.x) < epsilon || abs(m2d(D.x) - m2d(D_old.x)) < gamma)
            {
                Xopt = D;      // Zwrócenie optymalnego punktu d(i)
                Xopt.flag = 0; // Flaga sukcesu
                break;         // Zakończenie algorytmu
            }

            // Sprawdzenie, czy liczba wywołań funkcji celu przekroczyła Nmax
            if (solution::f_calls > Nmax)
            {
                Xopt = D;      // Zwrócenie optymalnego punktu d(i)
                Xopt.flag = 1; // Flaga przekroczenia liczby iteracji
                break;         // Zakończenie algorytmu
            }

            // Aktualizacja starego punktu d(i)
            D_old = D;

            // Inkrementacja iteratora
            i++;
        }

        Xopt.flag = 0; // Flaga sukcesu
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
		int n = get_dim(XB);
		for (int j = 0; j < n; j++)
		{
			// wektor jednostkowy e_j
			matrix ej(n, 1, 0.0);
			ej(j, 0) = 1.0;

			// Testowanie x + s * e_j
			solution temp_forward = XB;

			// w kierunku dodatnim
			temp_forward.x = XB.x + s * ej; 
			temp_forward.fit_fun(ff, ud1, ud2);

			if (temp_forward.y < XB.y)
			{
				XB = temp_forward;
			}
			else
			{
				// Testowanie x - s * e_j
				solution temp_backward = XB;

				// w kierunku ujemnym
				temp_backward.x = XB.x - s * ej;  
				temp_backward.fit_fun(ff, ud1, ud2);

				if (temp_backward.y < XB.y) 
				{
					XB = temp_backward; 
				}
			}
		}

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
