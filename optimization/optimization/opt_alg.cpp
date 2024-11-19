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
		Xopt.ud = trans(x0);
		solution XB, XBOld, X(x0);
		X.fit_fun(ff, ud1, ud2);
		while (true)
		{
			XB = X;
			X = HJ_trial(ff, XB, s, ud1, ud2);
			if (X.y < XB.y)
			{
				while (true)
				{
					XBOld = XB;
					XB = X;

					X.x = 2 * XB.x - XBOld.x;
					X.fit_fun(ff, ud1, ud2);
					X = HJ_trial(ff, X, s, ud1, ud2);

					if (solution::f_calls > Nmax)
					{
						Xopt = XB;
						Xopt.flag = 0;
						return Xopt;
					}
					if (X.y >= XB.y)
					{
						break;
					}
				}
				XB = X;
			}
			else
			{
				s = alpha * s;
			}

			if (solution::f_calls > Nmax)
			{
				Xopt = XB;
				Xopt.flag = 0;
				break;
			}

			if (s < epsilon)
			{
				Xopt = XB;
				Xopt.flag = 1;
				break;
			}
		}
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
		matrix e = ident_mat(n);
		solution X;
		for (int j = 0; j < n; j++)
		{
			X.x = XB.x + s * e[j]; // Przesunięcie w kierunku e_j
			X.fit_fun(ff, ud1, ud2);

			if (X.y < XB.y)
			{
				XB = X;
			}
			else
			{
				X.x = XB.x - s * e[j]; // Przesunięcie w przeciwnym kierunku e_j
				X.fit_fun(ff, ud1, ud2);

				if (X.y < XB.y)
				{
					XB = X;
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
		Xopt.ud = trans(x0);
		solution XB(x0), X;
		int n = get_dim(XB);
		matrix l(n,1), p(n,1), s(s0), D = ident_mat(n);
		XB.fit_fun(ff, ud1, ud2);
		do {
			for(int i = 0; i < n; i++) {
				X.x = XB.x +s(i)*D[i];
				X.fit_fun(ff, ud1, ud2);
				if(X.y<XB.y) {
					XB=X;
					l(i) += s(i);
					s(i) *= alpha;
				}
				else {
					s(i) = - beta * s(i);
					++p(i);
				}	
			}
			Xopt.ud.add_row(trans(XB.x));
			bool change = true;
			for(int i = 0; i < n; i++) {
				if(l(i) == 0 || p(i) ==0) {
					change = false;
					break;
				}
			}
			if(change) {
				matrix Q(n,n), v(n,1);
				for(int i = 0; i < n; i++) {
					for(int j = 0; j <= i; j++) {
						Q(i,j) = l(i);
					}
				}
				Q = Q*D;
				v = Q[0] / norm(Q[0]);
				D.set_col(v, 0);
				for (int i = 1; i < n; i++) {
					matrix tmp(n,1);
					for(int j = 0; j < i; j++) {
						tmp = tmp + trans(Q[i])*D[j]*D[j];
					}
					v = (Q[i]-tmp)/norm(Q[i]-tmp);
					D.set_col(v,i);
				}
				s = s0;
				l = matrix(n,1);
				p = matrix(n,1);
			}
			double max_s = abs(s(0));
			for(int i = 0; i < n; i++) {
				if(max_s<abs(s(i))) {
					max_s = abs(s(i));
				}
			}
			if(max_s < epsilon) {
				Xopt = XB;
				Xopt.flag = 0;
				break;
			}
			if(solution::f_calls > Nmax) {
				Xopt = XB;
				Xopt.flag = 1;
				break;
			}
		}while(true);
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
		solution Xopt(x0); 
		solution Xprev;   
		int iter = 0;      

		// Statyczne zmienne do przechowywania danych
		static matrix(*global_ff)(matrix, matrix, matrix) = nullptr;
		static double global_c = 0.0;
		static matrix global_ud1, global_ud2;

		// Funkcja celu z karą jako funkcja globalna
		auto penalized_ff = [](matrix x, matrix ud1, matrix ud2) -> matrix {
			// Obliczenie ograniczeń
			matrix g = global_ff(x, global_ud1, global_ud2);

			// Obliczanie kary zewnętrznej
			double penalty = 0.0;
			int constraints = get_len(g);
			for (int i = 0; i < constraints; ++i)
			{
				penalty += std::max(0.0, g(i));
			}

			// Obliczanie wartości funkcji celu z karą
			matrix f_value = global_ff(x, global_ud1, global_ud2);
			return f_value + global_c * penalty;
		};

		// Przypisanie wskaźników globalnych
		global_ff = ff;
		global_ud1 = ud1;
		global_ud2 = ud2;

		while (iter < Nmax)
		{
			// Aktualizacja współczynnika kary w globalnym kontekście
			global_c = c;

			// Wyznaczanie minimum funkcji z karą
			Xprev = Xopt;
			Xopt.fit_fun(penalized_ff, ud1, ud2);

			// Sprawdzenie kryterium zakończenia
			if (norm(Xopt.x - Xprev.x) < epsilon)
			{
				Xopt.flag = 0; 
				return Xopt;
			}
			c *= dc; // Aktualizacja współczynnika kary
			iter++;
			if (iter >= Nmax)
			{
				Xopt.flag = 1; 
				break;
			}
		}

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
		// Initialize the simplex
		int* size = get_size(x0);
		int n = size[0]; // Assuming x0 is a column vector
		delete[] size;

		std::vector<matrix> simplex(n + 1, x0); // Simplex vertices
		for (int i = 0; i < n; ++i) {
			matrix unit_vector(n, 1, 0.0);
			unit_vector(i, 0) = 1.0;
			simplex[i + 1] = x0 + s * unit_vector;
		}

		// Create a solution instance
		solution sol;

		// Iterative Nelder-Mead Process
		int f_calls = 0;
		int min_idx = 0;
		int max_idx = 0;
		matrix pmin, pmax;

		while (f_calls < Nmax) {
			// Evaluate function at simplex vertices
			std::vector<double> f_values;
			for (const auto& vertex : simplex) {
				f_values.push_back(m2d(sol.fit_fun(ff, vertex, ud1)));
			}

			// Find indices of pmin and pmax
			auto min_max = std::minmax_element(f_values.begin(), f_values.end());
			min_idx = std::distance(f_values.begin(), min_max.first); // Declare and initialize min_idx
			max_idx = std::distance(f_values.begin(), min_max.second); // Declare and initialize max_idx

			// Define and initialize pmin and pmax
			matrix pmin = simplex[min_idx]; // pmin initialized using min_idx
			matrix pmax = simplex[max_idx]; // pmax initialized using max_idx

			// Centroid calculation
			matrix centroid(n, 1, 0.0);
			for (int i = 0; i <= n; ++i) {
				if (i != max_idx) centroid = centroid + simplex[i];
			}
			centroid = centroid / n;

			// Reflection
			matrix podb = centroid + alpha * (centroid - pmax);
			double f_podb = m2d(sol.fit_fun(ff, podb, ud1));
			++f_calls;

			if (f_podb < f_values[min_idx]) {
				// Expansion
				matrix pe = centroid + gamma * (podb - centroid);
				double f_pe = m2d(sol.fit_fun(ff, pe, ud1));
				++f_calls;

				if (f_pe < f_podb) simplex[max_idx] = pe;
				else simplex[max_idx] = podb;

			}
			else if (f_podb < f_values[max_idx]) {
				simplex[max_idx] = podb;
			}
			else {
				matrix pz = centroid + beta * (pmax - centroid);
				double f_pz = m2d(sol.fit_fun(ff, pz, ud1));
				++f_calls;

				if (f_pz < f_values[max_idx]) simplex[max_idx] = pz;
				else {
					for (int i = 0; i <= n; ++i) {
						if (i != min_idx) simplex[i] = delta * (simplex[i] + pmin);
					}
				}
			}

			// Termination check
			double max_dist = 0.0;
			for (int i = 0; i <= n; ++i) {
				max_dist = std::max(max_dist, norm(simplex[i] - pmin));
			}
			if (max_dist < epsilon) break;
		}

		// Return solution
		solution Xopt;
		Xopt.x = simplex[min_idx];
		Xopt.y = sol.fit_fun(ff, pmin, ud1);
		Xopt.flag = (f_calls >= Nmax) ? 1 : 0;
		Xopt.f_calls = f_calls;
		return Xopt;
	}
	catch (std::string ex_info)
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
