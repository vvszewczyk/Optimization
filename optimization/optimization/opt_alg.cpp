﻿#include "opt_alg.h"

solution MC(matrix (*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1,
            matrix ud2)
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
        throw("solution MC(...):\n" + ex_info);
    }
}

double *expansion(matrix (*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1,
                  matrix ud2)
{
    try
    {
        auto *p = new double[2]{0, 0};
        auto i = 0;
        solution X0(x0), X1(x0 + d);
        X0.fit_fun(ff, ud1, ud2);
        X1.fit_fun(ff, ud1, ud2);
        if (X1.y == X0.y)
        {
            p[0] = m2d(X0.x);
            p[1] = m2d(X1.x);
            return p;
        }
        else if (X1.y > X0.y)
        {
            d = -d;
            X1.x = x0 + d;
            X1.fit_fun(ff, ud1, ud2);
            if (X1.y >= X0.y)
            {
                p[0] = m2d(X0.x - d);
                p[1] = m2d(X1.x);
                return p;
            }
        }
        solution X2;
        while (true)
        {
            ++i;
            X2.x = x0 + pow(alpha, i) * d;
            X2.fit_fun(ff, ud1, ud2);
            if (X2.y >= X1.y || solution::f_calls > Nmax)
                break;
            X0 = X1;
            X1 = X2;
        }
        if (d > 0)
        {
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
        throw("double* expansion(...):\n" + ex_info);
    }
}

solution fib(matrix (*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
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
        throw("solution fib(...):\n" + ex_info);
    }
}

solution lag(matrix (*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax,
             matrix ud1, matrix ud2)
{
    try
    {
        solution Xopt; // Zmienna przechowująca rozwiązanie

        // Krok 2: inicjalizacja a(0), b(0)
        solution A(a), B(b), C;
        C.x = (a + b) / 2; // c(0) = (a + b) / 2

        solution D, D_old(a); // d(i) i kopia pomocnicza

        A.fit_fun(ff, ud1, ud2); // Wywołanie funkcji dla A
        B.fit_fun(ff, ud1, ud2); // Wywołanie funkcji dla B
        C.fit_fun(ff, ud1, ud2); // Wywołanie funkcji dla C

        double l, m;
        int i = 0; // Inicjalizacja licznika iteracji
        std::cout << "i = " << i << ": Range = " << m2d(B.x) - m2d(A.x) << std::endl;

        // Pętla główna optymalizacji Lagrange'a
        while (i < Nmax) // repeat
        {
            // Krok 4: Obliczanie licznika l
            l = m2d(A.y * (pow(m2d(B.x), 2) - pow(m2d(C.x), 2)) + B.y * (pow(m2d(C.x), 2) - pow(m2d(A.x), 2)) +
                    C.y * (pow(m2d(A.x), 2) - pow(m2d(B.x), 2)));

            // Krok 5: Obliczanie mianownika m
            m = m2d(A.y * (m2d(B.x) - m2d(C.x)) + B.y * (m2d(C.x) - m2d(A.x)) + C.y * (m2d(A.x) - m2d(B.x)));

            // Sprawdzenie, czy m <= 0
            if (m <= 0)
            {
                Xopt = D_old;   // Zwrócenie poprzedniego rozwiązania
                Xopt.flag = -1; // Flaga błędu (error)
                return Xopt;    // Zakończenie algorytmu
            }

            // Krok 9: Obliczanie nowego punktu d(i)
            D.x = 0.5 * l / m;
            D.fit_fun(ff, ud1, ud2); // Wywołanie funkcji dla D

            // Krok 10-19: Aktualizacja przedziału w zależności od pozycji d(i)
            if (A.x < D.x && D.x < C.x) // Krok 10: a(i) < d(i) < c(i)
            {
                if (D.y < C.y) // Krok 11: f(d(i)) < f(c(i))
                {
                    // A = A    // krok 12: a(i+1) = a(i)
                    B = C; // krok 14: b(i+1) = c(i)
                    C = D; // krok 13: c(i+1) = d(i) odwrotnie ponieważ nadpisujemy C
                }
                else // Krok 15
                {
                    A = D; // krok 16: a(i+1) = d(i)
                           // C = C    // krok 17: c(i+1) = c(i)
                           // D = D    // krok 18: d(i+1) = d(i)
                }
            } // krok 19 & 20
            else if (C.x < D.x && D.x < B.x) // Krok 21: c(i) < d(i) < b(i)
            {
                if (D.y < C.y) // Krok 22: f(d(i)) < f(c(i))
                {
                    A = C; // krok 23: a(i+1) = c(i)
                    C = D; // krok 24: c(i+1) = d(i)
                           // B = B     // krok 25: b(i+1) = b(i)
                }
                else // Krok 26
                {
                    // A = A;   // krok 27: a(i+1) = a(i)
                    // C = C;   // krok 28: c(i+1) = c(i)
                    B = D; // krok 29: b(i+1) = d(i)
                }
            } // krok 30
            else // Krok 31: Jeśli d(i) nie mieści się w przedziale [a, b]
            {
                Xopt = D_old;   // Zwrócenie poprzedniego rozwiązania
                Xopt.flag = -1; // Krok 32: Flaga błędu (error)
                if (m2d(B.x) - m2d(A.x) < 1)
                    Xopt.flag = 0;
                return Xopt; // Zakończenie algorytmu
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
        return Xopt;   // Zwrócenie optymalnego rozwiązania
    }
    catch (string ex_info)
    {
        throw("solution lag(...):\n" + ex_info); // Obsługa błędów
    }
}

solution HJ(matrix (*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax,
            matrix ud1, matrix ud2)
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
        throw("solution HJ(...):\n" + ex_info);
    }
}

solution HJ_trial(matrix (*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
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
        throw("solution HJ_trial(...):\n" + ex_info);
    }
}

solution Rosen(matrix (*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon,
               int Nmax, matrix ud1, matrix ud2)
{
    try
    {
        solution Xopt;
        Xopt.ud = trans(x0);
        solution XB(x0), X;
        int n = get_dim(XB);
        matrix l(n, 1), p(n, 1), s(s0), D = ident_mat(n);
        XB.fit_fun(ff, ud1, ud2);
        do
        {
            for (int i = 0; i < n; i++)
            {
                X.x = XB.x + s(i) * D[i];
                X.fit_fun(ff, ud1, ud2);
                if (X.y < XB.y)
                {
                    XB = X;
                    l(i) += s(i);
                    s(i) *= alpha;
                }
                else
                {
                    s(i) = -beta * s(i);
                    ++p(i);
                }
            }
            Xopt.ud.add_row(trans(XB.x));
            bool change = true;
            for (int i = 0; i < n; i++)
            {
                if (l(i) == 0 || p(i) == 0)
                {
                    change = false;
                    break;
                }
            }
            if (change)
            {
                matrix Q(n, n), v(n, 1);
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j <= i; j++)
                    {
                        Q(i, j) = l(i);
                    }
                }
                Q = Q * D;
                v = Q[0] / norm(Q[0]);
                D.set_col(v, 0);
                for (int i = 1; i < n; i++)
                {
                    matrix tmp(n, 1);
                    for (int j = 0; j < i; j++)
                    {
                        tmp = tmp + trans(Q[i]) * D[j] * D[j];
                    }
                    v = (Q[i] - tmp) / norm(Q[i] - tmp);
                    D.set_col(v, i);
                }
                s = s0;
                l = matrix(n, 1);
                p = matrix(n, 1);
            }
            double max_s = abs(s(0));
            for (int i = 0; i < n; i++)
            {
                if (max_s < abs(s(i)))
                {
                    max_s = abs(s(i));
                }
            }
            if (max_s < epsilon)
            {
                Xopt = XB;
                Xopt.flag = 0;
                break;
            }
            if (solution::f_calls > Nmax)
            {
                Xopt = XB;
                Xopt.flag = 1;
                break;
            }
        } while (true);
        return Xopt;
    }
    catch (string ex_info)
    {
        throw("solution Rosen(...):\n" + ex_info);
    }
}

solution pen(matrix (*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1,
             matrix ud2)
{
    try
    {
        solution Xopt(x0), Xprev(0);
        matrix S(2, 1);
        S(0) = c;
        S(1) = ud1(0);

        /*std::cout << "Debug pen: Initialization\n";
        std::cout << "Xopt.x = " << Xopt.x << "\n";
        std::cout << "Xprev.x = " << Xprev.x << "\n\n";
        std::cout << solution::f_calls<< "\n\n";
        */
        bool change = true;
        double check = 0;
        while (change)
        {
            ud2(0, 0) = c; // Aktualizacja wartości c w ud2
            Xprev = Xopt;
            Xopt = sym_NM(ff, Xopt.x, ud1(1), ud1(2), ud1(3), ud1(4), ud1(5), ud1(6), Nmax, ud1, ud2);
            c = c * dc;
            check = norm(Xopt.x - Xprev.x);
            if (check < epsilon)
            {
                change = false;
            }
        }

        // Debugowanie: Wyświetlenie końcowego rozwiązania
        // std::cout << "Final penalty solution:\n";
        // std::cout << "x = " << Xopt.x << "\n";
        // std::cout << "y (function value): " << Xopt.y << "\n";
        // std::cout << "Flag: " << Xopt.flag << "\n";

        return Xopt;
    }
    catch (string ex_info)
    {
        throw("solution pen(...):\n" + ex_info);
    }
}

solution sym_NM(matrix (*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma,
                double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
    try
    {
        solution Xopt;
        int DIM = get_len(x0); // Dimensionality of the problem

        // Step 1: Initialize simplex
        matrix p(DIM, DIM + 1);    // Simplex points
        matrix e = ident_mat(DIM); // Identity matrix for direction vectors

        p.set_col(x0, 0); // Initial point x0 as the first vertex
        for (int i = 1; i <= DIM; ++i)
        {
            // Construct remaining vertices of the simplex
            p.set_col(x0 + e[i - 1] * s, i);
        }

        // Evaluate function at each vertex of the simplex
        matrix f_values(DIM + 1, 1);

        for (int i = 0; i <= DIM; ++i)
        {
            Xopt.x = p[i];                                 // Set the current point in Xopt
            f_values(i) = m2d(Xopt.fit_fun(ff, ud1, ud2)); // Automatic increment inside fit_fun
        }

        // cout << "Initial simplex:\n";
        // for (int i = 0; i <= DIM; ++i)
        //{
        //	cout << "Vertex " << i << ": x = " << p[i] << ", f(x) = " << f_values(i)
        //<< "\n";
        // }

        double max_norm;
        do
        {
            // Debugowanie: Wyświetlenie aktualnego stanu sympleksu
            // cout << "Iteration: " << solution::f_calls << "\n";
            // cout << "Current simplex:\n";
            // for (int i = 0; i <= DIM; ++i)
            //{
            //	cout << "Vertex " << i << ": x = " << p[i] << ", f(x) = " <<
            // f_values(i) << "\n";
            //}
            // Step 2: Order vertices by function value
            int p_min = 0, p_max = 0;
            for (int i = 1; i <= DIM; ++i)
            {
                if (f_values(i) < f_values(p_min))
                    p_min = i;
                if (f_values(i) > f_values(p_max))
                    p_max = i;
            }

            // Calculate the centroid of the simplex excluding the worst point
            matrix centroid = matrix(DIM, 1, 0.0);
            for (int i = 0; i <= DIM; ++i)
            {
                if (i == p_max)
                    continue;
                centroid = centroid + p[i];
            }
            centroid = centroid / DIM;

            // Step 3: Reflect the worst point
            Xopt.x = centroid + alpha * (centroid - p[p_max]);  // Set Xopt to reflection point
            double f_reflect = m2d(Xopt.fit_fun(ff, ud1, ud2)); // Automatic increment
            // cout << "Reflection point: x = " << Xopt.x << ", f(x) = " << f_reflect
            // << "\n";

            if (f_reflect < f_values(p_min))
            {
                // Step 4a: Expansion
                Xopt.x = centroid + gamma * (Xopt.x - centroid);   // Set Xopt to expansion point
                double f_expand = m2d(Xopt.fit_fun(ff, ud1, ud2)); // Automatic increment
                // cout << "Expansion point: x = " << Xopt.x << ", f(x) = " << f_expand
                // << "\n";

                if (f_expand < f_reflect)
                {
                    // Accept expansion point
                    p.set_col(Xopt.x, p_max);
                    f_values(p_max) = f_expand;
                }
                else
                {
                    // Accept reflection point
                    p.set_col(centroid + alpha * (centroid - p[p_max]), p_max);
                    f_values(p_max) = f_reflect;
                }
            }
            else if (f_reflect < f_values(p_max))
            {
                // Step 4b: Accept reflection
                p.set_col(Xopt.x, p_max);
                f_values(p_max) = f_reflect;
            }
            else
            {
                // Step 4c: Contraction
                Xopt.x = centroid + beta * (p[p_max] - centroid);    // Set Xopt to contraction point
                double f_contract = m2d(Xopt.fit_fun(ff, ud1, ud2)); // Automatic increment
                // cout << "Contraction point: x = " << Xopt.x << ", f(x) = " <<
                // f_contract << "\n";

                if (f_contract < f_values(p_max))
                {
                    // Accept contraction point
                    p.set_col(Xopt.x, p_max);
                    f_values(p_max) = f_contract;
                }
                else
                {
                    // Step 4d: Shrink the simplex
                    // cout << "Shrinking simplex:\n";
                    for (int i = 0; i <= DIM; ++i)
                    {
                        if (i == p_min)
                            continue;
                        p.set_col(p[p_min] + delta * (p[i] - p[p_min]), i);
                        Xopt.x = p[i];                                 // Set Xopt for shrinking
                        f_values(i) = m2d(Xopt.fit_fun(ff, ud1, ud2)); // Automatic increment
                        // cout << "New vertex " << i << ": x = " << p[i] << ", f(x) = " <<
                        // f_values(i) << "\n";
                    }
                }
            }

            // Check termination condition
            max_norm = 0.0;
            for (int i = 0; i <= DIM; ++i)
            {
                if (i == p_min)
                    continue;
                double norm_diff = norm(p[p_min] - p[i]);
                if (norm_diff > max_norm)
                    max_norm = norm_diff;
            }

            // Debugowanie: Wyświetlenie maksymalnej normy
            // cout << "Max norm: " << max_norm << "\n";

            if (solution::f_calls > Nmax)
            {
                Xopt.flag = -2; // Max function calls exceeded
                break;
            }

        } while (max_norm >= epsilon);

        // Return the optimal solution
        int p_min = 0;
        for (int i = 1; i <= DIM; ++i)
        {
            if (f_values(i) < f_values(p_min))
                p_min = i;
        }

        Xopt.x = p[p_min];
        Xopt.y = f_values(p_min);
        Xopt.flag = 1; // Success

        // Debugowanie: Wyświetlenie końcowego rozwiązania
        // cout << "Final solution:\n";
        // cout << "x = " << Xopt.x << "\n";
        // cout << "y (function value): " << Xopt.y << "\n";
        // cout << "Flag: " << Xopt.flag << "\n";
        return Xopt;
    }
    catch (string ex_info)
    {
        throw("solution sym_NM(...):\n" + ex_info);
    }
}

// (Steepest Descent)
solution SD(matrix (*ff)(matrix, matrix, matrix), matrix (*gf)(matrix, matrix, matrix), matrix x0, double h0,
            double epsilon, int Nmax, matrix ud1, matrix ud2)
{
    try
    {
        // cout << "xpp1" << endl;
        solution Xopt;
        Xopt.x = x0;
        Xopt.y = Xopt.fit_fun(ff, ud1,
                              ud2); // Wartość funkcji celu w punkcie początkowym
        // cout << "Xopt.y check 1" << endl;
        matrix x_prev;
        double norm_diff;

        do
        {
            // Gradient w punkcie x
            matrix grad = Xopt.grad(gf, ud1, ud2);
            // std::cout << "\nDEBUG: GRAD "<< grad << std::endl;
            //  Kierunek d = -gradient
            matrix d = -grad;
            // std::cout << "\nDEBUG: d " << d << std::endl;
            //  Wersja stałokrokowa
            double h = h0;
            // std::cout << "\nDEBUG: h " << h << std::endl;

            // Wersja zmiennokrokowa
            if (h0 < 0)
            {
                //       | x0(0)   d0(0) |
                // ud2 = |               |
                //       | x0(1)   d0(1) |

                matrix ud2_temp(2, 2);
                ud2_temp(0, 0) = Xopt.x(0);
                ud2_temp(1, 0) = Xopt.x(1);
                ud2_temp(0, 1) = d(0);
                ud2_temp(1, 1) = d(1);
                double *expansionResults = expansion(ff, 0.0, 0.05, 1.2, Nmax, ud1, ud2_temp);
                solution sol_Golden =
                    golden(ff4T, expansionResults[0], expansionResults[1], epsilon, Nmax, ud1, ud2_temp);
                h = m2d(sol_Golden.x);
                std::cout << "Aktualny krok h: " << h << std::endl;
            }
            // Zapis poprzedniego punktu
            x_prev = Xopt.x;
            // std::cout << "\nDEBUG: Xopt.x \n" << Xopt.x << std::endl;
            //  Aktualizacja punktu w x

            matrix hd = h * d;
            // std::cout << "\nDEBUG: hd \n" << hd << std::endl;
            Xopt.x = Xopt.x + h * d;
            // cout << "Xopt.x: "<< Xopt.x << endl;
            //  Obliczenie wartości funkcji celu w nowym punkcie
            Xopt.y = Xopt.fit_fun(ff, ud1, ud2);
            // cout << "Xopt.y check 2" << endl;
            //  Sprawdzenie warunku stopu
            norm_diff = norm(Xopt.x - x_prev);
            // Wyniki
            if (norm_diff < epsilon || solution::f_calls > Nmax || solution::g_calls > Nmax)
            {
                Xopt.flag = (norm_diff < epsilon) ? 0 : -1;
                break;
            }
        } while (true);

        return Xopt;
    }
    catch (string ex_info)
    {
        throw("solution SD(...):\n" + ex_info);
    }
}

// (Conjugate Gradient)
solution CG(matrix (*ff)(matrix, matrix, matrix), matrix (*gf)(matrix, matrix, matrix), matrix x0, double h0,
            double epsilon, int Nmax, matrix ud1, matrix ud2)
{
    try
    {
        solution Xopt;
        int i = 0;
        Xopt.x = x0;
        Xopt.y = Xopt.fit_fun(ff, ud1,
                              ud2); // Wartość funkcji celu w punkcie początkowym
        Xopt.flag = 0;
        matrix g_prev = Xopt.grad(gf, ud1, ud2); // Gradient początkowy
        matrix d_prev = -g_prev;                 // Kierunek początkowy przeciwny do gradientu
        matrix x_prev;
        double norm_diff;

        /*std::cout << "DEBUG BEFORE LOOP: " << i
                << ", x: (" << Xopt.x(0, 0) << ", " << Xopt.x(1, 0) << ", " <<
           Xopt.x(2,0)<<"), "
                << ", f(x): " << m2d(Xopt.y)
                << ", f_calls: " << solution::f_calls
                << ", g_calls: " << solution::g_calls
                << ", flag: " << Xopt.flag
                << std::endl;*/

        do
        {
            x_prev = Xopt.x;
            // std::cout << "DEBUG: x_prev: " << x_prev << std::endl;

            // Metoda stałokrokowa
            double h = h0;
            // std::cout << "DEBUG: h: " << h << std::endl;
            //  Metoda zmiennokrokowa (złoty podział)
            if (h0 < 0)
            {
                matrix ud2_temp(2, 2);
                ud2_temp(0, 0) = Xopt.x(0);
                ud2_temp(1, 0) = Xopt.x(1);
                ud2_temp(0, 1) = d_prev(0);
                ud2_temp(1, 1) = d_prev(1);
                double *expansionResults = expansion(ff, 0.0, 0.05, 1.2, Nmax, ud1, ud2_temp);
                solution sol_Golden =
                    golden(ff, expansionResults[0], expansionResults[1], epsilon, Nmax, ud1, ud2_temp);
                h = m2d(sol_Golden.x);
            }

            matrix scale_h0(1, 1);
            scale_h0(0, 0) = h;
            Xopt.x = Xopt.x + (scale_h0 * d_prev);

            Xopt.y = Xopt.fit_fun(ff, ud1, ud2);
            /*std::cout << "DEBUG: Xopt.y: " << Xopt.y << std::endl;
            std::cout << "DEBUG: i = " << i << std::endl;*/
            i++;

            matrix g_new = Xopt.grad(gf, ud1, ud2);

            // Obliczenie współczynnika beta
            double numerator = norm(g_new) * norm(g_new);
            // std::cout << "DEBUG: numerator: " << numerator << std::endl;
            double denominator = norm(g_prev) * norm(g_prev);
            // std::cout << "DEBUG: denominator: " << denominator << std::endl;
            double beta = numerator / denominator;
            /*if(denominator != 0)
            {
                    beta = numerator / denominator;
            }*/

            // Aktualizacja kierunku
            matrix scale_beta(1, 1);
            scale_beta(0, 0) = beta;
            matrix d_new = (-g_new) + (scale_beta * d_prev);

            // Obliczenie normy różnicy punktów
            norm_diff = norm(Xopt.x - x_prev);

            // Debug Print
            // if (i % 100 == 0 || norm_diff < epsilon)
            //{ // Print co 100 iteracji lub na koniec
            //	std::cout << "CG Iteracja: " << i
            //		<< ", x: (" << Xopt.x(0, 0) << ", " << Xopt.x(1, 0) <<
            //", " << Xopt.x(2, 0) << "), "
            //		<< ", f(x): " << m2d(Xopt.y)
            //		<< ", norm_diff: " << norm_diff
            //		<< ", f_calls: " << solution::f_calls
            //		<< ", g_calls: " << solution::g_calls
            //		<< ", flag: " << Xopt.flag
            //		<< std::endl;
            //}
            g_prev = g_new;
            d_prev = d_new;

            if (norm_diff < epsilon || solution::f_calls > Nmax || solution::g_calls > Nmax)
            {
                Xopt.flag = (norm_diff < epsilon) ? 0 : -1;
                break;
            }
        } while (true);

        return Xopt;
    }
    catch (std::string ex_info)
    {
        throw("solution CG(...):\n" + ex_info);
    }
}

solution Newton(matrix (*ff)(matrix, matrix, matrix), matrix (*gf)(matrix, matrix, matrix),
                matrix (*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1,
                matrix ud2)
{
    try
    {
        solution Xopt;
        Xopt.x = x0;
        Xopt.y = Xopt.fit_fun(ff, ud1,
                              ud2); // Wartość funkcji celu w punkcie początkowym

        matrix x_prev;
        double norm_diff;

        do
        {
            // Gradient w punkcie x
            matrix grad = Xopt.grad(gf, ud1, ud2);
            // Hessian
            matrix H_inv;
            try
            {
                matrix H = Xopt.hess(Hf, ud1, ud2);
                H_inv = inv(H); // H^-1 (macierz odwrotna)
            }
            catch (const std::string &ex)
            {
                std::cerr << "Error computing hessian or hessian inverse: " << ex << std::endl;
                Xopt.flag = -1; // Nie można odwrócić macierzy
                break;
            }
            // Kierunek d = -H^-1*grad
            matrix d = -1 * H_inv * grad;
            // Metoda stałokrokowa
            double h = h0;

            // Metoda zmiennokrokowa (złoty podział)
            if (h0 < 0)
            {
                matrix ud2_temp(2, 2);
                ud2_temp(0, 0) = Xopt.x(0);
                ud2_temp(1, 0) = Xopt.x(1);
                ud2_temp(0, 1) = d(0);
                ud2_temp(1, 1) = d(1);
                double *expansionResults = expansion(ff, 0.0, 0.05, 1.2, Nmax, ud1, ud2_temp);
                solution sol_Golden =
                    golden(ff, expansionResults[0], expansionResults[1], epsilon, Nmax, ud1, ud2_temp);
                h = m2d(sol_Golden.x);
            }
            // Zapis poprzedniego punktu
            x_prev = Xopt.x;
            // Aktualizacja punktu w x
            Xopt.x = Xopt.x + h * d;
            // Obliczenie wartości funkcji celu w nowym punkcie
            Xopt.y = Xopt.fit_fun(ff, ud1, ud2);

            // Sprawdzenie warunku stopu
            norm_diff = norm(Xopt.x - x_prev);

            if (norm_diff < epsilon || solution::f_calls > Nmax || solution::g_calls > Nmax || solution::H_calls > Nmax)
            {
                Xopt.flag = (norm_diff < epsilon) ? 0 : -1;
                break;
            }
        } while (true);

        return Xopt;
    }
    catch (string ex_info)
    {
        throw("solution Newton(...):\n" + ex_info);
    }
}

solution golden(matrix (*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1,
                matrix ud2)
{
    try
    {
        solution Xopt;

        // Inicjalizacja parametrów
        double alpha = (sqrt(5.0) - 1.0) / 2.0; // Stała złotego podziału
        double c = b - alpha * (b - a);
        double d = a + alpha * (b - a);

        matrix fc = ff(matrix(c), ud1,
                       ud2); // Obliczenie wartości funkcji celu w punkcie c
        matrix fd = ff(matrix(d), ud1,
                       ud2); // Obliczenie wartości funkcji celu w punkcie d

        // Pętla iteracyjna
        while ((b - a) > epsilon)
        {
            if (m2d(fc) < m2d(fd))
            {
                // Zawężenie przedziału z prawej strony
                // a = a;
                b = d;
                d = c;
                c = b - alpha * (b - a);
                fd = fc;                      // Aktualizacja fd
                fc = ff(matrix(c), ud1, ud2); // Aktualizacja fc
            }
            else
            {
                // Zawężenie przedziału z lewej strony
                a = c;
                // b = b;
                c = d;
                d = a + alpha * (b - a);
                fc = fd;                      // Aktualizacja fc
                fd = ff(matrix(d), ud1, ud2); // Aktualizacja fd
            }
            if (solution::f_calls > Nmax)
            {
                Xopt.flag = -1; // Przekroczono maksymalną liczbę wywołań funkcji celu
                break;
            }
            Xopt.flag = 0;
        }

        // Ustawienie optymalnego punktu w wyniku
        double x_opt = (a + b) / 2.0;
        Xopt.x = matrix(1, 1, x_opt);
        Xopt.y = ff(Xopt.x, ud1, ud2);

        return Xopt;
    }
    catch (string ex_info)
    {
        throw("solution golden(...):\n" + ex_info);
    }
}

solution Powell(matrix (*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
    try
    {
        solution Xopt;
        int n = get_len(x0);

        // Baza kierunków: d1 = [1,0]^T, d2 = [0,1]^T
        std::vector<matrix> d;

        for (int i = 0; i < n; i++)
        {
            matrix direction(n, 1, 0.0);
            direction(i, 0) = 1.0;
            d.push_back(direction);
        }

        matrix x = x0; // Punkt aktualny
        int i = 0;     // Licznik

        while (i < Nmax)
        {
            matrix p0 = x; // p0(i) = x(i)

            // Przejście po kierunkach
            for (int j = 0; j < n; j++)
            {
                // ud2 2x2: wiersz 0: punkt startowy, wiersz 1: kierunek
                matrix ud2_line(2, n);

                for (int k = 0; k < n; k++)
                {
                    ud2_line(0, k) = x(k, 0);
                    ud2_line(1, k) = d[j](k, 0);
                }

                // Ekspansja, aby wyznaczyć przedział
                double h0 = 0.0;
                double step = 1.0;
                double alpha = 1.2;
                std::unique_ptr<double[]> interval(expansion(ff, h0, step, alpha, Nmax, ud1, ud2_line));
                double a = interval[0];
                double b = interval[1];

                // Złoty podział
                solution sol_line = golden(ff, a, b, epsilon, Nmax, ud1, ud2_line);
                double h = m2d(sol_line.x);

                // Aktualizacja punktu
                for (int k = 0; k < n; k++)
                {
                    x(k, 0) += d[j](k, 0) * h;
                }
            }

            matrix pn = x; // pn(i) po przejściu wszystkich kierunków

            double norm_diff = 0.0;

            for (int k = 0; k < n; k++)
            {
                norm_diff += std::pow(pn(k, 0) - p0(k, 0), 2);
            }

            norm_diff = std::sqrt(norm_diff);

            if (norm_diff < epsilon)
            {
                Xopt.x = x;
                // Wywołujemy ff z ud2 = NaN, aby otrzymać f1 i f2
                matrix ud2_empty(1, 1, std::numeric_limits<double>::quiet_NaN());
                Xopt.y = ff(x, ud1, ud2_empty);
                Xopt.flag = 0;
                return Xopt;
            }

            if (solution::f_calls > Nmax)
            {
                Xopt.x = x;
                matrix ud2_empty(1, 1, std::numeric_limits<double>::quiet_NaN());
                Xopt.y = ff(x, ud1, ud2_empty);
                Xopt.flag = -1;
                return Xopt;
            }

            // Aktualizacja bazy kierunków
            if (n > 1)
            {
                for (int j = 0; j < n - 1; ++j)
                {
                    d[j] = d[j + 1];
                }
            }
            matrix d_new(n, 1);
            for (int k = 0; k < n; ++k)
            {
                d_new(k, 0) = pn(k, 0) - p0(k, 0);
            }
            d[n - 1] = d_new;

            // Wyznaczanie h_{n+1} wzdłuż nowego kierunku d_n
            matrix ud2_new_line(2, n);
            for (int k = 0; k < n; ++k)
            {
                ud2_new_line(0, k) = pn(k, 0);
                ud2_new_line(1, k) = d[n - 1](k, 0);
            }

            double h0_2 = 0.0;
            double step_2 = 1.0;
            double alpha_2 = 2.0;
            std::unique_ptr<double[]> interval2(expansion(ff, h0_2, step_2, alpha_2, Nmax, ud1, ud2_new_line));
            double a2 = interval2[0];
            double b2 = interval2[1];

            solution sol_line2 = golden(ff, a2, b2, epsilon, Nmax, ud1, ud2_new_line);
            double h2 = m2d(sol_line2.x);

            // Aktualizacja punktu
            for (int k = 0; k < n; ++k)
            {
                x(k, 0) += d[n - 1](k, 0) * h2;
            }

            i++;
        }

        // Teoretycznie nigdy tu nie dojdziemy, ale jakby co:
        Xopt.x = x;
        matrix ud2_empty(1, 1, std::numeric_limits<double>::quiet_NaN());
        Xopt.y = ff(x, ud1, ud2_empty);
        Xopt.flag = -1;
        return Xopt;
    }
    catch (string ex_info)
    {
        throw("solution Powell(...):\n" + ex_info);
    }
}

solution EA(matrix (*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0,
            double epsilon, int Nmax, matrix ud1, matrix ud2)
{
    try
    {
        // Parametry mutacji samoadaptacyjnej
        double alpha = std::sqrt((double)N);
        double beta = std::pow(2.0 * (double)N, -0.25);

        // Zerujemy liczniki (o ile chcemy śledzić f_calls)
        solution::clear_calls();

        // ========================
        // Funkcje do losowania (lokalnie, bez makr)
        // ========================
        auto randU = []() -> double {
            // U(0,1)
            return (double)rand() / (double)RAND_MAX;
        };

        auto randN = [&]() -> double {
            // N(0,1) metodą Box-Muller (wersja uproszczona)
            double u1 = randU();
            double u2 = randU();
            return sqrt(-2.0 * log(u1)) * cos(2.0 * 3.1415 * u2);
        };

        // ========================
        // 1) Inicjalizacja populacji
        // ========================
        vector<solution> population(mi);

        for (int i = 0; i < mi; i++)
        {
            // Tworzymy osobnika s z wektorem x (N×1)
            solution s(matrix(N, 1));
            // Rezerwujemy w s.ud miejsce na wektor sigma (N×1)
            s.ud = matrix(N, 1);

            // Losujemy każdą współrzędną x(d,0) w przedziale [ lb(d,0), ub(d,0) ]
            for (int d = 0; d < N; d++)
            {
                double r = randU();
                double xd = lb(d, 0) + r * (ub(d, 0) - lb(d, 0));
                s.x(d, 0) = xd;
                // Początkowa sigma
                s.ud(d, 0) = sigma0(d, 0);
            }

            // Obliczamy wartość funkcji celu:
            s.fit_fun(ff, ud1, ud2);
            // => s.y(0,0) będzie mieć f(x)

            population[i] = s;
        }

        // Pomocnicza lambda do znalezienia indeksu najlepszego (po s.y(0,0))
        auto bestIndex = [&](const vector<solution> &pop) {
            int bIdx = 0;
            double bVal = pop[0].y(0, 0);
            for (int j = 1; j < (int)pop.size(); j++)
            {
                double val = pop[j].y(0, 0);
                if (val < bVal)
                {
                    bVal = val;
                    bIdx = j;
                }
            }
            return bIdx;
        };

        // ========================
        // Główna pętla
        // ========================
        while (true)
        {
            // Sprawdzamy warunek stopu
            int bIdxPop = bestIndex(population);
            double fBest = population[bIdxPop].y(0, 0);

            // Czy f(x*) <= epsilon?
            if (fBest <= epsilon)
            {
                // Sukces
                population[bIdxPop].flag = 0;
                return population[bIdxPop];
            }

            // Czy przekroczyliśmy Nmax wywołań funkcji celu?
            if (solution::f_calls >= Nmax)
            {
                // Koniec – limit wyczerpany
                population[bIdxPop].flag = 1;
                return population[bIdxPop];
            }

            // ========================
            // 2) Koło ruletki: fi_j = 1 / f_j
            // ========================
            vector<double> fi(mi), q(mi + 1, 0.0);
            double sumFi = 0.0;
            for (int j = 0; j < mi; j++)
            {
                double fj = population[j].y(0, 0);
                fi[j] = 1.0 / fj;
                sumFi += fi[j];
            }
            for (int j = 1; j <= mi; j++)
            {
                q[j] = q[j - 1] + (fi[j - 1] / sumFi);
            }

            // ========================
            // 3) Generowanie potomstwa (lambda)
            // ========================
            vector<solution> offspring(lambda);

            for (int off = 0; off < lambda; off++)
            {
                // Wybór rodzica A
                double rA = randU();
                int idxA = 0;
                for (int j = 1; j <= mi; j++)
                {
                    if (rA <= q[j])
                    {
                        idxA = j - 1;
                        break;
                    }
                }
                // Wybór rodzica B
                double rB = randU();
                int idxB = 0;
                for (int j = 1; j <= mi; j++)
                {
                    if (rB <= q[j])
                    {
                        idxB = j - 1;
                        break;
                    }
                }

                // Krzyżowanie: child.x(d,0) = A.x(d,0) + (1-rC)*B.x(d,0)
                double rC = randU();
                solution child(matrix(N, 1));
                child.ud = matrix(N, 1);

                for (int d = 0; d < N; d++)
                {
                    double xA = population[idxA].x(d, 0);
                    double xB = population[idxB].x(d, 0);
                    double sA = population[idxA].ud(d, 0);
                    double sB = population[idxB].ud(d, 0);

                    child.x(d, 0) = xA + (1.0 - rC) * xB;
                    // Sigma np. średnia
                    child.ud(d, 0) = 0.5 * (sA + sB);
                }

                // Mutacja samoadaptacyjna
                double globalN = randN(); // globalne N(0,1)
                for (int d = 0; d < N; d++)
                {
                    double localN = randN();
                    double sigma_old = child.ud(d, 0);

                    // sigma_d = sigma_d * exp(alpha*gN + beta*lN)
                    double sigma_new = sigma_old * exp(alpha * globalN + beta * localN);
                    child.ud(d, 0) = sigma_new;

                    // x_d = x_d + sigma_d * N(0,1)
                    double step = sigma_new * randN();
                    child.x(d, 0) += step;

                    // Granice
                    if (child.x(d, 0) < lb(d, 0))
                        child.x(d, 0) = lb(d, 0);
                    if (child.x(d, 0) > ub(d, 0))
                        child.x(d, 0) = ub(d, 0);
                }

                // Obliczamy f(child.x)
                child.fit_fun(ff, ud1, ud2);

                offspring[off] = child;
            }

            // ========================
            // 4) Selekcja (μ+λ) -> μ
            // ========================
            vector<solution> combined;
            combined.reserve(mi + lambda);

            // stare osobniki:
            for (int i = 0; i < mi; i++)
                combined.push_back(population[i]);
            // potomki
            for (int i = 0; i < lambda; i++)
                combined.push_back(offspring[i]);

            // Sortowanie rosnąco po y(0,0)
            std::sort(combined.begin(), combined.end(),
                      [&](const solution &A, const solution &B) { return A.y(0, 0) < B.y(0, 0); });

            // Pierwsze mi -> nowa populacja
            for (int i = 0; i < mi; i++)
                population[i] = combined[i];
        }
    }
    catch (string ex_info)
    {
        throw("solution EA(...):\n" + ex_info);
    }
}
