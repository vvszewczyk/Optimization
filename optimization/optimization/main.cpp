/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

#include "opt_alg.h"
#include <iostream>
#include <random>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

bool isGlobalMinimum(const solution &sol, double tol)
{
    double x1 = sol.x(0, 0);
    double x2 = sol.x(1, 0);
    // Oczekiwane min: (3, -2) => y=0
    // Tu sprawdzamy, czy y < tol
    double yVal = sol.y(0, 0);
    return (yVal <= tol);
}

void save_to_file_lab_5(std::ofstream &file, double a_val, double w, const matrix &x0, const solution &sol,
                        double f1_val, double f2_val, long f_calls)
{
    file << a_val << ","       // wartość a
         << w << ","           // waga w
         << x0(0, 0) << ","    // początkowy x1
         << x0(1, 0) << ","    // początkowy x2
         << sol.x(0, 0) << "," // rozwiązanie x1
         << sol.x(1, 0) << "," // rozwiązanie x2
         << f1_val << ","      // wartość f1
         << f2_val << ","      // wartość f2
         << f_calls << ","     // liczba wywołań funkcji celu
         << sol.flag << "\n";  // flaga statusu
}

int main()
{
    try
    {
        lab6();
    }
    catch (string EX_INFO)
    {
        cerr << "ERROR:\n";
        cerr << EX_INFO << endl << endl;
    }
    system("pause");
    return 0;
}

void lab0()
{
    // Funkcja testowa
    double epsilon = 1e-2;
    int Nmax = 10000;
    matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);
    solution opt;
    a(0) = -1;
    a(1) = 2;
    opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
    cout << opt << endl << endl;
    solution::clear_calls();

    // Wahadlo
    Nmax = 1000;
    epsilon = 1e-2;
    lb = 0;
    ub = 5;
    double teta_opt = 1;
    opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
    cout << opt << endl << endl;
    solution::clear_calls();

    // Zapis symulacji do pliku csv
    matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{m2d(opt.x), 0.5});
    matrix *Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
    ofstream Sout("symulacja_lab0.csv");
    Sout << hcat(Y[0], Y[1]);
    Sout.close();
    Y[0].~matrix();
    Y[1].~matrix();
}

void lab1()
{
    // Generator liczb losowych
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(-100, 100);

    // Delimiter
    std::string delimiter = ",";

    double *expansionResults = new double[2];
    double d = 1.0;

    // Współczynniki ekspansji
    double alpha1 = 1.01;
    double alpha2 = 1.2;
    double alpha3 = 1.3;

    int Nmax = 10000;
    double epsilon = 1e-5;
    double gamma = 1e-180;

    ofstream expToFile("output/lab1/expansion.csv");
    ofstream fibToFile("output/lab1/fibonacci.csv");
    ofstream lagToFile("output/lab1/lagrange.csv");

    // TABELA 1, 2
    for (int i = 0; i < 100; i++)
    {
        int x0 = dis(gen);

        expansionResults = expansion(df1, x0, d, alpha3, Nmax, 0);
        expToFile << x0 << delimiter << expansionResults[0] << delimiter << expansionResults[1] << delimiter
                  << solution::f_calls << endl;
        solution::clear_calls();

        solution fibonacci = fib(df1, expansionResults[0], expansionResults[1], epsilon);
        double x_fib = m2d(fibonacci.x);
        double y_fib = m2d(fibonacci.y);
        std::string minimum_type_fib = (x_fib >= 62.72 && x_fib <= 62.73) ? "GLOBALNE" : "LOKALNE";
        fibToFile << x_fib << delimiter << y_fib << delimiter << solution::f_calls << delimiter << minimum_type_fib
                  << "\n";
        solution::clear_calls();

        solution lagrange = lag(df1, expansionResults[0], expansionResults[1], epsilon, gamma, Nmax);
        double x_lag = m2d(lagrange.x);
        double y_lag = m2d(lagrange.y);
        std::string minimum_type_lag = (x_lag >= 62.71 && x_lag <= 62.74) ? "GLOBALNE" : "LOKALNE";
        std::string minimum_type_lag2 = (lagrange.flag == -1) ? "ERROR" : minimum_type_lag;
        lagToFile << x_lag << delimiter << y_lag << delimiter << solution::f_calls << delimiter << minimum_type_lag2
                  << "\n";
        solution::clear_calls();
    }

    // Wykres + ostatni rekord z tabeli 1
    solution fibonacci2 = fib(df1, -100, 100, epsilon);

    double x_fib = m2d(fibonacci2.x);
    double y_fib = m2d(fibonacci2.y);
    std::string minimum_type_fib2 = (x_fib >= 62.72 && x_fib <= 62.73) ? "GLOBALNE" : "LOKALNE";
    std::cout << x_fib << delimiter << y_fib << delimiter << solution::f_calls << delimiter << minimum_type_fib2
              << "\n";
    std::cout << "Fibonacci" << std::endl;
    std::cout << fibonacci2;
    solution::clear_calls();

    solution lagrange2 = lag(df1, -100, 100, epsilon, gamma, Nmax);
    double x_lag = m2d(lagrange2.x);
    double y_lag = m2d(lagrange2.y);
    std::string minimum_type_lag2 = (x_lag >= 62.72 && x_lag <= 62.73) ? "GLOBALNE" : "LOKALNE";
    std::string minimum_type_lag3 = (lagrange2.flag == -1) ? "ERROR" : minimum_type_lag2;
    std::cout << x_lag << delimiter << y_lag << delimiter << solution::f_calls << delimiter << minimum_type_lag3
              << "\n";
    std::cout << "Lagrange\n";
    std::cout << lagrange2;
    solution::clear_calls();

    // Tabela 3

    std::cout << "\n\nTABELA 3\n\n";
    solution fibEx2 = fib(f2, 1e-4, 1e-2, epsilon);
    std::cout << "Fib" << std::endl << fibEx2 << endl;

    solution::clear_calls();

    solution lagEx2 = lag(f2, 1e-4, 1e-2, epsilon, gamma, 2000);
    std::cout << "Lag" << std::endl << lagEx2 << endl;

    // Symulacja
    std::cout << "\n\nSymulacja\n\n";
    ofstream simFib("output/lab1Symulacja_fibonacci.csv");
    ofstream simLag("output/lab1Symulacja_lagrange.csv");

    matrix symFib = matrix(3, new double[3]{5, 1, 20});
    matrix *solved1 = solve_ode(f1, 0, 1, 2000, symFib, NAN, fibEx2.x(0));
    simFib << solved1[1] << std::endl << std::endl;
    std::cout << solved1[1] << std::endl << std::endl;

    matrix symLag = matrix(3, new double[3]{5, 1, 20});
    matrix *solved2 = solve_ode(f1, 0, 1, 2000, symLag, NAN, lagEx2.x(0));
    simLag << solved2[1] << std::endl;
    std::cout << solved2[1] << std::endl;

    delete[] expansionResults;
}

void lab2()
{
    // -------- Funkcja Testująca (Tabele 1 i 2) --------- //
    double epsilon = 1e-06;
    int Nmax = 1000;
    double alphaHJ = 0.5;
    double alphaR = 2.0;
    double beta = 0.5;
    std::vector<double> steps = {0.1, 0.5, 1.0};
    std::string delimiter = ",";

    // Otwarcie dwóch osobnych plików CSV dla każdej z metod
    std::ofstream fileHJ("output/lab2/Hooke-Jeeves.csv");
    std::ofstream fileRosen("lab2/Rosenbrock.csv");

    if (!fileHJ.is_open() || !fileRosen.is_open())
    {
        std::cerr << "Nie mozna otworzyc pliku poczatkowych\n";
        return;
    }

    // Nagłówki kolumn z separatorem średnika
    fileHJ << "Metoda,Dlugosc kroku,x0(1),x0(2),x(1),x(2),y*,f_calls,Minimum Globalne\n";
    fileRosen << "Metoda,Dlugosc kroku,x0(1),x0(2),x(1),x(2),y*,f_calls,Minimum Globalne\n";

    // Generator liczb losowych
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1.0, 1.0);

    for (double step : steps)
    {
        for (int i = 0; i < 100; i++)
        {
            matrix x0(2, 1);
            x0(0, 0) = dis(gen);
            x0(1, 0) = dis(gen);

            // Testowanie metody Hooke'a-Jeevesa i zapis do pliku Hooke-Jeeves.csv
            solution y0HJ = HJ(ff2T, x0, step, alphaHJ, epsilon, Nmax);
            fileHJ << "Hooke-Jeeves" << delimiter << step << delimiter << x0(0, 0) << delimiter << x0(1, 0) << delimiter
                   << y0HJ.x(0, 0) << delimiter << y0HJ.x(1, 0) << delimiter << m2d(y0HJ.y) << delimiter
                   << solution::f_calls << delimiter << ((abs(m2d(y0HJ.y)) < epsilon) ? "TAK" : "NIE") << delimiter
                   << "\n";

            solution::clear_calls();

            // Testowanie metody Rosenbrocka i zapis do pliku Rosenbrock.csv
            matrix s0(2, 1, step);
            solution y0Rosen = Rosen(ff2T, x0, s0, alphaR, beta, epsilon, Nmax);
            fileRosen << "Rosenbrock" << delimiter << step << delimiter << x0(0, 0) << delimiter << x0(1, 0)
                      << delimiter << y0Rosen.x(0, 0) << delimiter << y0Rosen.x(1, 0) << delimiter << m2d(y0HJ.y)
                      << delimiter << solution::f_calls << delimiter
                      << ((abs(m2d(y0Rosen.y)) < epsilon) ? "TAK" : "NIE") << "\n";

            solution::clear_calls();
        }
    }

    // Zamknięcie plików
    fileHJ.close();
    fileRosen.close();

    std::cout << "Wyniki zapisane do plikow Hooke-Jeeves.csv i Rosenbrock.csv\n";

    // -------- Problem rzeczywisty (Tabela 3 i symulacja) --------- //
    // Otwarcie plików CSV dla wyników problemu rzeczywistego
    std::ofstream HookeSymulacja("output/lab2/HookeSymulacja.csv");
    std::ofstream RosenbrockSymulacja("output/lab2/RosenbrockSymulacja.csv");
    std::ofstream Symulacja("output/lab2/Symulacja.csv");

    if (!HookeSymulacja.is_open() || !RosenbrockSymulacja.is_open() || !Symulacja.is_open())
    {
        std::cerr << "Nie mozna otworzyc pliku symulacji\n";
        return;
    }

    double step = 0.8;
    double k_values[2] = {2.5, 7.5};
    matrix x0(2, k_values);

    solution::clear_calls();

    // Optymalizacja metodą Hooke-Jeeves
    solution HookeR = HJ(ff2R, x0, step, alphaHJ, epsilon, Nmax);
    HookeSymulacja << HookeR.x(0) << delimiter << HookeR.x(1) << delimiter << m2d(HookeR.y) << delimiter
                   << solution::f_calls << delimiter << HookeR.flag << "\n";
    solution::clear_calls();

    // Optymalizacja metodą Rosenbrocka
    solution RosenbrockR = Rosen(ff2R, x0, matrix(2, 1, step), alphaR, beta, epsilon, Nmax);
    RosenbrockSymulacja << RosenbrockR.x(0) << delimiter << RosenbrockR.x(1) << delimiter << m2d(RosenbrockR.y)
                        << delimiter << solution::f_calls << delimiter << RosenbrockR.flag << "\n";
    solution::clear_calls();

    // Zamknięcie plików wynikowych dla optymalizacji
    HookeSymulacja.close();
    RosenbrockSymulacja.close();

    // -------- Symulacja z optymalnymi parametrami -------- //
    // Parametry czasowe
    double t0 = 0;
    double dt = 0.1;
    double tend = 100;

    // Parametry przemieszczenia
    matrix y0(2, new double[2]{0.0, 0.0}); // początkowy

    // Symulacja dla wyników optymalnych z Hooke-Jeeves
    matrix *yz1 = solve_ode(df2, t0, dt, tend, y0, HookeR.x(0), HookeR.x(1));
    Symulacja << "Symulacja dla Hooke-Jeeves\n" << yz1[1] << "\n\n";
    delete[] yz1;
    solution::clear_calls();

    // Symulacja dla wyników optymalnych z Rosenbrocka
    matrix *yz2 = solve_ode(df2, t0, dt, tend, y0, RosenbrockR.x(0), RosenbrockR.x(1));
    Symulacja << "Symulacja dla Rosenbrock\n" << yz2[1] << "\n";
    delete[] yz2;
    solution::clear_calls();

    // Zamknięcie pliku symulacji
    Symulacja.close();

    std::cout << "Wyniki zapisane do plikow symulacji.\n";
}

void lab3()
{
    double epsilon = 1e-3;
    int Nmax = 1000;
    double s = 1.0;     // Rozmiar kroku
    double alpha = 1.0; // Współczynnik odbicia
    double beta = 0.5;  // Współczynnik kontrakcji
    double gamma = 2.0; // Współczynnik ekspansji
    double delta = 0.5; // Współczynnik zmniejszenia sympleksu
    double c = 1.0;     // Wstępna wartość kary w `pen`
    double dc = 2.0;    // Mnożnik zwiększający karę w `pen`

    // Lista wartości `a`
    std::vector<double> a_values = {4.0, 4.4934, 5.0};

    // Otwarcie pliku CSV
    std::ofstream results_file("output/lab3/results_table_combined.csv");

    if (!results_file.is_open())
    {
        std::cerr << "Nie można otworzyć pliku do zapisu wyników!\n";
        return;
    }

    // Nagłówki tabeli
    results_file << "x1(0),x2(0),"
                 << "x1*(pen),x2*(pen),r*(pen),y*(pen),Liczba wywołań funkcji celu (pen),"
                 << "x1*(sym_NM),x2*(sym_NM),r*(sym_NM),y*(sym_NM),Liczba wywołań funkcji celu (sym_NM)\n";

    // Generator punktów początkowych
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-5.0, 5.0); // Punkt początkowy w dopuszczalnym zakresie

    for (double a : a_values)
    {
        for (int i = 0; i < 100; i++)
        {
            matrix x0(2, 1);
            bool valid_point = false;
            while (!valid_point)
            {
                x0(0, 0) = dis(gen);
                x0(1, 0) = dis(gen);
                // Sprawdzenie ograniczeń
                double x1 = x0(0, 0);
                double x2 = x0(1, 0);
                double g1 = -x1 + 1;                           // g1(x1) <= 0
                double g2 = -x2 + 1;                           // g2(x2) <= 0
                double g3 = sqrt(pow(x1, 2) + pow(x2, 2)) - a; // g3(x1, x2) <= 0

                valid_point = (g1 <= 0) && (g2 <= 0) && (g3 <= 0);
            }
            // Przygotowanie `ud1` jako macierz 6x1 dla sym_NM i macierz 1x1 dla ograniczeń
            matrix ud1(7, 1);
            ud1(0, 0) = a;     // Parametr a dla ograniczeń
            ud1(1, 0) = s;     // Rozmiar kroku
            ud1(2, 0) = alpha; // Odbicie
            ud1(3, 0) = beta;  // Kontrakcja
            ud1(4, 0) = gamma; // Ekspansja
            ud1(5, 0) = delta; // Zmniejszenie sympleksu
            ud1(6, 0) = epsilon;

            matrix ud2(1, 1, c); // wstępna wartość kary

            try
            {
                // Wywołanie funkcji `pen`
                solution::clear_calls(); // Reset licznika funkcji celu
                solution opt_pen = pen(ff3T, x0, c, dc, epsilon, Nmax, ud1, ud2);
                double r_pen = std::sqrt(opt_pen.x(0, 0) * opt_pen.x(0, 0) + opt_pen.x(1, 0) * opt_pen.x(1, 0));
                cout << solution::f_calls << endl;
                // Wywołanie funkcji `sym_NM`
                matrix ud2_2(1, 1, r_pen); // wstępna wartość kary
                solution::clear_calls();   // Reset licznika funkcji celu
                solution opt_sym = sym_NM(ff3T, x0, ud1(1, 0), ud1(2, 0), ud1(3, 0), ud1(4, 0), ud1(5, 0), ud1(6, 0),
                                          Nmax, matrix(1, 1, a), ud2_2);
                double r_sym = std::sqrt(opt_sym.x(0, 0) * opt_sym.x(0, 0) + opt_sym.x(1, 0) * opt_sym.x(1, 0));

                // Zapis wyników do jednej linii
                results_file << x0(0, 0) << "," << x0(1, 0) << "," << opt_pen.x(0, 0) << "," << opt_pen.x(1, 0) << ","
                             << r_pen << "," << opt_pen.y(0, 0) << "," << solution::f_calls << "," << opt_sym.x(0, 0)
                             << "," << opt_sym.x(1, 0) << "," << r_sym << "," << opt_sym.y(0, 0) << ","
                             << solution::f_calls << "\n";
            }
            catch (const std::exception &e)
            {
                std::cerr << "Błąd podczas optymalizacji: " << e.what() << "\n";
            }
        }
    }

    results_file.close();
    std::cout << "Wyniki zapisane w pliku results_table_combined.csv\n";

    // Problem rzeczywisty
    epsilon = 1e-3; // Dokładność
    Nmax = 1000;    // Maksymalna liczba iteracji
    c = 10.0;       // Współczynnik kary początkowej
    dc = 10.0;      // Współczynnik zwiększenia kary

    // Ograniczenia dla v0x i omega
    matrix x0(2, new double[2]{5.0, .0}); // x0 = 5.0 omega = 0.0

    // Parametry do sympleksu przekazane przez ud1
    // s = 1.0, alpha = 1.0, beta = 0.5, gamma = 2.0, delta = 0.5, epsilon = 1e-3, Nmax = 1000
    matrix ud1(7, new double[7]{1.0, 1.0, 0.5, 2.0, 0.5, epsilon, static_cast<double>(Nmax)});

    // Współczynnik kary dla funkcji celu
    matrix ud2(1, 1, 100.0);

    // Debugowanie inicjalizacji
    // cout << "Debug: Initializing ud1 and ud2...\n";
    // cout << "ud1: " << ud1 << "\n";
    // cout << "ud2: " << ud2 << "\n";

    // Rozwiązanie problemu optymalizacji
    solution::clear_calls;
    solution opt = pen(ff3R, x0, c, dc, epsilon, Nmax, ud1, ud2);

    // Wyświetlenie wyników optymalizacji
    cout << "Optymalne wyniki:\n";
    cout << "v0x = " << opt.x(0) << " m/s, omega = " << opt.x(1) << " rad/s\n";
    cout << "Wartość funkcji celu (negatywne x_end): " << opt.y << "\n";
    cout << "Flaga zakończenia: " << opt.flag << "\n";
    cout << "Liczba wywołan funkcji celu: " << solution::f_calls << "\n";

    // Symulacja trajektorii dla optymalnych parametrów
    matrix Y0(4, new double[4]{0, opt.x(0), 100, 0}); // Warunki początkowe
    matrix *Y = solve_ode(df3, 0.0, 0.01, 7.0, Y0, matrix(), matrix(1, 1, opt.x(1)));

    // Zapis trajektorii do pliku CSV
    ofstream file("output/lab3/optimal_trajectory.csv");
    file << hcat(Y[0], Y[1]); // Łączenie czasu i współrzędnych w jeden plik
    file.close();

    // Wyświetlenie trajektorii
    // cout << "Trajektoria zapisana w pliku 'optimal_trajectory.csv'.\n";

    // Zwolnienie pamięci
    delete[] Y;
}

void lab4()
{
    // Parametry
    double epsilon = 1e-6;
    int Nmax = 1000;
    std::vector<double> h_values = {0.05, 0.12, -1};
    // std::vector<double> h_values = { -1.0 }; // for testing golden separetion
    std::string delimiter = ",";

    // Otwórz plik do zapisu wyników dla metody SD
    std::ofstream results_file_SD("output/lab4/lab4_SD_results.csv");
    results_file_SD << "Method,h0,x0(1),x0(2),x*(1),x*(2),y*,f_calls,g_calls,flag\n";

    // Otwórz plik do zapisu wyników dla metody CG
    std::ofstream results_file_CG("output/lab4/lab4_CG_results.csv");
    results_file_CG << "Method,h0,x0(1),x0(2),x*(1),x*(2),y*,f_calls,g_calls,flag\n";

    // Otwórz plik do zapisu wyników dla metody Newtona
    std::ofstream results_file_Newton("output/lab4/lab4_Newton_results.csv");
    results_file_Newton << "Method,h0,x0(1),x0(2),x*(1),x*(2),y*,f_calls,g_calls,H_calls,flag\n";

    // Generator liczb losowych
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-10.0, 10.0);

    for (double h0 : h_values)
    {
        for (int i = 0; i < 100; i++)
        {
            matrix x0(2, 1);
            x0(0, 0) = dis(gen);
            x0(1, 0) = dis(gen);
            double x0_expansion = dis(gen);

            // Metoda Gradientu Prostego (SD)
            solution::clear_calls();
            solution sol_SD = SD(ff4T, gf4T, x0, h0, epsilon, Nmax);

            results_file_SD << "SD" << delimiter << h0 << delimiter << x0(0, 0) << delimiter << x0(1, 0) << delimiter
                            << sol_SD.x(0, 0) << delimiter << sol_SD.x(1, 0) << delimiter << m2d(sol_SD.y) << delimiter
                            << solution::f_calls << delimiter << solution::g_calls << delimiter << sol_SD.flag << "\n";

            // Metoda Gradientów Sprzężonych (CG)
            solution::clear_calls();
            solution sol_CG = CG(ff4T, gf4T, x0, h0, epsilon, Nmax);

            results_file_CG << "CG" << delimiter << h0 << delimiter << x0(0, 0) << delimiter << x0(1, 0) << delimiter
                            << sol_CG.x(0, 0) << delimiter << sol_CG.x(1, 0) << delimiter << m2d(sol_CG.y) << delimiter
                            << solution::f_calls << delimiter << solution::g_calls << delimiter << sol_CG.flag << "\n";

            // Metoda Newtona
            solution::clear_calls();
            solution sol_Newton = Newton(ff4T, gf4T, hf4T, x0, h0, epsilon, Nmax);

            results_file_Newton << "Newton" << delimiter << h0 << delimiter << x0(0, 0) << delimiter << x0(1, 0)
                                << delimiter << sol_Newton.x(0, 0) << delimiter << sol_Newton.x(1, 0) << delimiter
                                << m2d(sol_Newton.y) << delimiter << solution::f_calls << delimiter << solution::g_calls
                                << delimiter << solution::H_calls << delimiter << sol_Newton.flag << "\n";
        }
    }

    // Zamknięcie plików
    results_file_SD.close();
    results_file_CG.close();
    results_file_Newton.close();

    std::cout << "Wyniki zapisane w output/lab4/ do plików lab4_SD_results.csv, lab4_CG_results.csv i "
                 "lab4_Newton_results.csv\n";

    ////PROBLEM RZECZYWISTU
    // int m = 3;  // Liczba przykładów
    // int n = 100; // Liczba cech w X (zakładając 80 cech)

    //// Wczytujemy dane
    // matrix X = loadDataToMatrix("XData.txt", m, n);
    // matrix Y = loadDataToMatrix1("YData.txt", n);  // Tylko 1 kolumna dla Y

    //// Parametry optymalizacji
    // double epsilon = 1e-6;  // Tolerancja
    // int Nmax = 1000;        // Maksymalna liczba iteracji
    // vector<double> step_sizes = { 0.01, 0.001, 0.0001 };  // Długości kroków

    //// Tworzymy macierz początkowych wartości (theta = [0, 0, 0])

    // matrix x0(3, 1, 0.0);

    // vector<vector<double>> results;  // Wyniki do zapisania w pliku CSV

    ///*double J_theta = m2d(logisticCostFunction(x0, X, Y));
    // cout << "J(theta) = " << J_theta << endl;

    // matrix grad = computeGradient(x0, X, Y);
    // cout << "Gradient J(theta):" << endl;
    // for (int i = 0; i < get_size(grad)[0]; ++i) {
    //	cout << "grad(" << i << ") = " << grad(i) << endl;
    // }*/

    //// Uruchamiamy optymalizację dla każdej długości kroku
    // for (double h0 : step_sizes) {
    //	// Spadek gradientu (SD)
    //	std::cout << "DEBUG: x0 " << x0 << std::endl;
    //	solution CG_result = CG(logisticCostFunction, computeGradient, x0, h0, epsilon, Nmax, X, Y);
    //	std::cout << "DEBUG: CG_result " << CG_result << std::endl;
    //	//std::cout << "DEBUG: CG_result: " << CG_result << std::endl;
    //	double accuracy_SD = computeAccuracy(X, Y, CG_result.x);  // Obliczanie dokładności

    //	// Zapisujemy wyniki do tabeli
    //	vector<double> row;  // Tworzymy pusty wektor typu std::vector<double>

    //	// Dodajemy elementy do wektora
    //	row.push_back(h0);  // Długość kroku
    //	row.push_back(m2d(CG_result.x(0, 0)));  // θ0
    //	row.push_back(m2d(CG_result.x(1, 0)));  // θ1
    //	row.push_back(m2d(CG_result.x(2, 0)));  // θ2
    //	row.push_back(m2d(logisticCostFunction(CG_result.x, X, Y)));  // J(θ*)
    //	row.push_back(accuracy_SD);  // P(θ*)
    //	row.push_back(static_cast<double>(solution::g_calls));  // g_calls jako double
    //	solution::clear_calls();
    //	// Dodajemy row (std::vector<double>) do results
    //	results.push_back(row);

    //}

    //// Zapisujemy wyniki do pliku CSV
    // writeResultsToCSV("optimization_results.csv", results);
    // cout << "Wyniki zapisane do pliku optimization_results.csv" << endl;
}

void lab5()
{
    //// Parametry
    // double epsilon = 1e-8;
    // int Nmax = 20000;  // Maksymalna liczba wywołań funkcji celu
    //// Wagi (w) od 0 do 1 z krokiem 0.01
    //// Parametry a = {1,10,100}

    // std::vector<double> a_values = { 1.0, 10.0, 100.0 };

    //// Przygotowanie pliku do zapisu wyników
    // std::ofstream results_file_a1("output/lab5/results_powell_test_a1.csv");
    // std::ofstream results_file_a10("output/lab5/results_powell_test_a10.csv");
    // std::ofstream results_file_a100("output/lab5/results_powell_test_a100.csv");
    // results_file_a1 << "a,w,x1_0,x2_0,sol_x1,sol_x2,f1,f2,f_calls,flag\n";
    // results_file_a10 << "a,w,x1_0,x2_0,sol_x1,sol_x2,f1,f2,f_calls,flag\n";
    // results_file_a100 << "a,w,x1_0,x2_0,sol_x1,sol_x2,f1,f2,f_calls,flag\n";

    //// Generator liczb losowych dla punktu startowego
    // std::random_device rd;
    // std::mt19937 gen(rd());
    // std::uniform_real_distribution<> dis(-10.0, 10.0);

    //// Pętla po wagach w
    // for (int i = 0; i <= 100; i++)
    //{
    //	double w = i * 0.01;

    //	// Losowy punkt startowy x0
    //	matrix x0(2, 1);
    //	x0(0, 0) = dis(gen);
    //	x0(1, 0) = dis(gen);

    //	// Pętla po wartościach a
    //	for (auto a_val : a_values)
    //	{
    //		// Ud1: tu przechowamy wagę w (ud1(0)) i parametr a (ud1(1))
    //		// Zgodnie z założeniami z transkryptu:
    //		// ud1(0) = waga, ud1(1) = a
    //		matrix ud1(2, 1);
    //		ud1(0) = w;
    //		ud1(1) = a_val;

    //		// ud2 pusty: ud2 będzie ustawiany w trakcie line search przez algorytm
    //		matrix ud2;

    //		// Czyścimy liczniki wywołań funkcji celu itp.
    //		solution::clear_calls();

    //		// Wywołujemy metodę Powella
    //		solution sol = Powell(ff5T, x0, epsilon, Nmax, ud1, ud2);

    //		// Utwórz macierz 1x1 wypełnioną NaN, aby uzyskać f1 i f2
    //		matrix ud2_empty(1, 1, std::numeric_limits<double>::quiet_NaN());
    //		matrix y_val = ff5T(sol.x, ud1, ud2_empty);

    //		double f1_val = y_val(0, 0);
    //		double f2_val = y_val(1, 0);

    //		/*std::cout <<"w\t"<< w << std::endl;
    //		std::cout <<"x01\t"<< x0(0, 0) << std::endl;
    //		std::cout <<"x02\t"<< x0(1, 0) << std::endl;
    //		std::cout << "solx1\t" << sol.x(0, 0) << std::endl;
    //		std::cout << "solx2\t" << sol.x(1, 0) << std::endl;
    //		std::cout << "f1\t" << f1_val << std::endl;
    //		std::cout << "f2\t" << f2_val << std::endl;
    //		std::cout << "fcals\t" << solution::f_calls << std::endl;
    //		std::cout << "flag\t" << sol.flag << std::endl<< std::endl;*/

    //		if (a_val == 1)
    //		{
    //			results_file_a1
    //				<< a_val << ","
    //				<< w << ","
    //				<< x0(0, 0) << ","
    //				<< x0(1, 0) << ","
    //				<< sol.x(0, 0) << ","
    //				<< sol.x(1, 0) << ","
    //				<< f1_val << ","
    //				<< f2_val << ","
    //				<< solution::f_calls << ","
    //				<< sol.flag << "\n";
    //		}
    //		else if (a_val == 10)
    //		{
    //			results_file_a10
    //				<< a_val << ","
    //				<< w << ","
    //				<< x0(0, 0) << ","
    //				<< x0(1, 0) << ","
    //				<< sol.x(0, 0) << ","
    //				<< sol.x(1, 0) << ","
    //				<< f1_val << ","
    //				<< f2_val << ","
    //				<< solution::f_calls << ","
    //				<< sol.flag << "\n";
    //		}
    //		else if (a_val == 100)
    //		{
    //			results_file_a100
    //				<< a_val << ","
    //				<< w << ","
    //				<< x0(0, 0) << ","
    //				<< x0(1, 0) << ","
    //				<< sol.x(0, 0) << ","
    //				<< sol.x(1, 0) << ","
    //				<< f1_val << ","
    //				<< f2_val << ","
    //				<< solution::f_calls << ","
    //				<< sol.flag << "\n";
    //		}

    //	}
    //}

    // results_file_a1.close();
    // results_file_a10.close();
    // results_file_a100.close();
    // std::cout << "Wyniki zapisane w pliku results_powell_test.csv\n";

    // PROBLEM RZECZYWISTY
    //  Tworzymy generator:
    std::random_device rd;
    std::mt19937 gen(rd());
    // Rozkłady jednostajne (uniform) w określonych przedziałach:
    //   - l ~ U(0.2, 1.0)
    //   - d ~ U(0.01, 0.05)
    std::uniform_real_distribution<double> dist_l(0.2, 1.0);
    std::uniform_real_distribution<double> dist_d(0.01, 0.05);
    // Otwieramy plik CSV, żeby zapisywać wyniki
    std::ofstream csv("wyniki_ff5R.csv");
    if (!csv.is_open())
    {
        std::cerr << "Nie udalo sie otworzyc pliku wyniki_ff5R.csv\n";
    }
    // Nagłówek:
    csv << "l(0) [mm],d(0) [mm],l* [mm],d* [mm],masa* [kg],ugiecie* [mm],naprężenie* [Pa],Liczba wywolan funkcji "
           "celu\n";
    // Parametry metody Powella
    double epsilon = 1e-6;
    int Nmax = 10000;
    // Pętla po w = 0..1 (krok 0.01 => 101 punktów)
    for (int i = 0; i <= 100; i++)
    {
        double w = i / 100.0;
        // Tworzymy ud1 = [w], 1x1:
        matrix ud1(1, 1);
        ud1(0, 0) = w;
        // Ud2 = [NaN], 1x1, by w ff5R(...) zwracać wektor (masa, ugięcie, naprężenie) przy wywołaniu z NaN
        matrix ud2NaN(1, 1, std::numeric_limits<double>::quiet_NaN());
        // Generujemy punkt startowy x0: (l0, d0)
        double l0 = dist_l(gen);
        double d0 = dist_d(gen);
        matrix x0(2, 1);
        x0(0, 0) = l0;
        x0(1, 0) = d0;
        // Zerujemy licznik wywołań funkcji celu
        solution::f_calls = 0;
        // Wywołujemy metodę Powella
        solution sol = Powell(ff5R, x0, epsilon, Nmax, ud1, ud2NaN);
        cout << "iter" << endl;
        // Odczytujemy rozwiązanie
        double l_star = sol.x(0, 0);
        double d_star = sol.x(1, 0);
        std::cout << sol << std::endl;
        // Aby uzyskać (masa, ugięcie, naprężenie), ponownie wywołujemy ff5R(...) z ud2=NaN

        double mass = sol.y(0, 0);       // Masa już obliczona przez Powell
        double deflection = sol.y(1, 0); // Masa już obliczona
        double stress = sol.y(2, 0);

        // Skalarna wartość funkcji celu
        double f_val = sol.y(0, 0);
        long calls = solution::f_calls;

        l0 = l0 * 1000;
        d0 = d0 * 1000;
        l_star = l_star * 1000;
        d_star = d_star * 1000;
        deflection = deflection * 1000;

        // Zapis do CSV
        csv << l0 << "," // l(0)
            << d0 << "," // d(0)
            << l_star << "," << d_star << "," << mass << "," << deflection << "," << stress << "," << calls << "\n";
        // Możesz też coś wypisać na ekran
    }
    csv.close();
    std::cout << "\nZakonczono iteracje. Wyniki w pliku: wyniki_ff5R.csv\n";

    // Przygotowanie pliku do zapisu wyników
    std::ofstream results_file_a1("output/lab5/results_powell_test_a1.csv");
    std::ofstream results_file_a10("output/lab5/results_powell_test_a10.csv");
    std::ofstream results_file_a100("output/lab5/results_powell_test_a100.csv");
    results_file_a1 << "a,w,x1_0,x2_0,sol_x1,sol_x2,f1,f2,f_calls,flag\n";
    results_file_a10 << "a,w,x1_0,x2_0,sol_x1,sol_x2,f1,f2,f_calls,flag\n";
    results_file_a100 << "a,w,x1_0,x2_0,sol_x1,sol_x2,f1,f2,f_calls,flag\n";

    // Generator liczb losowych dla punktu startowego
    // std::random_device rd;
    // std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-10.0, 10.0);

    // Parametry
    epsilon = 1e-8;
    Nmax = 20000;
    std::vector<double> a_values = {1.0, 10.0, 100.0};

    // Pętla po wagach w
    for (int i = 0; i <= 100; i++)
    {
        // Wagi (w) od 0 do 1 z krokiem 0.01
        double w = i * 0.01;

        // Losowy punkt startowy x0
        matrix x0(2, 1);
        x0(0, 0) = dis(gen);
        x0(1, 0) = dis(gen);

        // Pętla po wartościach a
        for (auto a_val : a_values)
        {
            // Ud1: tu przechowamy wagę w (ud1(0)) i parametr a (ud1(1))
            // Zgodnie z założeniami z transkryptu:
            // ud1(0) = waga, ud1(1) = a
            matrix ud1(2, 1);
            ud1(0) = w;
            ud1(1) = a_val;

            // ud2 pusty: ud2 będzie ustawiany w trakcie line search przez algorytm
            matrix ud2;

            // Czyścimy liczniki wywołań funkcji celu itp.
            solution::clear_calls();

            // Wywołujemy metodę Powella
            solution sol = Powell(ff5T, x0, epsilon, Nmax, ud1, ud2);

            // Utwórz macierz 1x1 wypełnioną NaN, aby uzyskać f1 i f2
            matrix ud2_empty(1, 1, std::numeric_limits<double>::quiet_NaN());
            matrix y_val = ff5T(sol.x, ud1, ud2_empty);

            double f1_val = y_val(0, 0);
            double f2_val = y_val(1, 0);

            if (a_val == 1)
            {
                save_to_file_lab_5(results_file_a1, a_val, w, x0, sol, f1_val, f2_val, solution::f_calls);
            }
            else if (a_val == 10)
            {
                save_to_file_lab_5(results_file_a10, a_val, w, x0, sol, f1_val, f2_val, solution::f_calls);
            }
            else if (a_val == 100)
            {
                save_to_file_lab_5(results_file_a100, a_val, w, x0, sol, f1_val, f2_val, solution::f_calls);
            }
        }
    }

    results_file_a1.close();
    results_file_a10.close();
    results_file_a100.close();
    std::cout << "Wyniki zapisane w pliku results_powell_test.csv\n";

    ////PROBLEM RZECZYWISTY
    //// Tworzymy generator:
    // std::random_device rd;
    // std::mt19937 gen(rd());
    //// Rozkłady jednostajne (uniform) w określonych przedziałach:
    ////   - l ~ U(0.2, 1.0)
    ////   - d ~ U(0.01, 0.05)
    // std::uniform_real_distribution<double> dist_l(0.2, 1.0);
    // std::uniform_real_distribution<double> dist_d(0.01, 0.05);
    //// Otwieramy plik CSV, żeby zapisywać wyniki
    // std::ofstream csv("wyniki_ff5R.csv");
    // if (!csv.is_open())
    //{
    //	std::cerr << "Nie udalo sie otworzyc pliku wyniki_ff5R.csv\n";
    // }
    //// Nagłówek:
    // csv << "l(0) [mm],d(0) [mm],l* [mm],d* [mm],masa* [kg],ugiecie* [mm],naprężenie* [Pa],Liczba wywolan funkcji
    // celu\n";
    //// Parametry metody Powella
    // double epsilon = 1e-6;
    // int Nmax = 10000;
    //// Pętla po w = 0..1 (krok 0.01 => 101 punktów)
    // for (int i = 0; i <= 100; i++)
    //{
    //	double w = i / 100.0;
    //	// Tworzymy ud1 = [w], 1x1:
    //	matrix ud1(1, 1);
    //	ud1(0, 0) = w;
    //	// Ud2 = [NaN], 1x1, by w ff5R(...) zwracać wektor (masa, ugięcie, naprężenie) przy wywołaniu z NaN
    //	matrix ud2NaN(1, 1, std::numeric_limits<double>::quiet_NaN());
    //	// Generujemy punkt startowy x0: (l0, d0)
    //	double l0 = dist_l(gen);
    //	double d0 = dist_d(gen);
    //	matrix x0(2, 1);
    //	x0(0, 0) = l0;
    //	x0(1, 0) = d0;
    //	// Zerujemy licznik wywołań funkcji celu
    //	solution::f_calls = 0;
    //	// Wywołujemy metodę Powella
    //	solution sol = Powell(ff5R, x0, epsilon, Nmax, ud1, ud2NaN);

    //	std::cout << sol << std::endl;
    //	// Odczytujemy rozwiązanie
    //	double l_star = sol.x(0, 0);
    //	double d_star = sol.x(1, 0);
    //	double mass = sol.y(0, 0);
    //	double deflection = sol.y(1, 0);
    //	double stress = sol.y(2, 0);
    //	long calls = solution::f_calls;

    //	l0 *=  1000;
    //	d0 *= 1000;
    //	l_star *= 1000;
    //	d_star *= 1000;
    //	deflection *= 1000;

    //	// Zapis do CSV
    //	csv << l0 << ","   // l(0)
    //		<< d0 << ","   // d(0)
    //		<< l_star << ","
    //		<< d_star << ","
    //		<< mass << ","
    //		<< deflection << ","
    //		<< stress << ","
    //		<< calls << "\n";
    //	// Możesz też coś wypisać na ekran
    //}
    // csv.close();
    // std::cout << "\nZakonczono iteracje. Wyniki w pliku: wyniki_ff5R.csv\n";
}

// Funkcja do wczytania danych z pliku
matrix loadExperimentalData(const std::string &filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Nie można otworzyć pliku: " << filename << std::endl;
        return matrix(); // Zwraca pustą macierz
    }

    std::vector<double> x1_vals;
    std::vector<double> x2_vals;
    std::string line;

    while (std::getline(file, line))
    {
        std::replace(line.begin(), line.end(), ',', '.'); // Zamiana ',' na '.'
        std::replace(line.begin(), line.end(), ';', ' '); // Zamiana ';' na ' '

        std::istringstream iss(line);
        double x1, x2;
        if (iss >> x1 >> x2)
        {
            x1_vals.push_back(x1);
            x2_vals.push_back(x2);
        }
    }

    file.close();

    // Tworzenie macierzy na podstawie wczytanych danych
    int rows = x1_vals.size();
    matrix data_matrix(rows, 2); // Tworzymy macierz o wymiarach rows x 2

    for (int i = 0; i < rows; ++i)
    {
        data_matrix(i, 0) = x1_vals[i]; // Pierwsza kolumna
        data_matrix(i, 1) = x2_vals[i]; // Druga kolumna
    }

    return data_matrix;
}

void lab6()
{
    int N = 2;
    matrix lb(N, 1), ub(N, 1);
    for (int i = 0; i < N; i++)
    {
        lb(i, 0) = -10.0;
        ub(i, 0) = 10.0;
    }
    // Parametry algorytmu
    int mi = 5;
    int lambda = 10;
    double eps = 1e-2;
    int Nmax = 10000;

    vector<double> sigmaValues = {0.01, 0.1, 1.0, 10.0, 100.0};
    matrix ud1, ud2; // puste w tym przykładzie
    const int NUM_RUNS = 100;

    // Tabela do zbierania wyników do statystyk
    // stats[i] będzie przechowywać 100 rozwiązań dla danej sigma.
    // Można też od razu uśredniać, ale tu zbieramy wszystkie, by ewentualnie przejrzeć.
    vector<vector<solution>> stats(sigmaValues.size());

    for (size_t sIdx = 0; sIdx < sigmaValues.size(); sIdx++)
    {
        double sVal = sigmaValues[sIdx];

        cout << "\nPoczatkowa wartosc zakresu mutacji: " << sVal << "\n";
        cout << "Lp., x1*, x2*, y*, Liczba wywolan funkcji celu, Minimum globalne[tak/nie]\n";

        for (int run = 1; run <= NUM_RUNS; run++)
        {
            // Ustawiamy poczatkowe sigma (2x1)
            matrix sigma0(N, 1);
            for (int d = 0; d < N; d++)
            {
                sigma0(d, 0) = sVal;
            }

            solution::clear_calls();

            solution bestSol = EA(ff6T, N, lb, ub, mi, lambda, sigma0, eps, Nmax, ud1, ud2);

            double x1 = bestSol.x(0, 0);
            double x2 = bestSol.x(1, 0);
            double fVal = bestSol.y(0, 0);
            int calls = solution::f_calls;
            bool isMinGlobal = isGlobalMinimum(bestSol, 1e-6);

            cout << run << ", " << x1 << ", " << x2 << ", " << fVal << ", " << calls << ", "
                 << (isMinGlobal ? "tak" : "nie") << "\n";

            stats[sIdx].push_back(bestSol);
        }
    }

    // ///////////////////////////////////////////////// //
    // -------------- PROBLEM RZECZYWISTY -------------- //
    // ///////////////////////////////////////////////// //

    N = 2; // Liczba zmiennych decyzyjnych: b1 i b2
    matrix lb_real(N, 1), ub_real(N, 1);
    lb_real(0, 0) = 0.1; // Dolna granica dla b1
    lb_real(1, 0) = 0.1; // Dolna granica dla b2
    ub_real(0, 0) = 3.0; // Górna granica dla b1
    ub_real(1, 0) = 3.0; // Górna granica dla b2

    // Parametry algorytmu
    mi = 5;
    lambda = 10;
    sigmaValues = {0.1};
    eps = 1e-1; // 1e-1 fast calculations and ok results, 1e-2 better result but 3 mins calculation
    Nmax = 10000;

    cout << "\nOptymalizacja problemu rzeczywistego:\n";
    // OPTYMALIZACJA
    matrix sigma0_real(N, 1, sigmaValues[0]);

    solution::clear_calls();

    matrix ud1_real;
    matrix ud2_real = loadExperimentalData("polozenia.txt"); // Experimantal data

    solution bestSolReal = EA(ff6R, N, lb_real, ub_real, mi, lambda, sigma0_real, eps, Nmax, ud1_real, ud2_real);

    double b1 = bestSolReal.x(0, 0);
    double b2 = bestSolReal.x(1, 0);
    double J = bestSolReal.y(0, 0);
    int calls = solution::f_calls;
    cout << "b1*, b2*, J(b*), Liczba wywolan funkcji celu\n";
    cout << b1 << ", " << b2 << ", " << J << ", " << calls << "\n";
}
