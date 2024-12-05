﻿/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

#include"opt_alg.h"
#include <random>

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


int main()
{
	try
	{
		lab4();
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
	//Funkcja testowa
	double epsilon = 1e-2;
	int Nmax = 10000;
	matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);
	solution opt;
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Wahadlo
	Nmax = 1000;
	epsilon = 1e-2;
	lb = 0;
	ub = 5;
	double teta_opt = 1;
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(opt.x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
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

	double* expansionResults = new double[2];
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
		expToFile << x0 << delimiter << expansionResults[0] << delimiter << expansionResults[1] << delimiter << solution::f_calls << endl;
		solution::clear_calls();

		solution fibonacci = fib(df1, expansionResults[0], expansionResults[1], epsilon);
		double x_fib = m2d(fibonacci.x);
		double y_fib = m2d(fibonacci.y);
		std::string minimum_type_fib = (x_fib >= 62.72 && x_fib <= 62.73) ? "GLOBALNE" : "LOKALNE";
		fibToFile << x_fib << delimiter << y_fib << delimiter << solution::f_calls << delimiter << minimum_type_fib << "\n";
		solution::clear_calls();

		solution lagrange = lag(df1, expansionResults[0], expansionResults[1], epsilon, gamma, Nmax);
		double x_lag = m2d(lagrange.x);
		double y_lag = m2d(lagrange.y);
		std::string minimum_type_lag = (x_lag >= 62.71 && x_lag <= 62.74) ? "GLOBALNE" : "LOKALNE";
		std::string minimum_type_lag2 = (lagrange.flag == -1) ? "ERROR" : minimum_type_lag;
		lagToFile << x_lag << delimiter << y_lag << delimiter << solution::f_calls << delimiter << minimum_type_lag2 << "\n";
		solution::clear_calls();

	}

	// Wykres + ostatni rekord z tabeli 1
	solution fibonacci2 = fib(df1, -100, 100, epsilon);

	double x_fib = m2d(fibonacci2.x);
	double y_fib = m2d(fibonacci2.y);
	std::string minimum_type_fib2 = (x_fib >= 62.72 && x_fib <= 62.73) ? "GLOBALNE" : "LOKALNE";
	std::cout << x_fib << delimiter << y_fib << delimiter << solution::f_calls << delimiter << minimum_type_fib2 << "\n";
	std::cout << "Fibonacci" << std::endl;
	std::cout << fibonacci2;
	solution::clear_calls();

	solution lagrange2 = lag(df1, -100, 100, epsilon, gamma, Nmax);
	double x_lag = m2d(lagrange2.x);
	double y_lag = m2d(lagrange2.y);
	std::string minimum_type_lag2 = (x_lag >= 62.72 && x_lag <= 62.73) ? "GLOBALNE" : "LOKALNE";
	std::string minimum_type_lag3 = (lagrange2.flag == -1) ? "ERROR" : minimum_type_lag2;
	std::cout << x_lag << delimiter << y_lag << delimiter << solution::f_calls << delimiter << minimum_type_lag3 << "\n";
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

	matrix symFib = matrix(3, new double[3] {5, 1, 20});
	matrix* solved1 = solve_ode(f1, 0, 1, 2000, symFib, NAN, fibEx2.x(0));
	simFib << solved1[1] << std::endl << std::endl;
	std::cout << solved1[1] << std::endl << std::endl;

	matrix symLag = matrix(3, new double[3] {5, 1, 20});
	matrix* solved2 = solve_ode(f1, 0, 1, 2000, symLag, NAN, lagEx2.x(0));
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
	std::vector<double> steps = { 0.1, 0.5, 1.0 };
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
				<< y0HJ.x(0, 0) << delimiter << y0HJ.x(1, 0) << delimiter << m2d(y0HJ.y) << delimiter << solution::f_calls
				<<delimiter<< ((abs(m2d(y0HJ.y)) < epsilon) ? "TAK" : "NIE")<< delimiter << "\n";

			solution::clear_calls();

			// Testowanie metody Rosenbrocka i zapis do pliku Rosenbrock.csv
			matrix s0(2, 1, step);
			solution y0Rosen = Rosen(ff2T, x0, s0, alphaR, beta, epsilon, Nmax);
			fileRosen << "Rosenbrock" << delimiter << step << delimiter << x0(0, 0) << delimiter << x0(1, 0) << delimiter
				 << y0Rosen.x(0, 0) << delimiter << y0Rosen.x(1, 0) << delimiter <<m2d(y0HJ.y) <<delimiter<< solution::f_calls << delimiter
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
	double k_values[2] = { 2.5, 7.5 };
	matrix x0(2, k_values);

	solution::clear_calls();

	// Optymalizacja metodą Hooke-Jeeves
	solution HookeR = HJ(ff2R, x0, step, alphaHJ, epsilon, Nmax);
	HookeSymulacja << HookeR.x(0) << delimiter << HookeR.x(1) << delimiter
		<< m2d(HookeR.y) << delimiter << solution::f_calls << delimiter << HookeR.flag << "\n";
	solution::clear_calls();

	// Optymalizacja metodą Rosenbrocka
	solution RosenbrockR = Rosen(ff2R, x0, matrix(2, 1, step), alphaR, beta, epsilon, Nmax);
	RosenbrockSymulacja << RosenbrockR.x(0) << delimiter << RosenbrockR.x(1) << delimiter
		<< m2d(RosenbrockR.y) << delimiter << solution::f_calls << delimiter << RosenbrockR.flag << "\n";
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
	matrix y0(2, new double[2] {0.0, 0.0}); // początkowy

	// Symulacja dla wyników optymalnych z Hooke-Jeeves
	matrix* yz1 = solve_ode(df2, t0, dt, tend, y0, HookeR.x(0), HookeR.x(1));
	Symulacja << "Symulacja dla Hooke-Jeeves\n" << yz1[1] << "\n\n";
	delete[] yz1;
	solution::clear_calls();


	// Symulacja dla wyników optymalnych z Rosenbrocka
	matrix* yz2 = solve_ode(df2, t0, dt, tend, y0, RosenbrockR.x(0), RosenbrockR.x(1));
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
	double s = 1.0;      // Rozmiar kroku
	double alpha = 1.0;  // Współczynnik odbicia
	double beta = 0.5;   // Współczynnik kontrakcji
	double gamma = 2.0;  // Współczynnik ekspansji
	double delta = 0.5;  // Współczynnik zmniejszenia sympleksu
	double c = 1.0;      // Wstępna wartość kary w `pen`
	double dc = 2.0;     // Mnożnik zwiększający karę w `pen`

	// Lista wartości `a`
	std::vector<double> a_values = { 4.0, 4.4934, 5.0 };

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

	for (double a : a_values) {
		for (int i = 0; i < 100; i++) {
			matrix x0(2, 1);
			bool valid_point = false;
			while (!valid_point) {
				x0(0, 0) = dis(gen);
				x0(1, 0) = dis(gen);
				// Sprawdzenie ograniczeń
				double x1 = x0(0, 0);
				double x2 = x0(1, 0);
				double g1 = -x1 + 1;                          // g1(x1) <= 0
				double g2 = -x2 + 1;                          // g2(x2) <= 0
				double g3 = sqrt(pow(x1, 2) + pow(x2, 2)) - a; // g3(x1, x2) <= 0

				valid_point = (g1 <= 0) && (g2 <= 0) && (g3 <= 0);
			}
			// Przygotowanie `ud1` jako macierz 6x1 dla sym_NM i macierz 1x1 dla ograniczeń
			matrix ud1(7, 1);
			ud1(0, 0) = a;       // Parametr a dla ograniczeń
			ud1(1, 0) = s;       // Rozmiar kroku
			ud1(2, 0) = alpha;   // Odbicie
			ud1(3, 0) = beta;    // Kontrakcja
			ud1(4, 0) = gamma;   // Ekspansja
			ud1(5, 0) = delta;   // Zmniejszenie sympleksu
			ud1(6, 0) = epsilon;

			matrix ud2(1, 1, c); // wstępna wartość kary

			try {
				// Wywołanie funkcji `pen`
				solution::clear_calls(); // Reset licznika funkcji celu
				solution opt_pen = pen(ff3T, x0, c, dc, epsilon, Nmax, ud1, ud2);
				double r_pen = std::sqrt(opt_pen.x(0, 0) * opt_pen.x(0, 0) + opt_pen.x(1, 0) * opt_pen.x(1, 0));
				cout << solution::f_calls << endl;
				// Wywołanie funkcji `sym_NM`
				matrix ud2_2(1, 1, r_pen); // wstępna wartość kary
				solution::clear_calls(); // Reset licznika funkcji celu
				solution opt_sym = sym_NM(ff3T, x0, ud1(1, 0), ud1(2, 0), ud1(3, 0), ud1(4, 0), ud1(5, 0), ud1(6, 0), Nmax, matrix(1, 1, a), ud2_2);
				double r_sym = std::sqrt(opt_sym.x(0, 0) * opt_sym.x(0, 0) + opt_sym.x(1, 0) * opt_sym.x(1, 0));

				// Zapis wyników do jednej linii
				results_file << x0(0, 0) << "," << x0(1, 0) << ","
					<< opt_pen.x(0, 0) << "," << opt_pen.x(1, 0) << "," << r_pen << "," << opt_pen.y(0, 0) << "," << solution::f_calls << ","
					<< opt_sym.x(0, 0) << "," << opt_sym.x(1, 0) << "," << r_sym << "," << opt_sym.y(0, 0) << "," << solution::f_calls << "\n";
			}
			catch (const std::exception& e) {
				std::cerr << "Błąd podczas optymalizacji: " << e.what() << "\n";
			}
		}
	}

	results_file.close();
	std::cout << "Wyniki zapisane w pliku results_table_combined.csv\n";

	// Problem rzeczywisty
	epsilon = 1e-3;   // Dokładność
	Nmax = 1000;         // Maksymalna liczba iteracji
	c = 10.0;          // Współczynnik kary początkowej
	dc = 10.0;        // Współczynnik zwiększenia kary

	// Ograniczenia dla v0x i omega
	matrix x0(2, new double[2] {5.0, .0});     // x0 = 5.0 omega = 0.0

	// Parametry do sympleksu przekazane przez ud1
	// s = 1.0, alpha = 1.0, beta = 0.5, gamma = 2.0, delta = 0.5, epsilon = 1e-3, Nmax = 1000
	matrix ud1(7, new double[7] {1.0, 1.0, 0.5, 2.0, 0.5, epsilon, static_cast<double>(Nmax)});

	// Współczynnik kary dla funkcji celu
	matrix ud2(1, 1, 100.0);

	// Debugowanie inicjalizacji
	//cout << "Debug: Initializing ud1 and ud2...\n";
	//cout << "ud1: " << ud1 << "\n";
	//cout << "ud2: " << ud2 << "\n";

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
	matrix Y0(4, new double[4] {0, opt.x(0), 100, 0}); // Warunki początkowe
	matrix* Y = solve_ode(df3, 0.0, 0.01, 7.0, Y0, matrix(), matrix(1, 1, opt.x(1)));

	// Zapis trajektorii do pliku CSV
	ofstream file("output/lab3/optimal_trajectory.csv");
	file << hcat(Y[0], Y[1]); // Łączenie czasu i współrzędnych w jeden plik
	file.close();

	// Wyświetlenie trajektorii
	//cout << "Trajektoria zapisana w pliku 'optimal_trajectory.csv'.\n";

	// Zwolnienie pamięci
	delete[] Y;
}

void lab4()
{
	// Parametry
	double epsilon = 1e-6;
	int Nmax = 1000;
	std::vector<double> h_values = { 0.05, 0.12 };
	std::string delimiter = ",";
	// Dla ekspancji
	double* expansionResults = new double[2];
	double d = 1.0; // Kierunek
	//double alfa = 1.1;
	double alfa = 2.0; // Współczynnik ekspansji


	// Otwórz plik do zapisu wyników dla metody SD
	std::ofstream results_file_SD("output/lab4/lab4_SD_results.csv");
	results_file_SD << "Method,h0,x0(1),x0(2),x*(1),x*(2),y*,f_calls,g_calls,flag\n";

	// Otwórz plik do zapisu wyników dla metody CG
	std::ofstream results_file_CG("output/lab4/lab4_CG_results.csv");
	results_file_CG << "Method,h0,x0(1),x0(2),x*(1),x*(2),y*,f_calls,g_calls,flag\n";

	// Otwórz plik do zapisu wyników dla metody Newtona
	std::ofstream results_file_Newton("output/lab4/lab4_Newton_results.csv");
	results_file_Newton << "Method,h0,x0(1),x0(2),x*(1),x*(2),y*,f_calls,g_calls,H_calls,flag\n";

	// Otwórz plik do zapisu wyników dla metody Złotego podziału
	std::ofstream results_file_Golden("output/lab4/lab4_Golden_results.csv");
	results_file_Golden << "Method,x0,x*(0),y*,f_calls,flag\n";

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

			results_file_SD << "SD" << delimiter
				<< h0 << delimiter
				<< x0(0, 0) << delimiter << x0(1, 0) << delimiter
				<< sol_SD.x(0, 0) << delimiter << sol_SD.x(1, 0) << delimiter
				<< m2d(sol_SD.y) << delimiter
				<< solution::f_calls << delimiter << solution::g_calls << delimiter
				<< sol_SD.flag << "\n";

			// Metoda Gradientów Sprzężonych (CG)
			solution::clear_calls();
			solution sol_CG = CG(ff4T, gf4T, x0, h0, epsilon, Nmax);

			results_file_CG << "CG" << delimiter
				<< h0 << delimiter
				<< x0(0, 0) << delimiter << x0(1, 0) << delimiter
				<< sol_CG.x(0, 0) << delimiter << sol_CG.x(1, 0) << delimiter
				<< m2d(sol_CG.y) << delimiter
				<< solution::f_calls << delimiter << solution::g_calls << delimiter
				<< sol_CG.flag << "\n";

			// Metoda Newtona
			solution::clear_calls();
			solution sol_Newton = Newton(ff4T, gf4T, hf4T, x0, h0, epsilon, Nmax);

			results_file_Newton << "Newton" << delimiter
				<< h0 << delimiter
				<< x0(0, 0) << delimiter << x0(1, 0) << delimiter
				<< sol_Newton.x(0, 0) << delimiter << sol_Newton.x(1, 0) << delimiter
				<< m2d(sol_Newton.y) << delimiter
				<< solution::f_calls << delimiter << solution::g_calls << delimiter
				<< solution::H_calls << delimiter << sol_Newton.flag << "\n";

			// Złoty podział
			// Inicjalizacja ud1 i ud2
			matrix ud1;
			ud1(0, 0) = NAN;

			//       | x0(0)   d0(0) |
			// ud2 = |               |
			//       | x0(1)   d0(1) |

			matrix ud2(2, 2);
			ud2(0, 0) = x0(0);
			ud2(1, 0) = x0(1);
			ud2(0, 1) = 0.0;
			ud2(1, 1) = 0.0;

			// Na początek za pomocą metody ekspacji obliczamy przedział na którym znajduję się minimum
			// Przez to, że x1 i x2 mogą być w różnej odległości od siebie to wybieramy też losowe x0
			d = x0(1, 0) - x0(0, 0);
			expansionResults = expansion(ff4T, x0_expansion, d, alfa, Nmax, ud1, ud2); // which x0?

			solution sol_Golden = golden(ff4T, expansionResults[0], expansionResults[1], epsilon, Nmax, ud1, ud2);

			results_file_Golden << "Golden" << delimiter
				<< x0_expansion << delimiter << sol_Golden.x(0, 0)
				<< delimiter << m2d(sol_Golden.y) << delimiter
				<< solution::f_calls << delimiter << sol_Golden.flag << "\n";
		}
	}

	// Zamknięcie plików
	results_file_SD.close();
	results_file_CG.close();
	results_file_Newton.close();
	results_file_Golden.close();

	std::cout << "Wyniki zapisane w output/lab4/ do plików lab4_SD_results.csv, lab4_CG_results.csv, lab4_Newton_results.csv i lab4_Golden_results.csv\n";
}




void lab5()
{

}
  
void lab6()
{

}
