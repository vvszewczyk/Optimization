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
		lab3();
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

	ofstream expToFile("./expansion.csv");
	ofstream fibToFile("./fibonacci.csv");
	ofstream lagToFile("./lagrange.csv");

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
	ofstream simFib("Symulacja_fibonacci.csv");
	ofstream simLag("Symulacja_lagrange.csv");

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
	std::ofstream fileHJ("lab2/Hooke-Jeeves.csv");
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
	std::ofstream HookeSymulacja("lab2/HookeSymulacja.csv");
	std::ofstream RosenbrockSymulacja("lab2/RosenbrockSymulacja.csv");
	std::ofstream Symulacja("lab2/Symulacja.csv");

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

}

void lab4()
{

}

void lab5()
{

}
  
void lab6()
{

}
