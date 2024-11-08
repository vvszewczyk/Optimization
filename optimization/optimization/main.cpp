/*********************************************
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
		lab2();
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
	ofstream simFib("simulation_fibonacci.csv");
	ofstream simLag("simulation_lagrange.csv");

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

	// -------- Funkcja Testująca --------- //
	double epsilon = 1e-06;
	int Nmax = 1000;
	double alphaHJ = 0.5;
	double alphaR = 0.5;
	double beta = 0.5;
	std::vector<double> steps = { 0.1, 0.5, 1.0 };
	std::string delimiter = ";";

	// Otwarcie dwóch osobnych plików CSV dla każdej z metod
	std::ofstream fileHJ("Hooke-Jeeves.csv");
	std::ofstream fileRosen("Rosenbrock.csv");

	if (!fileHJ.is_open() || !fileRosen.is_open())
	{
		std::cerr << "Cannot open file.\n";
		return;
	}

	// Nagłówki kolumn z separatorem średnika
	fileHJ << "Metoda;Długość kroku;x0(1);x0(2);x(1);x(2);f_calls;Minimum Globalne;Sukces\n";
	fileRosen << "Metoda;Długość kroku;x0(1);x0(2);x(1);x(2);f_calls;Minimum Globalne;Sukces\n";

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
			solution optHJ = HJ(df2, x0, step, alphaHJ, epsilon, Nmax);
			fileHJ << "Hooke-Jeeves;" << step << delimiter << x0(0, 0) << delimiter << x0(1, 0) << delimiter
				<< optHJ.x(0, 0) << delimiter << optHJ.x(1, 0) << delimiter << optHJ.f_calls << delimiter
				<< optHJ.y << delimiter << (optHJ.flag == 1 ? "TAK" : "NIE") << "\n";

			solution::clear_calls();

			// Testowanie metody Rosenbrocka i zapis do pliku Rosenbrock.csv
			matrix s0(2, 1, step);
			solution optRosen = Rosen(df2, x0, s0, alphaR, beta, epsilon, Nmax);
			fileRosen << "Rosenbrock;" << step << delimiter << x0(0, 0) << delimiter << x0(1, 0) << delimiter
				<< optRosen.x(0, 0) << delimiter << optRosen.x(1, 0) << delimiter << optRosen.f_calls << delimiter
				<< optRosen.y << delimiter << (optRosen.flag == 1 ? "TAK" : "NIE") << "\n";

			solution::clear_calls(); 
		}
	}

	// Zamknięcie plików
	fileHJ.close();
	fileRosen.close();

	std::cout << "Wyniki zapisane do plików Hooke-Jeeves.csv i Rosenbrock.csv\n";

	// -------- Problem rzeczywisty --------- //
	
	double s = 0.1;
	matrix s0 = matrix(2, 1, s);

	matrix x0 = matrix(2, 1, 5);
	cout << x0 << "\n\n";

	// Pliki CSV 
	std::ofstream HookeSimulation("HookeSimulation.csv");
	std::ofstream RosenbrockSimulation("RosenbrockSimulation.csv");
	std::ofstream Simulation("Simulation.csv");

	solution HookeR = HJ(df2, x0, s, alphaHJ, epsilon, Nmax);
	int a = solution::f_calls;
	HookeSimulation << HookeR.x << ":" << HookeR.y << delimiter << a << delimiter << endl;
	solution::clear_calls();

	solution RosenbrockR = Rosen(df2, x0, s0, alphaR, beta, epsilon, Nmax);
	int b = solution::f_calls;
	RosenbrockSimulation << RosenbrockR.x << ":" << RosenbrockR.y << delimiter << b << delimiter << endl;
	solution::clear_calls();
	HookeSimulation.close();

	double t0 = 0;
	double td = 0.1;
	double tend = 100;

	// Symulacja
	matrix y0_1 = matrix(2, 1);
	matrix ud1(2, new double[2] {M_PI, 0});
	matrix* yz1 = solve_ode(Q, t0, td, tend, y0_1, ud1, HookeR.x);
	cout << yz1[1] << endl;
	Simulation << yz1[1] << endl;

	Simulation << "\n\n\n";

	matrix y0_2 = matrix(2, 1);
	matrix* yz2 = solve_ode(Q, t0, td, tend, y0_2, ud1, RosenbrockR.x);
	cout << yz2[1] << endl;
	Simulation << yz2[1] << endl;
	Simulation.close();
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
