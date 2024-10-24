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
	std::random_device rd;  // Ziarno losowe
	std::mt19937 gen(rd()); // Generator Mersenne Twister
	std::uniform_int_distribution<> dis(-100, 100); // Definiujemy zakres od -100 do 100

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
		expToFile << x0 << "," << expansionResults[0] << "," << expansionResults[1] << "," << solution::f_calls << endl;
		solution::clear_calls();

		solution fibonacci = fib(df1, expansionResults[0], expansionResults[1], epsilon);
		double x_fib = m2d(fibonacci.x);
		double y_fib = m2d(fibonacci.y);
		std::string minimum_type_fib = (x_fib >= 62.72 && x_fib <= 62.73) ? "GLOBALNE" : "LOKALNE";
		fibToFile << x_fib << "," << y_fib << "," << solution::f_calls << "," << minimum_type_fib << "\n";
		solution::clear_calls();

		solution lagrange = lag(df1, expansionResults[0], expansionResults[1], epsilon, gamma, Nmax);
		double x_lag = m2d(lagrange.x);
		double y_lag = m2d(lagrange.y);
		std::string minimum_type_lag = (x_lag >= 62.71 && x_lag <= 62.74) ? "GLOBALNE" : "LOKALNE";
		std::string minimum_type_lag2 = (lagrange.flag == -1) ? "ERROR" : minimum_type_lag;
		lagToFile << x_lag << "," << y_lag << "," << solution::f_calls << "," << minimum_type_lag2 << "\n";
		solution::clear_calls();

	}

	// WYKRES + ostatni rekord z tabeli 1
	solution fibonacci2 = fib(df1, -100, 100, epsilon);


	double x_fib = m2d(fibonacci2.x);
	double y_fib = m2d(fibonacci2.y);
	std::string minimum_type_fib2 = (x_fib >= 62.72 && x_fib <= 62.73) ? "GLOBALNE" : "LOKALNE";
	std::cout << x_fib << "," << y_fib << "," << solution::f_calls << "," << minimum_type_fib2 << "\n";
	std::cout <<"Fibonacci" << std::endl;
	std::cout << fibonacci2;
	solution::clear_calls();

	solution lagrange2 = lag(df1, -100, 100, epsilon, gamma, Nmax);
	double x_lag = m2d(lagrange2.x);
	double y_lag = m2d(lagrange2.y);
	std::string minimum_type_lag2 = (x_lag >= 62.72 && x_lag <= 62.73) ? "GLOBALNE" : "LOKALNE";
	std::string minimum_type_lag3 = (lagrange2.flag == -1) ? "ERROR" : minimum_type_lag2;
	std::cout << x_lag << "," << y_lag << "," << solution::f_calls << "," << minimum_type_lag3 << "\n";
	std::cout << "Lagrange" << std::endl;
	std::cout << lagrange2;
	solution::clear_calls();
	

	// TABELA 3

	std::cout << "\n\nTABELA 3\n\n";
	solution fibEx2 = fib(f2, 1e-4, 1e-2, epsilon);
	std::cout << "Fib" << std::endl << fibEx2 << endl;

	solution::clear_calls();

	solution lagEx2 = lag(f2, 1e-4, 1e-2, epsilon, gamma, 2000);
	std::cout << "Lag" << std::endl << lagEx2 << endl;


	// SYMULACJA

	std::cout << "\n\nSYMULACJA\n\n";
	ofstream simFib("simulation_fibonacci.csv");
	ofstream simLag("simulation_lagrange.csv");

	matrix symFib = matrix(3, new double[3]{5, 1, 20});
	matrix* solved1 = solve_ode(f1, 0, 1, 2000, symFib, NAN, fibEx2.x(0));
	simFib << solved1[1] << std::endl << std::endl;
	std::cout << solved1[1] << std::endl << std::endl;

	matrix symLag = matrix(3, new double[3] {5, 1, 20});
	matrix* solved2 = solve_ode(f1, 0, 1, 2000, symLag, NAN, lagEx2.x(0));
	simLag << solved2[1] << std::endl;
	std::cout << solved2[1] << std::endl;

	delete [] expansionResults;
}

void lab2()
{
	double epsilon = 1e-2;
	int Nmax = 100;
	matrix x0(2, 1, 0.0);  // Początkowy punkt
	double s = 0.1;        // Rozmiar kroku

	solution opt = HJ(df2, x0, s, 0.9, epsilon, Nmax);  // Testowanie optymalizacji Hooke-Jeeves
	cout << "Optimized solution: " << opt << endl;
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
