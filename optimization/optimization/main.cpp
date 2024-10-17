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
		lab1();
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
	std::cout << "Global minumum x = 62.7, y =  " << m2d(df1(62.7)) << "\n\n";

	// Random starting point
	std::random_device rd;  // Ziarno losowe
	std::mt19937 gen(rd()); // Generator Mersenne Twister
	std::uniform_int_distribution<> dis(-100, 100); // Definiujemy zakres od -100 do 100

	double* expansionResults = new double[2];
	double d = 1.0;
	double alpha = 2.0;
	int Nmax = 10000;
	double epsilon = 1e-2;
	double gamma = 1e-200;

	for (int i = 0; i < 100; i++)
	{
		int x0 = dis(gen);
		expansionResults = expansion(df1, x0, d, alpha, Nmax);
		//std::cout << "expansionResults[0] = " << expansionResults[0] << "  expansionResults[1] = " << expansionResults[1] << "\n";
		solution::clear_calls();

		solution fibonacci = fib(df1, expansionResults[0], expansionResults[1], epsilon);
		cout << "\nfib(): " << m2d(fibonacci.x) << "\n";
		solution::clear_calls();

		solution lagrange = lag(df1, expansionResults[0], expansionResults[1], epsilon, gamma, Nmax);
		//cout << "\nSolution result flag = " << lagrange.flag << "\n";
		cout <<"x0 = " << x0 << ": m2d(lagrange.x) = " << m2d(lagrange.x) << ", m2d(lagrange.y) = " << m2d(lagrange.y) << ", solution::f_calls = " << solution::f_calls << "\n";
		solution::clear_calls();
		x0++;
	}

	delete [] expansionResults;
}

void lab2()
{

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
