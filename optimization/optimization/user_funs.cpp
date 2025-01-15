#include "user_funs.h"

#include "solution.h"
#include <cassert>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

matrix ff0T(matrix x, matrix ud1, matrix ud2)
{
    matrix y;
    y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);
    return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)
{
    matrix y;
    matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{m2d(x), 0.5});
    matrix *Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);
    int n = get_len(Y[0]);
    double teta_max = Y[1](0, 0);
    for (int i = 1; i < n; ++i)
        if (teta_max < Y[1](i, 0))
            teta_max = Y[1](i, 0);
    y = abs(teta_max - m2d(ud1));
    Y[0].~matrix();
    Y[1].~matrix();
    return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
    matrix dY(2, 1);
    double m = 1, l = 0.5, b = 0.5, g = 9.81;
    double I = m * pow(l, 2);
    dY(0) = Y(1);
    dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;
    return dY;
}

// LAB1
matrix df1(matrix x, matrix ud1, matrix ud2)
{
    matrix y;
    y = -cos(0.1 * x(0)) * exp(-pow(0.1 * x(0) - (2 * M_PI), 2)) + 0.002 * pow(0.1 * x(0), 2);
    return y;
}

matrix f1(double t, matrix Y, matrix ud1, matrix ud2)
{
    const double PA = 0.5; // Pole podstawy zbiornika A [m^2]
    const double TA = 90;  // Temperatura w zbiorniku A [°C]

    const double a = 0.98; // Wspó³czynnik lepkoœci cieczy
    const double b = 0.63; // Wspó³czynnik zwê¿enia strumienia
    const double g = 9.81; // Przyspieszenie ziemskie [m/s^2]

    const double PB = 1.0;        // Pole podstawy zbiornika B [m^2]
    const double DB = 0.00365665; // Pole przekroju otworu w zbiorniku B [m^2]
    const double Fin = 0.01;      // Szybkoœæ nap³ywu wody do zbiornika B [m^3/s]
    const double Tin = 20;        // Temperatura nap³ywaj¹cej wody [°C]

    matrix dY(3, 1);

    const double FAout = (Y(0) > 0) ? (a * b * m2d(ud2) * sqrt(2 * g * Y(0) / PA)) : 0;
    const double FBout = (Y(1) > 1) ? (a * b * DB * sqrt(2 * g * Y(1) / PB)) : 0;

    dY(0) = -FAout;
    dY(1) = FAout + Fin - FBout;
    dY(2) = (Fin / Y(1)) * (Tin - Y(2)) + (FAout / Y(1)) * (TA - Y(2));

    return dY;
}

matrix f2(matrix x, matrix ud1, matrix ud2)
{
    matrix y;
    matrix Y0 = matrix(3, new double[3]{5, 1, 20});
    matrix *Y = solve_ode(f1, 0, 1, 2000, Y0, ud1, x);
    int n = get_len(Y[0]);
    double max = Y[1](0, 2);
    for (size_t i = 0; i < n; i++)
    {
        if (max < Y[1](i, 2))
        {
            max = Y[1](i, 2);
        }
        y = abs(max - 50);
    }

    return y;
}

// LAB2
matrix ff2T(matrix x, matrix ud2, matrix ud1)
{
    double x1 = x(0, 0);
    double x2 = x(1, 0);
    matrix result(1, 1);
    result(0, 0) = pow(x1, 2) + pow(x2, 2) - cos(2.5 * M_PI * x1) - cos(2.5 * M_PI * x2) + 2;
    return result;
}

// Moment kontorlny
double M(double alpha, double omega, double k1, double k2)
{
    const double alphaRef = M_PI; // Docelowy kąt [π rad]
    const double omegaRef = 0.0;  // Docelowa prędkość kątowa [0 rad/s]

    double alphaDiff = alphaRef - alpha;
    double omegaDiff = omegaRef - omega;
    return k1 * alphaDiff + k2 * omegaDiff;
}

matrix df2(double t, matrix Y, matrix ud1, matrix ud2)
{
    // Parametry fizyczne
    const double l = 1.0;                         // Długość ramienia [m]
    const double mr = 1.0;                        // Masa ramienia [kg]
    const double mc = 5.0;                        // Masa ciężarka [kg]
    const double b = 0.5;                         // Współczynnik tarcia [Nms]
    const double I = ((mr / 3) + mc) * pow(l, 2); // Moment bezwładności [kg * m^2]

    // Współczynniki sterowania k1 i k2
    double k1 = m2d(ud1);
    double k2 = m2d(ud2);

    // Obecna pozycja i prędkość
    double alpha = m2d(Y(0, 0));
    double omega = m2d(Y(1, 0));

    // Obliczenie momentu kontrolnego M(t)
    double controlMomentum = M(alpha, omega, k1, k2);

    // Obliczenie pochodnych
    matrix dY(2, 1);
    dY(0, 0) = omega;
    dY(1, 0) = (controlMomentum - b * omega) / I;

    return dY;
}

matrix ff2R(matrix x, matrix ud1, matrix ud2)
{
    // Optymalne wartości współczynników sterowania k1 i k2 z wektora x
    double k1 = m2d(x(0));
    double k2 = m2d(x(1));

    // Ustawienie warunków czasowych
    double t0 = 0.0;
    double td = 0.1;
    double tend = 100.0;

    // początkowy kąt i prędkość
    matrix Y0(2, new double[2]{0.0, 0.0});

    // Symulacja ruchu ramienia
    matrix *results = solve_ode(df2, t0, td, tend, Y0, matrix(1, 1, k1), matrix(1, 1, k2));

    // Obliczenie jakości funkcji Q jako sumy kwadratów różnic od wartości
    // docelowych
    double Q = 0.0;
    for (int i = 0; i < get_len(results[0]); i++)
    {
        double alphaDiff = results[1](i, 0) - M_PI; // Różnica pozycji od wartości docelowej
        double omegaDiff = results[1](i, 1);        // Różnica prędkości od zera

        double controlMomentum = M(m2d(results[1](i, 0)), m2d(results[1](i, 1)), k1, k2);

        Q += (10 * pow(alphaDiff, 2) + pow(omegaDiff, 2) + pow(controlMomentum, 2)) * td;
    }

    delete[] results;
    return matrix(Q);
}

// LAB3
matrix ff3T(matrix x, matrix ud1, matrix ud2)
{
    double x1 = x(0, 0);
    double x2 = x(1, 0);
    double a = ud1(0, 0); // Parametr a
    double c = ud2(0, 0); // Współczynnik kary

    matrix result(1, 1);
    double denominator = M_PI * sqrt(pow(x1 / M_PI, 2) + pow(x2 / M_PI, 2));
    if (denominator == 0)
    {
        std::cerr << "ff3T - dividing by 0\n";
        result(0, 0) = 0;
    }
    else
    {
        result(0, 0) = sin(M_PI * sqrt(pow(x1 / M_PI, 2) + pow(x2 / M_PI, 2))) / denominator;
    }

    // Obliczenie kar za naruszenie ograniczeń
    double g1 = -x1 + 1;                           // g1(x1) <= 0
    double g2 = -x2 + 1;                           // g2(x2) <= 0
    double g3 = sqrt(pow(x1, 2) + pow(x2, 2)) - a; // g3(x1, x2) <= 0

    double penalty = 0.0;

    // Dodaj karę za każde naruszenie ograniczenia
    if (g1 > 0)
        penalty += pow(g1, 2);
    if (g2 > 0)
        penalty += pow(g2, 2);
    if (g3 > 0)
        penalty += pow(g3, 2);

    // Dodaj sumę kar do funkcji celu, skalowaną przez współczynnik kary c
    result(0, 0) += c * penalty;

    return result;
}

matrix df3(double t, matrix Y, matrix ud1, matrix ud2)
{
    double m = 0.6;          // masa piłki w kg
    double r = 0.12;         // promień piłki w metrach
    double g = 9.81;         // przyspieszenie grawitacyjne
    double C = 0.47;         // współczynnik oporu powietrza
    double rho = 1.2;        // gęstość powietrza
    double S = M_PI * r * r; // powierzchnia przekroju piłki

    // Wektory prędkości
    double vx = Y(1);                   // prędkość pozioma
    double vy = Y(3);                   // prędkość pionowa
    double v = sqrt(vx * vx + vy * vy); // całkowita prędkość

    // Siły
    double Dx = 0.5 * C * rho * S * v * vx;                // siła oporu pozioma
    double Dy = 0.5 * C * rho * S * v * vy;                // siła oporu pionowa
    double FMx = rho * vy * ud2(0) * 2 * M_PI * pow(r, 3); // siła Magnusa pozioma
    double FMy = rho * vx * ud2(0) * 2 * M_PI * pow(r, 3); // siła Magnusa pionowa

    // Debugowanie sił
    // cout << "Debug df3:\n";
    // cout << "vx = " << vx << ", vy = " << vy << ", v = " << v << "\n";
    // cout << "Dx = " << Dx << ", Dy = " << Dy << "\n";
    // cout << "FMx = " << FMx << ", FMy = " << FMy << "\n";

    // Równania różniczkowe
    matrix dY(4, 1);
    dY(0) = vx;                      // dx/dt
    dY(1) = (-Dx - FMx) / m;         // d^2x/dt^2
    dY(2) = vy;                      // dy/dt
    dY(3) = (-m * g - Dy - FMy) / m; // d^2y/dt^2
    return dY;
}

matrix ff3R(matrix x, matrix ud1, matrix ud2)
{
    double v0x = x(0);   // początkowa prędkość pozioma
    double omega = x(1); // początkowa rotacja

    matrix Y0(4, new double[4]{0, v0x, 100, 0}); // warunki początkowe
    matrix *Y = solve_ode(df3, 0.0, 0.01, 7.0, Y0, ud1, matrix(1, 1, omega));

    // cout << "Debug ff3R: Results from solve_ode\n";
    // cout << "Time vector:\n" << Y[0] << "\n";
    // cout << "State matrix:\n" << Y[1] << "\n";

    int n = get_len(Y[0]);
    int i0 = 0, i50 = 0;

    // Znajdowanie indeksów dla y = 0 i y = 50
    for (int i = 0; i < n; i++)
    {
        if (abs(Y[1](i, 2) - 50) < abs(Y[1](i50, 2) - 50))
        {
            i50 = i;
        }
        if (abs(Y[1](i, 2)) < abs(Y[1](i0, 2)))
        {
            i0 = i;
        }
    }

    double x_end = Y[1](i0, 0); // X na końcu
    double x50 = Y[1](i50, 0);  // X na wysokości 50

    // Dodanie debugowania trajektorii
    // cout << "Debugging trajectory for current x:\n";
    // cout << "v0x = " << v0x << ", omega = " << omega << "\n";
    // cout << "x_end = " << Y[1](i0, 0) << ", x50 = " << Y[1](i50, 0) << "\n";

    // Debuguj konkretną trajektorię dla tego przypadku
    // double debug_v0x = 4.0;
    // double debug_omega = 1.5;
    // matrix debug_Y0(4, new double[4] {0, debug_v0x, 100, 0});
    // matrix* debug_Y = solve_ode(df3, 0.0, 0.01, 7.0, debug_Y0, ud1, matrix(1,
    // 1, debug_omega)); cout << "Debug x50 for (v0x = 4.0, omega = 1.5): " <<
    // debug_Y[1](i50, 0) << "\n";

    // Debugowanie wartości
    // cout << "Debug ff3R:\n";
    // cout << "v0x = " << v0x << ", omega = " << omega << "\n";
    // cout << "x_end = " << x_end << ", x50 = " << x50 << "\n";

    // Ograniczenia
    double penalty = 0.0;
    if (abs(v0x) > 10)
    {
        penalty += ud1(0) * pow(abs(v0x) - 10, 2);
        // cout << "Penalty for v0x: " << penalty << "\n";
    }
    if (abs(omega) > 15)
    {
        penalty += ud1(0) * pow(abs(omega) - 15, 2);
        // cout << "Penalty for omega: " << penalty << "\n";
    }
    if (abs(x50 - 5) <= 0.5)
    {
        // cout << "Success: x50 = " << x50 << " is within the acceptable range
        // [4.5, 5.5].\n";
    }
    else
    {
        penalty += ud1(0) * pow(abs(x50 - 5) - 0.5, 2);
        // cout << "Warning: x50 = " << x50 << " is out of range [4.5, 5.5]. Penalty
        // applied: " << penalty << "\n";
    }

    // Wynik
    double result = -x_end + penalty;
    // cout << "Result (function value): " << result << "\n";

    Y[0].~matrix();
    Y[1].~matrix();
    return matrix(1, 1, result);
}

// LAB4
matrix ff4T(matrix x, matrix ud1, matrix ud2) // funkcja celu
{
    matrix y;
    if (isnan(ud2(0, 0)))
    {
        double x1 = x(0, 0);
        double x2 = x(1, 0);

        y = pow(x1 + (2 * x2) - 7, 2) + pow((2 * x1) + x2 - 5, 2);
    }
    else
    {
        //       | x0(0)   d0(0) |
        // ud2 = |               |
        //       | x0(1)   d0(1) |
        //
        //           x0 + h * d
        y = ff4T(ud2[0] + x * ud2[1], NAN, NAN);
    }

    return y;
}

matrix gf4T(matrix x, matrix ud1, matrix ud2) // gradient funkcji celu
{
    double x1 = x(0, 0);
    double x2 = x(1, 0);
    double u = x1 + 2 * x2 - 7;
    double v = 2 * x1 + x2 - 5;

    matrix grad(2, 1);
    grad(0, 0) = 2 * u + 4 * v; // Pochodna po x1
    grad(1, 0) = 4 * u + 2 * v; // Pochodna po x2

    return grad;
}

matrix hf4T(matrix _x, matrix _ud1, matrix _ud2) // hessian funkcji celu
{
    // Hessian na podstawie funkcji ff4T
    // Przez to, że ff4T jest funkcją kwadratową to hessian
    // jest zawsze macierzą skalarów o tych samych wartościach.

    //     | ∂²f / ∂x₁² = 10       ∂²f / ∂x₁∂x₂ = 8 |
    // H = |                                        |
    //     | ∂²f / ∂x₂∂x₁ = 8     ∂²f / ∂x₂² = 10   |

    matrix hessian(2, 2);
    hessian(0, 0) = 10.0;
    hessian(0, 1) = 8.0;
    hessian(1, 0) = 8.0;
    hessian(1, 1) = 10.0;

    return hessian;
}

double sigmoid(matrix X, matrix theta)
{
    double z = 0.0;
    int *size = get_size(X); // Pobierz rozmiar macierzy X
    int n = size[1];         // Liczba kolumn w macierzy X (liczba cech)

    // Obliczamy iloczyn skalarny X * theta
    for (int i = 0; i < n; ++i)
    {
        z += X(0, i) * theta(i, 0); // Iloczyn skalarny
    }

    // Zwracamy wynik funkcji sigmoidalnej
    return 1.0 / (1.0 + std::exp(-z)); // Sigmoid
}

matrix logisticCostFunction(matrix x, matrix ud1, matrix ud2)
{
    // Używamy get_size() do pobrania rozmiarów macierzy
    matrix cost(1, 1, 0.0); // Zainicjalizuj koszt na 0.0
    // std::cout << "DEBUG: cost " << cost << std::endl;

    for (int i = 0; i < 100; ++i)
    {
        matrix xi(3, 1);
        for (int j = 0; j < 3; ++j)
        {
            xi(j, 0) = ud1(j, i);
        }
        matrix yi(1, 1, 0.0);
        yi(0, 0) = ud2(i, 0);

        double z = 0.0;
        for (int j = 0; j < 3; ++j)
        {
            z += xi(j, 0) * x(j, 0);
        }
        // std::cout << "DEBUG: z \n" << z << std::endl;

        double sigm = 1.0 / (1.0 + std::exp(-z));

        if (sigm < 1e-15)
        {
            sigm = 1e-15;
        }
        else if (sigm > 1 - 1e-15)
        {
            sigm = 1 - 1e-15;
        }
        // std::cout << "DEBUG: yi \n" << yi(0,0) << std::endl;
        /*std::cout << "DEBUG: std::log(sigm) \n" << std::log(sigm) << std::endl;
        std::cout << "DEBUG: 1 - yi(0, 0) \n" << 1 - yi(0, 0) << std::endl;
        std::cout << "DEBUG: std::log(1 - sigm) \n" << std::log(1 - sigm) <<
        std::endl;*/
        double cost_value = (yi(0, 0) * std::log(sigm)) + ((1 - yi(0, 0)) * std::log(1 - sigm));
        // std::cout << "DEBUG: cost_value \n" << cost_value << std::endl;
        cost(0, 0) += cost_value;
    }

    cost = -cost / 100.0;

    return cost;
}

matrix computeGradient(matrix theta, matrix X, matrix Y)
{
    matrix grad(3, 1, 0.0); // Inicjalizacja gradientu (3 x 1)
    int m = 100;            // Liczba przykładów (100)

    // Iteracja po wszystkich próbkach (po wierszach w X)
    for (int i = 0; i < m; ++i)
    {
        matrix xi(3, 1, 0.0); // Wektor cech i-tego przykładu (3 x 1)

        // Pobranie cech i-tego przykładu (wiersz z X)
        for (int j = 0; j < 3; ++j)
        {                       // j - indeks cechy
            xi(j, 0) = X(j, i); // Wypełniamy xi cechami z i-tego wiersza X (X(j, i))
        }

        // Obliczanie wartości funkcji sigmoidalnej dla tej próbki
        double z = 0.0;
        for (int j = 0; j < 3; ++j)
        {                                // j - indeks cechy
            z += xi(j, 0) * theta(j, 0); // Iloczyn skalarny (z = theta^T * X)
        }

        // Obliczanie funkcji sigmoidalnej (hipotezy) h(x_i)
        double sigm = 1.0 / (1.0 + std::exp(-z)); // Sigmoid (hipoteza h(x_i))

        // Prawdziwa etykieta y_i z Y
        double yi = Y(i, 0);

        // Dodajemy składnik do gradientu dla każdej cechy
        for (int j = 0; j < 3; ++j)
        {                                // j - indeks cechy
            double sklad = sigm - yi;    // Różnica między hipotezą a etykietą
            grad(j) += sklad * xi(j, 0); // Gradient dla każdej cechy
        }
    }

    // Dzielenie przez liczbę przykładów (średni gradient)
    grad = grad / m;
    return grad;
}

double computeAccuracy(matrix X, matrix Y, matrix theta)
{
    int correct = 0;

    // Iteracja po wszystkich przykładach
    for (int i = 0; i < 100; ++i)
    {
        matrix xi(3, 1, 0.0);
        double z = 0.0;
        for (int j = 0; j < 3; ++j)
        {                                // j - indeks cechy
            z += xi(j, 0) * theta(j, 0); // Iloczyn skalarny (z = theta^T * X)
        }
        double sigm = 1.0 / (1.0 + std::exp(-z)); // Sigmoid (hipoteza h(x_i))

        int predicted = (sigm >= 0.5) ? 1 : 0;
        int actual = Y(i, 0); // Rzeczywista etykieta y_i

        // Jeśli przewidywanie jest poprawne, zwiększamy licznik
        if (predicted == actual)
        {
            correct++;
        }
    }

    return (double)correct / 100 * 100; // Procent poprawnie sklasyfikowanych przypadków
}

void writeResultsToCSV(const string &filename, const vector<vector<double>> &results)
{
    ofstream outFile(filename);

    // Nagłówki kolumn
    outFile << "Długość kroku,θ0*,θ1*,θ2*,J(θ*),P(θ*),g_calls" << endl;

    // Zapis wyników
    for (const auto &result : results)
    {
        outFile << result[0] << ","   // Długość kroku
                << result[1] << ","   // θ0
                << result[2] << ","   // θ1
                << result[3] << ","   // θ2
                << result[4] << ","   // J(θ*)
                << result[5] << ","   // P(θ*)
                << result[6] << endl; // g_calls
    }

    outFile.close();
    cout << "Wyniki zapisane do pliku: " << filename << endl;
}

string trim(const string &value)
{
    size_t start = value.find_first_not_of(" \t\n\r");
    size_t end = value.find_last_not_of(" \t\n\r");
    return (start == string::npos) ? "" : value.substr(start, end - start + 1);
}

matrix loadDataToMatrix(const string &filename, int rows, int cols)
{
    matrix m(rows, cols, 0.0); // Tworzymy macierz o wymiarach rows x cols
    ifstream file(filename);

    if (!file.is_open())
    {
        cerr << "Nie mozna otworzyc pliku: " << filename << endl;
        return m; // Zwraca pustą macierz, jeśli nie udało się otworzyć pliku
    }

    string line;
    int row = 0;

    while (getline(file, line) && row < rows)
    {
        stringstream ss(line);
        string value;
        int col = 0;
        while (getline(ss, value, ';') && col < cols)
        {
            value = trim(value); // Usuwamy ewentualne białe znaki
            try
            {
                if (!value.empty())
                {
                    m(row, col) = stod(value); // Konwertujemy na double i zapisujemy w macierzy
                    col++;
                }
            }
            catch (const std::invalid_argument &e)
            {
                cerr << "Niepoprawna liczba: " << value << endl;
            }
            catch (const std::out_of_range &e)
            {
                cerr << "Liczba poza zakresem: " << value << endl;
            }
        }
        row++;
    }

    file.close();
    return m;
}

// LAB5
//                     ud1 = [a]   ud2 = [x1, d1]
//							[w]
//[x2, d2]
matrix ff5T(matrix x, matrix ud1, matrix ud2)
{
    matrix y;
    // Sprawdzamy scenariusz:
    if (std::isnan(ud2(0, 0)))
    {
        // Poprawny odczyt a z macierzy ud1 (2x1)
        double a = ud1(1, 0); // Upewnij się, że to prawidłowy indeks

        double x1 = x(0, 0);
        double x2 = x(1, 0);

        double f1 = a * (pow(x1 - 2, 2) + pow(x2 - 2, 2));
        double f2 = (1.0 / a) * (pow(x1 + 2, 2) + pow(x2 + 2, 2));

        y = matrix(2, 1);
        y(0, 0) = f1;
        y(1, 0) = f2;
    }
    else
    {
        // Przypadek przeszukiwania linii
        double w = ud1(0, 0); // Upewnij się, że to prawidłowy indeks
        double a = ud1(1, 0); // Upewnij się, że to prawidłowy indeks
        double step = x(0, 0);

        // Odczytujemy p_start i d z ud2 (2x2)
        matrix p_start(2, 1), d(2, 1);
        p_start(0, 0) = ud2(0, 0);
        p_start(1, 0) = ud2(0, 1);

        d(0, 0) = ud2(1, 0);
        d(1, 0) = ud2(1, 1);

        matrix new_point = p_start + d * step;

        // Wywołujemy ff5T ponownie, ale teraz z ud2 = NaN (by nie wchodzić w line
        // search)
        matrix ud2_empty(2, 1, std::numeric_limits<double>::quiet_NaN());
        matrix yt = ff5T(new_point, ud1, ud2_empty);

        double f1 = yt(0, 0);
        double f2 = yt(1, 0);
        double F = w * f1 + (1.0 - w) * f2;

        y = matrix(1, 1, F);
    }
    return y;
}

void test_ff5T()
{
    const double epsilon = 1e-6; // Tolerancja dla porównań

    matrix ud1(2, 1);
    ud1(0, 0) = 0.5; // w = 0.5
    ud1(1, 0) = 1.0; // a = 1
    matrix ud2_empty(2, 1, std::numeric_limits<double>::quiet_NaN());

    // Test punktu (0, 0) dla w = 0.5
    matrix x(2, 1);
    x(0, 0) = 0.0;
    x(1, 0) = 0.0;
    matrix y = ff5T(x, ud1, ud2_empty);

    // Oczekujemy f1 = 8, f2 = 8
    bool f1_correct = std::abs(y(0, 0) - 8.0) < epsilon;
    bool f2_correct = std::abs(y(1, 0) - 8.0) < epsilon;

    if (!f1_correct || !f2_correct)
    {
        std::cerr << "Test failed:\n";
        std::cerr << "Expected f1 = 8, got " << y(0, 0) << "\n";
        std::cerr << "Expected f2 = 8, got " << y(1, 0) << "\n";
    }

    assert(f1_correct);
    assert(f2_correct);

    std::cout << "All ff5T() tests passed successfully." << std::endl;
}

matrix ff5R(matrix x, matrix ud1, matrix ud2)
{
    // Rozróżniamy dwa przypadki:

    // 1) Gdy ud2(0,0) jest NaN => zwracamy [masa; ugiecie; naprezenie]
    if (std::isnan(ud2(0, 0)))
    {
        // Tutaj x ma wymiar (2x1): x(0,0) = l, x(1,0) = d
        matrix y(3, 1);

        double ro = 7800.0; // kg/m^3
        double P = 1e3;     // N
        double E = 207e9;   // Pa

        double l = x(0, 0); // [m]
        double d = x(1, 0); // [m]

        // Masa:
        double masa = ro * l * M_PI * std::pow(d, 2) / 4.0;

        // Ugięcie:
        double ugiecie = 64.0 * P * std::pow(l, 3) / (3.0 * E * M_PI * std::pow(d, 4));

        // Naprężenie:
        double sigma = 32.0 * P * l / (M_PI * std::pow(d, 3));

        y(0, 0) = masa;
        y(1, 0) = ugiecie;
        y(2, 0) = sigma;

        return y;
    }
    else
    {
        // 2) Gdy ud2(0,0) != NaN => obliczamy skalar funkcji celu (1x1)
        //    do zminimalizowania w line search.

        // Uwaga: w line search x(0,0) to krok h
        //        w ud2 mamy p_start i direction

        // Odczytujemy wage w:
        double w = ud1(0, 0);

        // Odczytujemy h = x(0,0)
        double h = x(0, 0);

        // Odczytujemy p_start i direction z ud2 (2x2)
        // wiersz 0 => p_start = (p_start_l, p_start_d)
        // wiersz 1 => direction = (dir_l, dir_d)
        double p0_l = ud2(0, 0);
        double p0_d = ud2(0, 1);

        double dir_l = ud2(1, 0);
        double dir_d = ud2(1, 1);

        // Wyliczamy aktualny punkt xt:
        double l = p0_l + h * dir_l;
        double d = p0_d + h * dir_d;

        // Teraz obliczamy masę i ugięcie jak w objective function:
        double ro = 7800.0;
        double P = 1e3;
        double E = 207e9;

        double masa = ro * l * M_PI * std::pow(d, 2) / 4.0;
        double ugiecie = 64.0 * P * std::pow(l, 3) / (3.0 * E * M_PI * std::pow(d, 4));
        double sigma = 32.0 * P * l / (M_PI * std::pow(d, 3));

        // Funkcja celu = w*masa + (1-w)*ugiecie
        double f = w * masa + (1.0 - w) * ugiecie;

        // Dodajemy kary:
        double c = 1e20; // duza kara

        // Ograniczenia:
        if (l < 0.2)
            f += c * std::pow(0.2 - l, 2);
        if (l > 1.0)
            f += c * std::pow(l - 1.0, 2);

        if (d < 0.01)
            f += c * std::pow(0.01 - d, 2);
        if (d > 0.05)
            f += c * std::pow(d - 0.05, 2);

        if (ugiecie > 0.005)
        {
            f += c * std::pow(ugiecie - 0.005, 2);
        }
        if (sigma > 300e6)
        {
            f += c * std::pow(sigma - 300e6, 2);
        }

        // Zwracamy (1x1) macierz z f
        matrix f_mat(1, 1);
        f_mat(0, 0) = f;

        return f_mat;
    }
}

matrix ff6T(matrix x, matrix ud2, matrix ud1)
{
    double x1 = x(0, 0);
    double x2 = x(1, 0);
    matrix result(1, 1);
    result(0, 0) = pow(x1, 2) + pow(x2, 2) - cos(2.5 * M_PI * x1) - cos(2.5 * M_PI * x2) + 2;
    return result;
}


