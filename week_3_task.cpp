#include <fstream>
#include <iostream>
#include <memory>
#define _USE_MATH_DEFINES
#include <cmath>

//function calculation (считает правильно, сверено с калькулятором) 
double func(double x) {
    return (pow(sin(x / 3.0), 3)) * (atan(x));
}

//функция расчета узлов 
void knots(double a, double b, double *x) {
    const int n = 10;
    double h = (b - a) / n;
    for (int i = 0; i < n; i++) {
        x[i] = a + h * i;
    }
}

//полином лагранжа (done somehow)
double Lagrange(int n, double x, double* x_array, double* f_array) {
    double sum = 0.0;
    double product = 1.0;

    double *p_sum = &sum;
    double* p_product = &product;

    for (int i = 0; i < n; i++) {
        *p_product = 1.0; 

        for (int j = 0; j < n; j++) {
            if (j != i) {
                *p_product = product * (x - x_array[j]) / (x_array[i] - x_array[j]);
            }
        }

        *p_sum += f_array[i] * (*p_product);
    }
    return sum;
}

// полином ньютона 
// calculating diggeerence 
double div_dif(double* x_array, double* f_array, int n) {
    double sum = 0.0;
    double product = 1.0;
    double* p_product = &product;
    double* p_sum = &sum;

    for (int i = 0; i <= n; i++) {
        *p_product = 1;
        for (int j = 0; j <= n; j++) {
            if (j != i) { //znamenatel
                *p_product = product * (x_array[i] - x_array[j]);
            }
        }
        *p_sum = sum + f_array[i] / (product);
    }
    //std::cout << sum << "\n";
    return sum;
}

//Calc Newton's polinom 
double Newton(double *x_array, double *f_array, double x, int n ) {
    
    double sum = div_dif(x_array, f_array, n-1);
    double *p_sum = &sum;
    for (int i = n-2; i>=0; i--) {
        *p_sum = sum * (x - x_array[i]) + div_dif(x_array, f_array, i);
    }
    return sum;
}


// узлы чебышева
double Chebyshev_knots(double a, double b, int n, int i, double PI) {
    return ((a + b) / 2 + (b - a) / 2 * cos((2 * i + 1) * PI / (2 * n)));
}

int main()
{
    double const PI = acos(-1.0);
    double a = 0.0;
    double b = 10.0;
    const int n = 11;
    //std::cout << "type in the number of nodes: \n";
    //std::cin >> n;
    
    double h = (b - a) / n;
    
    double x[n];
    double f[n];
    double Lag[n];
    double New[n];

    //x array init and f array init (common)
    for (int i = 0; i < n; i++) {
        x[i]  = h * i;
        f[i] = func(x[i]);
    }
    std::cout << "Common knots" << "\n";
    std::cout << "x " << "f(x) " << "Lagrange " << "Newton" << "\n";
    for (int i = 0; i < n; i++) {
        Lag[i] = Lagrange(n,x[i], x, f);
        New[i] = Newton(x, f, x[i], n);
        std::cout << x[i] << " " << f[i] << " " << Lag[i] <<" "<< New[i] << "\n";
    }

    double x_ch[n];
    double f_ch[n];
    double Lag_ch[n];
    double New_ch[n];

    //x array init and f array init (chebyshev)
    for (int i = 0; i < n; i++) {
        x_ch[i] = Chebyshev_knots(a, b, n, i, PI);
        f_ch[i] = func(x_ch[i]);
    }
    std::cout << "Chebyshev knots" << "\n";
    std::cout << "x " << "f(x) " << "Lagrange " << "Newton" << "\n";
    for (int i = 0; i < n; i++) {
        Lag_ch[i] = Lagrange(n, x_ch[i], x_ch, f_ch);
        New_ch[i] = Newton(x_ch, f_ch, x_ch[i], n);
        std::cout << x_ch[i] << " " << f_ch[i] << " " << Lag_ch[i] << " " << New_ch[i] << "\n";
    }

    //additional knotes for the plot
    const int n_add = 100;
    double h_add = (b - a) / n_add;
    double x_add[n_add];
    double f_add[n_add];
    double L_add[n_add];
    double N_add[n_add];

    const char csv_file_name1[64] = "data_fixed.csv";
    std::ofstream csv_file;
    csv_file.open(csv_file_name1);
    csv_file << "x,f,Lagrange,Newton\n";

    for (int i = 0; i < n_add; i++) {
        x_add[i] = a + h_add * i;
        f_add[i] = func(x_add[i]);
    }
    for (int i = 0; i < n_add; i++) {
        L_add[i] = Lagrange(n, x_add[i], x, f);
        N_add[i] = Newton(x, f, x_add[i], n);
        csv_file << x_add[i] << "," << f_add[i] << "," << L_add[i] << "," << N_add[i] << "\n";
    }
    csv_file.close();

    double L_ch_add[n_add];
    double N_ch_add[n_add];

    const char csv_file_name2[64] = "data_Chebyshev.csv";
    csv_file.open(csv_file_name2);
    csv_file << "x,f,Lagrange,Newton\n";


    for (int i = 0; i < n_add; i++) {
        L_ch_add[i] = Lagrange(n, x_add[i], x_ch, f_ch);
        N_ch_add[i] = Newton(x_ch, f_ch, x_add[i], n);
        csv_file << x_add[i] << "," << f_add[i] << "," << L_ch_add[i] << "," << N_ch_add[i] << "\n";
    }
    csv_file.close();

    return 0;
}