#include <cstdlib>
#include <iostream>
#include "Matrix.h"
#include <random>
#include <fstream>
#include <chrono>
#include <vector>
#include <string>
#include <eigen3/Eigen/Eigen>
#include <iomanip>

double generate_values(unsigned int i, unsigned int n);
template <typename MultiplyFunc>
void benchmark_method(std::ofstream& csv_file, unsigned int n, Matrix* a, Matrix* b, MultiplyFunc multiply_func);
void run_matrix_benchmarks(unsigned int max_n);
void verify_with_eigen(unsigned int n);

int main() {
    srand(4);
    std::cout << "Test dla macierzy 4x4:\n";

    std::cout << "\n-------------------------------------------------\n\n";

    std::cout << "Test dla macierzy 17x17 (złożony peeling!):\n";
    verify_with_eigen(100);

    return 0;
}

double generate_values(unsigned int i, unsigned int n) {
    return (double)(1e-8
        + (1 - 1e-8) * rand() / RAND_MAX);
}

void verify_with_eigen(unsigned int n) {

    Matrix A(n, n);
    Matrix B(n, n);
    A.fillMatrix(generate_values);
    B.fillMatrix(generate_values);
    Eigen::MatrixXd eigen_A(n, n);
    Eigen::MatrixXd eigen_B(n, n);
    for (unsigned int i = 0; i < n; ++i) {
        for (unsigned int j = 0; j < n; ++j) {
            eigen_A(i, j) = A(i, j);
            eigen_B(i, j) = B(i, j);
        }
    }


    Eigen::MatrixXd eigen_C = eigen_A * eigen_B;


    Matrix C_binet    = Matrix::binet_multiply(&A, &B);
    Matrix C_strassen = Matrix::strassen_multiply(&A, &B);
    Matrix C_laderman = Matrix::laderman_multiply(&A, &B);


    auto calculate_max_error = [&](const Matrix& custom_C) {
        double max_err = 0.0;
        for (unsigned int i = 0; i < n; ++i) {
            for (unsigned int j = 0; j < n; ++j) {
                double diff = std::abs(custom_C(i, j) - eigen_C(i, j));
                if (diff > max_err) {
                    max_err = diff;
                }
            }
        }
        return max_err;
    };

    double err_binet    = calculate_max_error(C_binet);
    double err_strassen = calculate_max_error(C_strassen);
    double err_laderman = calculate_max_error(C_laderman);


    std::cout << "--- RAPORT BLEDOW (Max Absolute Error vs Eigen) ---\n";
    std::cout << "Rozmiar macierzy: " << n << " x " << n << "\n\n";


    std::cout << std::scientific << std::setprecision(6);
    std::cout << "Binet:    " << err_binet << "\n";
    std::cout << "Strassen: " << err_strassen << "\n";
    std::cout << "Laderman: " << err_laderman << "\n\n";

    double tolerance = 1e-9;
    double max_all = std::max({err_binet, err_strassen, err_laderman});

    if (max_all > tolerance) {
        std::cout << "OSTRZEZENIE: Przynajmniej jeden algorytm wykazuje znaczacy blad numeryczny!\n";
    } else {
        std::cout << "SUKCES: Wszystkie algorytmy licza idealnie (w granicach bledu zmiennoprzecinkowego).\n";
    }
}


template <typename MultiplyFunc>
void benchmark_method(std::ofstream& csv_file, unsigned int n, Matrix* a, Matrix* b, MultiplyFunc multiply_func) {
    const int NUM_RUNS = 10;
    double total_time_ms = 0.0;
    long long flops_count = 0;

    for (int i = 0; i < NUM_RUNS; ++i) {
        auto start = std::chrono::high_resolution_clock::now();

        Matrix c = multiply_func(a, b);

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = end - start;
        total_time_ms += elapsed.count();


        if (i == NUM_RUNS - 1) {
            flops_count = c.flops;
        }
    }

    double avg_time_ms = total_time_ms / NUM_RUNS;

    csv_file << n << "," << avg_time_ms << "," << flops_count << "\n";
}


void run_matrix_benchmarks(unsigned int max_n) {
    std::ofstream file_binet("binet_results.csv");
    std::ofstream file_strassen("strassen_results.csv");
    std::ofstream file_laderman("laderman_results.csv");

    if (!file_binet.is_open() || !file_strassen.is_open() || !file_laderman.is_open()) {
        std::cerr << "Blad przy otwieraniu plikow CSV!" << std::endl;
        return;
    }


    std::string header = "Size,AvgTime_ms,Flops\n";
    file_binet << header;
    file_strassen << header;
    file_laderman << header;

    std::cout << "Rozpoczynam testy macierzy od N=10 do N=" << max_n << "...\n";

    for (unsigned int n = 1000; n <= max_n; ++n) {
        if (n % 10 == 0) {
            std::cout << "Testowanie wymiaru: " << n << "x" << n << "..." << std::endl << std::flush;
        }

        Matrix a(n, n);
        Matrix b(n, n);
        a.fillMatrix(generate_values);
        b.fillMatrix(generate_values);

        benchmark_method(file_binet, n, &a, &b, Matrix::binet_multiply);
        benchmark_method(file_strassen, n, &a, &b, Matrix::strassen_multiply);
        benchmark_method(file_laderman, n, &a, &b, Matrix::laderman_multiply);
    }

    file_binet.close();
    file_strassen.close();
    file_laderman.close();

    std::cout << "\nTesty zakonczone sukcesem! Wyniki zapisano do plikow .csv\n";
}