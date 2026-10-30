//
// Created by valleron on 10/15/26.
//
#pragma once
#ifndef ALGORYTMY_MACIERZOWE_MATRIX_H
#define ALGORYTMY_MACIERZOWE_MATRIX_H
#include <vector>
class MatrixView;
class ConstMatrixView;
class Matrix {
public:
    Matrix(unsigned int rows, unsigned int columns);
    unsigned int n_rows() const;
    unsigned int n_columns() const;
    void fillMatrix(double (*f)(unsigned int, unsigned int));
    static Matrix binet_multiply(Matrix* a, Matrix* b);
    static Matrix strassen_multiply(Matrix* a, Matrix* b);
    static Matrix laderman_multiply(Matrix* a, Matrix* b);

    double& operator()(unsigned int r, unsigned int c);
    double operator()(unsigned int r, unsigned int c) const;
    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix & other) const;
    Matrix& operator+=(const ConstMatrixView &rhs);
    Matrix& operator-=(const ConstMatrixView &rhs);
    void print();
    unsigned int flops;
private:
    std::vector<double> data;
    unsigned int rows;
    unsigned int columns;



    static void binet_recursion(ConstMatrixView a, ConstMatrixView b, MatrixView c);
    static void strassen_recursion(ConstMatrixView a, ConstMatrixView b, MatrixView c);
    static void laderman_recursion(const ConstMatrixView &a, const ConstMatrixView &b, MatrixView c);
    static Matrix compute_helpers(const ConstMatrixView &left, const ConstMatrixView &right, unsigned n);
    static Matrix compute_helpers_laderman(const ConstMatrixView &left, const ConstMatrixView &right, unsigned n);
};

#endif //ALGORYTMY_MACIERZOWE_MATRIX_H