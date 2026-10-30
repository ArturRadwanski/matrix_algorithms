//
// Created by valleron on 10/29/26.
//

#include "Matrix.h"

#include <iostream>
#include <iomanip>
#include "MatrixView.h"
#include "ConstMatrixView.h"

Matrix & Matrix::operator+=(const ConstMatrixView &rhs) {
    for (unsigned int i = 0; i < rows; ++i) {
        for (unsigned int j = 0; j < columns; ++j) {
            (*this)(i, j) += rhs(i, j);
        }
    }
    return *this;
}

Matrix & Matrix::operator-=(const ConstMatrixView &rhs) {
    for (unsigned int i = 0; i < rows; ++i) {
        for (unsigned int j = 0; j < columns; ++j) {
            (*this)(i, j) -= rhs(i, j);
        }
    }
    return *this;
}


inline Matrix operator+(Matrix&& lhs, const ConstMatrixView &rhs) {
    lhs += rhs;
    return std::move(lhs);
}

inline Matrix operator-(Matrix&& lhs, const ConstMatrixView &rhs) {
    lhs -= rhs;
    return std::move(lhs);
}

double & Matrix::operator()(unsigned int r, unsigned int c) {
    return data[r * columns + c];
}

double Matrix::operator()(unsigned int r, unsigned int c) const {
    return data[r * columns + c];
}

Matrix Matrix::operator+(const Matrix &other) const {
    if (rows != other.rows || columns != other.columns) {
        throw std::invalid_argument("Wymiary macierzy do dodawania musza byc identyczne");
    }
    Matrix result(rows, columns);
    for (size_t i = 0; i < data.size(); ++i) {
        result.data[i] = data[i] + other.data[i];
    }
    return result;
}

Matrix Matrix::operator-(const Matrix &other) const {
    if (rows != other.rows || columns != other.columns) {
        throw std::invalid_argument("Wymiary macierzy do dodawania musza byc identyczne");
    }
    Matrix result(rows, columns);
    for (size_t i = 0; i < data.size(); ++i) {
        result.data[i] = data[i] - other.data[i];
    }
    return result;
}

inline Matrix operator+(Matrix&& lhs, const Matrix& rhs) {
    lhs += ConstMatrixView(rhs);
    return std::move(lhs);
}

inline Matrix operator-(Matrix&& lhs, const Matrix& rhs) {
    lhs -= ConstMatrixView(rhs);
    return std::move(lhs);
}