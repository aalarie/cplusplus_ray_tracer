//
//  matrix.h
//
//  Author: Antoine Alarie
//
//  Library providing matrix arithmetic
//

#ifndef _matrix_h_
#define _matrix_h_

class Matrix
{
private:
    int num_rows, num_cols;
public:
    double** member;
    
    Matrix(){};
    Matrix(int rows, int cols);
    void free_members();
    int rows();
    int cols();
    void set_rows(int rows);
    void set_cols(int cols);
    double norm();
    void print();
    double dot(Matrix A);
    void normalize();
    void transpose();
    Matrix operator+(Matrix A);
    Matrix operator-(Matrix A);
    Matrix operator-();
    Matrix operator*(Matrix A);
    Matrix operator*(double a);
    bool operator==(Matrix A);
};

Matrix operator*(double a, Matrix& matrix);
Matrix cross_product(Matrix A, Matrix B);
double determinant(Matrix A);
Matrix cofactor(Matrix A);
Matrix mat_inverse(Matrix A);

Matrix mat_rotation_x(double theta);
Matrix mat_rotation_y(double theta);
Matrix mat_rotation_z(double theta);
Matrix mat_identity(int size);

#endif
