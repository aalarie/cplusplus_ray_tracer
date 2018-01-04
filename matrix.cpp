//
//  matrix.cpp
//
//  Author: Antoine Alarie
//
//  Library providing matrix arithmetic
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "matrix.h"

////////////Matrix class definition

//  Constructs a rows x cols matrix
//
Matrix::Matrix(int rows, int cols)
{
    num_rows = rows;
    num_cols = cols;
    
    member = (double**) malloc(num_rows * sizeof(double*));
    for (int i = 0; i < rows; ++i)
        member[i] = (double*) malloc(num_cols * sizeof(double));
}

void Matrix::free_members()
{
    for (int i = 0; i < num_rows; ++i)
        free(member[i]);
    free(member);
}

//  Returns the number of rows in this matrix
//
int Matrix::rows()
{
    return num_rows;
}

//  Returns the number of columns in this matrix
//
int Matrix::cols()
{
    return num_cols;
}

//  Set number of rows
//
void Matrix::set_rows(int rows)
{
    int diff = rows - num_rows;
    num_rows = rows;
    
    member = (double**) realloc(member, num_rows * sizeof(double*));
    
    if (diff > 0)
    {
        for (int i = num_rows - diff; i < num_rows; ++i)
            member[i] = (double*) malloc(num_cols * sizeof(double));
    }
}

//  Set number of columns
//
void Matrix::set_cols(int cols)
{
    num_cols = cols;
    
    for (int i = 0; i < num_rows; ++i)
        member[i] = (double*) realloc(member[i], num_cols * sizeof(double));
}

//  Prints the contents of this matrix to the screen
//
void Matrix::print()
{
    for (int i = 0; i < num_rows; ++i)
    {
        printf("|\t");
        for (int j = 0; j < num_cols; ++j)
        {
            if (j+1 == num_cols)
                printf("%.2f", member[i][j]);
            else
                printf("%.2f\t\t\t", member[i][j]);
        }
        printf("\t|\n");
    }
}

/////Operations

//  Adds matrix A to this matrix
//
Matrix Matrix::operator+(Matrix A)
{
    if (num_rows != A.rows() || num_cols != A.cols())
    {
        printf("Matrix Error: Incompatible Matrix Sizes\n");
        exit(EXIT_FAILURE);
    }
    
    Matrix C = Matrix(num_rows, num_cols);
    
    for (int i = 0; i < num_rows; ++i)
        for (int j = 0; j < num_cols; ++j)
            C.member[i][j] = member[i][j] + A.member[i][j];
    
    return C;
}

// Subtracts this matrix by matrix A
//
Matrix Matrix::operator-(Matrix A)
{
    if (num_rows != A.rows() || num_cols != A.cols())
    {
        printf("Matrix Error: Incompatible Matrix Sizes\n");
        exit(EXIT_FAILURE);
    }
    
    Matrix C = Matrix(num_rows, num_cols);
    
    for (int i = 0; i < num_rows; ++i)
        for (int j = 0; j < num_cols; ++j)
            C.member[i][j] = member[i][j] - A.member[i][j];
    
    return C;
}

Matrix Matrix::operator-()
{
    Matrix C = Matrix(num_rows, num_cols);
    
    for (int i = 0; i < num_rows; ++i)
        for (int j = 0; j < num_cols; ++j)
            C.member[i][j] = -1 * member[i][j];
    
    return C;
}

Matrix Matrix::operator*(Matrix A)
{
    if (num_cols != A.rows())
    {
        printf("Matrix Error: Incompatible Matrix Sizes\n");
        exit(EXIT_FAILURE);
    }
    
    Matrix C = Matrix(num_rows, A.cols());
    double sum = 0;
    
    for (int i = 0; i < C.rows(); ++i) {
        for (int j = 0; j < C.cols(); ++j) {
            for (int k = 0; k < num_cols; ++k)
                sum += member[i][k] * A.member[k][j];
            C.member[i][j] = sum;
            sum = 0;
        }
    }
    
    return C;
}

Matrix Matrix::operator*(double a)
{
    Matrix C = Matrix(num_rows, num_cols);
    
    for (int i = 0; i < num_rows; ++i)
        for (int j = 0; j < num_cols; ++j)
            C.member[i][j] = a * member[i][j];
    
    return C;
}

bool Matrix::operator==(Matrix A)
{
    if (num_cols != A.rows())
        return false;
    
    for (int i = 0; i < num_rows; ++i) {
        for (int j = 0; j < num_cols; ++j) {
            if (A.member[i][j] != member[i][j])
                return false;
        }
    }
    
    return true;
}


double Matrix::dot(Matrix A)
{
    if (num_rows != A.rows() && num_cols != A.cols())
    {
        printf("Matrix Error: Incompatible Matrix Sizes\n");
        exit(EXIT_FAILURE);
    }
    
    double product = 0;
    
    for (int i = 0; i < num_rows; ++i)
        for (int j = 0; j < num_cols; ++j)
            product += member[i][j] * A.member[i][j];
    
    return product;
}

//  Returns the euclidean norm of this matrix if is a vector
//  Otherwise returns -1
//
double Matrix::norm()
{
    double sum = 0;
    
    if (num_rows == 1)
    {
        for (int i = 0; i < num_cols; ++i)
            sum += pow(member[0][i], 2);
    }
    else if (num_cols == 1)
    {
        for (int i = 0; i < num_rows; ++i)
            sum += pow(member[i][0], 2);
    }
    else
        return -1;
    
    return sqrt(sum);
}

//  Normalizes this matrix if it is a vector
//  Otherwise prints a warning and does nothing
//
void Matrix::normalize()
{
    double n;
    n = norm();
    
    if (n == -1)
    {
        printf("Matrix Warning: Not a Vector\n");
        return;
    }
    
    for (int i = 0; i < num_rows; ++i)
        for (int j = 0; j < num_cols; ++j)
            member[i][j] /= n ;
    
}

//  Transposes this matrix
//
void Matrix::transpose()
{
    Matrix B = Matrix(num_cols, num_rows);
    
    for (int i = 0; i < num_rows; ++i)
        for (int j = 0; j < num_cols; ++j)
            B.member[j][i] = member[i][j] ;
    
    free_members();
    member=B.member;
    num_rows = B.num_rows;
    num_cols = B.num_cols;
}

////////////Matrix operations definition

// Returns the result of the cross product A x B
//
Matrix cross_product(Matrix A, Matrix B)
{
    Matrix product;
    
    if (((A.rows() == 1 && A.cols() == 3) || (A.rows() == 3 && A.cols() == 1)) &&
        ((B.rows() == 1 && B.cols() == 3) || (B.rows() == 3 && B.cols() == 1)))
    {
        product = Matrix(3,1);
        
        if (A.rows() == 1 && A.cols() == 3 && B.rows() == 1 && B.cols() == 3)
        {
            product.member[0][0] = A.member[0][1] * B.member[0][2] - A.member[0][2] * B.member[0][1];
            product.member[1][0] = -(A.member[0][0] * B.member[0][2] - A.member[0][2] * B.member[0][0]);
            product.member[2][0] = A.member[0][0] * B.member[0][1] - A.member[0][1] * B.member[0][0];
        }
        else if (A.rows() == 1 && A.cols() == 3)
        {
            product.member[0][0] = A.member[0][1] * B.member[2][0] - A.member[0][2] * B.member[1][0];
            product.member[1][0] = -(A.member[0][0] * B.member[2][0] - A.member[0][2] * B.member[0][0]);
            product.member[2][0] = A.member[0][0] * B.member[1][0] - A.member[0][1] * B.member[0][0];
        }
        else if (B.rows() == 1 && B.cols() == 3)
        {
            product.member[0][0] = A.member[1][0] * B.member[0][2] - A.member[2][0] * B.member[0][1];
            product.member[1][0] = -(A.member[0][0] * B.member[0][2] - A.member[2][0] * B.member[0][0]);
            product.member[2][0] = A.member[0][0] * B.member[0][1] - A.member[1][0] * B.member[0][0];
        }
        else
        {
            product.member[0][0] = A.member[1][0] * B.member[2][0] - A.member[2][0] * B.member[1][0];
            product.member[1][0] = -(A.member[0][0] * B.member[2][0] - A.member[2][0] * B.member[0][0]);
            product.member[2][0] = A.member[0][0] * B.member[1][0] - A.member[1][0] * B.member[0][0];
        }
    }
    else
    {
        printf("Matrix Error: Incompatible Matrix Sizes\n");
        exit(EXIT_FAILURE);
    }
    
    return product;
}


double determinant(Matrix A)
{
    if (A.rows() != A.cols())
    {
        printf("Matrix Error: Not a Square Matrix\n");
        exit(EXIT_FAILURE);
    }
    
    if (A.rows() == 1)
        return A.member[0][0];
    
    if (A.rows() == 2)
        return A.member[0][0] * A.member[1][1] - A.member[1][0] * A.member[0][1];
    
    double det = 0;
    int j2;
    Matrix M;
    for (int p = 0; p < A.rows(); ++p)
    {
        int h = 0;
        int k = 0;
        M = Matrix(A.rows()-1, A.cols()-1);
        for (int i = 1; i < A.rows(); ++i)
        {
            j2 = 1;
            for (int j = 0; j < A.rows(); ++j)
            {
                if (j == p)
                    continue;
                M.member[h][k] = A.member[i][j];
                k++;
                if (k == A.rows()-1)
                {
                    h++;
                    k = 0;
                }
            }
        }
        det += pow(-1, p) * A.member[0][p] * determinant(M);
        M.free_members();
    }
    
    return det;
}


Matrix cofactor(Matrix A)
{
    double det;
    Matrix B = Matrix(A.rows(), A.cols());
    Matrix C = Matrix(A.rows()-1, A.cols()-1);
    
    for (int j = 0; j < A.rows(); ++j) {
        for (int i = 0; i < A.rows(); ++i) {
            int i1 = 0;
            for (int ii = 0; ii < A.rows(); ++ii) {
                if (ii == i)
                    continue;
                int j1 = 0;
                for (int jj = 0; jj < A.rows(); ++jj) {
                    if (jj == j)
                        continue;
                    C.member[i1][j1] = A.member[ii][jj];
                    j1++;
                }
                i1++;
            }
            det = determinant(C);
            B.member[i][j] = pow(-1, i+j+2) * det;
        }
    }
    
    C.free_members();
    return B;
}


Matrix mat_inverse(Matrix A)
{
    double det;
    Matrix B = Matrix(A.rows(), A.cols());
    
    det = determinant(A);
    B = cofactor(A);
    B.transpose();
    for (int i = 0; i < B.rows(); ++i)
        for (int j = 0; j < B.cols(); ++j)
            B.member[i][j] /= det;
    
    return B;
}

Matrix operator*(double a, Matrix& B)
{
    Matrix C = Matrix(B.rows(), B.cols());
    
    for (int i = 0; i < B.rows(); ++i)
        for (int j = 0; j < B.cols(); ++j)
            C.member[i][j] = a * B.member[i][j];
    
    return C;
}

//////Rotation Matrices
Matrix mat_rotation_x(double theta)
{
    Matrix R = Matrix(3,3);
    theta = theta * M_PI/180; //Convert to radians
    R.member[0][0] = 1;    R.member[0][1] = 0;             R.member[0][2] = 0;
    R.member[1][0] = 0;    R.member[1][1] = cos(theta);    R.member[1][2] = -sin(theta);
    R.member[2][0] = 0;    R.member[2][1] = sin(theta);    R.member[2][2] = cos(theta);
    
    return R;
}


Matrix mat_rotation_y(double theta)
{
    Matrix R = Matrix(3,3);
    theta = theta * M_PI/180; //Convert to radians
    R.member[0][0] = cos(theta);    R.member[0][1] = 0;    R.member[0][2] = sin(theta);
    R.member[1][0] = 0;             R.member[1][1] = 1;    R.member[1][2] = 0;
    R.member[2][0] = -sin(theta);   R.member[2][1] = 0;    R.member[2][2] = cos(theta);
    
    return R;
}


Matrix mat_rotation_z(double theta)
{
    Matrix R = Matrix(3,3);
    theta = theta * M_PI/180; //Convert to radians
    R.member[0][0] = cos(theta);    R.member[0][1] = -sin(theta);   R.member[0][2] = 0;
    R.member[1][0] = sin(theta);    R.member[1][1] = cos(theta);    R.member[1][2] = 0;
    R.member[2][0] = 0;             R.member[2][1] = 0;             R.member[2][2] = 1;
    
    return R;
}

////Identity Matrix
Matrix mat_identity(int size)
{
    Matrix A = Matrix(size, size);
    
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
        {
            if (i == j)
                A.member[i][j] = 1;
            else
                A.member[i][j] = 0;
        }
    
    return A;
}

/*int main()
{
    Matrix a = Matrix(3,1);
    Matrix b = Matrix(3,1);
    double result;
    
    a.member[0][0] = 5;
    a.member[1][0] = 1;
    a.member[2][0] = 4;
    
    b.member[0][0] = 1; //b.member[0][1] = 2; b.member[0][2] = 0;
    b.member[1][0] = -1; //b.member[1][1] = 1; b.member[1][2] = 1;
    b.member[2][0] = 1; //b.member[2][1] = 2; b.member[2][2] = 3;
    
    result = a.dot(b);
    
    //a.print(); printf("\n");
    b.print(); printf("\n");
    printf("%.2f\n", result);
}*/

