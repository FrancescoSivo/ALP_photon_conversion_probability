#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>
#include <math.h>
#include <complex>
using namespace std;

//define the type real vector of long double
typedef vector<long double> vec_real;
//define the type real matrix of long double
typedef vector<vec_real> mat_real;
//define the type complex vector of long double
typedef vector<complex<long double>> vec_complex;
//define the type complex matrix of long double
typedef vector<vec_complex> mat_complex;

//function to print a real vector
void print_vec(vec_real &v){
    //! print a real vector
    //@param v: the vector to be printed
    //@return none
    cout<<"[";
    for(unsigned int i=0; i<v.size(); i++){
        cout << v[i] << " , ";
    }
    cout << "]" << endl;
}
//function to print a real matrix
void print_mat(mat_real &m){
    //! print a real matrix
    //@param m: the matrix to be printed
    //@return none
    cout<<"|";
    for(unsigned int i=0; i<m.size(); i++){
        for(unsigned int j=0; j<m[i].size(); j++){
            cout << m[i][j] << " | ";
        }
        cout << endl;
        cout << "|";
    }
}
//function to print a complex vector
void print_vec(vec_complex &v){
    //! print a complex vector
    //@param v: the vector to be printed
    //@return none
    cout<<"[";
    for(unsigned int i=0; i<v.size(); i++){
        cout << v[i] << " , ";
    }
    cout << "]" << endl;
}
//function to print a complex matrix
void print_mat(mat_complex &m){
    //! print a complex matrix
    //@param m: the matrix to be printed
    //@return none
    cout<<"|";
    for(unsigned int i=0; i<m.size(); i++){
        for(unsigned int j=0; j<m[i].size(); j++){
            cout << m[i][j] << " | ";
        }
        cout << endl;
        cout << "|";
    }
}
//function that initializes a real matrix of size n x m initialized with 0
mat_real init_mat_real(int n, int m){
    //! initialize a real matrix of size n x m initialized with 0
    //@param n: number of rows
    //@param m: number of columns
    //@return: a real matrix of size n x m initialized with 0
    mat_real mat(n);
    for(int i = 0; i < n; i++){
        mat[i].resize(m);
        for(int j = 0; j < m; j++){
            mat[i][j] = 0.0;
        }
    }
    return mat;
}
//function that initializes a complex matrix of size n x m initialized with 0
mat_complex init_mat_complex(int n, int m){
    //! initialize a complex matrix of size n x m initialized with 0
    //@param n: number of rows
    //@param m: number of columns
    //@return: a complex matrix of size n x m initialized with 0
    mat_complex mat(n);
    for(int i = 0; i < n; i++){
        mat[i].resize(m);
        for(int j = 0; j < m; j++){
            mat[i][j] = 0.0;
        }
    }
    return mat;
}
//function that takes a real matrix and returns a real vector containing the diagonal elements of the matrix and checks if the matrix is square
vec_real diag_mat(mat_real mat){
    //! takes a real matrix and returns a real vector containing the diagonal elements of the matrix and checks if the matrix is square
    //@param mat: a real matrix
    //@return: a real vector containing the diagonal elements of the matrix and checks if the matrix is square
    int n = mat.size();
    int m = mat[0].size();
    if(n != m){
        cout << "The matrix is not square" << endl;
        exit(1);
    }
    vec_real diag(n);
    for(int i = 0; i < n; i++){
        diag[i] = mat[i][i];
    }
    return diag;
}
//function that takes a complex matrix and returns a complex vector containing the diagonal elements of the matrix and checks if the matrix is square
vec_complex diag_mat(mat_complex mat){
    //! takes a complex matrix and returns a complex vector containing the diagonal elements of the matrix and checks if the matrix is square
    //@param mat: a complex matrix
    //@return: a complex vector containing the diagonal elements of the matrix and checks if the matrix is square
    int n = mat.size();
    int m = mat[0].size();
    if(n != m){
        cout << "The matrix is not square" << endl;
        exit(1);
    }
    vec_complex diag(n);
    for(int i = 0; i < n; i++){
        diag[i] = mat[i][i];
    }
    return diag;
}
//funtion that takes a real matrix and returns a real matrix which is the transpose of the input matrix
mat_real transp_mat(mat_real mat){
    //! takes a real matrix and returns a real matrix which is the transpose of the input matrix
    //@param mat: a real matrix
    //@return: a real matrix which is the transpose of the input matrix
    int n = mat.size();
    int m = mat[0].size();
    mat_real transp(m);
    for(int i = 0; i < m; i++){
        transp[i].resize(n);
        for(int j = 0; j < n; j++){
            transp[i][j] = mat[j][i];
        }
    }
    return transp;
}
//funtion that takes a complex matrix and returns a complex matrix which is the transpose of the input matrix
mat_complex transp_mat(mat_complex mat){
    //! takes a complex matrix and returns a complex matrix which is the transpose of the input matrix
    //@param mat: a complex matrix
    //@return: a complex matrix which is the transpose of the input matrix
    int n = mat.size();
    int m = mat[0].size();
    mat_complex transp(m);
    for(int i = 0; i < m; i++){
        transp[i].resize(n);
        for(int j = 0; j < n; j++){
            transp[i][j] = mat[j][i];
        }
    }
    return transp;
}
//function that returns the trace of a real matrix and checks if the matrix is square
long double trace_mat(mat_real mat){
    //! returns the trace of a real matrix and checks if the matrix is square
    //@param mat: a real matrix
    //@return: the trace of a real matrix and checks if the matrix is square
    int n = mat.size();
    int m = mat[0].size();
    if(n != m){
        cout << "The matrix is not square" << endl;
        exit(1);
    }
    long double trace = 0;
    for(int i = 0; i < n; i++){
        trace += mat[i][i];
    }
    return trace;
}
//function that returns the trace of a complex matrix and checks if the matrix is square
complex<long double> trace_mat(mat_complex mat){
    //! returns the trace of a complex matrix and checks if the matrix is square
    //@param mat: a complex matrix
    //@return: the trace of a complex matrix and checks if the matrix is square
    int n = mat.size();
    int m = mat[0].size();
    if(n != m){
        cout << "The matrix is not square" << endl;
        exit(1);
    }
    complex<long double> trace = 0;
    for(int i = 0; i < n; i++){
        trace += mat[i][i];
    }
    return trace;
}
//function that sum a real matrix with a real matrix and returns the result
mat_real sum_mat(mat_real mat1, mat_real mat2){
    //! sum a real matrix with a real matrix and returns the result
    //@param mat1: a real matrix
    //@param mat2: a real matrix
    //@return: the result of the sum of the two matrices
    int n = mat1.size();
    int m = mat1[0].size();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            mat1[i][j] = mat1[i][j] + mat2[i][j];
        }
    }
    return mat1;
}
//function that sum a complex matrix with a complex matrix and returns the result
mat_complex sum_mat(mat_complex mat1, mat_complex mat2){
    //! sum a complex matrix with a complex matrix and returns the result
    //@param mat1: a complex matrix
    //@param mat2: a complex matrix
    //@return: the result of the sum of the two matrices
    int n = mat1.size();
    int m = mat1[0].size();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            mat1[i][j] = mat1[i][j] + mat2[i][j];
        }
    }
    return mat1;
}
//function that sum a real matrix with a complex matrix and returns the result
mat_complex sum_mat(mat_real mat1, mat_complex mat2){
    //! sum a real matrix with a complex matrix and returns the result
    //@param mat1: a real matrix
    //@param mat2: a complex matrix
    //@return: the result of the sum of the two matrices
    int n = mat1.size();
    int m = mat1[0].size();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            mat2[i][j] = mat1[i][j] + mat2[i][j];
        }
    }
    return mat2;
}
//function that sum a complex matrix with a real matrix and returns the result
mat_complex sum_mat(mat_complex mat1, mat_real mat2){
    //! sum a complex matrix with a real matrix and returns the result
    //@param mat1: a complex matrix
    //@param mat2: a real matrix
    //@return: the result of the sum of the two matrices
    int n = mat1.size();
    int m = mat1[0].size();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            mat1[i][j] = mat1[i][j] + mat2[i][j];
        }
    }
    return mat1;
}
//function that subtracts a real matrix with a real matrix and returns the result
mat_real sub_mat(mat_real mat1, mat_real mat2){
    //! subtracts a real matrix with a real matrix and returns the result
    //@param mat1: a real matrix
    //@param mat2: a real matrix
    //@return: the result of the subtraction of the two matrices
    int n = mat1.size();
    int m = mat1[0].size();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            mat1[i][j] = mat1[i][j] - mat2[i][j];
        }
    }
    return mat1;
}
//function that subtracts a complex matrix with a complex matrix and returns the result
mat_complex sub_mat(mat_complex mat1, mat_complex mat2){
    //! subtracts a complex matrix with a complex matrix and returns the result
    //@param mat1: a complex matrix
    //@param mat2: a complex matrix
    //@return: the result of the subtraction of the two matrices
    int n = mat1.size();
    int m = mat1[0].size();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            mat1[i][j] = mat1[i][j] - mat2[i][j];
        }
    }
    return mat1;
}
//function that subtracts a real matrix with a complex matrix and returns the result
mat_complex sub_mat(mat_real mat1, mat_complex mat2){
    //! subtracts a real matrix with a complex matrix and returns the result
    //@param mat1: a real matrix
    //@param mat2: a complex matrix
    //@return: the result of the subtraction of the two matrices
    int n = mat1.size();
    int m = mat1[0].size();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            mat2[i][j] = mat1[i][j] - mat2[i][j];
        }
    }
    return mat2;
}
//function that subtracts a complex matrix with a real matrix and returns the result
mat_complex sub_mat(mat_complex mat1, mat_real mat2){
    //! subtracts a complex matrix with a real matrix and returns the result
    //@param mat1: a complex matrix
    //@param mat2: a real matrix
    //@return: the result of the subtraction of the two matrices
    int n = mat1.size();
    int m = mat1[0].size();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            mat1[i][j] = mat1[i][j] - mat2[i][j];
        }
    }
    return mat1;
}
//function that multiplies a real matrix with a real scalar and returns the result
mat_real mult_mat_scalar(mat_real mat, long double scalar){
    //! multiplies a real matrix with a real scalar and returns the result
    //@param mat: a real matrix
    //@param scalar: a real scalar
    //@return: the result of the multiplication of the matrix with the scalar
    int n = mat.size();
    int m = mat[0].size();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            mat[i][j] *= scalar;
        }
    }
    return mat;
}
//function that multiplies a complex matrix with a real scalar and returns the result
mat_complex mult_mat_scalar(mat_complex mat, long double scalar){
    //! multiplies a complex matrix with a real scalar and returns the result
    //@param mat: a complex matrix
    //@param scalar: a real scalar
    //@return: the result of the multiplication of the matrix with the scalar
    int n = mat.size();
    int m = mat[0].size();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            mat[i][j] *= scalar;
        }
    }
    return mat;
}
//function that multiplies a real matrix with a complex scalar and returns the result
mat_complex mult_mat_scalar(mat_real mat, complex<long double> scalar){
    //! multiplies a real matrix with a complex scalar and returns the result
    //@param mat: a real matrix
    //@param scalar: a complex scalar
    //@return: the result of the multiplication of the matrix with the scalar
    int n = mat.size();
    int m = mat[0].size();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            mat[i][j] *= scalar;
        }
    }
    return mat;
}
//function that multiplies a complex matrix with a complex scalar and returns the result
mat_complex mult_mat_scalar(mat_complex mat, complex<long double> scalar){
    //! multiplies a complex matrix with a complex scalar and returns the result
    //@param mat: a complex matrix
    //@param scalar: a complex scalar
    //@return: the result of the multiplication of the matrix with the scalar
    int n = mat.size();
    int m = mat[0].size();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            mat[i][j] *= scalar;
        }
    }
    return mat;
}
//function that multiplies a real scalar with a real matrix and returns the result
mat_real mult_mat_scalar(long double scalar, mat_real mat){
    //! multiplies a real scalar with a real matrix and returns the result
    //@param scalar: a real scalar
    //@param mat: a real matrix
    //@return: the result of the multiplication of the scalar with the matrix
    int n = mat.size();
    int m = mat[0].size();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            mat[i][j] *= scalar;
        }
    }
    return mat;
}
//function that multiplies a complex scalar with a complex matrix and returns the result
mat_complex mult_mat_scalar(complex<long double> scalar, mat_complex mat){
    //! multiplies a complex scalar with a complex matrix and returns the result
    //@param scalar: a complex scalar
    //@param mat: a complex matrix
    //@return: the result of the multiplication of the scalar with the matrix
    int n = mat.size();
    int m = mat[0].size();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            mat[i][j] *= scalar;
        }
    }
    return mat;
}
//function that multiplies a real matrix with a complex scalar and returns the result
mat_complex mult_mat_scalar(complex<long double> scalar, mat_real mat){
    //! multiplies a real matrix with a complex scalar and returns the result
    //@param scalar: a complex scalar
    //@param mat: a real matrix
    //@return: the result of the multiplication of the scalar with the matrix
    int n = mat.size();
    int m = mat[0].size();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            mat[i][j] *= scalar;
        }
    }
    return mat;
}
//function that multiplies a complex matrix with a real scalar and returns the result
mat_complex mult_mat_scalar(long double scalar, mat_complex mat){
    //! multiplies a complex matrix with a real scalar and returns the result
    //@param scalar: a real scalar
    //@param mat: a complex matrix
    //@return: the result of the multiplication of the scalar with the matrix
    int n = mat.size();
    int m = mat[0].size();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            mat[i][j] *= scalar;
        }
    }
    return mat;
}
//function that returns the product of two real matrices
mat_real mult_mat_mat(mat_real mat1, mat_real mat2){
    //! returns the product of two real matrices
    //@param mat1: a real matrix
    //@param mat2: a real matrix
    //@return: the product of the two matrices
    int n = mat1.size();
    int m = mat2[0].size();
    int p = mat2.size();
    mat_real mult(n);
    for(int i = 0; i < n; i++){
        mult[i].resize(m);
        for(int j = 0; j < m; j++){
            mult[i][j] = 0;
            for(int k = 0; k < p; k++){
                mult[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
    return mult;
}
//function that returns the product of two complex matrices
mat_complex mult_mat_mat(mat_complex mat1, mat_complex mat2){
    //! returns the product of two complex matrices
    //@param mat1: a complex matrix
    //@param mat2: a complex matrix
    //@return: the product of the two matrices
    int n = mat1.size();
    int m = mat2[0].size();
    int p = mat2.size();
    mat_complex mult(n);
    for(int i = 0; i < n; i++){
        mult[i].resize(m);
        for(int j = 0; j < m; j++){
            mult[i][j] = 0;
            for(int k = 0; k < p; k++){
                mult[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
    return mult;
}
//function that returns the product of a real matrix and a complex matrix
mat_complex mult_mat_mat(mat_real mat1, mat_complex mat2){
    //! returns the product of a real matrix and a complex matrix
    //@param mat1: a real matrix
    //@param mat2: a complex matrix
    //@return: the product of the two matrices
    int n = mat1.size();
    int m = mat2[0].size();
    int p = mat2.size();
    mat_complex mult(n);
    for(int i = 0; i < n; i++){
        mult[i].resize(m);
        for(int j = 0; j < m; j++){
            mult[i][j] = 0;
            for(int k = 0; k < p; k++){
                mult[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
    return mult;
}
//function that returns the product of a complex matrix and a real matrix
mat_complex mult_mat_mat(mat_complex mat1, mat_real mat2){
    //! returns the product of a complex matrix and a real matrix
    //@param mat1: a complex matrix
    //@param mat2: a real matrix
    //@return: the product of the two matrices
    int n = mat1.size();
    int m = mat2[0].size();
    int p = mat2.size();
    mat_complex mult(n);
    for(int i = 0; i < n; i++){
        mult[i].resize(m);
        for(int j = 0; j < m; j++){
            mult[i][j] = 0;
            for(int k = 0; k < p; k++){
                mult[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
    return mult;
}
//function that returns the matrix composed by e absolute values of the elements of a complex matrix
mat_real matrix_abs(mat_complex mat){
    //! returns the matrix composed by e absolute values of the elements of a complex matrix
    //@param mat: a complex matrix
    //@return: the matrix composed by e absolute values of the elements of the matrix
    int n = mat.size();
    int m = mat[0].size();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            mat[i][j] = abs(mat[i][j]);
        }
    }
    return mat;
}
//function that returns the sum of two real vectors
vec_real sum_vect(vec_real vect1, vec_real vect2){
    //! returns the sum of two real vectors
    //@param vect1: a real vector
    //@param vect2: a real vector
    //@return: the sum of the two vectors
    int n = vect1.size();
    for(int i = 0; i < n; i++){
        vect1[i] = vect1[i] + vect2[i];
    }
    return vect1;
}
//function that returns the sum of two complex vectors
vec_complex sum_vect(vec_complex vect1, vec_complex vect2){
    //! returns the sum of two complex vectors
    //@param vect1: a complex vector
    //@param vect2: a complex vector
    //@return: the sum of the two vectors
    int n = vect1.size();
    for(int i = 0; i < n; i++){
        vect1[i] = vect1[i] + vect2[i];
    }
    return vect1;
}
//function that returns the sum of a real vector and a complex vector
vec_complex sum_vect(vec_real vect1, vec_complex vect2){
    //! returns the sum of a real vector and a complex vector
    //@param vect1: a real vector
    //@param vect2: a complex vector
    //@return: the sum of the two vectors
    int n = vect1.size();
    for(int i = 0; i < n; i++){
        vect2[i] = vect1[i] + vect2[i];
    }
    return vect2;
}
//function that returns the sum of a complex vector and a real vector
vec_complex sum_vect(vec_complex vect1, vec_real vect2){
    //! returns the sum of a complex vector and a real vector
    //@param vect1: a complex vector
    //@param vect2: a real vector
    //@return: the sum of the two vectors
    int n = vect1.size();
    for(int i = 0; i < n; i++){
        vect1[i] = vect1[i] + vect2[i];
    }
    return vect1;
}
//function that returns the commutator of two real matrices
mat_real commutator(mat_real mat1, mat_real mat2){
    //! returns the commutator of two real matrices
    //@param mat1: a real matrix
    //@param mat2: a real matrix
    //@return: the commutator of the two matrices
    return sub_mat(mult_mat_mat(mat1, mat2),mult_mat_mat(mat2, mat1));
}
//function that returns the commutator of two complex matrices
mat_complex commutator(mat_complex mat1, mat_complex mat2){
    //! returns the commutator of two complex matrices
    //@param mat1: a complex matrix
    //@param mat2: a complex matrix
    //@return: the commutator of the two matrices
    return sub_mat(mult_mat_mat(mat1, mat2),mult_mat_mat(mat2, mat1));
}
//function that returns the commutator of a real matrix and a complex matrix
mat_complex commutator(mat_real mat1, mat_complex mat2){
    //! returns the commutator of a real matrix and a complex matrix
    //@param mat1: a real matrix
    //@param mat2: a complex matrix
    //@return: the commutator of the two matrices
    return sub_mat(mult_mat_mat(mat1, mat2),mult_mat_mat(mat2, mat1));
}
//function that returns the commutator of a complex matrix and a real matrix
mat_complex commutator(mat_complex mat1, mat_real mat2){
    //! returns the commutator of a complex matrix and a real matrix
    //@param mat1: a complex matrix
    //@param mat2: a real matrix
    //@return: the commutator of the two matrices
    return sub_mat(mult_mat_mat(mat1, mat2),mult_mat_mat(mat2, mat1));
}
//function that returns the anti-commutator of 2 real matrices
mat_real anti_commutator(mat_real mat1, mat_real mat2){
    //! returns the anti-commutator of 2 real matrices
    //@param mat1: a real matrix
    //@param mat2: a real matrix
    //@return: the anti-commutator of the two matrices
    return sub_mat(mult_mat_mat(mat1, mat2),mult_mat_mat(mat2, mat1));
}
//function that returns the anti-commutator of 2 complex matrices
mat_complex anti_commutator(mat_complex mat1, mat_complex mat2){
    //! returns the anti-commutator of 2 complex matrices
    //@param mat1: a complex matrix
    //@param mat2: a complex matrix
    //@return: the anti-commutator of the two matrices
    return sub_mat(mult_mat_mat(mat1, mat2),mult_mat_mat(mat2, mat1));
}
//function that returns the anti-commutator of a real matrix and a complex matrix
mat_complex anti_commutator(mat_real mat1, mat_complex mat2){
    //! returns the anti-commutator of a real matrix and a complex matrix
    //@param mat1: a real matrix
    //@param mat2: a complex matrix
    //@return: the anti-commutator of the two matrices
    return sub_mat(mult_mat_mat(mat1, mat2),mult_mat_mat(mat2, mat1));
}
//function that returns the anti-commutator of a complex matrix and a real matrix
mat_complex anti_commutator(mat_complex mat1, mat_real mat2){
    //! returns the anti-commutator of a complex matrix and a real matrix
    //@param mat1: a complex matrix
    //@param mat2: a real matrix
    //@return: the anti-commutator of the two matrices
    return sub_mat(mult_mat_mat(mat1, mat2),mult_mat_mat(mat2, mat1));
}
//function that returns the maximum value of a real matrix
long double max_value(mat_real mat){
    //! returns the maximum value of a real matrix
    //@param mat: a real matrix
    //@return: the maximum value of the matrix
    long double max = mat[0][0];
    for(unsigned int i = 0; i < mat.size(); i++){
        for(unsigned int j = 0; j < mat[0].size(); j++){
            if(mat[i][j] > max){
                max = mat[i][j];
            }
        }
    }
    return max;
}
//function that returns the minimum value of a real matrix
long double min_value(mat_real mat){
    //! returns the minimum value of a real matrix
    //@param mat: a real matrix
    //@return: the minimum value of the matrix
    long double min = mat[0][0];
    for(unsigned int i = 0; i < mat.size(); i++){
        for(unsigned int j = 0; j < mat[0].size(); j++){
            if(mat[i][j] < min){
                min = mat[i][j];
            }
        }
    }
    return min;
}
#endif /* MATRIX_H */
//end of matrix.h