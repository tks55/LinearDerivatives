#include "linalg.hpp"
#include <stdexcept>
#include <cmath>

/*Parameters: const Matrix& a, const Vector& b; Desc: Solves the matrix equation, Ax = b, for
the value of x, given an invertible matrix.*/
Vector LinAlg::Solve(const Matrix& a, const Vector& b) {
    Matrix a_inv = a.Inv();
    Vector result = a_inv * b;
    return (result);
}

/*Parameters: const Matrix& matrix; Desc: Returns the PLU Decomposition of the given matrix, via the
usage of the PLUDecomp function in the Matrix class.*/
std::tuple<Matrix, Matrix, Matrix> LinAlg::PLUDecomp(const Matrix& matrix) {
    std::tuple<Matrix, Matrix, Matrix> plu_decomp = matrix.PLUDecomp();
    return (plu_decomp);
}

/*Parameters: const Matrix& matrix, size_t max_iter; Desc: Returns the eigenvector corresponding to absolute
largest eigenvalue, if one exists, via the method of power iteration.*/
std::pair<double, Vector> LinAlg::PowerIter(const Matrix& matrix, size_t max_iter) {
    std::pair<size_t, size_t> dims = matrix.GetDims();
    if (dims.first != dims.second) {
        throw std::invalid_argument("CANNOT ITERATE OVER NON-SQUARE MATRIX!");
    }
    std::pair<double, Vector> power_pair;
    Vector iter_vector = Vector(dims.first, 1.0);
    iter_vector = iter_vector.Norm();
    Vector last_vector;
    size_t current_iter = 0;
    while (current_iter < max_iter) {
        last_vector = iter_vector;
        iter_vector = matrix * iter_vector;
        Vector normal_iter_vector = iter_vector.Norm();
        if (normal_iter_vector == last_vector) {
            for (size_t row = 0; row < dims.first; row++) {
                if (std::abs(last_vector.GetEntry(row)) > Vector::GetPrecisionDelta()) {
                    double eigenvalue = (iter_vector.GetEntry(row) / last_vector.GetEntry(row));
                    power_pair.first = eigenvalue;
                    power_pair.second = iter_vector.Norm();
                    return (power_pair);
                }
            }
            break;
        }
        iter_vector = normal_iter_vector;
        current_iter++;
    }
    throw std::invalid_argument("UNABLE TO DETERMINE DOMINANT EIGENVECTOR/VALUE!");
    return (power_pair);
}

/*Parameters: const Matrix& matrix; Desc: Converts a Matrix into an array of (column) Vectors,
returning the array.*/
Vector* LinAlg::ToColumnVectors(const Matrix& matrix) {
    std::pair<size_t, size_t> dims = matrix.GetDims();
    Vector* column_vectors = new Vector[dims.second];
    for (size_t curr_col = 0; curr_col < dims.second; curr_col++) {
        double* curr_vect = new double[dims.first];
        for (size_t curr_row = 0; curr_row < dims.first; curr_row++) {
            curr_vect[curr_row] = matrix.GetEntry(curr_row, curr_col);
        }
        column_vectors[curr_col] = Vector(dims.first, curr_vect);
    }
    return (column_vectors);
}

/*Parameters: Vector* column_vectors, size_t num_vectors; Desc: Converts an array of same-sized
column Vectors of size num_vectors to a Matrix, returning the Matrix.*/
Matrix LinAlg::ToMatrix(Vector* column_vectors, size_t num_vectors) {
    if (num_vectors <= 0) {
        throw std::invalid_argument("UNABLE TO CONVERT TO MATRIX--NON-POSITIVE NUMBER OF COLUMN VECTORS!");
    }
    size_t num_rows = column_vectors[0].size();
    Matrix return_matrix = Matrix(num_rows, num_vectors);
    for (size_t curr_col = 0; curr_col < num_vectors; curr_col++) {
        Vector curr_vect = column_vectors[curr_col];
        if (curr_vect.size() != num_rows) {
            throw std::invalid_argument("UNABLE TO CONVERT TO MATRIX--INVALID COLUMN VECTORS/NUMBER OF COLUMN VECTORS PASSED!");
        }
        for (size_t curr_row = 0; curr_row < num_rows; curr_row++) {
            return_matrix.ChangeEntry(curr_row, curr_col, curr_vect.GetEntry(curr_row));
        }
    }
    delete[] column_vectors;
    return (return_matrix);
}

/*Parameters: const Matrix& matrix, std::function<double(double)> function; Desc: Applies a given function to all of the values
of an input Matrix, returning the modified/transformed Matrix.*/
Matrix LinAlg::ApplyFunction(const Matrix& matrix, std::function<double(double)> function) {
    Matrix output_matrix = Matrix(matrix);
    std::pair<size_t, size_t> dims = matrix.GetDims();
    for (size_t curr_row = 0; curr_row < dims.first; curr_row++) {
        for (size_t curr_col = 0; curr_col < dims.second; curr_col++) {
            double curr_val = output_matrix.GetEntry(curr_row, curr_col);
            double modified_val = function(curr_val);
            output_matrix.ChangeEntry(curr_row, curr_col, modified_val);
        }
    }
    return (output_matrix);
}