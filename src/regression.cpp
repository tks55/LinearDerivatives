#include "regression.hpp"
#include "linalg.hpp"

/*Parameters: const Matrix& input_values, const Vector& target_vector; Desc: Fits a regression line, minimizing
the total squared error between true target values and target values predicted by the given input values, via the
usage of normal equations (A^T * Ax = A^T * b) to solve for the coefficents of a hypothesized/user-defined line of
best fit. Returns the coefficents for each predictor variable.*/
Vector Regression::LinearRegression(const Matrix& input_values, const Vector& target_vector) {
    Matrix norm_matrix = input_values.T() * input_values;
    Vector coefs = LinAlg::Solve(norm_matrix, input_values.T() * target_vector);
    return (coefs);
}