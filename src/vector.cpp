#include "vector.hpp"
#include "stdexcept"
#include "cmath"

Vector::Vector(size_t size, double starting_val): Matrix(size, 1, starting_val, false){}

Vector::Vector(size_t size, double* init_array): Matrix(size, 1, 0, false){
    for (size_t curr_row = 0; curr_row < size; curr_row++)  {
        this->ChangeEntry(curr_row, 0, init_array[curr_row]);
    }
}

Vector::Vector(const Matrix& rhs): Matrix(rhs) {
    if (this->GetDims().second != 1) {
        throw (std::invalid_argument("INVALID MATRIX PASSED INTO VECTOR CONVERSION!"));
    }
}

double Vector::Dot(const Vector& rhs) {
    double product = Dot(*this, rhs);
    return (product);
}

double Vector::Dot(const Vector& lhs, const Vector& rhs) {
    size_t size = lhs.GetDims().first;
    if (size != rhs.GetDims().first) {
        throw (std::invalid_argument("INVALID VECTOR DOT PRODUCT DIMENSIONS!"));
    }
    double product = 0.0;
    for (size_t curr_row = 0; curr_row < size; curr_row++) {
        product += (lhs.GetEntry(curr_row, 0) * rhs.GetEntry(curr_row, 0));
    }
    return (product);
}

Vector Vector::Norm() {
    Vector normal_vect = Norm(*this);
    return (normal_vect);
}

Vector Vector::Norm(const Vector& vect) {
    double mag = std::sqrt(Dot(vect, vect));
    if (mag == 0) {
        throw std::invalid_argument("ZERO VECTOR--UNABLE TO NORMALIZE VECTOR!");
    }
    Vector normal_vect = (vect / mag);
    return (normal_vect);
}