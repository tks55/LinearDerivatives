#include "vector.hpp"
#include "stdexcept"
#include "cmath"

Vector::Vector(size_t size, double starting_val): Matrix(size, 1, starting_val, false){}

Vector::Vector(size_t size, double* init_array): Matrix(size, 1, 0, false){
    for (size_t curr_row = 0; curr_row < size; curr_row++)  {
        this->ChangeEntry(curr_row, init_array[curr_row]);
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
        product += (lhs.GetEntry(curr_row) * rhs.GetEntry(curr_row));
    }
    return (product);
}

Vector Vector::Norm() {
    Vector normal_vect = Norm(*this);
    return (normal_vect);
}

Vector Vector::Norm(const Vector& vect) {
    double mag = Magnitude(vect);
    if (mag == 0) {
        throw std::invalid_argument("ZERO VECTOR--UNABLE TO NORMALIZE VECTOR!");
    }
    Vector normal_vect = (vect / mag);
    return (normal_vect);
}

Vector Vector::ProjectOnto(const Vector& rhs) {
    Vector projection_vector = ProjectOnto(*this, rhs);
    return (projection_vector);
}

Vector Vector::ProjectOnto(const Vector& lhs, const Vector& rhs) {
    double product = Dot(lhs, rhs);
    product /= Dot(rhs, rhs);
    Vector projection_vector = product * rhs;
    return (projection_vector);
}

double Vector::Magnitude() {
    double mag = Magnitude(*this);
    return (mag);
}

double Vector::Magnitude(const Vector& vect) {
    double mag = std::sqrt(Dot(vect, vect));
    return (mag);
}

Vector Vector::RandomVectorFloat(size_t seed, size_t size, double lower_bound, double upper_bound) {
    Vector new_vect = RandomMatrixFloat(seed, size, 1, lower_bound, upper_bound);
    return (new_vect);
}

Vector Vector::RandomVectorInt(size_t seed, size_t size, int lower_bound, int upper_bound) {
    Vector new_vect = RandomMatrixInt(seed, size, 1, lower_bound, upper_bound);
    return (new_vect);
}

void Vector::ChangeEntry(size_t entry, double val) {
    this->Matrix::ChangeEntry(entry, 0, val);
}

double Vector::GetEntry(size_t entry) const {
    double val = this->Matrix::GetEntry(entry, 0);
    return (val);
}