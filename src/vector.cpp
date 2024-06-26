#include "vector.hpp"
#include "stdexcept"
#include "cmath"

/*Constructor (by default/starting parameters)*/
Vector::Vector(size_t size, double starting_val): Matrix(size, 1, starting_val, false){}

/*Constructor (by free store array)*/
Vector::Vector(size_t size, double* init_array): Matrix(size, 1, 0, false){
    for (size_t curr_row = 0; curr_row < size; curr_row++)  {
        this->ChangeEntry(curr_row, init_array[curr_row]);
    }
    delete[] init_array;
}

/*Cast Constructor (by matrix input)*/
Vector::Vector(const Matrix& rhs): Matrix(rhs) {
    if (this->GetDims().second != 1) {
        throw (std::invalid_argument("INVALID MATRIX PASSED INTO VECTOR CONVERSION!"));
    }
}

/*Parameters: const Vector& rhs; Desc: Computes the dot product of the current vector with another vector.*/
double Vector::Dot(const Vector& rhs) {
    double product = Dot(*this, rhs);
    return (product);
}

/*Parameters: const Vector& lhs, const Vector& rhs; Desc: Computes the dot product of two vectors.*/
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

/*Parameters: None; Desc: Returns the normalized (unit) vector of the current vector.*/
Vector Vector::Norm() {
    Vector normal_vect = Norm(*this);
    return (normal_vect);
}

/*Parameters: const Vector& vect; Desc: Returns the normalized (unit) vector of the input vector.*/
Vector Vector::Norm(const Vector& vect) {
    double mag = Magnitude(vect);
    if (mag == 0) {
        throw std::invalid_argument("ZERO VECTOR--UNABLE TO NORMALIZE VECTOR!");
    }
    Vector normal_vect = (vect / mag);
    return (normal_vect);
}

/*Parameters: const Vector& rhs; Desc: Returns the projection of the current vector onto the rhs vector.*/
Vector Vector::ProjectOnto(const Vector& rhs) {
    Vector projection_vector = ProjectOnto(*this, rhs);
    return (projection_vector);
}

/*Parameters: const Vector& lhs, const Vector& rhs; Desc: Returns the projection of the lhs vector 
onto the rhs vector.*/
Vector Vector::ProjectOnto(const Vector& lhs, const Vector& rhs) {
    double product = Dot(lhs, rhs);
    product /= Dot(rhs, rhs);
    Vector projection_vector = product * rhs;
    return (projection_vector);
}

/*Parameters: None; Desc: Returns the magnitude of the current vector.*/
double Vector::Magnitude() {
    double mag = Magnitude(*this);
    return (mag);
}

/*Parameters: const Vector& vect; Desc: Returns the magnitude of the input vector.*/
double Vector::Magnitude(const Vector& vect) {
    double mag = std::sqrt(Dot(vect, vect));
    return (mag);
}

/*Parameters: size_t seed, size_t size, double lower_bound, double upper_bound;
Desc: Creates a Random Vector of Floats of size size using a given seed. Uses the 64-Bit
implementation of the Mersenne Twister PRNG, applied to a uniform real distribution
(upper bound exclusive).*/
Vector Vector::RandomVectorFloat(size_t seed, size_t size, double lower_bound, double upper_bound) {
    Vector new_vect = RandomMatrixFloat(seed, size, 1, lower_bound, upper_bound);
    return (new_vect);
}

/*Parameters: size_t seed, size_t size, int lower_bound, int upper_bound;
Desc: Creates a Random Vector of Integers of size size using a given seed. Uses the 64-Bit
implementation of the Mersenne Twister PRNG, applied to a uniform integer distribution
(upper bound inclusive).*/
Vector Vector::RandomVectorInt(size_t seed, size_t size, int lower_bound, int upper_bound) {
    Vector new_vect = RandomMatrixInt(seed, size, 1, lower_bound, upper_bound);
    return (new_vect);
}

/*Parameters: size_t entry, double val; Desc: Changes the value of a given vector at entry
to val. Zero-Indexed.*/
void Vector::ChangeEntry(size_t entry, double val) {
    this->Matrix::ChangeEntry(entry, 0, val);
}

/*Parameters: size_t row, size_t col; Desc: Returns the value of a given vector at entry. 
Zero-Indexed.*/
double Vector::GetEntry(size_t entry) const {
    double val = this->Matrix::GetEntry(entry, 0);
    return (val);
}