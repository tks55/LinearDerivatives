#include "vector.hpp"
#include "stdexcept"

Vector::Vector(size_t size, double starting_val): Matrix(size, 1, starting_val, false){}

Vector::Vector(size_t size, double* init_array): Matrix(size, 1, 0, false){
    for (size_t curr_row = 0; curr_row < size; curr_row++)  {
        this->ChangeEntry(curr_row, 0, init_array[curr_row]);
    }
}

double Vector::Dot(const Vector& rhs) {
    size_t size = this->GetDims().first;
    if (size != rhs.GetDims().first) {
        throw (std::invalid_argument("INVALID VECTOR DOT PRODUCT DIMENSIONS!"));
    }
    double product = 0.0;
    for (size_t curr_row = 0; curr_row < size; curr_row++) {
        product += (this->GetEntry(curr_row, 0) * rhs.GetEntry(curr_row, 0));
    }
    return (product);
}
