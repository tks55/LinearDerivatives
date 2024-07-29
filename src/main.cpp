#include "regression.hpp"
#include "matrix.hpp"
#include "vector.hpp"
#include "filereader.hpp"
#include <iostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

int main() {
    std::tuple<std::string, std::vector<std::string>, std::pair<std::string, Vector>, Matrix> input_tuple = FileReader::ReadParams("./params.txt");
    std::vector<std::string> predictor_names = std::get<1>(input_tuple);
    Vector predicted_variable = std::get<2>(input_tuple).second;
    Matrix regression_matrix = std::get<3>(input_tuple);
    Vector coefs = Regression::LinearRegression(regression_matrix, predicted_variable);
    for (size_t i = 0; i < predictor_names.size(); i++) {
        std::cout << (predictor_names.at(i) + ": " + std::to_string(coefs.GetEntry(i))) << std::endl;
    }
}