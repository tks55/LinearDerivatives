#include "filereader.hpp"
#include "function_list.hpp"
#include "linalg.hpp"
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <functional>


/*Parameters: std::string file_path; Desc: Reads in a given .csv file, storing the contents in a tuple containing
a map (from column names to column indices), an array of strings (containing row indices) and a Matrix (containing
the numerical data stored). Returns this tuple for further usage, either by the user or other functions.*/
std::tuple<std::map<std::string, size_t>, std::string*, Matrix> FileReader::ReadCSV(std::string file_path) {
    if (!CheckFileExtension(file_path, ".csv")) {
        throw std::invalid_argument("INVALID CSV FILE TYPE!");
    }
    std::ifstream file(file_path);
    if (!file.is_open()) {
        throw std::invalid_argument("INVALID CSV FILE PATH!");
    }
    size_t num_lines = GetLines(file_path);
    if (num_lines == 0) {
        throw std::invalid_argument("CSV FILE IS EMPTY!");
    }
    std::string first_line;
    std::getline(file, first_line);
    std::map<std::string, size_t> column_names = ColumnToIndex(first_line);
    size_t num_valid_lines = num_lines - 1;
    size_t num_valid_columns = column_names.size() - 1;
    std::string* indices = new std::string[num_valid_lines];
    Matrix csv_matrix = Matrix(num_valid_lines, num_valid_columns, 0, false);
    std::string curr_line;
    size_t curr_row = 0;
    while (std::getline(file, curr_line)) {
        std::stringstream curr_line_stream(curr_line);
        std::string curr_entry;
        size_t curr_col = 0;
        while (std::getline(curr_line_stream, curr_entry, ',')) {
            if (curr_col == 0) {
                indices[curr_row] = curr_entry;
            } else {
                double curr_val;
                try {
                    curr_val = std::stod(curr_entry);
                } catch (...) {
                    throw std::invalid_argument("INVALID FILE FORMATTING!");
                }
                csv_matrix.ChangeEntry(curr_row, curr_col - 1, curr_val);
            }
            curr_col++;
        }
        curr_row++;
    }
    std::tuple<std::map<std::string, size_t>, std::string*, Matrix> return_tuple{column_names, indices, csv_matrix};
    return (return_tuple);
}

/*PRIVATE/HELPER FUNCTION; Parameters: std::string file_path; Desc: Counts and returns the number of lines in an open file.*/
size_t FileReader::GetLines(std::string file_path) {
    std::ifstream file_lines(file_path);
    if (!file_lines.is_open()) {
        throw std::invalid_argument("INVALID FILE PATH!");
    }
    size_t num_lines = 0;
    std::string line_count;
    while (std::getline(file_lines, line_count)) {
        num_lines++;
    }
    return (num_lines);
}

/*PRIVATE/HELPER FUNCTION; Parameters: std::string file_path, std::string file_extension; Desc: Compares the file extension
of the given file_path and compares it with the desired file extension, returning true if they are the same, and false if
they differ.*/
bool FileReader::CheckFileExtension(std::string file_path, std::string file_extension) {
    size_t extension_size = file_extension.size();
    std::string given_extension = file_path.substr(file_path.size() - extension_size);
    bool is_same = (given_extension == file_extension);
    return (is_same);
}

/*PRIVATE/HELPER FUNCTION; Parameters: std::string first_line; Desc: Extracts the column names from the first line of a .csv file,
mapping them to their respective column indices. Returns the map of column names to indices.*/
std::map<std::string, size_t> FileReader::ColumnToIndex(std::string first_line) {
    std::map<std::string, size_t> column_names;
    std::stringstream first_line_stream(first_line);
    size_t curr_val = 0;
    std::string col_name;
    while (std::getline(first_line_stream, col_name, ',')) {
        if (col_name.back() == '\r') {
            col_name = col_name.substr(0, col_name.size() - 1);
        }
        column_names[col_name] = curr_val;
        curr_val++;
    }
    return (column_names);
}

/*Parameters: std::string file_path; Desc: Reads in the provided params.txt file, processing the mode and given parameters to create a
design matrix based on user-defined properties/functions. Returns a tuple containing a string representation of the desired operating mode,
a vector of strings containing the arranged column names of the user-defined design matrix, a pair containing the name and values of the 
desired prediction vector, and the Design Matrix.*/
std::tuple<std::string, std::vector<std::string>, std::pair<std::string, Vector>, Matrix> FileReader::ReadParams(std::string file_path) {
    if (!CheckFileExtension(file_path, ".txt")) {
        throw std::invalid_argument("INVALID PARAMS FILE TYPE!");
    }
    std::ifstream file(file_path);
    if (!file.is_open()) {
        throw std::invalid_argument("INVALID PARAMS FILE PATH!");
    }
    std::pair<std::string, std::string> preprocessed_data = ParamsPreprocessor(file);
    std::string mode = preprocessed_data.first;
    std::string given_fp = preprocessed_data.second;
    std::tuple<std::map<std::string, size_t>, std::string*, Matrix> csv_data = ReadCSV(given_fp);
    std::map<std::string, size_t> column_names = std::get<0>(csv_data);
    //Default behavior--may change in future
    delete[] std::get<1>(csv_data);
    Vector* column_vector_list = LinAlg::ToColumnVectors(std::get<2>(csv_data));
    FunctionList function_list = FunctionList();
    std::map<std::string, Vector> key_transform;
    std::string predict_column;
    Vector predict_vector;
    std::string function_line;
    std::pair<std::string, Vector> predicted_info;
    while (std::getline(file, function_line)) {
        if ((function_line.size() < 3) || (function_line.at(0) == '/')) {
            continue;
        } else if (function_line.substr(0, 7) == "PREDICT") {
            predict_column = function_line.substr(8, function_line.size() - 9);
            size_t predict_vector_position = (column_names.at(predict_column) - 1);
            predict_vector = column_vector_list[predict_vector_position];
            predicted_info = {predict_column, predict_vector};
            continue;
        }
        size_t pos = function_line.find("=");
        if ((pos == std::string::npos) || (pos == (function_line.size() - 1))) {
            throw std::invalid_argument("INVALID FUNCTION SETUP!");
        }
        std::string current_key = function_line.substr(0, pos);
        std::function<double(double)> current_function = function_list.GetFunction(current_key);
        std::string variables = function_line.substr(pos + 1);
        std::stringstream variables_stream(variables);
        std::string variable;
        while (std::getline(variables_stream, variable, ',')) {
            if (variable.back() == '\r') {
                variable = variable.substr(0, variable.size() - 1);
            }
            std::string transformed_variable_index = (current_key + "(" + variable + ")");
            size_t current_vect_position = (column_names.at(variable) - 1);
            Vector current_vect = column_vector_list[current_vect_position];
            key_transform[transformed_variable_index] = LinAlg::ApplyFunction(current_vect, current_function);
        }
    }
    delete[] column_vector_list;
    size_t curr_index = 0;
    Vector* transformed_vectors = new Vector[key_transform.size()];
    std::vector<std::string> keys;
    for (auto const& [key, transformed_vector] : key_transform) {
        transformed_vectors[curr_index] = transformed_vector;
        keys.push_back(key);
        curr_index++;
    }
    Matrix transformed_matrix = LinAlg::ToMatrix(transformed_vectors, key_transform.size());
    std::tuple<std::string, std::vector<std::string>, std::pair<std::string, Vector>, Matrix> return_tuple{mode, keys, predicted_info, transformed_matrix};
    return (return_tuple);
}

/*PRIVATE/HELPER FUNCTION; Parameters: std::ifstream& file; Desc: Processes the initalization variables (mode and csv file path) from the params.txt file,
returning the values in a std::pair.*/
std::pair<std::string, std::string> FileReader::ParamsPreprocessor(std::ifstream& file) {
    std::string given_fp;
    std::string mode = "DEFAULT";
    std::string setup_lines;
    while (std::getline(file, setup_lines)) {
        if (setup_lines.size() < 10) {
            continue;
        } else if (setup_lines.substr(0, 8) == "SET_MODE") {
            mode = setup_lines.substr(9);
        } else if (setup_lines.substr(0, 8) == "INPUT_FP") {
            if (setup_lines.at(9) != '\"') {
                throw std::invalid_argument("INVALID FILE PATH FORMATTING!");
            }
            size_t find_end_quote = setup_lines.find("\"", 11);
            given_fp = setup_lines.substr(10, find_end_quote - 10);
            break;
        } else {
            continue;
        }
    }
    std::pair<std::string, std::string> preprocessed_data{mode, given_fp};
    return (preprocessed_data);
}