#include "function_list.hpp"
#include <cmath>
#include <stdexcept>

/*Parameters: double power; Desc: Returns a function that inputs a double value and returns
the value of the input to a given power.*/
std::function<double(double)> FunctionList::ToPower(double power) {
    std::function<double(double)> to_power = [power] (double input) -> double {
        return pow(input, power);
    };
    return (to_power);
}

/*NON-MEMBER FUNCTION; Returns the double value of 1.0 for all double inputs.*/
double Const(double input) {
    if (input == 0) {
        return (1);
    }
    return ((FunctionList::ToPower(0))(input));
}

/*Constructor (Default)*/
FunctionList::FunctionList() {
    //Constant
    functions["CONST"] = Const;

    //Polynomials
    functions["LINEAR"] = ToPower(1);
    functions["QUADRATIC"] = ToPower(2);
    functions["CUBIC"] = ToPower(3);
    functions["QUARTIC"] = ToPower(4);
    functions["QUINTIC"] = ToPower(5);

    //Inverse
    functions["INVERSE"] = ToPower(-1);

    //Root
    functions["SQRT"] = sqrt;
    functions["CBRT"] = cbrt;

    //Trigonmetric
    functions["SIN"] = sin;
    functions["COS"] = cos;
    functions["TAN"] = tan;

    //Logarithmic/Exponential
    functions["EXP"] = exp;
    functions["LOGE"] = log;
    functions["LOG2"] = log2;
    functions["LOG10"] = log10;
}

/*Parameters: std::string name, std::function<double(double)> function; Desc: Adds a function of input
and output double to the instance variable map, functions, with a key of the provided name.*/
void FunctionList::AddFunction(std::string name, std::function<double(double)> function) {
    functions[name] = function;
}

/*Parameters: std::string name; Desc: Retrieves a function from the instance variable map, functions, given
a string key of the provided name.*/
std::function<double(double)> FunctionList::GetFunction(std::string name) {
    if (functions.count(name) == 0) {
        throw std::invalid_argument("INVALID FUNCTION NAME!");
    }
    std::function<double(double)> function = functions[name];
    return (function);
}

/*Parameters: std::string name; Desc: Deletes a function from the instance variable map, functions, given
a string key of the provided name.*/
void FunctionList::DeleteFunction(std::string name) {
    size_t return_val = functions.erase(name);
    if (return_val == 0) {
        throw std::invalid_argument("INVALID FUNCTION NAME!");
    }
}

/*Parameters: const FunctionList& fl; Desc: Returns a string listing all the keys stored in the map of the given
FunctionList.*/
std::string FunctionList::ListFunctions(const FunctionList& fl) {
    std::string output_string = "";
    for (auto const& [name, function] : fl.functions) {
        output_string += (name + "\n");
    }
    return (output_string);
}

/*Parameters: none Desc: Returns a string listing all the keys stored in the map of the current FunctionList.*/
std::string FunctionList::ListFunctions() {
    std::string output_string = ListFunctions(*this);
    return (output_string);
}