#ifndef FUNCTION_LIST_HPP
#define FUNCTION_LIST_HPP

#include <string>
#include <unordered_map>
#include <functional>

class FunctionList {
    public:

        //Constructor
        FunctionList();

        //Interface Functions 
        void AddFunction(std::string name, std::function<double(double)> function);
        std::function<double(double)> GetFunction(std::string name);
        static std::function<double(double)> ToPower(double power);
        void DeleteFunction(std::string name);
        static std::string ListFunctions(const FunctionList& fl);
        std::string ListFunctions();
        
    private:

        //Instance Variables
        std::unordered_map<std::string, std::function<double(double)>> functions;
        
};

#endif