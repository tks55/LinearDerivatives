#ifndef FILEREADER_HPP
#define FILEREADER_HPP

#include "matrix.hpp"
#include "vector.hpp"
#include <tuple>
#include <map>
#include <string>
#include <vector>

class FileReader {

    public:

        //Class Functions
        static std::tuple<std::map<std::string, size_t>, std::string*, Matrix> ReadCSV(std::string file_path);
        static std::tuple<std::string, std::vector<std::string>, std::pair<std::string, Vector>, Matrix> ReadParams(std::string file_path);

    private:

        //Private/Helper Functions
        static size_t GetLines(std::string file_path);
        static bool CheckFileExtension(std::string file_path, std::string file_extension);
        static std::map<std::string, size_t> ColumnToIndex(std::string first_line);
        static std::pair<std::string, std::string> ParamsPreprocessor(std::ifstream& file);

};

#endif