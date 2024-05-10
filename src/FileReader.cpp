#include "FileReader.h"

FileReader::FileReader() {
}

std::vector<std::vector<std::string>> FileReader::getData() {
    std::string line;
    // Clear data_ before reading new data
    data_.clear();

    // Read and process the first line
    if (ignoreFirstLine) {
        if (std::getline(file_, line)) {
        }
    }
    // Read and process the rest of the lines
    while (std::getline(file_, line)) {
        std::istringstream iss(line);
        std::string value;
        std::vector<std::string> v;
        while (std::getline(iss, value, ',')) {
            v.push_back(value);
        }
        data_.push_back(v);
    }
    return data_;
}

bool FileReader::setFile(std::string file, bool ignore) {
    file_.open(file);
    if (file_.fail()) {
        return false;
    }
    ignoreFirstLine = ignore;
    return true;
}

std::pair<double, double> FileReader::getCoordinates(unsigned id) {
    file_.seekg(0, std::ios::beg);
    std::string line;

    std::string idStr = std::to_string(id);

    std::pair<double, double> res;

    while (std::getline(file_, line)) {
        if (line.compare(0, idStr.length(), idStr) == 0) {
            std::istringstream iss(line);
            std::string value;
            std::vector<std::string> v;
            while (std::getline(iss, value, ',')) {
                v.push_back(value);
            }

            res.first = std::stod(v.at(1));
            res.second = std::stod(v.at(2));

            return res;
        }
    }

}
