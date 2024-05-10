#include "FileReader.h"

FileReader::FileReader() {
}

std::vector<std::vector<std::string>> FileReader::getData() {
    std::string line;
    // Clear data_ before reading new data
    data_.clear();

    // Read and process the first line
    if (std::getline(file_, line)) {
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

bool FileReader::setFile(std::string file) {
    file_.open(file);
    if (file_.fail()) {
        return false;
    }
    return true;
}
