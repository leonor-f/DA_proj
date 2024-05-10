//
// Created by leono on 01/05/2024.
//

#include "Application.h"

#include "FileReader.h"
#include <iostream>
#include <unordered_map>

#include "NetworkPoint.h"

using namespace std;

Application::Application(DATASET mode) {
    FileReader file;
    switch (mode) {
        case SHIPPING:
            file.setFile("../../dataset/Toy-Graphs/shipping.csv");
            break;
        case STADIUMS:
            file.setFile("../../dataset/Toy-Graphs/shipping.csv");
            break;
        case TOURISM:
            file.setFile("../../dataset/Toy-Graphs/shipping.csv");
            break;
        default:
            return;
    }

    auto data = file.getData();
    for (const auto &d : data) {
        for (const auto &s : d) {
            std::cout << s << " ";
        }
        std::cout <<std::endl;
    }

}
