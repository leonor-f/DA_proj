//
// Created by leono on 01/05/2024.
//

#include "Application.h"

#include "FileReader.h"
#include <iostream>

using namespace std;

Application::Application(DATASET mode) {
    FileReader file;
    switch (mode) {
        case SHIPPING:
            file.setFile("../dataset/Toy-Graphs/shipping.csv");
            break;
        case STADIUMS:
            file.setFile("../dataset/Toy-Graphs/shipping.csv");
            break;
        case TOURISM:
            file.setFile("../dataset/Toy-Graphs/shipping.csv");
            break;
        default:
            return;
    }

    auto data = file.getData();


}
