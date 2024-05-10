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

    for (const auto &line : data) {
        NetworkPoint source(std::stoi(line.at(0)));
        NetworkPoint dest(std::stoi(line.at(1)));

        network_.addVertex(source);
        network_.addVertex(dest);
        network_.addEdge(source, dest, std::stod(line.at(2)));
    }

    /*
    for (const auto &vertex : network_.getVertexSet()) {
        std::cout << vertex.second->getInfo().getId() << std::endl;

        for (const auto &edge : vertex.second->getAdj()) {
            std::cout << edge->getDest()->getInfo().getId() << " " << edge->getWeight() << std::endl;
        }

        std::cout << std::endl;
    }
    */
}
