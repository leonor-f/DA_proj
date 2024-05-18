#include "Application.h"

#include <iostream>
#include <unordered_map>
#include <cmath>


#include "NetworkPoint.h"

using namespace std;

/*double convertToRadians(double degree) {
    return degree * M_PI / 180.0;
}

double  haversine(double lat1, double lon1, double lat2, double lon2) {
    double rad_lat1 = convertToRadians(lat1);
    double rad_lon1 = convertToRadians(lon1);
    double rad_lat2 = convertToRadians(lat2);
    double rad_lon2 = convertToRadians(lon2);

    double delta_lat = rad_lat2 - rad_lat1;
    double delta_lon = rad_lon2 - rad_lon1;

    double aux = pow(sin(delta_lat / 2.0), 2) + cos(rad_lat1) * cos(rad_lat2) * pow(sin(delta_lon / 2.0), 2);
    double c = 2.0 * atan2(sqrt(aux), sqrt(1.0 - aux));

    const double earthradius = 6371000.0; // in meters

    return earthradius * c;
}*/

void Application::fullyConnectGraph() {
    for (const auto &p: network_->getVertexSet()) {
        auto point1 = p.second->getInfo();
        for (const auto &pp: network_->getVertexSet()) {
            auto point2 = pp.second->getInfo();
            if (point1 == point2) continue;

            bool edgeExists = false;

            // Check if there is already an edge between point1 and point2
            for (const auto &edge: p.second->getAdj()) {
                if (edge->getDest()->getInfo() == point2) {
                    edgeExists = true;
                    break;
                }
            }

            if (!edgeExists) {
                // Calculate distance between point1 and point2
                auto point1Coords = nodes.getCoordinates(point1.getId());
                auto point2Coords = nodes.getCoordinates(point2.getId());
                double distance = network_->haversine(point1Coords.second, point1Coords.first, point2Coords.second,
                                            point2Coords.first);

                // Add edge between point1 and point2
                network_->addEdge(point1, point2, distance);
            }
        }
    }
}

Application::Application() {
    network_ = nullptr;
    isFileRead = false;
    menu();
}

void Application::menu() {
    int choice;

    do {
        cout << endl << "+-----------------------------------------------------------------------------+" << endl;
        cout << "|  Welcome to the Routing Algorithm for Ocean Shipping and Urban Deliveries!  |" << endl;
        cout << "|                             We hope to be useful!                           |" << endl;
        cout << "+-----------------------------------------------------------------------------+" << endl << endl;
        cout << "1. Load and parse data" << endl;
        cout << "2. Backtracking Algorithm" << endl;
        cout << "3. Triangular Approximation Heuristic" << endl;
        cout << "4. Other Heuristics" << endl;
        cout << "5. TSP in the Real World" << endl << endl;

        cout << "0. Exit" << endl << endl;
        cout << "Enter your choice and press ENTER: ";
        cin >> choice;

        string sel;
        switch (choice) {
            case 1:
                loadData();
                break;
            case 2:
                backtracking();
                break;
            case 3:
                triangular();
                break;
            case 4:
                other();
                break;
            case 5:
                realWorld();
                break;
            case 0:
                cout << "Exiting program. Goodbye!" << endl;
                exit(EXIT_SUCCESS);
            default:
                cout << "Invalid choice. Please try again." << endl;
        }
    } while (choice != 0);
}

void Application::loadData() {
    int choice;

    /*if (network_ != nullptr) {
        delete network_;
        network_ = nullptr;
    }*/

    network_ = new Graph;

    do {
        cout << endl;
        cout << "1. Toy Graph" << endl;
        cout << "2. Extra Medium-Size Graphs" << endl;
        cout << "3. Real-world-Graphs" << endl << endl;

        cout << "0. Go back" << endl << endl;
        cout << "Enter your choice and press ENTER: ";
        cin >> choice;

        string sel;
        switch (choice) {
            case 1:
                loadToyGraph();
                break;
            case 2:
                loadMediumGraph();
                break;
            case 3:
                loadRealGraph();
                break;
            case 0:
                menu();
                break;
            default:
                cout << "Invalid choice. Please try again." << endl;
        }
    } while (choice != 0);
}

void Application::backtracking() {
    if (!isFileRead) {
        cout << "\nPlease select a graph to read:\n";
        loadData();
    }

    clock_t start = clock(); // começar a contar o tempo de execução

    // calcular backtracking
    vector<unsigned int> path;
    double minDist = network_->tspBT(path);

    clock_t end = clock(); // terminar a contagem do tempo de execução

    double duration =
            static_cast<double>(end - start) / CLOCKS_PER_SEC * 1000; // milisegundos

    cout << endl << "Path:" << endl;
    for (auto i = 0; i < path.size() - 1; i++)
        cout << path[i] << " -> ";
    cout << path[path.size() - 1];

    cout << endl << "Minimum cost: " << minDist << endl;
    cout << "Execution time: " << duration << " milliseconds" << endl;

    goBack();
}

void Application::triangular() {
    if (!isFileRead) {
        cout << "\nPlease select a graph to read:\n";
        loadData();
    }
    cout << "Fully connecting graph..." << endl;
    if (needToConnect) fullyConnectGraph();
    cout << "Done! Continuing..." << endl;
    auto g = network_->aproxTSP();

    double total = 0.0;

    for (size_t i = 0; i < g.size() - 1; i++) {
        cout << "From " << g.at(i).getId() << " to " << g.at(i + 1).getId() << " with weight of "
             << network_->getEdgeWeight(g.at(i), g.at(i + 1)) << endl;
        total += network_->getEdgeWeight(g.at(i), g.at(i + 1));
    }

    cout << "Total weight: " << total << endl;

    //do not forget to compare with backtracking for the small graphs
    goBack();
}


void Application::other() {
    if (!isFileRead) {
        cout << "\nPlease select a graph to read:\n";
        loadData();
    }

    clock_t start = clock();

    // calcular heurística
    vector<unsigned int> path;
    double heuristicDist = network_->tspHeuristic(path);

    clock_t end = clock();
    double duration = static_cast<double>(end - start) / CLOCKS_PER_SEC * 1000; // milisegundos

    cout << endl << "Heuristic Path:" << endl;
    for (auto i = 0; i < path.size() - 1; i++)
        cout << path[i] << " -> ";
    cout << path[path.size() - 1];

    cout << endl << "Total cost: " << heuristicDist << endl;
    cout << "Execution time: " << duration << " milliseconds" << endl;

    goBack();
}

void Application::realWorld() {
    if (!isFileRead) {
        cout << "\nPlease select a graph to read:\n";
        loadData();
    }

    string origin;
    cout << "Enter the **origin** of the Travelling Salesperson Problem: ";
    cin >> origin;

    // cout << "Path does not exist! << endl;
    goBack();
}

void Application::goBack() {
    int choice;
    do {
        cout << endl << "+--------------------+" << endl << "| Press 0 to go back |" << endl
             << "+--------------------+" << endl;
        cin >> choice;
    } while (choice != 0);

    menu();
}

void Application::loadToyGraph() {
    if (network_ != nullptr) {
        int choice;

        vector<vector<string>> data;

        needToConnect = false;

        do {
            cout << endl;
            cout << "1. Shipping" << endl;
            cout << "2. Stadiums" << endl;
            cout << "3. Tourism" << endl << endl;

            cout << "0. Go back" << endl << endl;
            cout << "Enter your choice and press ENTER: ";
            cin >> choice;

            string sel;
            switch (choice) {
                case 1:
                    file.setFile("../../dataset/Toy-Graphs/shipping.csv", true);
                    data = file.getData();

                    for (const auto &line: data) {
                        NetworkPoint source(std::stoi(line.at(0)));
                        NetworkPoint dest(std::stoi(line.at(1)));

                        network_->addVertex(source);
                        network_->addVertex(dest);
                        network_->addBidirectionalEdge(source, dest, std::stod(line.at(2)));
                    }
                    isFileRead = true;
                    menu();
                    break;
                case 2:
                    file.setFile("../../dataset/Toy-Graphs/stadiums.csv", true);
                    data = file.getData();

                    for (const auto &line: data) {
                        NetworkPoint source(std::stoi(line.at(0)));
                        NetworkPoint dest(std::stoi(line.at(1)));

                        network_->addVertex(source);
                        network_->addVertex(dest);
                        network_->addBidirectionalEdge(source, dest, std::stod(line.at(2)));
                    }
                    isFileRead = true;
                    menu();
                    break;
                case 3:
                    file.setFile("../../dataset/Toy-Graphs/tourism.csv", true);
                    data = file.getData();

                    for (const auto &line: data) {
                        NetworkPoint source(std::stoi(line.at(0)));
                        NetworkPoint dest(std::stoi(line.at(1)));

                        network_->addVertex(source);
                        network_->addVertex(dest);
                        network_->addBidirectionalEdge(source, dest, std::stod(line.at(2)));
                    }
                    isFileRead = true;
                    menu();
                    break;
                case 0:
                    menu();
                    break;
                default:
                    cout << "Invalid choice. Please try again." << endl;
            }
        } while (choice != 0);
    }
}

void Application::loadMediumGraph() {
    if (network_ != nullptr) {
        int choice;
        vector<vector<string>> data;
        needToConnect = true;
        cout << endl;
        cout << "Enter the number of edges {25,50,75,100,200,300,400,500,600,700,800,900}:";
        cin >> choice;

        if (!file.setFile("../../dataset/Extra_Fully_Connected_Graphs/edges_" + to_string(choice) + ".csv", false)) {
            cout << "File does not exist.\n";
            return;
        }
        nodes.setFile("../../dataset/Extra_Fully_Connected_Graphs/nodes.csv", true);
        data = file.getData();

        for (const auto &line: data) {
            NetworkPoint source(std::stoi(line.at(0)), nodes.getCoordinates(stoi(line.at(0))).second,
                                nodes.getCoordinates(stoi(line.at(0))).first);
            NetworkPoint dest(std::stoi(line.at(1)), nodes.getCoordinates(stoi(line.at(0))).second,
                              nodes.getCoordinates(stoi(line.at(0))).first);
            network_->addVertex(source);
            network_->addVertex(dest);
            network_->addBidirectionalEdge(source, dest, std::stod(line.at(2)));
        }
        isFileRead = true;
        menu();
    }
}

void Application::loadRealGraph() {
    if (network_ != nullptr) {
        int choice;

        vector<vector<string>> data;
        needToConnect = true;
        cout << endl;
        cout << "Select the graph {1,2,3}:";
        cin >> choice;
        if (!file.setFile("../../dataset/Real-world-Graphs/graph" + to_string(choice) + "/edges.csv", true)) {
            cout << "File does not exist.\n";
            return;
        }

        nodes.setFile("../../dataset/Real-world-Graphs/graph" + to_string(choice) + "/nodes.csv", true);

        data = file.getData();

        for (const auto &line: data) {
            NetworkPoint source(stoi(line.at(0)), nodes.getCoordinates(stoi(line.at(0))).second,
                                nodes.getCoordinates(stoi(line.at(0))).first);
            NetworkPoint dest(stoi(line.at(1)), nodes.getCoordinates(stoi(line.at(0))).second,
                              nodes.getCoordinates(stoi(line.at(0))).first);
            network_->addVertex(source);
            network_->addVertex(dest);
            network_->addBidirectionalEdge(source, dest, stod(line.at(2)));
        }
        isFileRead = true;
        menu();
    }
}
