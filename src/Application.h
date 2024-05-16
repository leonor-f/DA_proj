#ifndef PROJ2_APPLICATION_H
#define PROJ2_APPLICATION_H

#include <string>

#include "FileReader.h"
#include "Graph.h"

class Application {
    //Graph *network_ = new Graph;
    Graph *network_;
    FileReader file;
    FileReader nodes;
    void fullyConnectGraph();

public:
    Application();
    void menu();
    void loadData();
    void backtracking();
    void triangular();
    void other();
    void realWorld(std::string c);
    void goBack();
    void loadToyGraph();
    void loadMediumGraph();
    void loadRealGraph();
    bool isFileRead;
};

#endif //PROJ2_APPLICATION_H