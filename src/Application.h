//
// Created by leono on 01/05/2024.
//

#ifndef PROJ2_APPLICATION_H
#define PROJ2_APPLICATION_H

#include <string>
#include "Graph.h"


class Application {


public:
    Application();
    void menu();
    void loadData();
    void backtracking();
    void triangular();
    void other();
    void realWorld(std::string c);
    void goBack();
};

#endif //PROJ2_APPLICATION_H
