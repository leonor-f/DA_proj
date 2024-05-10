//
// Created by leono on 01/05/2024.
//

#ifndef PROJ2_APPLICATION_H
#define PROJ2_APPLICATION_H

#include <string>

#include "Graph.h"

enum DATASET {
    SHIPPING,
    STADIUMS,
    TOURISM,
};

class Application {
    Graph network_;

public:
    Application(DATASET mode);
};

#endif //PROJ2_APPLICATION_H
