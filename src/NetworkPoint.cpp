//
// Created by alexandre on 03-05-2024.
//

#include "NetworkPoint.h"

NetworkPoint::NetworkPoint(unsigned int id) {
    id_ = id;
}

unsigned NetworkPoint::getId() const {
    return id_;
}
