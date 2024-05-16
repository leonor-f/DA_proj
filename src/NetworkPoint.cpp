#include "NetworkPoint.h"

NetworkPoint::NetworkPoint(unsigned int id) {
    id_ = id;
}

unsigned NetworkPoint::getId() const {
    return id_;
}

bool NetworkPoint::operator==(const NetworkPoint &other) {
    return this->id_ == other.id_;
}
