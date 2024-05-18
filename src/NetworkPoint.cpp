#include "NetworkPoint.h"

NetworkPoint::NetworkPoint(unsigned int id) {
    id_ = id;
}

NetworkPoint::NetworkPoint(unsigned int id, double lat, double lon) {
    id_ = id;
    lat_ = lat;
    lon_ = lon;
}

unsigned NetworkPoint::getId() const {
    return id_;
}

double NetworkPoint::getLat() const {
    return lat_;
}

double NetworkPoint::getLon() const {
    return lon_;
}

bool NetworkPoint::operator==(const NetworkPoint &other) {
    return this->id_ == other.id_;
}
