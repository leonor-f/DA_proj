#ifndef SRC_NETWORKPOINT_H
#define SRC_NETWORKPOINT_H


class NetworkPoint {
    unsigned id_;
    double lat_; // latitude
    double lon_; // longitude
public:
    NetworkPoint(unsigned id);
    NetworkPoint(unsigned int id, double lat, double lon);

    [[nodiscard]] unsigned getId() const;
    [[nodiscard]] double getLat() const;
    [[nodiscard]] double getLon() const;

    bool operator==(const NetworkPoint &other);
};


#endif //SRC_NETWORKPOINT_H
