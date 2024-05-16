#ifndef SRC_NETWORKPOINT_H
#define SRC_NETWORKPOINT_H


class NetworkPoint {
    unsigned id_;
public:
    NetworkPoint(unsigned id);
    unsigned getId() const;
    bool operator==(const NetworkPoint &other);

};


#endif //SRC_NETWORKPOINT_H
