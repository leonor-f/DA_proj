#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <unordered_map>
#include <float.h>
#include "NetworkPoint.h"
#include "Application.cpp"
#include <unordered_set>

template <class T>
class Edge;

#define INF std::numeric_limits<double>::max()

/************************* Vertex  **************************/

template <class T>
class Vertex {
public:
    Vertex(T in);

    T getInfo() const;
    std::vector<Edge<T> *> getAdj() const;
    bool isVisited() const;
    bool isProcessing() const;
    unsigned int getIndegree() const;
    double getDist() const;
    Edge<T> *getPath() const;
    std::vector<Edge<T> *> getIncoming() const;

    void setInfo(T info);
    void setVisited(bool visited);
    void setProcesssing(bool processing);
    void setIndegree(unsigned int indegree);
    void setDist(double dist);
    void setPath(Edge<T> *path);
    Edge<T> * addEdge(Vertex<T> *dest, double w);
    bool removeEdge(T in);
    void removeOutgoingEdges();

    const Edge<T> *getEdge(const Vertex<T> *dest) const;

protected:
    T info;                // info node
    std::vector<Edge<T> *> adj;  // outgoing edges

    // auxiliary fields
    bool visited = false; // used by DFS, BFS, Prim ...
    bool processing = false; // used by isDAG (in addition to the visited attribute)
    unsigned int indegree; // used by topsort
    double dist = 0;
    Edge<T> *path = nullptr;

    std::vector<Edge<T> *> incoming; // incoming edges

    void deleteEdge(Edge<T> *edge);

};

/********************** Edge  ****************************/

template <class T>
class Edge {
public:
    Edge(Vertex<T> *orig, Vertex<T> *dest, double w);

    Vertex<T> * getDest() const;
    double getWeight() const;
    bool isSelected() const;
    Vertex<T> * getOrig() const;
    Edge<T> *getReverse() const;
    double getFlow() const;

    void setSelected(bool selected);
    void setReverse(Edge<T> *reverse);
    void setFlow(double flow);
protected:
    Vertex<T> * dest; // destination vertex
    double weight; // edge weight, can also be used for capacity

    // auxiliary fields
    bool selected = false;

    // used for bidirectional edges
    Vertex<T> *orig;
    Edge<T> *reverse = nullptr;

    double flow; // for flow-related problems
};

/********************** Graph  ****************************/

class Graph {
public:
    ~Graph();
    /*
    * Auxiliary function to find a vertex with a given the content.
    */
    Vertex<NetworkPoint> *findVertex(const NetworkPoint &in) const;
    /*
     *  Adds a vertex with a given content or info (in) to a graph (this).
     *  Returns true if successful, and false if a vertex with that content already exists.
     */
    bool addVertex(const NetworkPoint &in);
    bool removeVertex(const NetworkPoint &in);

    /*
     * Adds an edge to a graph (this), given the contents of the source and
     * destination vertices and the edge weight (w).
     * Returns true if successful, and false if the source or destination vertex does not exist.
     */
    bool addEdge(const NetworkPoint &sourc, const NetworkPoint &dest, double w);
    bool removeEdge(const NetworkPoint &source, const NetworkPoint &dest);
    bool addBidirectionalEdge(const NetworkPoint &sourc, const NetworkPoint &dest, double w);

    Edge<NetworkPoint> *getNearestNeighbor(Vertex<NetworkPoint> *v) const;

    int getNumVertex() const;
    std::unordered_map<unsigned, Vertex<NetworkPoint> *> getVertexSet() const;

    std:: vector<NetworkPoint> dfs() const;
    std:: vector<NetworkPoint> dfs(const NetworkPoint & source) const;
    void dfsVisit(Vertex<NetworkPoint> *v,  std::vector<NetworkPoint> & res) const;
    std::vector<NetworkPoint> bfs(const NetworkPoint & source) const;

    bool isDAG() const;
    bool dfsIsDAG(Vertex<NetworkPoint> *v) const;
    std::vector<NetworkPoint> topsort() const;
    /**
     * @brief Auxiliary function to copy a graph.
     * @Complexity O(V + E)
     * @return A pointer to the new graph.
     */
    Graph * copyGraph();
    std::vector<NetworkPoint> aproxTSP();
    Graph computeMST(Vertex<NetworkPoint> *root);

    double getEdgeWeight(const NetworkPoint &a, const NetworkPoint &b) const;

    /**
     * @brief Recursive helper function for solving the Traveling Salesman Problem using backtracking.
     *
     * This function attempts to find the shortest path that visits all vertices in the graph exactly once
     * and returns to the starting vertex. It explores all possible paths using backtracking and updates
     * the minimum distance and corresponding path whenever a shorter path is found.
     *
     * @param curIndex The current index in the path being constructed.
     * @param curDist The current total distance of the path being constructed.
     * @param curPath The current path being constructed.
     * @param minDist Reference to the minimum distance found so far.
     * @param path Reference to the path corresponding to the minimum distance found so far.
     * @complexity O(V!)
     */
    void tspBTRec(unsigned int curIndex, double curDist, std::vector<unsigned int> &curPath, double &minDist,
                  std::vector<unsigned int> &path) const;

    /**
     * @brief Solves the Traveling Salesman Problem using backtracking.
     *
     * This function initializes the necessary variables and calls the recursive helper function `tspBTRec`
     * to compute the shortest possible route that visits each vertex exactly once and returns to the starting point.
     * It returns the minimum distance of such a path and stores the corresponding path in the provided vector.
     *
     * @param path Reference to a vector where the path of the minimum distance will be stored.
     * @return The minimum distance of the Traveling Salesman Problem solution.
     * @complexity O(V!)
     */
    double tspBT(std::vector<unsigned int> &path) const;

protected:
    std::unordered_map<unsigned, Vertex<NetworkPoint> *> vertexSet;    // vertex set
};

void deleteMatrix(int **m, int n);
void deleteMatrix(double **m, int n);


#endif //GRAPH_H