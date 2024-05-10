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
     * @brief Calculates the maximum flow to the specified city.
     * @Complexity:  O(E) if Edmonds-Karp has already been called.
     * If not: O(V * E^2)
     * @param city The city vertex to calculate the maximum flow for.
     * @return The maximum flow to the specified city.
     */
    double getMaxFlow(NetworkPoint city);

    /**
     * @brief Checks water supply for each city in the graph.
     * @Complexity:  O(V + E) if Edmonds-Karp has already been called.
     * If not: O(V * E^2)
     * @return A vector of pairs, where each pair contains the city code, demand, and maximum flow deficit.
     */
    std::vector<std::pair<std::string, std::pair<double, double>>> checkWaterSupply();

    /**
     * @brief Prints the values returned from the 'checkWaterSupply()' function.
     * @Complexity O(n)
     * @param supply The vector returned from the 'checkWaterSupply()' function.
     */
    void printWaterSupply(std::vector<std::pair<std::string, std::pair<double, double>>> supply);

    /**
    * @brief Calculates average difference, variance, and standard deviation of the flow to all cities
    * @Complexity:  O(V + E) if Edmonds-Karp has already been called.
    * If not: O(V * E^2)
    * @return A vector containing the calculated metrics: average difference, variance, and standard deviation.
    */
    std::vector<double> calculateMetrics();

    /**
     * @brief Balances the load by redistributing flow in the graph.
     * @Complexity O(V + E)
     * @param averageDifference The average difference used for load balancing.
     */
    void balanceLoad(double averageDifference);

    /**
     * @brief Prints the results obtained from the 'calculateMetrics()' function
     * @Complexity O(n)
     * @param metric A vector containing the calculated metrics.
     */
    void printMetrics(std::vector<double> metric);

    /**
     * @brief Auxiliary function to copy a graph.
     * @Complexity O(V + E)
     * @return A pointer to the new graph.
     */
    Graph * copyGraph();

protected:
    std::unordered_map<unsigned, Vertex<NetworkPoint> *> vertexSet;    // vertex set
    /**
     * @brief Calculates the minimal residual capacity along the augmenting path.
     * @Complexity: O(V)
     * @param s Pointer to the source vertex.
     * @param t Pointer to the sink vertex.
     * @return The minimal residual capacity along the augmenting path.
     */
    unsigned findMinimalResidualAlongPath(Vertex<NetworkPoint> *s, Vertex<NetworkPoint> *p) const;
    /**
    * @brief Finds an augmenting path from the super source vertex to the super sink vertex using BFS.
    * @Complexity: O(V + E)
    * @param s Pointer to the source vertex.
    * @param t Pointer to the sink vertex.
    * @return True if an augmenting path is found; otherwise, false.
    */
    bool findAugPath(Vertex<NetworkPoint> *s, Vertex<NetworkPoint> *t) ;
    /**
     * @brief Augments flow along the path from sink to source.
     * @Complexity: O(V)
     * @param s Pointer to the source vertex.
     * @param t Pointer to the sink vertex.
     * @param f Amount of flow to be augmented along the path.
     */
    void augment(Vertex<NetworkPoint> *s, Vertex<NetworkPoint> *t, double f);
    double ** distMatrix = nullptr;   // dist matrix for Floyd-Warshall
    int **pathMatrix = nullptr;   // path matrix for Floyd-Warshall
    /**
     * @brief Adds a super source and super sink to the graph.
     * @Complexity: O(V)
     */
    void addSuperSourceAndSink();

    /**
     * serves to store if the maximum flow edmond karps algorithm has already been ran
     */
    bool maxFlowRan = false;

    /**
     * @brief Tests and visits a vertex for the 'findAugPath' function.
     * Time Complexity: O(1)
     * @param q Reference to the queue.
     * @param e Pointer to the edge connecting to the vertex.
     * @param w Pointer to the vertex to be tested and visited.
     * @param residual Residual capacity of the edge connecting to the vertex.
     */
    void testAndVisit(std::queue<Vertex<NetworkPoint> *> &q, Edge<NetworkPoint> *e, Vertex<NetworkPoint> *w, double residual);

    /**
    * @brief Sets the flow of all edges to 0 and runs the Edmonds-Karp algorithm
    * @Complexity: O(V * E^2)
    */
    void edmondsKarp();
};

void deleteMatrix(int **m, int n);
void deleteMatrix(double **m, int n);


#endif //GRAPH_H