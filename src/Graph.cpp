#include <valarray>
#include "Graph.h"

/************************* Vertex  **************************/

template<class T>
Vertex<T>::Vertex(T in): info(in) {}

/*
 * Auxiliary function to add an outgoing edge to a vertex (this),
 * with a given destination vertex (d) and edge weight (w).
 */
template<class T>
Edge<T> *Vertex<T>::addEdge(Vertex<T> *d, double w) {
    auto newEdge = new Edge<T>(this, d, w);
    adj.push_back(newEdge);
    d->incoming.push_back(newEdge);
    return newEdge;
}

/*
 * Auxiliary function to remove an outgoing edge (with a given destination (d))
 * from a vertex (this).
 * Returns true if successful, and false if such edge does not exist.
 */
template<class T>
bool Vertex<T>::removeEdge(T in) {
    bool removedEdge = false;
    auto it = adj.begin();
    while (it != adj.end()) {
        Edge<T> *edge = *it;
        Vertex<T> *dest = edge->getDest();
        if (dest->getInfo() == in) {
            it = adj.erase(it);
            deleteEdge(edge);
            removedEdge = true; // allows for multiple edges to connect the same pair of vertices (multigraph)
        } else {
            it++;
        }
    }
    return removedEdge;
}

/*
 * Auxiliary function to remove an outgoing edge of a vertex.
 */
template<class T>
void Vertex<T>::removeOutgoingEdges() {
    auto it = adj.begin();
    while (it != adj.end()) {
        Edge<T> *edge = *it;
        it = adj.erase(it);
        deleteEdge(edge);
    }
}

template<class T>
T Vertex<T>::getInfo() const {
    return this->info;
}

template<class T>
std::vector<Edge<T> *> Vertex<T>::getAdj() const {
    return this->adj;
}

template<class T>
const Edge<T> *Vertex<T>::getEdge(const Vertex<T> *dest) const {
    for (const auto &edge: adj) {
        if (edge->getDest() == dest) {
            return edge;
        }
    }
    return nullptr;
}

template<class T>
bool Vertex<T>::isVisited() const {
    return this->visited;
}

template<class T>
bool Vertex<T>::isProcessing() const {
    return this->processing;
}

template<class T>
unsigned int Vertex<T>::getIndegree() const {
    return this->indegree;
}

template<class T>
double Vertex<T>::getDist() const {
    return this->dist;
}

template<class T>
Edge<T> *Vertex<T>::getPath() const {
    return this->path;
}

template<class T>
std::vector<Edge<T> *> Vertex<T>::getIncoming() const {
    return this->incoming;
}

template<class T>
void Vertex<T>::setInfo(T in) {
    this->info = in;
}

template<class T>
void Vertex<T>::setVisited(bool visited) {
    this->visited = visited;
}

template<class T>
void Vertex<T>::setProcesssing(bool processing) {
    this->processing = processing;
}

template<class T>
void Vertex<T>::setIndegree(unsigned int indegree) {
    this->indegree = indegree;
}

template<class T>
void Vertex<T>::setDist(double dist) {
    this->dist = dist;
}

template<class T>
void Vertex<T>::setPath(Edge<T> *path) {
    this->path = path;
}

template<class T>
void Vertex<T>::deleteEdge(Edge<T> *edge) {
    Vertex<T> *dest = edge->getDest();
    // Remove the corresponding edge from the incoming list
    auto it = dest->incoming.begin();
    while (it != dest->incoming.end()) {
        if ((*it)->getOrig()->getInfo() == info) {
            it = dest->incoming.erase(it);
        } else {
            it++;
        }
    }
    delete edge;
}

/********************** Edge  ****************************/

template<class T>
Edge<T>::Edge(Vertex<T> *orig, Vertex<T> *dest, double w): orig(orig), dest(dest), weight(w) {}

template<class T>
Vertex<T> *Edge<T>::getDest() const {
    return this->dest;
}

template<class T>
double Edge<T>::getWeight() const {
    return this->weight;
}

template<class T>
Vertex<T> *Edge<T>::getOrig() const {
    return this->orig;
}

template<class T>
Edge<T> *Edge<T>::getReverse() const {
    return this->reverse;
}

template<class T>
bool Edge<T>::isSelected() const {
    return this->selected;
}

template<class T>
double Edge<T>::getFlow() const {
    return flow;
}

template<class T>
void Edge<T>::setSelected(bool selected) {
    this->selected = selected;
}

template<class T>
void Edge<T>::setReverse(Edge<T> *reverse) {
    this->reverse = reverse;
}

template<class T>
void Edge<T>::setFlow(double flow) {
    this->flow = flow;
}

/********************** Graph  ****************************/

int Graph::getNumVertex() const {
    return vertexSet.size();
}

std::unordered_map<unsigned, Vertex<NetworkPoint> *> Graph::getVertexSet() const {
    return vertexSet;
}

/*
 * Auxiliary function to find a vertex with a given content.
 */
Vertex<NetworkPoint> *Graph::findVertex(const NetworkPoint &in) const {
    auto a = vertexSet.find(in.getId());
    if (a != vertexSet.end()) return a->second;
    return nullptr;
}

/*
 *  Adds a vertex with a given content or info (in) to a graph (this).
 *  Returns true if successful, and false if a vertex with that content already exists.
 */
bool Graph::addVertex(const NetworkPoint &in) {
    if (findVertex(in) != nullptr)
        return false;
    vertexSet.emplace(in.getId(), new Vertex<NetworkPoint>(in));
    return true;
}

/*
 *  Removes a vertex with a given content (in) from a graph (this), and
 *  all outgoing and incoming edges.
 *  Returns true if successful, and false if such vertex does not exist.
 */
bool Graph::removeVertex(const NetworkPoint &in) {
    auto it = vertexSet.find(in.getId());

    if (it == vertexSet.end()) return false;

    it->second->removeOutgoingEdges();

    auto itt = vertexSet.begin();

    while (itt != vertexSet.end()) {
        itt->second->removeEdge(it->second->getInfo());
        ++itt;
    }

    vertexSet.erase(it);

    delete it->second;
    return true;

    /*for (auto it = vertexSet.begin(); it != vertexSet.end(); it++) {
        if ((*it)->getInfo() == in) {
            auto v = *it;
            v->removeOutgoingEdges();
            for (auto u : vertexSet) {
                u->removeEdge(v->getInfo());
            }
            vertexSet.erase(it);
            delete v;
            return true;
        }
    }`*/
}

/*
 * Adds an edge to a graph (this), given the contents of the source and
 * destination vertices and the edge weight (w).
 * Returns true if successful, and false if the source or destination vertex does not exist.
 */
bool Graph::addEdge(const NetworkPoint &sourc, const NetworkPoint &dest, double w) {
    auto v1 = findVertex(sourc);
    auto v2 = findVertex(dest);
    if (v1 == nullptr || v2 == nullptr) {
        std::cout << "nullptr\n";
        return false;
    }
    v1->addEdge(v2, w);
    return true;
}

/*
 * Removes an edge from a graph (this).
 * The edge is identified by the source (sourc) and destination (dest) contents.
 * Returns true if successful, and false if such edge does not exist.
 */
bool Graph::removeEdge(const NetworkPoint &sourc, const NetworkPoint &dest) {
    Vertex<NetworkPoint> *srcVertex = findVertex(sourc);
    if (srcVertex == nullptr) {
        return false;
    }
    return srcVertex->removeEdge(dest);
}

bool Graph::addBidirectionalEdge(const NetworkPoint &sourc, const NetworkPoint &dest, double w) {
    auto v1 = findVertex(sourc);
    auto v2 = findVertex(dest);
    if (v1 == nullptr || v2 == nullptr)
        return false;
    auto e1 = v1->addEdge(v2, w);
    auto e2 = v2->addEdge(v1, w);
    e1->setReverse(e2);
    e2->setReverse(e1);
    return true;
}

Edge<NetworkPoint> *Graph::getNearestNeighbor(Vertex<NetworkPoint> *v) const {
    for (auto e: v->getAdj())
        if (!e->getDest()->isVisited())
            return e;
    return nullptr;
}

/****************** DFS ********************/

/*
 * Performs a depth-first search (dfs) traversal in a graph (this).
 * Returns a vector with the contents of the vertices by dfs order.
 */
std::vector<NetworkPoint> Graph::dfs() const {
    std::vector<NetworkPoint> res;
    for (auto v: vertexSet)
        v.second->setVisited(false);
    for (auto v: vertexSet)
        if (!v.second->isVisited())
            dfsVisit(v.second, res);
    return res;
}

/*
 * Performs a depth-first search (dfs) in a graph (this) from the source node.
 * Returns a vector with the contents of the vertices by dfs order.
 */
std::vector<NetworkPoint> Graph::dfs(const NetworkPoint &source) const {
    std::vector<NetworkPoint> res;
    // Get the source vertex
    auto s = findVertex(source);
    if (s == nullptr) {
        return res;
    }
    // Set that no vertex has been visited yet
    for (auto v: vertexSet) {
        v.second->setVisited(false);
    }
    // Perform the actual DFS using recursion
    dfsVisit(s, res);

    return res;
}

/*
 * Auxiliary function that visits a vertex (v) and its adjacent, recursively.
 * Updates a parameter with the list of visited node contents.
 */
void Graph::dfsVisit(Vertex<NetworkPoint> *v, std::vector<NetworkPoint> &res) const {
    v->setVisited(true);
    res.push_back(v->getInfo());
    for (auto &e: v->getAdj()) {
        auto w = e->getDest();
        if (!w->isVisited()) {
            dfsVisit(w, res);
        }
    }
}

/****************** BFS ********************/
/*
 * Performs a breadth-first search (bfs) in a graph (this), starting
 * from the vertex with the given source contents (source).
 * Returns a vector with the contents of the vertices by bfs order.
 */
std::vector<NetworkPoint> Graph::bfs(const NetworkPoint &source) const {
    std::vector<NetworkPoint> res;
    // Get the source vertex
    auto s = findVertex(source);
    if (s == nullptr) {
        return res;
    }

    // Set that no vertex has been visited yet
    for (auto v: vertexSet) {
        v.second->setVisited(false);
    }

    // Perform the actual BFS using a queue
    std::queue<Vertex<NetworkPoint> *> q;
    q.push(s);
    s->setVisited(true);
    while (!q.empty()) {
        auto v = q.front();
        q.pop();
        res.push_back(v->getInfo());
        for (auto &e: v->getAdj()) {
            auto w = e->getDest();
            if (!w->isVisited()) {
                q.push(w);
                w->setVisited(true);
            }
        }
    }
    return res;
}

/****************** isDAG  ********************/
/*
 * Performs a depth-first search in a graph (this), to determine if the graph
 * is acyclic (acyclic directed graph or DAG).
 * During the search, a cycle is found if an edge connects to a vertex
 * that is being processed in the stack of recursive calls (see theoretical classes).
 * Returns true if the graph is acyclic, and false otherwise.
 */

bool Graph::isDAG() const {
    for (auto v: vertexSet) {
        v.second->setVisited(false);
        v.second->setProcesssing(false);
    }
    for (auto v: vertexSet) {
        if (!v.second->isVisited()) {
            if (!dfsIsDAG(v.second)) return false;
        }
    }
    return true;
}

/**
 * Auxiliary function that visits a vertex (v) and its adjacent, recursively.
 * Returns false (not acyclic) if an edge to a vertex in the stack is found.
 */
bool Graph::dfsIsDAG(Vertex<NetworkPoint> *v) const {
    v->setVisited(true);
    v->setProcesssing(true);
    for (auto e: v->getAdj()) {
        auto w = e->getDest();
        if (w->isProcessing()) return false;
        if (!w->isVisited()) {
            if (!dfsIsDAG(w)) return false;
        }
    }
    v->setProcesssing(false);
    return true;
}

/****************** toposort ********************/
/*
 * Performs a topological sorting of the vertices of a graph (this).
 * Returns a vector with the contents of the vertices by topological order.
 * If the graph has cycles, returns an empty vector.
 * Follows the algorithm described in theoretical classes.
 */

std::vector<NetworkPoint> Graph::topsort() const {
    std::vector<NetworkPoint> res;

    for (auto v: vertexSet) {
        v.second->setIndegree(0);
    }
    for (auto v: vertexSet) {
        for (auto e: v.second->getAdj()) {
            unsigned int indegree = e->getDest()->getIndegree();
            e->getDest()->setIndegree(indegree + 1);
        }
    }

    std::queue<Vertex<NetworkPoint> *> q;
    for (auto v: vertexSet) {
        if (v.second->getIndegree() == 0) {
            q.push(v.second);
        }
    }

    while (!q.empty()) {
        Vertex<NetworkPoint> *v = q.front();
        q.pop();
        res.push_back(v->getInfo());
        for (auto e: v->getAdj()) {
            auto w = e->getDest();
            w->setIndegree(w->getIndegree() - 1);
            if (w->getIndegree() == 0) {
                q.push(w);
            }
        }
    }

    if (res.size() != vertexSet.size()) {
        //std::cout << "Impossible topological ordering!" << std::endl;
        res.clear();
        return res;
    }

    return res;
}

void preOrderTraversal(Graph &MST, Vertex<NetworkPoint> *current, std::vector<NetworkPoint> &L) {
    if (current == nullptr)
        return;

    L.push_back(current->getInfo());
    for (const auto &edge: current->getAdj()) {
        if (L.size() == 1 || L[L.size() - 2].getId() != edge->getDest()->getInfo().getId()) {
            preOrderTraversal(MST, edge->getDest(), L);
        }
    }
}

std::vector<NetworkPoint> Graph::aproxTSP() {
    auto root = vertexSet.begin()->second;
    std::unordered_set<unsigned> visited;
    std::vector<NetworkPoint> tour;

    Graph mst = computeMST(root);

    auto mst_root = mst.findVertex(root->getInfo());

    std::vector<NetworkPoint> l;

    preOrderTraversal(mst, mst_root, l);

    for (const auto &vertex: l) {
        if (visited.find(vertex.getId()) == visited.end()) {
            tour.push_back(vertex);
            visited.insert(vertex.getId());
        }
    }

    return tour;

}

struct EdgeComparator {
    bool operator()(const Edge<NetworkPoint> *lhs, const Edge<NetworkPoint> *rhs) const {
        // Compare edges based on their weights
        return lhs->getWeight() > rhs->getWeight();
    }
};

//this function is working as intended
Graph Graph::computeMST(Vertex<NetworkPoint> *root) {
    Graph MST;

    if (root == nullptr)
        return MST;

    std::priority_queue<Edge<NetworkPoint> *, std::vector<Edge<NetworkPoint> *>, EdgeComparator> pq;

    // Set to track visited vertices
    std::unordered_set<Vertex<NetworkPoint> *> visited;
    visited.insert(root);

    // Add all edges incident to the root vertex to the priority queue
    for (auto edge: root->getAdj()) {
        pq.push(edge);
    }

    // Add the root vertex to the MST
    MST.addVertex(root->getInfo());

    // Main loop of Prim's algorithm
    while (!pq.empty()) {
        // Get the edge with the smallest weight
        Edge<NetworkPoint> *minEdge = pq.top();
        pq.pop();

        // Get the destination vertex of the edge
        Vertex<NetworkPoint> *destVertex = minEdge->getDest();

        // If the destination vertex is already visited, skip this edge
        if (visited.find(destVertex) != visited.end()) {
            continue;
        }

        // Add the destination vertex to the MST
        MST.addVertex(destVertex->getInfo());

        // Add the edge to the MST
        MST.addEdge(minEdge->getOrig()->getInfo(), destVertex->getInfo(), minEdge->getWeight());

        // Mark the destination vertex as visited
        visited.insert(destVertex);

        // Add all edges incident to the destination vertex to the priority queue
        for (auto edge: destVertex->getAdj()) {
            if (visited.find(edge->getDest()) == visited.end()) {
                pq.push(edge);
            }
        }
    }

    return MST;
}

double Graph::getEdgeWeight(const NetworkPoint &a, const NetworkPoint &b) const {
    auto vertex = findVertex(a);
    for (const auto &edge: vertex->getAdj()) {
        if (edge->getDest()->getInfo() == b) {
            return edge->getWeight();
        }
    }
    return 0;
}


inline void deleteMatrix(int **m, int n) {
    if (m != nullptr) {
        for (int i = 0; i < n; i++)
            if (m[i] != nullptr)
                delete[] m[i];
        delete[] m;
    }
}

inline void deleteMatrix(double **m, int n) {
    if (m != nullptr) {
        for (int i = 0; i < n; i++)
            if (m[i] != nullptr)
                delete[] m[i];
        delete[] m;
    }
}

Graph::~Graph() {
    //deleteMatrix(distanceMatrix_, vertexSet.size());
    /*if (!vertexSet.empty()) {
        for (auto &v: vertexSet) {
            removeVertex(v.second->getInfo());
        }
        vertexSet.clear();
    }*/
}

Graph *Graph::copyGraph() {
    Graph *newGraph = new Graph();

    // Copy vertices
    for (const auto &p: vertexSet) {
        auto v = p.second;

        newGraph->addVertex(v->getInfo());
    }

    // Copy edges
    for (const auto &p: vertexSet) {
        auto v = p.second;

        for (auto e: v->getAdj()) {
            auto src = e->getOrig()->getInfo();
            auto dest = e->getDest()->getInfo();

            newGraph->addEdge(src, dest, e->getWeight());
        }
    }
    // newGraph->updateAllVerticesFlow();
    return newGraph;
}

double convertToRadians(double degree) {
    return degree * M_PI / 180.0;
}

double Graph::haversine(double lat1, double lon1, double lat2, double lon2) const {
    double rad_lat1 = convertToRadians(lat1);
    double rad_lon1 = convertToRadians(lon1);
    double rad_lat2 = convertToRadians(lat2);
    double rad_lon2 = convertToRadians(lon2);

    double delta_lat = rad_lat2 - rad_lat1;
    double delta_lon = rad_lon2 - rad_lon1;

    double aux = pow(sin(delta_lat / 2.0), 2) + cos(rad_lat1) * cos(rad_lat2) * pow(sin(delta_lon / 2.0), 2);
    double c = 2.0 * atan2(sqrt(aux), sqrt(1.0 - aux));

    const double earthradius = 6371000.0; // in meters

    return earthradius * c;
}

void Graph::tspBTRec(unsigned int curIndex, double curDist, std::vector<unsigned int> &curPath, double &minDist,
                     std::vector<unsigned int> &path) const {
    unsigned int n = getNumVertex();
    const Edge<NetworkPoint> *edge;

    if (curIndex == n && (edge = findVertex(curPath[n - 1])->getEdge(findVertex(0)))) {
        curDist += edge->getWeight();
        if (curDist < minDist) {
            minDist = curDist;
            path = curPath;
            path.push_back(0);
        }
        return;
    }

    for (unsigned int i = 1; i < n; i++) {
        edge = findVertex(curPath[curIndex - 1])->getEdge(findVertex(i));
        if (edge && (curDist + edge->getWeight() < minDist)) {
            bool isNewVertex = true;
            for (unsigned int j = 1; j < curIndex; j++) {
                if (curPath[j] == i) {
                    isNewVertex = false;
                    break;
                }
            }
            if (isNewVertex) {
                curPath[curIndex] = i;
                tspBTRec(curIndex + 1, curDist + edge->getWeight(), curPath, minDist, path);
            }
        }
    }
}

double Graph::tspBT(std::vector<unsigned int> &path) const {
    path.clear();

    std::vector<unsigned int> curPath(getNumVertex());
    double minDist = std::numeric_limits<double>::max();

    // path starts at node 0
    curPath[0] = 0;

    // in the first recursive call curIndex starts at 1 rather than 0
    tspBTRec(1, 0.0, curPath, minDist, path);

    return minDist;
}


// CLUSTERING

/*struct Cluster {
    std::vector<Vertex<NetworkPoint> *> vertices;
};

// Function to compute the distance between two clusters
double compute_cluster_distance(Cluster &cluster1, Cluster &cluster2) {
    double min_distance = std::numeric_limits<double>::infinity();
    for (const auto &v1: cluster1.vertices) {
        for (const auto &v2: cluster2.vertices) {
            double distance = haversine(v1->getInfo().getLat(), v1->getInfo().getLon(), v2->getInfo().getLat(),
                                        v2->getInfo().getLon());
            if (distance < min_distance) {
                min_distance = distance;
            }
        }
    }
    return min_distance;
}

// Function to perform hierarchical clustering on a graph
std::vector<Cluster> hierarchical_clustering(const Graph &graph, int num_clusters) {
    // Initialize clusters with each vertex as a separate cluster
    std::vector<Cluster> clusters;
    for (const auto &entry: graph.getVertexSet()) {
        Cluster cluster;
        cluster.vertices.push_back(entry.second);
        clusters.push_back(cluster);
    }

    // Perform hierarchical clustering until the desired number of clusters is reached
    while (clusters.size() > num_clusters) {
        // Find the closest pair of clusters
        double min_distance = std::numeric_limits<double>::infinity();
        int closest_cluster1 = -1, closest_cluster2 = -1;
        for (int i = 0; i < clusters.size(); ++i) {
            for (int j = i + 1; j < clusters.size(); ++j) {
                // Compute the distance between the centroids of the clusters
                double distance = compute_cluster_distance(clusters[i], clusters[j]); // Implement this function

                // Update the closest pair of clusters
                if (distance < min_distance) {
                    min_distance = distance;
                    closest_cluster1 = i;
                    closest_cluster2 = j;
                }
            }
        }

        // Merge the closest pair of clusters
        clusters[closest_cluster1].vertices.insert(clusters[closest_cluster1].vertices.end(),
                                                   clusters[closest_cluster2].vertices.begin(),
                                                   clusters[closest_cluster2].vertices.end());
        clusters.erase(clusters.begin() + closest_cluster2);
    }

    return clusters;
}

// Function to perform k-means clustering on a graph with dynamic number of clusters based on maximum number of vertices per cluster
[[maybe_unused]] vector<Cluster> k_means_clustering(const Graph &graph, int max_clusters) {
    // Initialize clusters
    vector<Cluster> clusters;

    // Initialize centroids with random vertices from the graph
    vector<Vertex<NetworkPoint> *> centroids;
    // Add code to randomly select initial centroids
    // ...

    // Initialize the maximum number of vertices per cluster
    const int max_vertices_per_cluster = ceil(graph.getNumVertex() / static_cast<double>(max_clusters));

    // Iterate until all vertices are assigned to clusters
    while (!centroids.empty()) {
        // Create a new cluster
        Cluster cluster;

        // Randomly select a centroid
        Vertex<NetworkPoint> *centroid = centroids.back();
        centroids.pop_back();
        cluster.vertices.push_back(centroid);

        // Assign vertices to the current cluster until it reaches the maximum number of vertices
        while (cluster.vertices.size() < max_vertices_per_cluster) {
            // Find the nearest unassigned vertex to the current centroid
            double min_distance = numeric_limits<double>::infinity();
            Vertex<NetworkPoint> *nearest_vertex = nullptr;
            for (const auto &entry: graph.getVertexSet()) {
                Vertex<NetworkPoint> *vertex = entry.second;
                // Check if the vertex is unassigned
                if (find(cluster.vertices.begin(), cluster.vertices.end(), vertex) == cluster.vertices.end()) {
                    auto v1 = centroid->getInfo();
                    auto v2 = vertex->getInfo();
                    double distance = haversine(v1.getLat(), v1.getLon(),
                                                v2.getLat(), v2.getLon());
                    if (distance < min_distance) {
                        min_distance = distance;
                        nearest_vertex = vertex;
                    }
                }
            }

            // Add the nearest unassigned vertex to the current cluster
            if (nearest_vertex != nullptr) {
                cluster.vertices.push_back(nearest_vertex);
            }
        }

        // Add the current cluster to the list of clusters
        clusters.push_back(cluster);
    }

    // Assign remaining vertices to the nearest clusters
    for (const auto &entry: graph.getVertexSet()) {
        Vertex<NetworkPoint> *vertex = entry.second;
        // Check if the vertex is unassigned
        if (none_of(clusters.begin(), clusters.end(), [&](const Cluster &cluster) {
            return find(cluster.vertices.begin(), cluster.vertices.end(), vertex) != cluster.vertices.end();
        })) {
            // Find the cluster whose centroid is closest to the vertex
            double min_distance = numeric_limits<double>::infinity();
            Cluster *nearest_cluster = nullptr;
            for (auto &cluster: clusters) {
                auto v1 = cluster.vertices[0]->getInfo();
                auto v2 = vertex->getInfo();
                double distance = haversine(v1.getLat(), v1.getLon(),
                                            v2.getLat(), v2.getLon()); // Assuming the first vertex is the centroid
                if (distance < min_distance) {
                    min_distance = distance;
                    nearest_cluster = &cluster;
                }
            }
            // Add the vertex to the nearest cluster
            if (nearest_cluster != nullptr) {
                nearest_cluster->vertices.push_back(vertex);
            }
        }
    }

    return clusters;
}*/


double Graph::tspHeuristic(std::vector<unsigned int> &path) const {
    path.clear();
    // nº of clusters
    unsigned int k = (int) round(std::sqrt(getNumVertex()));

    std::vector<std::vector<unsigned int>> clusters;
    clustering(clusters, k);

    std::vector<unsigned int> finalPath;
    finalPath.push_back(0);

    std::unordered_set<unsigned int> visited;
    visited.insert(0);

    /*for (const auto &cluster: clusters) {
        std::vector<unsigned int> clusterPath;
        double clusterDist = solveClusterTSP(cluster, clusterPath);
        for (unsigned int node: clusterPath) {
            if (visited.find(node) == visited.end()) {
                finalPath.push_back(node);
                visited.insert(node);
            }
        }
    }

    finalPath.push_back(0); // Ending at node 0
    path = finalPath;

    // calculate total distance of the final path
    double totalDist = 0.0;
    for (unsigned int i = 0; i < finalPath.size() - 1; ++i) {
        totalDist += findVertex(finalPath[i])->getEdge(findVertex(finalPath[i + 1]))->getWeight();
    }

    return totalDist;*/

    double totalDist = 0.0;
    unsigned int current = 0;

    for (const auto &cluster: clusters) {
        if (cluster.empty()) continue;

        std::vector<unsigned int> clusterPath;
        double clusterDist = solveClusterTSP(cluster, clusterPath);

        // find nearest neighbour of the current node
        unsigned int nearestInCluster = clusterPath[0];
        double minDist = std::numeric_limits<double>::max();
        for (unsigned int node: clusterPath) {
            if (auto edge = findVertex(current)->getEdge(findVertex(node))) {
                double dist = edge->getWeight();
                if (dist < minDist) {
                    minDist = dist;
                    nearestInCluster = node;
                }
            }
        }

        if (minDist == std::numeric_limits<double>::max()) continue; // No valid edge found

        // move to the nearest node in the current cluster
        totalDist += minDist;
        current = nearestInCluster;
        finalPath.push_back(current);
        visited.insert(current);

        // traverse the cluster path
        for (unsigned int node: clusterPath) {
            if (visited.find(node) == visited.end()) {
                totalDist += findVertex(current)->getEdge(findVertex(node))->getWeight();
                current = node;
                finalPath.push_back(current);
                visited.insert(current);
            }
        }
        totalDist += clusterDist;
    }

    if (auto edge = findVertex(current)->getEdge(findVertex(0))) {
        totalDist += edge->getWeight();
        finalPath.push_back(0);
    } else {
        std::cout << "Warning: No edge found returning to the start" << std::endl;
    }

    // return to the starting node 0
    //totalDist += findVertex(current)->getEdge(findVertex(0))->getWeight();
    //finalPath.push_back(0);

    path = finalPath;
    return totalDist;
}

void Graph::clustering(std::vector<std::vector<unsigned int>> &clusters, unsigned int k) const {
    std::vector<std::pair<double, double>> centroids(k);

    for (unsigned int i = 0; i < k; ++i) {
        centroids[i].first = findVertex(i)->getInfo().getLat();
        centroids[i].second = findVertex(i)->getInfo().getLon();
    }

    bool changed = true;

    while (changed) {
        changed = false;
        clusters.clear();
        clusters.resize(k);

        // assign vertices to the nearest centroid
        for (unsigned int i = 0; i < getNumVertex(); ++i) {
            NetworkPoint point = findVertex(i)->getInfo();
            unsigned int bestCluster = 0;
            double bestDist = std::numeric_limits<double>::max();

            for (unsigned int j = 0; j < k; ++j) {
                double dist = haversine(point.getLat(), point.getLon(), centroids[j].first,
                                        centroids[j].second);
                if (dist < bestDist) {
                    bestDist = dist;
                    bestCluster = j;
                }
            }
            clusters[bestCluster].push_back(i);
        }

        // update centroids
        for (unsigned int j = 0; j < k; ++j) {
            double lat = 0.0, lon = 0.0;
            for (unsigned int vertex: clusters[j]) {
                NetworkPoint point = findVertex(vertex)->getInfo();
                lat += point.getLat();
                lon += point.getLon();
            }

            lat /= (double) clusters[j].size();
            lon /= (double) clusters[j].size();

            if (haversine(lat, lon, centroids[j].first,
                          centroids[j].second) > 1e-5) {
                centroids[j].first = lat;
                centroids[j].second = lon;
                changed = true;
            }
        }
    }
}

double Graph::solveClusterTSP(const std::vector<unsigned int> &cluster, std::vector<unsigned int> &clusterPath) const {
    if (cluster.size() <= 1) {
        if (!cluster.empty())
            clusterPath.push_back(cluster[0]);
        return 0.0;
    }

    unsigned int current = cluster[0];
    std::unordered_set<unsigned int> unvisited(cluster.begin() + 1, cluster.end());
    //unvisited.erase(current);
    clusterPath.push_back(current);

    double totalDist = 0.0;

    while (!unvisited.empty()) {
        unsigned int next = getNearestVertex(current, unvisited);
        //std::cout << "Current: " << current << ", Next: " << next << std::endl;

        if (auto edge = findVertex(current)->getEdge(findVertex(next))) {
            totalDist += edge->getWeight();
            current = next;
            unvisited.erase(current);
            clusterPath.push_back(current);
        } else {
            std::cerr << "Warning: No edge found between vertices " << current << " and " << next << std::endl;
            unvisited.erase(current);
        }
    }

    return totalDist;
}

unsigned int Graph::getNearestVertex(unsigned int from, const std::unordered_set<unsigned int> &candidates) const {
    unsigned int nearest = *candidates.begin();
    double bestDist = std::numeric_limits<double>::max();

    for (unsigned int candidate: candidates) {
        if (auto edge = findVertex(from)->getEdge(findVertex(candidate))) {
            double dist = edge->getWeight();
            if (dist < bestDist) {
                bestDist = dist;
                nearest = candidate;
            }
        }
    }

    return nearest;
}