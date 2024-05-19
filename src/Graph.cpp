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
void Edge<T>::setReverse(Edge<T> *reverse) {
    this->reverse = reverse;
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

Graph::~Graph() {
    //deleteMatrix(distanceMatrix_, vertexSet.size());
    /*if (!vertexSet.empty()) {
        for (auto &v: vertexSet) {
            removeVertex(v.second->getInfo());
        }
        vertexSet.clear();
    }*/
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

bool Graph::isConnected() const {
    if (vertexSet.empty()) return true;

    std::unordered_set<unsigned> visited;
    std::queue<unsigned> q;

    auto startVertex = vertexSet.begin()->first;
    q.push(startVertex);
    visited.insert(startVertex);

    while (!q.empty()) {
        auto u = q.front();
        q.pop();
        for (const auto &e: vertexSet.at(u)->getAdj()) {
            auto v = e->getDest()->getInfo().getId();
            if (visited.find(v) == visited.end()) {
                visited.insert(v);
                q.push(v);
            }
        }
    }

    return visited.size() == vertexSet.size();
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

double Graph::getEdgeWeight(const NetworkPoint &a, const NetworkPoint &b) const {
    auto vertex = findVertex(a);
    for (const auto &edge: vertex->getAdj()) {
        if (edge->getDest()->getInfo() == b) {
            return edge->getWeight();
        }
    }
    return 0;
}

struct EdgeComparator {
    bool operator()(const Edge<NetworkPoint> *lhs, const Edge<NetworkPoint> *rhs) const {
        return lhs->getWeight() > rhs->getWeight();
    }
};

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

std::vector<NetworkPoint> Graph::aproxTSP() {
    auto root = findVertex(NetworkPoint(0));
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

    tour.push_back(tour.at(0));

    return tour;
}

double Graph::calculateTriangular(std::vector<NetworkPoint> g) {
    double total = 0.0;

    for (auto i = 0; i < g.size() - 1; i++) {
        total += getEdgeWeight(g.at(i), g.at(i + 1));
    }
    total += getEdgeWeight(g.at(g.size() - 1), g.at(0));

    return total;
}

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

double Graph::tspHeuristic(std::vector<unsigned int> &path) const {
    std::unordered_set<unsigned int> visited;

    unsigned int current = 0;
    path.push_back(current);
    visited.insert(current);
    double totalDist = 0.0;

    while (path.size() < vertexSet.size()) {
        double minDist = std::numeric_limits<double>::max();
        unsigned int next = 0;
        for (const auto &edge: vertexSet.at(current)->getAdj()) {
            auto neighbor = edge->getDest()->getInfo().getId();
            if (visited.find(neighbor) == visited.end() && edge->getWeight() < minDist) {
                minDist = edge->getWeight();
                next = neighbor;
            }
        }
        path.push_back(next);
        visited.insert(next);
        totalDist += minDist;
        current = next;
    }

    totalDist += vertexSet.at(current)->getEdge(vertexSet.at(path[0]))->getWeight();
    path.push_back(path[0]);

    return totalDist;
    /*
    path.clear();
    // nº de clusters
    unsigned int k = (int) round(std::sqrt(getNumVertex()));

    std::vector<std::vector<unsigned int>> clusters;
    clustering(clusters, k);

    std::vector<unsigned int> finalPath;
    finalPath.push_back(0);

    std::unordered_set<unsigned int> visited;
    visited.insert(0);

    for (const auto &cluster: clusters) {
        std::vector<unsigned int> clusterPath;
        solveClusterTSP(cluster, clusterPath);
        for (unsigned int node: clusterPath) {
            if (visited.find(node) == visited.end()) {
                finalPath.push_back(node);
                visited.insert(node);
            }
        }
    }

    finalPath.push_back(0); // terminar no node 0 (inicial)
    path = finalPath;

    double totalDist = 0.0;
    for (unsigned int i = 0; i < finalPath.size() - 1; i++) {
        totalDist += findVertex(finalPath[i])->getEdge(findVertex(finalPath[i + 1]))->getWeight();
    }

    return totalDist;*/
}

void Graph::clustering(std::vector<std::vector<unsigned int>> &clusters, unsigned int k) const {
    std::vector<std::pair<double, double>> centroids(k);

    for (unsigned int i = 0; i < k; i++) {
        centroids[i].first = findVertex(i)->getInfo().getLat();
        centroids[i].second = findVertex(i)->getInfo().getLon();
    }

    bool changed = true;

    while (changed) {
        changed = false;
        clusters.clear();
        clusters.resize(k);

        // assignar vértices ao centróide mais perto
        for (unsigned int i = 0; i < getNumVertex(); i++) {
            NetworkPoint point = findVertex(i)->getInfo();
            unsigned int bestCluster = 0;
            double bestDist = std::numeric_limits<double>::max();

            for (unsigned int j = 0; j < k; j++) {
                double dist = haversine(point.getLat(), point.getLon(), centroids[j].first,
                                        centroids[j].second);
                if (dist < bestDist) {
                    bestDist = dist;
                    bestCluster = j;
                }
            }
            clusters[bestCluster].push_back(i);
        }

        // atualizar centróides
        for (unsigned int j = 0; j < k; j++) {
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
    clusterPath.push_back(current);

    double totalDist = 0.0;

    while (!unvisited.empty()) {
        unsigned int next = getNearestVertex(current, unvisited);

        if (auto edge = findVertex(current)->getEdge(findVertex(next))) {
            totalDist += edge->getWeight();
            current = next;
            unvisited.erase(current);
            clusterPath.push_back(current);
        } else {
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

double Graph::tspRealWorld(std::vector<unsigned int> &path) const {
    if (!isConnected()) {
        // grafo não é conectado
        return -1;
    }

    std::unordered_set<unsigned int> visited;

    //unsigned int current = vertexSet.begin()->first;
    unsigned int current = path.at(0);
    //path.push_back(current);
    visited.insert(current);
    double totalDist = 0.0;

    while (path.size() < vertexSet.size()) {
        double minDist = std::numeric_limits<double>::max();
        unsigned int next = 0;
        for (const auto &edge: vertexSet.at(current)->getAdj()) {
            auto neighbor = edge->getDest()->getInfo().getId();
            if (visited.find(neighbor) == visited.end() && edge->getWeight() < minDist) {
                minDist = edge->getWeight();
                next = neighbor;
            }
        }
        if (next == 0) break; // não há mais vértices para visitar (grafo não é conectado)
        path.push_back(next);
        visited.insert(next);
        totalDist += minDist;
        current = next;
    }

    auto returnEdge = vertexSet.at(current)->getEdge(vertexSet.at(path[0]));
    if (returnEdge) {
        totalDist += returnEdge->getWeight();
        path.push_back(path[0]);
    } else {
        // não há nenhum path que retorna ao vértice inicial
        return -1;
    }

    return totalDist;
}