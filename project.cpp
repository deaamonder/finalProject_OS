#include <iostream>
#include <vector>
#include <list>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <string>
#include <sstream>
#include <thread>
#include <future>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <atomic>
#include <climits>
using namespace std;



class Graph {
private:
    unordered_map<int, list<pair<int, int>>> adjList;

public:
    void addEdge(int u, int v, int weight, bool isDirected = false) {
        adjList[u].push_back({v, weight});
        if (!isDirected) {
            adjList[v].push_back({u, weight});
        }
    }

    void printGraph() {
        for (auto& vertex : adjList) {
            cout << "Vertex " << vertex.first << ":";
            for (auto& neighbor : vertex.second) {
                cout << " -> (" << neighbor.first << ", " << neighbor.second << ")";
            }
            cout << endl;
        }
    }

    const unordered_map<int, list<pair<int, int>>>& getAdjList() const {
        return adjList;
    }
};

class MSTAlgorithm {
public:
    virtual void solve(Graph& graph) = 0;
    virtual ~MSTAlgorithm() = default;
};

class PrimMST : public MSTAlgorithm {
public:
    void solve(Graph& graph) override {
        cout << "Running Prim's Algorithm for MST...\n";

        const unordered_map<int, list<pair<int, int>>>& adjList = graph.getAdjList();

        if (adjList.empty()) {
            cout << "Graph is empty, no MST possible." << endl;
            return;
        }

        unordered_map<int, bool> inMST;
        unordered_map<int, int> key;
        unordered_map<int, int> parent;

        for (const auto& vertex : adjList) {
            key[vertex.first] = INT_MAX;
            inMST[vertex.first] = false;
        }

        priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;

        int startVertex = adjList.begin()->first;
        key[startVertex] = 0;
        parent[startVertex] = -1;
        pq.push({0, startVertex});

        while (!pq.empty()) {
            int u = pq.top().second;
            pq.pop();

            inMST[u] = true;

            for (const auto& neighbor : adjList.at(u)) {
                int v = neighbor.first;
                int weight = neighbor.second;

                if (!inMST[v] && key[v] > weight) {
                    key[v] = weight;
                    parent[v] = u;
                    pq.push({key[v], v});
                }
            }
        }

        int totalWeight = 0;
        cout << "Minimum Spanning Tree edges:\n";
        for (const auto& vertex : adjList) {
            int v = vertex.first;
            if (parent[v] != -1) {
                cout << "Edge: (" << parent[v] << " - " << v << ") with weight " << key[v] << "\n";
                totalWeight += key[v];
            }
        }
        cout << "Total weight of MST: " << totalWeight << endl;
    }
};

class KruskalMST : public MSTAlgorithm {
private:
    int findParent(int vertex, vector<int>& parent) {
        if (parent[vertex] != vertex) {
            parent[vertex] = findParent(parent[vertex], parent);
        }
        return parent[vertex];
    }

    void unionSets(int u, int v, vector<int>& parent, vector<int>& rank) {
        int rootU = findParent(u, parent);
        int rootV = findParent(v, parent);

        if (rootU != rootV) {
            if (rank[rootU] > rank[rootV]) {
                parent[rootV] = rootU;
            } else if (rank[rootU] < rank[rootV]) {
                parent[rootU] = rootV;
            } else {
                parent[rootV] = rootU;
                rank[rootU]++;
            }
        }
    }

public:
    void solve(Graph& graph) override {
        cout << "Running Kruskal's Algorithm for MST...\n";

        const unordered_map<int, list<pair<int, int>>>& adjList = graph.getAdjList();
        int numVertices = adjList.size();
        if (numVertices == 0) {
            cout << "Graph is empty, no MST possible." << endl;
            return;
        }

        vector<pair<int, pair<int, int>>> edges;
        for (const auto& vertex : adjList) {
            int u = vertex.first;
            for (const auto& neighbor : vertex.second) {
                int v = neighbor.first;
                int weight = neighbor.second;
                if (u < v) {
                    edges.push_back({weight, {u, v}});
                }
            }
        }

        sort(edges.begin(), edges.end());
        vector<int> parent(numVertices);
        vector<int> rank(numVertices, 0);
        for (int i = 0; i < numVertices; ++i) {
            parent[i] = i;
        }

        vector<pair<int, int>> mstEdges;
        int totalWeight = 0;

        for (const auto& edge : edges) {
            int weight = edge.first;
            int u = edge.second.first;
            int v = edge.second.second;

            int rootU = findParent(u, parent);
            int rootV = findParent(v, parent);

            if (rootU != rootV) {
                mstEdges.push_back({u, v});
                totalWeight += weight;
                unionSets(rootU, rootV, parent, rank);
            }
        }

        cout << "Minimum Spanning Tree edges:\n";
        for (const auto& edge : mstEdges) {
            cout << "Edge: (" << edge.first << " - " << edge.second << ")\n";
        }
        cout << "Total weight of MST: " << totalWeight << endl;
    }
};


class BoruvkaMST : public MSTAlgorithm {
private:
    int findParent(int vertex, vector<int>& parent) {
        if (parent[vertex] != vertex) {
            parent[vertex] = findParent(parent[vertex], parent);
        }
        return parent[vertex];
    }

    void unionSets(int u, int v, vector<int>& parent, vector<int>& rank) {
        int rootU = findParent(u, parent);
        int rootV = findParent(v, parent);

        if (rootU != rootV) {
            if (rank[rootU] > rank[rootV]) {
                parent[rootV] = rootU;
            } else if (rank[rootU] < rank[rootV]) {
                parent[rootU] = rootV;
            } else {
                parent[rootV] = rootU;
                rank[rootU]++;
            }
        }
    }

public:
    void solve(Graph& graph) override {
        cout << "Running Boruvka's Algorithm for MST...\n";

        const unordered_map<int, list<pair<int, int>>>& adjList = graph.getAdjList();
        int numVertices = adjList.size();

        if (numVertices == 0) {
            cout << "Graph is empty, no MST possible." << endl;
            return;
        }

        vector<int> parent(numVertices);
        vector<int> rank(numVertices, 0);

        for (int i = 0; i < numVertices; ++i) {
            parent[i] = i;
        }

        vector<pair<int, pair<int, int>>> mstEdges;
        int totalWeight = 0;

        int numComponents = numVertices;

        while (numComponents > 1) {
            vector<pair<int, int>> cheapest(numVertices, {-1, -1});

            for (const auto& vertex : adjList) {
                int u = vertex.first;
                for (const auto& neighbor : vertex.second) {
                    int v = neighbor.first;
                    int weight = neighbor.second;

                    int setU = findParent(u, parent);
                    int setV = findParent(v, parent);

                    if (setU != setV) {
                        if (cheapest[setU].first == -1 || cheapest[setU].first > weight) {
                            cheapest[setU] = {weight, v};
                        }
                        if (cheapest[setV].first == -1 || cheapest[setV].first > weight) {
                            cheapest[setV] = {weight, u};
                        }
                    }
                }
            }

            for (int u = 0; u < numVertices; ++u) {
                if (cheapest[u].first != -1) {
                    int v = cheapest[u].second;
                    int weight = cheapest[u].first;

                    int setU = findParent(u, parent);
                    int setV = findParent(v, parent);

                    if (setU != setV) {
                        unionSets(setU, setV, parent, rank);
                        mstEdges.push_back({u, {v, weight}});
                        totalWeight += weight;
                        numComponents--; 
                    }
                }
            }
        }

        cout << "Minimum Spanning Tree edges:\n";
        for (const auto& edge : mstEdges) {
            cout << "Edge: (" << edge.first << " - " << edge.second.first << ") with weight " << edge.second.second << "\n";
        }
        cout << "Total weight of MST: " << totalWeight << endl;
    }
};



class TarjanMST : public MSTAlgorithm {
private:
    void dfs(int v, int parent, vector<int>& low, vector<int>& disc, vector<bool>& inMST, vector<pair<int, int>>& mstEdges, int& time, const unordered_map<int, list<pair<int, int>>>& adjList, int& totalWeight) {
        disc[v] = low[v] = ++time;
        inMST[v] = true;

        for (const auto& neighbor : adjList.at(v)) {
            int u = neighbor.first;
            int weight = neighbor.second;

            if (!inMST[u]) { 
                mstEdges.push_back({v, u});
                totalWeight += weight;
                dfs(u, v, low, disc, inMST, mstEdges, time, adjList, totalWeight);
                low[v] = min(low[v], low[u]);
            } else if (u != parent) {
                low[v] = min(low[v], disc[u]);
            }
        }
    }

public:
    void solve(Graph& graph) override {
        cout << "Running Tarjan's Algorithm for MST...\n";

        const unordered_map<int, list<pair<int, int>>>& adjList = graph.getAdjList();
        int numVertices = adjList.size();

        if (numVertices == 0) {
            cout << "Graph is empty, no MST possible." << endl;
            return;
        }

        vector<int> low(numVertices, -1);
        vector<int> disc(numVertices, -1); 
        vector<bool> inMST(numVertices, false); 
        vector<pair<int, int>> mstEdges; 
        int time = 0;
        int totalWeight = 0; 

   
        for (int i = 0; i < numVertices; ++i) {
            if (!inMST[i]) {
                dfs(i, -1, low, disc, inMST, mstEdges, time, adjList, totalWeight);
            }
        }


        cout << "Minimum Spanning Tree edges:\n";
        for (const auto& edge : mstEdges) {
            cout << "Edge: (" << edge.first << " - " << edge.second << ")\n";
        }
        cout << "Total weight of MST: " << totalWeight << endl;
    }
};


class IntegerMST : public MSTAlgorithm {
public:
    void solve(Graph& graph) override {
        cout << "Running Integer MST Algorithm...\n";

        const unordered_map<int, list<pair<int, int>>>& adjList = graph.getAdjList();
        int numVertices = adjList.size();

        if (numVertices == 0) {
            cout << "Graph is empty, no MST possible." << endl;
            return;
        }


        unordered_map<int, bool> inMST;
        unordered_map<int, int> key;
        unordered_map<int, int> parent;
        int totalWeight = 0; 

        for (const auto& vertex : adjList) {
            key[vertex.first] = INT_MAX;
            inMST[vertex.first] = false;
        }


        priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;

        int startVertex = adjList.begin()->first;
        key[startVertex] = 0;
        parent[startVertex] = -1;
        pq.push({0, startVertex});

        while (!pq.empty()) {
            int u = pq.top().second;
            pq.pop();

            inMST[u] = true;

            for (const auto& neighbor : adjList.at(u)) {
                int v = neighbor.first;
                int weight = neighbor.second;

                if (!inMST[v] && key[v] > weight) {
                    key[v] = weight;
                    parent[v] = u;
                    pq.push({key[v], v});
                }
            }
        }


        cout << "Minimum Spanning Tree edges:\n";
        for (const auto& vertex : adjList) {
            int v = vertex.first;
            if (parent[v] != -1) {
                cout << "Edge: (" << parent[v] << " - " << v << ")\n";
                totalWeight += key[v];
            }
        }

        cout << "Total weight of MST: " << totalWeight << endl;
    }
};


class MSTFactory {
public:
    static MSTAlgorithm* getMSTSolver(const string& algorithm) {
        if (algorithm == "Prim") {
            return new PrimMST();
        } else if (algorithm == "Kruskal") {
            return new KruskalMST();
        } else if (algorithm == "Boruvka") {
            return new BoruvkaMST();
        } else if (algorithm == "Tarjan") {
            return new TarjanMST();
        } else if (algorithm == "Integer") {
            return new IntegerMST();
        } else {
            cout << "Unknown algorithm! Returning Prim's Algorithm by default.\n";
            return new PrimMST();
        }
    }
};


class MSTRequest {
public:
    string command;
    Graph graph;

    MSTRequest(const string& cmd) : command(cmd) {}
};


class ThreadPool {
private:
    vector<thread> workers;
    queue<function<void()>> tasks;
    mutex queueMutex;
    condition_variable condition;
    atomic<bool> stop;

public:
    ThreadPool(size_t numThreads) : stop(false) {
        for (size_t i = 0; i < numThreads; ++i) {
            workers.emplace_back([this] {
                while (true) {
                    function<void()> task;
                    {
                        unique_lock<mutex> lock(this->queueMutex);
                        this->condition.wait(lock, [this] {
                            return this->stop || !this->tasks.empty();
                        });
                        if (this->stop && this->tasks.empty()) return;
                        task = move(this->tasks.front());
                        this->tasks.pop();
                    }
                    task();
                }
            });
        }
    }

    ~ThreadPool() {
        {
            unique_lock<mutex> lock(queueMutex);
            stop = true;
        }
        condition.notify_all();
        for (thread &worker: workers) {
            worker.join();
        }
    }

    template<class F>
    void enqueue(F&& f) {
        {
            unique_lock<mutex> lock(queueMutex);
            tasks.emplace(forward<F>(f));
        }
        condition.notify_one();
    }
};


class MSTServer {
private:
    Graph graph;
    ThreadPool pool;

public:
    MSTServer(size_t numThreads) : pool(numThreads) {}

    void processRequest(const string& request) {
        stringstream ss(request);
        string command;
        ss >> command;

        if (command == "addEdge") {
            int u, v, weight;
            ss >> u >> v >> weight;
            graph.addEdge(u, v, weight);
            cout << "Added edge (" << u << ", " << v << ") with weight " << weight << endl;
        } else if (command == "printGraph") {
            graph.printGraph();
        } else if (command == "solveMST") {
            string algorithm;
            ss >> algorithm;
            pool.enqueue([this, algorithm]() {
                MSTAlgorithm* solver = MSTFactory::getMSTSolver(algorithm);
                solver->solve(graph);
                delete solver; 
            });
        } else {
            //cout << "Unknown command: " << command << endl;
        }
    }
};


int main() {
    MSTServer server(4); 


    string command;

    while (true) {
        cout << "\nEnter command (addEdge, printGraph, solveMST <algorithm>, or exit): ";
        getline(cin, command);

        if (command == "exit") {
            break;
        }

        server.processRequest(command);
    }

    return 0;
}