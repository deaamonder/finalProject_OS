#include <iostream>
#include <vector>
#include <queue>
#include <mutex>
#include <thread>
#include <condition_variable>
#include <functional>
#include <cmath>
#include <algorithm>

class Edge {
public:
    int src, dest, weight;
    Edge(int s, int d, int w) : src(s), dest(d), weight(w) {}
};

class Graph {
private:
    int vertices;
    std::vector<std::vector<Edge>> adjList;

public:
    Graph(int v) : vertices(v), adjList(v) {}

    void addEdge(int src, int dest, int weight) {
        adjList[src].emplace_back(src, dest, weight);
        adjList[dest].emplace_back(dest, src, weight);
    }

    const std::vector<Edge>& getEdges(int vertex) const {
        return adjList[vertex];
    }

    int numVertices() const {
        return vertices;
    }
};

class MSTMeasurements {
public:
    static int totalWeight(const std::vector<Edge>& mstEdges) {
        int weight = 0;
        for (const auto& edge : mstEdges) {
            weight += edge.weight;
        }
        return weight;
    }

    static int longestDistance(const Graph& mst) {
        int maxDistance = 0;
        for (int i = 0; i < mst.numVertices(); ++i) {
            maxDistance = std::max(maxDistance, bfsLongestPath(mst, i));
        }
        return maxDistance;
    }

    static double averageDistance(const Graph& mst) {
        int totalDistance = 0;
        int pairCount = 0;
        for (int i = 0; i < mst.numVertices(); ++i) {
            for (int j = i + 1; j < mst.numVertices(); ++j) {
                totalDistance += bfsShortestPath(mst, i, j);
                ++pairCount;
            }
        }
        return pairCount > 0 ? (double)totalDistance / pairCount : 0.0;
    }

    static int shortestDistance(const Graph& mst, int src, int dest) {
        return bfsShortestPath(mst, src, dest);
    }

private:
    static int bfsLongestPath(const Graph& mst, int start) {
        std::queue<std::pair<int, int>> q;
        std::vector<bool> visited(mst.numVertices(), false);
        q.push({start, 0});
        int maxDistance = 0;

        while (!q.empty()) {
            auto [vertex, dist] = q.front();
            q.pop();
            if (visited[vertex]) continue;
            visited[vertex] = true;
            maxDistance = std::max(maxDistance, dist);

            for (const auto& edge : mst.getEdges(vertex)) {
                if (!visited[edge.dest]) {
                    q.push({edge.dest, dist + edge.weight});
                }
            }
        }
        return maxDistance;
    }

    static int bfsShortestPath(const Graph& mst, int src, int dest) {
        std::queue<std::pair<int, int>> q;
        std::vector<bool> visited(mst.numVertices(), false);
        q.push({src, 0});

        while (!q.empty()) {
            auto [vertex, dist] = q.front();
            q.pop();
            if (vertex == dest) return dist;
            if (visited[vertex]) continue;
            visited[vertex] = true;

            for (const auto& edge : mst.getEdges(vertex)) {
                if (!visited[edge.dest]) {
                    q.push({edge.dest, dist + edge.weight});
                }
            }
        }
        return -1; // No path found
    }
};

class ActiveObject {
public:
    ActiveObject() : running(true), worker(&ActiveObject::process, this) {}

    ~ActiveObject() {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            running = false;
        }
        condVar.notify_all();
        if (worker.joinable()) {
            worker.join();
        }
    }

    void addTask(std::function<void()> task) {
        std::lock_guard<std::mutex> lock(mutex_);
        tasks.push(task);
        condVar.notify_one();
    }

private:
    void process() {
        while (true) {
            std::function<void()> task;
            {
                std::unique_lock<std::mutex> lock(mutex_);
                condVar.wait(lock, [this] { return !tasks.empty() || !running; });
                if (!running && tasks.empty()) {
                    break;
                }
                task = tasks.front();
                tasks.pop();
            }
            task();
        }
    }

    std::queue<std::function<void()>> tasks;
    std::mutex mutex_;
    std::condition_variable condVar;
    bool running;
    std::thread worker;
};

class ThreadPool {
private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    std::mutex mutex_;
    std::condition_variable condVar;
    bool shutdown = false;

public:
    ThreadPool(size_t threads) {
        for (size_t i = 0; i < threads; ++i) {
            workers.emplace_back([this] {
                while (true) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(mutex_);
                        condVar.wait(lock, [this] { return shutdown || !tasks.empty(); });
                        if (shutdown && tasks.empty()) return;
                        task = std::move(tasks.front());
                        tasks.pop();
                    }
                    task();
                }
            });
        }
    }

    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(mutex_);
            shutdown = true;
        }
        condVar.notify_all();
        for (std::thread& worker : workers) {
            if (worker.joinable()) {
                worker.join();
            }
        }
    }

    void enqueue(std::function<void()> task) {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            tasks.push(std::move(task));
        }
        condVar.notify_one();
    }

    

};

int main() {
    int vertices = 5; 
    Graph graph(vertices);

    while (true) {
        std::cout << "Enter command (addEdge, printGraph, solveMST <algorithm>, or exit): ";
        std::string command;
        std::cin >> command;

        if (command == "addEdge") {
            int src, dest, weight;
            std::cin >> src >> dest >> weight;
            graph.addEdge(src, dest, weight);
            std::cout << "Edge added.\n";
        } else if (command == "printGraph") {
            for (int i = 0; i < graph.numVertices(); ++i) {
                std::cout << "Edges from vertex " << i << ":\n";
                for (const auto& edge : graph.getEdges(i)) {
                    std::cout << "  to " << edge.dest << " with weight " << edge.weight << '\n';
                }
            }
        } else if (command == "solveMST") {
            std::string algorithm;
            std::cin >> algorithm;

            
            std::vector<Edge> mstEdges = {/* קשתות MST כתוצאה מאלגוריתם */};
            int totalWeight = MSTMeasurements::totalWeight(mstEdges);
            std::cout << "Total Weight of MST: " << totalWeight << '\n';

            
        } else if (command == "exit") {
            break;
        } else {
            std::cout << "Invalid command.\n";
        }
    }

    return 0;
}
