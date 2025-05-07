#include <iostream>
#include <vector>
#include <queue>
#include <list>
#include <stack>
#include <random> // For better random number generation
#include <limits> // For numeric_limits
#include <utility> // For std::pair

using namespace std;

// Node class representing a router
class Node {
public:
    std::vector<std::pair<int, int>> connections; // Connections to other routers (routerID, latency)

    Node() : connections() {} // Default constructor - initialize connections

    explicit Node(int /*routerId*/) : connections() {} // Constructor - initialize connections, routerId not used

    // [[nodiscard]] int getRouterID() const { return -1; } // Dummy implementation as it's not used

    void addConnection(int routerID, int latency) {
        connections.emplace_back(routerID, latency); // Use emplace_back
    }
};

// Graph class representing the entire network
class Graph {
public:
    int V;  // Number of routers (nodes)
    vector<Node> nodes;
private:
    std::mt19937 generator;
    std::uniform_real_distribution<double> distribution;

public:
    Graph(int V) : V(V), nodes(V), generator(std::random_device{}()), distribution(0.0, 1.0) {
        for (int i = 0; i < V; ++i) {
            nodes[i] = Node(i);  // Initialize nodes
        }
    }

    void addEdge(int u, int v, int latency) {
        nodes[u].addConnection(v, latency);  // Add bidirectional edge
        nodes[v].addConnection(u, latency);
    }

    // Dijkstra’s algorithm to find the shortest path from a source router
    vector<int> dijkstra(int src) {
        vector<int> dist(V, numeric_limits<int>::max()); // Initialize with max value
        dist[src] = 0;
        priority_queue<pair<int, int>, vector<pair<int, int>>, greater<>> pq; // Prefer transparent comparator
        pq.push({0, src});

        while (!pq.empty()) {
            int u = pq.top().second;
            pq.pop();

            for (const auto& neighbor : nodes[u].connections) { // Use const auto&
                int v = neighbor.first;
                int weight = neighbor.second;

                if (dist[u] + weight < dist[v]) {
                    dist[v] = dist[u] + weight;
                    pq.push({dist[v], v});
                }
            }
        }

        return dist; // Return distances to all routers from the source
    }

    [[nodiscard]] bool simulatePacketLoss(double lossProbability) {
        return distribution(generator) < lossProbability;
    }
};

// Linked list to handle packet retransmission (FIFO)
class PacketList {
public:
    list<int> packetList;

    void addPacket(int packetID) {
        packetList.emplace_back(packetID);
    }

    void printList() const {
        for (int packet : packetList) {
            cout << "Packet ID: " << packet << " ";
        }
        cout << endl;
    }
};

// Queue class to manage packet transmission (FIFO)
class PacketQueue {
public:
    queue<int> packetQueue;

    void enqueue(int packetID) {
        packetQueue.push(packetID);
    }

    void dequeue() {
        if (!packetQueue.empty()) {
            packetQueue.pop();
        }
    }

    void printQueue() const {
        queue<int> tempQueue = packetQueue;
        while (!tempQueue.empty()) {
            cout << "Packet ID: " << tempQueue.front() << " ";
            tempQueue.pop();
        }
        cout << endl;
    }
};

// Stack class to manage packet processing (LIFO)
class PacketStack {
public:
    stack<int> packetStack;

    void push(int packetID) {
        packetStack.push(packetID);
    }

    void pop() {
        if (!packetStack.empty()) {
            packetStack.pop();
        }
    }

    void printStack() const {
        stack<int> tempStack = packetStack;
        while (!tempStack.empty()) {
            cout << "Packet ID: " << tempStack.top() << " ";
            tempStack.pop();
        }
        cout << endl;
    }
};

// Selection Sort algorithm
void selectionSort(vector<int>& arr) {
    for (size_t i = 0; i < arr.size() - 1; ++i) {
        size_t minIdx = i;
        for (size_t j = i + 1; j < arr.size(); ++j) {
            if (arr[j] < arr[minIdx]) {
                minIdx = j;
            }
        }
        swap(arr[i], arr[minIdx]);
    }
}

// Bubble Sort algorithm
void bubbleSort(vector<int>& arr) {
    size_t n = arr.size();
    for (size_t i = 0; i < n - 1; ++i) {
        for (size_t j = 0; j < n - i - 1; ++j) {
            if (arr[j] > arr[j + 1]) {
                swap(arr[j], arr[j + 1]);
            }
        }
    }
}

int main() {
    // Random number generation is now handled by <random> in the Graph class

    // Create a graph representing the network of routers
    int V = 6;  // Number of routers
    Graph g(V);

    // Add connections (edges) between routers with latency
    g.addEdge(0, 1, 7);
    g.addEdge(0, 2, 9);
    g.addEdge(0, 5, 14);
    g.addEdge(1, 2, 10);
    g.addEdge(1, 4, 15);
    g.addEdge(2, 3, 11);
    g.addEdge(2, 5, 2);
    g.addEdge(3, 4, 6);
    g.addEdge(4, 5, 9);

    // Simulate packet transmission with a 10% loss probability
    PacketQueue packetQueue;
    PacketStack packetStack;
    PacketList retransmissionList;  // List to store lost packets

    // Enqueue some packets for transmission
    packetQueue.enqueue(1);
    packetQueue.enqueue(2);
    packetQueue.enqueue(3);

    // Stack-based packet processing (LIFO)
    packetStack.push(4);
    packetStack.push(5);

    cout << "Packet Queue: ";
    packetQueue.printQueue();

    cout << "Packet Stack (LIFO): ";
    packetStack.printStack();

    // Packet loss simulation with 10% probability
    double lossProbability = 0.1;

    while (!packetQueue.packetQueue.empty()) {
        int packet = packetQueue.packetQueue.front();
        packetQueue.dequeue();

        cout << "Processing Packet ID: " << packet << endl;

        // Simulate packet loss
        if (g.simulatePacketLoss(lossProbability)) {
            cout << "Packet ID: " << packet << " lost. Retransmitting..." << endl;
            retransmissionList.addPacket(packet);  // Add lost packet to retransmission list
        } else {
            cout << "Packet ID: " << packet << " transmitted successfully." << endl;
        }
    }

    // Print retransmitted packets
    cout << "\nRetransmitted Packets:" << endl;
    retransmissionList.printList();

    // Apply Dijkstra’s algorithm to find the shortest path from router 0
    cout << "\nShortest Paths from Router 0 (Dijkstra's):" << endl;
    vector<int> distDijkstra = g.dijkstra(0);
    for (size_t i = 0; i < distDijkstra.size(); ++i) {
        cout << "To Router " << i << ": Distance = " << distDijkstra[i] << endl;
    }

    // Sorting router distances (for demonstration)
    vector<int> routerDistances = distDijkstra;
    cout << "\nSorted Router Distances (Selection Sort):" << endl;
    selectionSort(routerDistances);
    for (int dist : routerDistances) {
        cout << dist << " ";
    }
    cout << endl;

    cout << "\nSorted Router Distances (Bubble Sort):" << endl;
    bubbleSort(routerDistances);
    for (int dist : routerDistances) {
        cout << dist << " ";
    }
    cout << endl;

    return 0;
}
