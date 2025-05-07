#include <iostream>
#include <vector>
#include <queue>
#include <list>
#include <stack>
#include <random>
#include <limits>
#include <utility>
#include <algorithm>
#include <chrono>
#include <stdexcept> // For exceptions


using namespace std;

// Node class to represent a router in the network
class Node {
public:
    vector<pair<int, int>> connections; // (neighbor ID, latency)

    Node() : connections() {}

    void addConnection(int routerID, int latency) {
        if (routerID < 0 || latency < 0) {
            throw invalid_argument("Router ID and latency must be non-negative.");
        }
        connections.emplace_back(routerID, latency);
    }
};

// Graph class to represent the network of routers
class Graph {
public:
    int V;       // Number of vertices (routers)
    vector<Node> nodes; // Array of nodes

private:
    mt19937 generator;             // Random number generator
    uniform_real_distribution<> distribution; // Distribution for probabilities

public:
    Graph(int V) : V(V), nodes(V), generator(random_device{}()), distribution(0.0, 1.0) {
        if (V <= 0) {
            throw invalid_argument("Number of vertices must be positive.");
        }
        for (int i = 0; i < V; ++i) {
            nodes[i] = Node(); // Initialize each node
        }
    }

    void addEdge(int u, int v, int latency) {
        if (u < 0 || u >= V || v < 0 || v >= V) {
            throw out_of_range("Vertex indices u and v are out of range.");
        }
        if (latency < 0) {
            throw invalid_argument("Latency must be non-negative.");
        }
        nodes[u].addConnection(v, latency);
        nodes[v].addConnection(u, latency); // Undirected graph
    }

    // Dijkstra's algorithm to find short paths from a source node
    vector<int> dijkstra(int src) {
        if (src < 0 || src >= V) {
            throw out_of_range("Source vertex src is out of range.");
        }

        vector<int> dist(V, numeric_limits<int>::max()); // Initialize distances to infinity
        dist[src] = 0;                                   // Distance from source to itself is 0

        // Priority queue to store (distance, node ID) pairs.  Use greater<> for min-heap.
        priority_queue<pair<int, int>, vector<pair<int, int>>, greater<>> pq;
        pq.push({0, src}); // Start with the source node

        while (!pq.empty()) {
            int u = pq.top().second;
            pq.pop();

            if (u >= nodes.size()) {
                throw runtime_error("Node index out of bounds in Dijkstra's algorithm.");
            }
            for (const auto& neighbor : nodes[u].connections) {
                int v = neighbor.first;
                int weight = neighbor.second;

                if (dist[u] + weight < dist[v]) {
                    dist[v] = dist[u] + weight;
                    pq.push({dist[v], v});
                }
            }
        }
        return dist;
    }

    // Simulate packet loss with a given probability
    bool simulatePacketLoss(double lossProbability) {
        if (lossProbability < 0.0 || lossProbability > 1.0) {
            throw invalid_argument("Loss probability must be between 0 and 1.");
        }
        return distribution(generator) < lossProbability;
    }
};

// Linked list for packet retransmission
class PacketList {
public:
    list<int> packetList;

    void addPacket(int packetID) {
        if (packetID < 0) {
            throw invalid_argument("Packet ID must be non-negative.");
        }
        packetList.push_back(packetID);
    }

    void printList() const {
        if (packetList.empty()) {
            cout << "Retransmission list is empty." << endl;
            return;
        }
        for (int packet : packetList) {
            cout << "Packet ID: " << packet << " ";
        }
        cout << endl;
    }
};

// Queue for packet transmission
class PacketQueue {
public:
    queue<int> packetQueue;

    void enqueue(int packetID) {
        if (packetID < 0) {
            throw invalid_argument("Packet ID must be non-negative.");
        }
        packetQueue.push(packetID);
    }

    void dequeue() {
        if (packetQueue.empty()) {
            // Consider throwing an exception or returning a special value
            // throw runtime_error("Cannot dequeue from an empty queue.");
            return; // Or return a special value like -1 to indicate an error
        }
        packetQueue.pop();
    }

    void printQueue() const {
        if (packetQueue.empty()) {
            cout << "Packet queue is empty." << endl;
            return;
        }
        queue<int> tempQueue = packetQueue; //copy the queue
        while (!tempQueue.empty()) {
            cout << "Packet ID: " << tempQueue.front() << " ";
            tempQueue.pop();
        }
        cout << endl;
    }
};

// Node structure for the binary tree
struct PacketTreeNode {
    int packetID;
    string sourceAddress;
    string destinationAddress;
    chrono::time_point<chrono::system_clock> timestamp;
    PacketTreeNode* left;
    PacketTreeNode* right;

    PacketTreeNode(int packetID, string src, string dest)
            : packetID(packetID), sourceAddress(src), destinationAddress(dest),
              timestamp(chrono::system_clock::now()), left(nullptr), right(nullptr) {
        if(packetID < 0) {
            throw invalid_argument("Packet ID must be non-negative.");
        }
    }
};

// Class for the binary tree, now specifically for network packets
class PacketTree {
public:
    PacketTreeNode* root;

    PacketTree() : root(nullptr) {}

    // Function to insert a packet into the binary tree
    void insertPacket(int packetID, string src, string dest) {
        if (packetID < 0) {
            throw invalid_argument("Packet ID must be non-negative.");
        }
        root = insertPacketNode(root, packetID, src, dest);
    }

    // Function to delete a packet from the binary tree
    void deletePacket(int packetID) {
        if (packetID < 0) {
            throw invalid_argument("Packet ID must be non-negative.");
        }
        root = deletePacketNodeRec(root, packetID);
    }

    // Function to find the minimum value node in a subtree
    PacketTreeNode* findMin(PacketTreeNode* node) {
        if (node == nullptr) {
            return nullptr;
        }
        while (node->left != nullptr)
            node = node->left;
        return node;
    }

    // Function to insert a packet node recursively
    PacketTreeNode* insertPacketNode(PacketTreeNode* node, int packetID, string src, string dest) {
        if (node == nullptr) {
            return new PacketTreeNode(packetID, src, dest);
        }

        if (packetID < node->packetID) {
            node->left = insertPacketNode(node->left, packetID, src, dest);
        } else if (packetID > node->packetID) {
            node->right = insertPacketNode(node->right, packetID, src, dest);
        } // else:  Duplicate value, do nothing

        return node;
    }

    // Function to delete a packet node recursively
    PacketTreeNode* deletePacketNodeRec(PacketTreeNode* node, int packetID) {
        if (node == nullptr) return nullptr;

        if (packetID < node->packetID)
            node->left = deletePacketNodeRec(node->left, packetID);
        else if (packetID > node->packetID)
            node->right = deletePacketNodeRec(node->right, packetID);
        else {
            // Node with only one child or no child
            if (node->left == nullptr) {
                PacketTreeNode* temp = node->right;
                delete node;
                return temp;
            } else if (node->right == nullptr) {
                PacketTreeNode* temp = node->left;
                delete node;
                return temp;
            }
            // Node with two children: Get the inorder successor (smallest
            // in the right subtree)
            PacketTreeNode* temp = findMin(node->right);

            if (!temp) {
                return nullptr;
            }
            // Copy the inorder successor's data to this node
            node->packetID = temp->packetID;
            node->sourceAddress = temp->sourceAddress;
            node->destinationAddress = temp->destinationAddress;
            node->timestamp = temp->timestamp;

            // Delete the inorder successor
            node->right = deletePacketNodeRec(node->right, temp->packetID);
        }
        return node;
    }


