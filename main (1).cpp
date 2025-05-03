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

    // Function to perform an inorder traversal of the packet tree
    void inorderTraversal() const {
        inorderTraversalRecursive(root);
        cout << endl;
    }

private:
    // Recursive function for inorder traversal
    void inorderTraversalRecursive(PacketTreeNode* node) const {
        if (node != nullptr) {
            inorderTraversalRecursive(node->left);
            cout << "Packet ID: " << node->packetID
                 << ", Source: " << node->sourceAddress
                 << ", Dest: " << node->destinationAddress
                 << ", Time: " << chrono::system_clock::to_time_t(node->timestamp) << " ";
            inorderTraversalRecursive(node->right);
        }
    }
};

// Stack for network packets
class PacketStack {
public:
    stack<PacketTreeNode*> packetStack;

    void pushPacket(PacketTreeNode* packet) {
        if (!packet) {
            throw invalid_argument("Packet cannot be null.");
        }
        packetStack.push(packet);
    }

    PacketTreeNode* popPacket() {
        if (packetStack.empty()) {
            return nullptr;
        }
        PacketTreeNode* topPacket = packetStack.top();
        packetStack.pop();
        return topPacket;
    }

    void printStack() const {
        if (packetStack.empty()) {
            cout << "Packet stack is empty." << endl;
            return;
        }
        stack<PacketTreeNode*> tempStack = packetStack;
        while (!tempStack.empty()) {
            PacketTreeNode* packet = tempStack.top();
            cout << "Packet ID: " << packet->packetID
                 << ", Source: " << packet->sourceAddress
                 << ", Dest: " << packet->destinationAddress
                 << ", Time: " << chrono::system_clock::to_time_t(packet->timestamp) << " ";
            tempStack.pop();
        }
        cout << endl;
    }
};

// Selection Sort algorithm
void selectionSort(vector<int>& arr) {
    if (arr.empty()) return;
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
    if (arr.empty()) return;
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
    try {
        // Create a graph representing the network of routers
        int V = 6;
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

        // Simulate packet transmission
        PacketQueue packetQueue;
        PacketStack packetStack;
        PacketList retransmissionList;

        packetQueue.enqueue(1);
        packetQueue.enqueue(2);
        packetQueue.enqueue(3);

        packetStack.pushPacket(new PacketTreeNode(4, "192.168.1.10", "192.168.1.20"));
        packetStack.pushPacket(new PacketTreeNode(5, "192.168.1.30", "192.168.1.40"));

        cout << "Packet Queue: ";
        packetQueue.printQueue();
        cout << "Packet Stack (LIFO): ";
        packetStack.printStack();

        double lossProbability = 0.1; // 10% loss probability

        while (!packetQueue.packetQueue.empty()) {
            int packet = packetQueue.packetQueue.front();
            packetQueue.dequeue();

            cout << "Processing Packet ID: " << packet << endl;
            if (g.simulatePacketLoss(lossProbability)) {
                cout << "Packet ID: " << packet << " lost. Retransmitting..." << endl;
                retransmissionList.addPacket(packet);
            } else {
                cout << "Packet ID: " << packet << " transmitted successfully." << endl;
            }
        }

        cout << "\nRetransmitted Packets: " << endl;
        retransmissionList.printList();

        // Dijkstra's algorithm
        cout << "\nShortest Paths from Router 0 (Dijkstra's):" << endl;
        vector<int> distDijkstra = g.dijkstra(0);
        for (int i = 0; i < V; ++i) {
            cout << "To Router " << i << ": Distance = " << distDijkstra[i] << endl;
        }

        // Performance Evaluation of Sorting Algorithms
        vector<int> sortingData = {5, 2, 9, 1, 5, 6}; // Example data
        vector<int> selectionSortData = sortingData;  //copy the data.
        vector<int> bubbleSortData = sortingData;     //copy the data.

        // Time Selection Sort
        auto start = chrono::high_resolution_clock::now();
        selectionSort(selectionSortData);
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> selectionSortTime = end - start;

        // Time Bubble Sort
        start = chrono::high_resolution_clock::now();
        bubbleSort(bubbleSortData);
        end = chrono::high_resolution_clock::now();
        chrono::duration<double> bubbleSortTime = end - start;

        cout << "\nSorted Data (Selection Sort): ";
        for (int val : selectionSortData) {
            cout << val << " ";
        }
        cout << endl;
        cout << "Selection Sort Time: " << selectionSortTime.count() << " seconds" << endl;

        cout << "\nSorted Data (Bubble Sort): ";
        for (int val : bubbleSortData) {
            cout << val << " ";
        }
        cout << endl;
        cout << "Bubble Sort Time: " << bubbleSortTime.count() << " seconds" << endl;

        // Demonstrate Binary Tree operations for network packets
        cout << "\nBinary Tree Operations (Network Packets):" << endl;
        PacketTree packetTree;
        packetTree.insertPacket(5, "192.168.1.1", "192.168.2.1");
        packetTree.insertPacket(3, "192.168.1.2", "192.168.2.2");
        packetTree.insertPacket(7, "192.168.1.3", "192.168.2.3");
        packetTree.insertPacket(2, "192.168.1.4", "192.168.2.4");
        packetTree.insertPacket(4, "192.168.1.5", "192.168.2.5");
        packetTree.insertPacket(6, "192.168.1.6", "192.168.2.6");
        packetTree.insertPacket(8, "192.168.1.7", "192.168.2.7");

        cout << "Inorder Traversal of Packet Tree: ";
        packetTree.inorderTraversal();

        packetTree.deletePacket(5);
        cout << "Inorder Traversal after deleting packet 5: ";
        packetTree.inorderTraversal();

        cout << "\nPacket Stack Operations:\n";
        PacketStack pStack;
        pStack.pushPacket(new PacketTreeNode(101, "10.0.0.1", "10.0.0.2"));
        pStack.pushPacket(new PacketTreeNode(102, "10.0.0.3", "10.0.0.4"));
        pStack.pushPacket(new PacketTreeNode(103, "10.0.0.5", "10.0.0.6"));

        cout << "Packets in stack:\n";
        pStack.printStack();

        PacketTreeNode* poppedPacket = pStack.popPacket();
        if (poppedPacket) {
            cout << "\nPopped Packet ID: " << poppedPacket->packetID << endl;
            delete poppedPacket; //remember to delete.
        }
        cout << "Packets in stack after pop:\n";
        pStack.printStack();
    } catch (const exception& e) {
        cerr << "An error occurred: " << e.what() << endl;
        return 1; // Return a non-zero value to indicate failure
    }

    return 0;
}
