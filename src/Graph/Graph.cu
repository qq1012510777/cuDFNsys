#include "Graph/Graph.cuh"

// ====================================================
// NAME:        Graph
// DESCRIPTION: constructor,
//              create adjacent list
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
template <typename T>
cuDFNsys::Graph<T>::Graph(const size_t &NumOfFractures_1,
                          std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<T>, cuDFNsys::Vector3<T>>> Intersection_map)
{
    Adj = new list<int>[NumOfFractures_1];
    V = NumOfFractures_1;
    // std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<T>,
    //                                     cuDFNsys::Vector3<T>>>::iterator i;
    for (auto i : Intersection_map)
    {
        addEdge(i.first.first, i.first.second);
        addEdge(i.first.second, i.first.first);
    }
}; // Graph
template cuDFNsys::Graph<double>::Graph(const size_t &NumOfFractures_1,
                                        std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<double>, cuDFNsys::Vector3<double>>> Intersection_map);
template cuDFNsys::Graph<float>::Graph(const size_t &NumOfFractures_1,
                                       std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<float>, cuDFNsys::Vector3<float>>> Intersection_map);

// ====================================================
// NAME:        DFS
// DESCRIPTION: DFS
//
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
template <typename T>
void cuDFNsys::Graph<T>::DFS(std::vector<std::vector<size_t>> &ListOfClusters)
{
    // Mark all the vertices as not visited
    vector<bool> visited(V, false);

    for (int i = 0; i < (int)V; i++)
    {
        vector<size_t> onecluster;
        onecluster.reserve(V);

        if (!visited[i])
        {
            DFSUtil(i, visited, onecluster);
        }

        onecluster.shrink_to_fit();

        if (onecluster.size() > 0)
        {
            ListOfClusters.push_back(onecluster);
        }
    }
}; // DFS
template void cuDFNsys::Graph<double>::DFS(std::vector<std::vector<size_t>> &ListOfClusters);
template void cuDFNsys::Graph<float>::DFS(std::vector<std::vector<size_t>> &ListOfClusters);

// ====================================================
// NAME:        DFSUtil
// DESCRIPTION: DFSUtil
//
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
template <typename T>
void cuDFNsys::Graph<T>::DFSUtil(int s, vector<bool> &visited, vector<size_t> &onecluster)
{
    // Create a stack for DFS
    stack<int> stack;

    // Push the current source node.
    stack.push(s);

    while (!stack.empty())
    {
        // Pop a vertex from stack and print it
        s = stack.top();
        stack.pop();

        // Stack may contain same vertex twice. So
        // we need to print the popped item only
        // if it is not visited.
        if (!visited[s])
        {
            //cout << s << " ";
            onecluster.push_back((size_t)s);
            visited[s] = true;
        }

        // Get all adjacent vertices of the popped vertex s
        // If a adjacent has not been visited, then push it
        // to the stack.
        for (auto i = Adj[s].begin(); i != Adj[s].end(); ++i)
        {
            if (!visited[*i])
            {
                stack.push(*i);
            };
        }
    }
}; // DFSUtil
template void cuDFNsys::Graph<double>::DFSUtil(int s, vector<bool> &visited, vector<size_t> &onecluster);
template void cuDFNsys::Graph<float>::DFSUtil(int s, vector<bool> &visited, vector<size_t> &onecluster);

// ====================================================
// NAME:        UseDFS
// DESCRIPTION: UseDFS
//
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
template <typename T>
inline void cuDFNsys::Graph<T>::UseDFS(std::vector<std::vector<size_t>> &S)
{
    S.clear();
    DFS(S);
    S.shrink_to_fit();
}; // UseDFS
template void cuDFNsys::Graph<double>::UseDFS(std::vector<std::vector<size_t>> &S);
template void cuDFNsys::Graph<float>::UseDFS(std::vector<std::vector<size_t>> &S);

// ====================================================
// NAME:        addEdge
// DESCRIPTION: create adjacent list
//
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
template <typename T>
inline void cuDFNsys::Graph<T>::addEdge(int v, int w)
{
    // Add w to v’s list.
    Adj[v].push_back(w);
}; // addEdge
template void cuDFNsys::Graph<double>::addEdge(int v, int w);
template void cuDFNsys::Graph<float>::addEdge(int v, int w);