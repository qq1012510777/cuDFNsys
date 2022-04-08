#include "Graph/Graph.cuh"

// ====================================================
// NAME:        Graph
// DESCRIPTION: constructor,
//              create adjacent list
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
cuDFNsys::Graph::Graph(const size_t &NumOfFractures_1,
                       std::map<pair<size_t, size_t>, pair<float3, float3>> Intersection_map)
{
    Adj = new list<int>[NumOfFractures_1];
    V = NumOfFractures_1;

    for (std::map<pair<size_t, size_t>, pair<float3, float3>>::iterator i =
             Intersection_map.begin();
         i != Intersection_map.end(); ++i)
    {
        addEdge(i->first.first, i->first.second);
        addEdge(i->first.second, i->first.first);
    }
}; // Graph

// ====================================================
// NAME:        DFS
// DESCRIPTION: DFS
//
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
void cuDFNsys::Graph::DFS(std::vector<std::vector<size_t>> &ListOfClusters)
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

// ====================================================
// NAME:        DFSUtil
// DESCRIPTION: DFSUtil
//
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
void cuDFNsys::Graph::DFSUtil(int s, vector<bool> &visited, vector<size_t> &onecluster)
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