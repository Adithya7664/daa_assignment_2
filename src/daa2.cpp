#include <bits/stdc++.h>
using namespace std;
using namespace std::chrono;

steady_clock::time_point START_TIME;
const double TIME_LIMIT = 3 * 60 * 60;

void check_time()
{
    auto now = steady_clock::now();
    double elapsed = duration_cast<seconds>(now - START_TIME).count();

    if (elapsed >= TIME_LIMIT)
    {
        cout << "\n[TLE] Time limit of 3 hours exceeded.\n";
        cout.flush();
        exit(0);
    }
}

struct Graph
{
    int n;
    vector<vector<int>> adj;

    Graph(int n) : n(n), adj(n) {}

    void add_edge(int u, int v)
    {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
};

double compute_density(Graph &G, vector<bool> &alive)
{
    int nodes = 0, edges = 0;

    for (int i = 0; i < G.n; i++)
    {
        if (alive[i])
        {
            nodes++;
            for (int v : G.adj[i])
                if (alive[v])
                    edges++;
        }
    }
    edges /= 2;
    return (nodes == 0) ? 0 : (double)edges / nodes;
}

pair<vector<int>, double> charikar(Graph &G)
{
    cout << "Running Charikar...\n";

    int n = G.n;
    vector<int> degree(n);
    vector<char> alive(n, 1);

    int remaining_nodes = n;
    int remaining_edges = 0;

    for (int i = 0; i < n; i++)
    {
        degree[i] = G.adj[i].size();
        remaining_edges += degree[i];
    }
    remaining_edges /= 2;

    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<>> pq;

    for (int i = 0; i < n; i++)
        pq.push({degree[i], i});

    double best_density = (double)remaining_edges / remaining_nodes;
    vector<char> best_alive = alive;

    while (!pq.empty())
    {
        pair<int, int> p = pq.top();
        int deg = p.first;
        int u = p.second;
        pq.pop();

        if (!alive[u] || deg != degree[u])
            continue;

        alive[u] = false;
        remaining_nodes--;

        for (int v : G.adj[u])
        {
            if (alive[v])
            {
                degree[v]--;
                remaining_edges--;
                pq.push({degree[v], v});
            }
        }

        if (remaining_nodes > 0)
        {
            double d = (double)remaining_edges / remaining_nodes;

            if (d > best_density)
            {
                best_density = d;
                best_alive = alive;
            }
        }
    }

    vector<int> best_nodes;
    for (int i = 0; i < n; i++)
        if (best_alive[i])
            best_nodes.push_back(i);

    cout << "Charikar done\n";
    return {best_nodes, best_density};
}

pair<vector<int>, double> greedy_pp(Graph &G, int T)
{
    cout << "Running Greedy++...\n";

    int n = G.n;
    vector<int> load(n, 0);

    double best_density = 0.0;
    vector<int> best_nodes;
    vector<char> best_alive(n, 1);
    for (int iter = 0; iter < T; iter++)
    {
        vector<int> degree(n);
        vector<char> alive(n, 1);

        int remaining_nodes = n;
        int remaining_edges = 0;

        for (int i = 0; i < n; i++)
        {
            degree[i] = G.adj[i].size();
            remaining_edges += degree[i];
        }
        remaining_edges /= 2;

        priority_queue<pair<int, int>, vector<pair<int, int>>, greater<>> pq;

        for (int i = 0; i < n; i++)
            pq.push({load[i] + degree[i], i});

        while (!pq.empty())
        {
            auto top = pq.top();
            pq.pop();

            int val = top.first;
            int u = top.second;

            if (!alive[u] || val != load[u] + degree[u])
                continue;
            // cout << "Removing node: " << u << " with degree: " << degree[u] << endl;
            // REAL removal
            int deg_u = degree[u];
            load[u] += deg_u;

            alive[u] = false;
            remaining_nodes--;

            for (int v : G.adj[u])
            {
                if (alive[v])
                {
                    degree[v]--;
                    remaining_edges--;

                    pq.push({load[v] + degree[v], v});
                }
            }

            if (remaining_nodes > 0)
            {
                double d = (double)remaining_edges / remaining_nodes;

                if (d > best_density)
                {
                    best_density = d;
                    best_alive = alive;
                }
            }
        }
    }
    for (int i = 0; i < n; i++)
        if (best_alive[i])
            best_nodes.push_back(i);

    cout << "Greedy++ done\n";
    return {best_nodes, best_density};
}

// Dinic Flow Implementation
struct Edge
{
    int v;
    double cap;
    int rev;
};

struct Dinic
{
    struct Edge
    {
        int to;
        double cap;
        int rev;
    };

    int N;
    vector<vector<Edge>> G;
    vector<int> level, ptr;

    Dinic(int n) : N(n), G(n), level(n), ptr(n) {}

    void add_edge(int u, int v, double cap)
    {
        Edge a = {v, cap, (int)G[v].size()};
        Edge b = {u, 0, (int)G[u].size()};
        G[u].push_back(a);
        G[v].push_back(b);
    }

    bool bfs(int s, int t)
    {
        fill(level.begin(), level.end(), -1);
        queue<int> q;
        q.push(s);
        level[s] = 0;

        while (!q.empty())
        {
            int u = q.front();
            q.pop();

            for (auto &e : G[u])
            {
                if (e.cap > 1e-9 && level[e.to] == -1)
                {
                    level[e.to] = level[u] + 1;
                    q.push(e.to);
                }
            }
        }

        return level[t] != -1;
    }

    double dfs(int u, int t, double pushed)
    {
        if (pushed < 1e-9)
            return 0;
        if (u == t)
            return pushed;

        for (int &cid = ptr[u]; cid < (int)G[u].size(); cid++)
        {
            Edge &e = G[u][cid];

            if (level[e.to] != level[u] + 1 || e.cap < 1e-9)
                continue;

            double tr = dfs(e.to, t, min(pushed, e.cap));
            if (tr < 1e-9)
                continue;

            e.cap -= tr;
            G[e.to][e.rev].cap += tr;
            return tr;
        }

        return 0;
    }

    double max_flow(int s, int t)
    {
        double flow = 0;

        while (bfs(s, t))
        {
            fill(ptr.begin(), ptr.end(), 0);

            while (double pushed = dfs(s, t, 1e18))
                flow += pushed;
        }

        return flow;
    }

    vector<bool> min_cut(int s)
    {
        vector<bool> vis(N, false);
        queue<int> q;

        q.push(s);
        vis[s] = true;

        while (!q.empty())
        {
            int u = q.front();
            q.pop();

            for (auto &e : G[u])
            {
                if (e.cap > 1e-9 && !vis[e.to])
                {
                    vis[e.to] = true;
                    q.push(e.to);
                }
            }
        }

        return vis;
    }
};

pair<vector<int>, double> fang_exact_with_nodes(Graph &G)
{
    cout << "Running Fang exact..." << endl;

    int n = G.n;
    int m = 0;
    for (int i = 0; i < n; i++)
        m += G.adj[i].size();
    m /= 2;

    double low = 0, high = m;
    vector<int> best_nodes;
    double best_density = 0;

    for (int it = 0; it < 30; it++)
    {
        double mid = (low + high) / 2;

        int S = n, T = n + 1;
        Dinic D(n + 2);

        for (int i = 0; i < n; i++)
        {
            int deg = G.adj[i].size();

            D.add_edge(S, i, m);

            double cap = m + 2 * mid - deg;
            if (cap < 0)
                cap = 0;

            D.add_edge(i, T, cap);
        }

        for (int u = 0; u < n; u++)
        {
            for (int v : G.adj[u])
            {
                if (u < v)
                {
                    D.add_edge(u, v, 1);
                    D.add_edge(v, u, 1);
                }
            }
        }

        double flow = D.max_flow(S, T);

        if (flow < (double)m * n)
        {
            low = mid;

            vector<bool> reachable = D.min_cut(S);
            vector<int> current_nodes;

            for (int i = 0; i < n; i++)
                if (reachable[i])
                    current_nodes.push_back(i);

            if (mid > best_density && !current_nodes.empty())
            {
                best_density = mid;
                best_nodes = current_nodes;
            }
        }
        else
        {
            high = mid;
        }
        cout << "Iteration " << it + 1 << ": Density = " << mid << ", Flow = " << flow << ", Nodes in cut = " << best_nodes.size() << endl;
    }

    cout << "Done\n";
    return {best_nodes, best_density};
}

vector<int> core_decomp(Graph &G)
{
    int n = G.n;
    vector<int> deg(n);

    int max_deg = 0;
    for (int i = 0; i < n; i++)
    {
        deg[i] = G.adj[i].size();
        max_deg = max(max_deg, deg[i]);
    }

    vector<vector<int>> bucket(max_deg + 1);

    for (int i = 0; i < n; i++)
        bucket[deg[i]].push_back(i);

    vector<int> core(n);
    vector<bool> removed(n, false);

    int curr_deg = 0;

    for (int k = 0; k <= max_deg; k++)
    {
        while (!bucket[k].empty())
        {
            int u = bucket[k].back();
            bucket[k].pop_back();

            if (removed[u])
                continue;

            removed[u] = true;
            core[u] = k;

            for (int v : G.adj[u])
            {
                if (!removed[v])
                {
                    int d = deg[v];
                    deg[v]--;

                    bucket[d - 1].push_back(v);
                }
            }
        }
    }

    return core;
}

pair<vector<int>, double> fang_core_flow_with_nodes(Graph &G)
{
    cout << "Running Fang Core+Flow..." << endl;

    vector<int> core = core_decomp(G);
    int max_core = *max_element(core.begin(), core.end());

    vector<int> mapping;
    vector<int> reverse_map(G.n, -1);
    cout << "Core decomposition done. Max core number: " << max_core << endl;
    vector<int> best_nodes;
    double best_density = 0;

    for (int k = max_core; k >= max_core - 5; k--)
    {
        vector<int> mapping;
        vector<int> reverse_map(G.n, -1);

        for (int i = 0; i < G.n; i++)
        {
            if (core[i] >= k)
            {
                reverse_map[i] = mapping.size();
                mapping.push_back(i);
            }
        }

        cout << "Trying core >= " << k << ", size: " << mapping.size() << endl;

        if (mapping.empty())
            continue;

        Graph H(mapping.size());

        for (int i = 0; i < mapping.size(); i++)
        {
            int u = mapping[i];
            for (int v : G.adj[u])
            {
                if (reverse_map[v] != -1 && i < reverse_map[v])
                {
                    H.add_edge(i, reverse_map[v]);
                }
            }
        }

        auto res = fang_exact_with_nodes(H);

        if (res.second > best_density)
        {
            best_density = res.second;
            best_nodes.clear();

            for (int u : res.first)
                best_nodes.push_back(mapping[u]);
        }
    }

    if (mapping.empty())
    {
        cout << "Core reduction empty, falling back to full graph\n"
             << endl;
        return fang_exact_with_nodes(G);
    }

    return {best_nodes, best_density};
}

int main(int argc, char *argv[])
{
    START_TIME = steady_clock::now();
    cout << "Program running" << endl;

    ifstream in(argv[1]);
    if (!in)
    {
        cout << "FILE NOT OPENED\n";
        return 1;
    }
    cout << "file opened" << endl;
    auto t_read_start = steady_clock::now();
    unordered_map<int, int> mmap;
    vector<pair<int, int>> edges;
    vector<int> revmap;

    string line;
    int idx = 0;
    while (getline(in, line))
    {
        if (line.empty() || line[0] == '#')
            continue;

        int u, v;
        if (sscanf(line.c_str(), "%d %d", &u, &v) != 2)
            continue;
        // cout << u << " " << v << endl;
        if (!mmap.count(u))
        {
            mmap[u] = idx++;
            revmap.push_back(u);
        }
        if (!mmap.count(v))
        {
            mmap[v] = idx++;
            revmap.push_back(v);
        }
        int um = mmap[u];
        int vm = mmap[v];
        edges.push_back({um, vm});
    }
    auto t_read_end = steady_clock::now();
    int n = mmap.size();
    // vector<vector<int>> adj(n);
    // for (auto &e : edges)
    // {
    //     int u = e.first;
    //     int v = e.second;
    //     if (u >= n || v >= n)
    //     {
    //         cout << "CRASH INDEX: " << u << " " << v << endl;
    //     }
    //     adj[u].push_back(v);
    //     adj[v].push_back(u);
    // }
    cout << "n: " << n << endl;
    // cout << "adj size: " << adj.size() << endl;
    Graph G(n);
    cout << "Adding edges to graph..." << endl;
    unordered_set<long long> seen;
    for (auto &e : edges)
    {
        int u = e.first;
        int v = e.second;

        if (u > v)
            swap(u, v);

        long long key = (long long)u * n + v;

        if (seen.count(key))
            continue;

        seen.insert(key);
        G.add_edge(u, v);
    }
    cout << "Edges added to graph." << endl;
    long long total_edges = 0;
    for (int i = 0; i < n; i++)
        total_edges += G.adj[i].size();

    cout << "Edges: " << total_edges / 2 << endl;
    cout << "-> RESULTS \n\n"
         << endl;

    auto t1_start = steady_clock::now();
    auto res1 = charikar(G);
    auto t1_end = steady_clock::now();

    cout << "[Charikar]" << endl;
    cout << "Density: " << res1.second << endl;
    cout << "Nodes (" << res1.first.size() << "): ";
    for (int u : res1.first)
        cout << u << " ";
    cout << "\n"
         << endl;

    int T = 5;

    auto t2_start = steady_clock::now();
    auto res2 = greedy_pp(G, T);
    auto t2_end = steady_clock::now();

    cout << "[Greedy++] (T = " << T << ")" << endl;
    cout << "Density: " << res2.second << endl;
    cout << "Nodes (" << res2.first.size() << "): ";
    for (int u : res2.first)
        cout << u << " ";
    cout << "\n"
         << endl;
    if (n > 1000000 || total_edges > 1000000)
    {
        auto total_end = steady_clock::now();

        cout << "\n[INFO] Large dataset detected (Skitter). Skipping flow algorithms.\n";

        cout << "\n TIMINGS \n";

        cout << "Reading time: "
             << (double)duration_cast<milliseconds>(t_read_end - t_read_start).count() / 1000
             << " s\n";

        cout << "Charikar time: "
             << (double)duration_cast<milliseconds>(t1_end - t1_start).count() / 1000
             << " s\n";

        cout << "Greedy++ time: "
             << (double)duration_cast<milliseconds>(t2_end - t2_start).count() / 1000
             << " s\n";

        cout << "Total execution time: "
             << (double)duration_cast<milliseconds>(total_end - START_TIME).count() / 1000
             << " s\n";

        cout << "\n DONE \n";
        return 0;
    }
    auto t3_start = steady_clock::now();
    auto res3 = fang_exact_with_nodes(G);
    auto t3_end = steady_clock::now();

    cout << "[Fang Exact - Flow]" << endl;
    cout << "Density: " << res3.second << endl;
    cout << "Nodes (" << res3.first.size() << "): ";
    for (int u : res3.first)
        cout << u << " ";
    cout << "\n"
         << endl;

    auto t4_start = steady_clock::now();
    auto res4 = fang_core_flow_with_nodes(G);
    auto t4_end = steady_clock::now();

    cout << "[Fang Core + Flow]" << endl;
    cout << "Density: " << res4.second << endl;
    cout << "Nodes (" << res4.first.size() << "): ";
    for (int u : res4.first)
        cout << u << " ";
    cout << "\n"
         << endl;

    auto total_end = steady_clock::now();

    cout << "\n TIMINGS \n";

    cout << "Reading time: "
         << (double)duration_cast<milliseconds>(t_read_end - t_read_start).count() / 1000
         << " s\n";

    cout << "Charikar time: "
         << (double)duration_cast<milliseconds>(t1_end - t1_start).count() / 1000
         << " s\n";

    cout << "Greedy++ time: "
         << (double)duration_cast<milliseconds>(t2_end - t2_start).count() / 1000
         << " s\n";

    cout << "Fang Exact time: "
         << (double)duration_cast<milliseconds>(t3_end - t3_start).count() / 1000
         << " s\n";

    cout << "Core+Flow time: "
         << (double)duration_cast<milliseconds>(t4_end - t4_start).count() / 1000
         << " s\n";

    cout << "Total execution time: "
         << (double)duration_cast<milliseconds>(total_end - START_TIME).count() / 1000
         << " s\n";
    cout << "DONE \n";
    return 0;
}