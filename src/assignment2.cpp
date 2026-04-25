#include <bits/stdc++.h>
using namespace std;

using ll = long long;
static constexpr double EPS      = 1e-9;
static constexpr double INF_CAP  = 1e18;

static chrono::steady_clock::time_point PROG_START;
static constexpr double TIME_LIMIT_SEC = 3.0 * 3600.0;

static inline uint64_t edge_key(int u, int v) {
    if (u > v) swap(u, v);
    return (static_cast<uint64_t>(static_cast<uint32_t>(u)) << 32) |
            static_cast<uint32_t>(v);
}

static inline double elapsed_sec() {
    auto now = chrono::steady_clock::now();
    return chrono::duration<double>(now - PROG_START).count();
}

static inline void check_global_time() {
    if (elapsed_sec() >= TIME_LIMIT_SEC) {
        cout << "\n[TIMEOUT] 3-hour wall-clock limit reached. Exiting.\n";
        exit(0);
    }
}

static size_t peak_rss_bytes() {
    ifstream f("/proc/self/status");
    string line;
    while (getline(f, line)) {
        if (line.rfind("VmRSS:", 0) == 0) {
            size_t val = 0;
            sscanf(line.c_str(), "VmRSS: %zu", &val);
            return val * 1024ULL;
        }
    }
    return 0;
}

static string bytes_to_human(size_t b) {
    if (b < 1024) return to_string(b) + " B";
    if (b < 1024*1024) return to_string(b/1024) + " KB";
    if (b < (size_t)1024*1024*1024) return to_string(b/(1024*1024)) + " MB";
    return to_string(b/(1024*1024*1024)) + " GB";
}

struct Graph {
    int n = 0;
    vector<vector<int>> adj;
    vector<pair<int,int>> edges;  

    Graph() = default;
    explicit Graph(int n_) : n(n_), adj(n_) {}

    void add_edge(int u, int v) {
        if (u == v) return;
        adj[u].push_back(v);
        adj[v].push_back(u);
        if (u > v) swap(u, v);
        edges.push_back({u, v});
    }

    void finalize() {
        for (auto &nb : adj) {
            sort(nb.begin(), nb.end());
            nb.erase(unique(nb.begin(), nb.end()), nb.end());
        }
        sort(edges.begin(), edges.end());
        edges.erase(unique(edges.begin(), edges.end()), edges.end());
    }
};

struct DensityResult {
    vector<int> nodes;
    double density = 0.0;
    bool timed_out = false;
};

struct FlowEdge {
    int   to, rev;
    double cap;
};

struct Dinic {
    int N;
    vector<vector<FlowEdge>> G;
    vector<int> level, it;

    explicit Dinic(int n) : N(n), G(n), level(n), it(n) {}

    void add_edge(int u, int v, double cap) {
        G[u].push_back({v, (int)G[v].size(),      cap });
        G[v].push_back({u, (int)G[u].size() - 1,  0.0});
    }

    bool bfs(int s, int t) {
        fill(level.begin(), level.end(), -1);
        queue<int> q;
        level[s] = 0;
        q.push(s);
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (auto &e : G[u])
                if (e.cap > EPS && level[e.to] < 0) {
                    level[e.to] = level[u] + 1;
                    q.push(e.to);
                }
        }
        return level[t] >= 0;
    }

    double dfs(int u, int t, double f) {
        if (u == t || f <= EPS) return f;
        for (int &i = it[u]; i < (int)G[u].size(); ++i) {
            FlowEdge &e = G[u][i];
            if (e.cap <= EPS || level[e.to] != level[u] + 1) continue;
            double got = dfs(e.to, t, min(f, e.cap));
            if (got > EPS) {
                e.cap -= got;
                G[e.to][e.rev].cap += got;
                return got;
            }
        }
        return 0.0;
    }

    double max_flow(int s, int t) {
        double flow = 0.0;
        while (bfs(s, t)) {
            fill(it.begin(), it.end(), 0);
            for (double d; (d = dfs(s, t, INF_CAP)) > EPS; )
                flow += d;
        }
        return flow;
    }

    vector<char> reachable_from(int s) const {
        vector<char> vis(N, 0);
        queue<int> q;
        vis[s] = 1; q.push(s);
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (const auto &e : G[u])
                if (e.cap > EPS && !vis[e.to]) {
                    vis[e.to] = 1;
                    q.push(e.to);
                }
        }
        return vis;
    }
};

DensityResult charikar(const Graph &G) {
    int n = G.n;
    vector<int>  deg(n);
    vector<char> alive(n, 1);

    priority_queue<pair<int,int>,
                   vector<pair<int,int>>,
                   greater<pair<int,int>>> pq;

    ll twice_edges = 0;
    for (int i = 0; i < n; ++i) {
        deg[i] = (int)G.adj[i].size();
        twice_edges += deg[i];
        pq.push({deg[i], i});
    }
    ll rem_nodes = n;
    double best = n ? double(twice_edges / 2) / n : 0.0;
    vector<char> best_alive = alive;

    while (!pq.empty()) {
        check_global_time();
        auto [d, u] = pq.top(); pq.pop();
        if (!alive[u] || d != deg[u]) continue;

        alive[u] = 0;
        --rem_nodes;
        for (int v : G.adj[u]) {
            if (alive[v]) {
                --deg[v];
                twice_edges -= 2;
                pq.push({deg[v], v});
            }
        }
        if (rem_nodes > 0) {
            double cur = double(twice_edges / 2) / double(rem_nodes);
            if (cur > best + EPS) { best = cur; best_alive = alive; }
        }
    }

    DensityResult res;
    res.density = best;
    for (int i = 0; i < n; ++i) if (best_alive[i]) res.nodes.push_back(i);
    return res;
}

DensityResult greedy_pp(const Graph &G, int T) {
    int n = G.n;
    vector<ll>   load(n, 0);    
    vector<char> best_alive(n, 1);
    double best_density = (n > 0 && !G.edges.empty())
                            ? double(G.edges.size()) / n : 0.0;

    for (int pass = 0; pass < T; ++pass) {
        check_global_time();
        vector<int>  deg(n);
        vector<char> alive(n, 1);

        priority_queue<pair<ll,int>,
                       vector<pair<ll,int>>,
                       greater<pair<ll,int>>> pq;

        ll twice_edges = 0;
        for (int i = 0; i < n; ++i) {
            deg[i] = (int)G.adj[i].size();
            twice_edges += deg[i];
            pq.push({load[i] + deg[i], i});
        }
        ll rem_nodes = n;
        double pass_best = n ? double(twice_edges / 2) / n : 0.0;
        vector<char> pass_alive = alive;

        while (!pq.empty()) {
            check_global_time();
            auto [val, u] = pq.top(); pq.pop();
            if (!alive[u] || val != load[u] + deg[u]) continue;

            load[u] += deg[u];        // update load before removal
            alive[u] = 0;
            --rem_nodes;

            for (int v : G.adj[u]) {
                if (alive[v]) {
                    --deg[v];
                    twice_edges -= 2;
                    pq.push({load[v] + deg[v], v});
                }
            }
            if (rem_nodes > 0) {
                double cur = double(twice_edges / 2) / double(rem_nodes);
                if (cur > pass_best + EPS) { pass_best = cur; pass_alive = alive; }
            }
        }
        if (pass_best > best_density + EPS) {
            best_density = pass_best;
            best_alive   = pass_alive;
        }
    }

    DensityResult res;
    res.density = best_density;
    for (int i = 0; i < n; ++i) if (best_alive[i]) res.nodes.push_back(i);
    return res;
}

struct LocalSubgraph {
    vector<int>          global_of_local;
    vector<int>          local_of_global;
    vector<vector<int>>  adj;               
};

static LocalSubgraph build_local(const Graph &G, const vector<int> &verts) {
    LocalSubgraph S;
    int k = (int)verts.size();
    S.global_of_local = verts;
    S.local_of_global.assign(G.n, -1);
    for (int i = 0; i < k; ++i) S.local_of_global[verts[i]] = i;
    S.adj.assign(k, {});
    for (int i = 0; i < k; ++i) {
        for (int gv : G.adj[verts[i]]) {
            int j = S.local_of_global[gv];
            if (j != -1) S.adj[i].push_back(j);
        }
        sort(S.adj[i].begin(), S.adj[i].end());
        S.adj[i].erase(unique(S.adj[i].begin(), S.adj[i].end()), S.adj[i].end());
    }
    return S;
}

struct Triangle { int a, b, c; };

static vector<Triangle> enumerate_triangles(const LocalSubgraph &S) {
    int n = (int)S.adj.size();
    vector<int> deg(n);
    for (int i = 0; i < n; ++i) deg[i] = (int)S.adj[i].size();

    vector<int> order(n); iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int x, int y) {
        return deg[x] != deg[y] ? deg[x] < deg[y]
                                : S.global_of_local[x] < S.global_of_local[y];
    });
    vector<int> rnk(n);
    for (int i = 0; i < n; ++i) rnk[order[i]] = i;

    vector<vector<int>> fwd(n);
    for (int u = 0; u < n; ++u)
        for (int v : S.adj[u])
            if (rnk[u] < rnk[v]) fwd[u].push_back(v);

    vector<char> mark(n, 0);
    vector<Triangle> tris;
    tris.reserve(512);

    for (int u : order) {
        check_global_time();
        for (int v : fwd[u]) mark[v] = 1;
        for (int v : fwd[u])
            for (int w : fwd[v])
                if (mark[w]) tris.push_back({u, v, w});
        for (int v : fwd[u]) mark[v] = 0;
    }
    return tris;
}

struct CoreResult {
    vector<int> core;   
    int         kmax = 0;
};

static CoreResult triangle_core_decomp(const LocalSubgraph &S) {
    int n = (int)S.adj.size();
    vector<Triangle> tris = enumerate_triangles(S);

    vector<vector<int>> inc(n);
    vector<int> cnt(n, 0);
    for (int tid = 0; tid < (int)tris.size(); ++tid) {
        auto &t = tris[tid];
        inc[t.a].push_back(tid); ++cnt[t.a];
        inc[t.b].push_back(tid); ++cnt[t.b];
        inc[t.c].push_back(tid); ++cnt[t.c];
    }

    CoreResult R;
    R.core.assign(n, 0);
    vector<char> removed(n, 0);
    vector<char> tri_alive(tris.size(), 1);

    priority_queue<pair<int,int>,
                   vector<pair<int,int>>,
                   greater<pair<int,int>>> pq;
    for (int i = 0; i < n; ++i) pq.push({cnt[i], i});

    while (!pq.empty()) {
        check_global_time();
        auto [c, v] = pq.top(); pq.pop();
        if (removed[v] || c != cnt[v]) continue;

        removed[v]  = 1;
        R.core[v]   = c;
        R.kmax      = max(R.kmax, c);

        for (int tid : inc[v]) {
            if (!tri_alive[tid]) continue;
            tri_alive[tid] = 0;
            auto &t = tris[tid];
            for (int u : {t.a, t.b, t.c}) {
                if (u == v || removed[u]) continue;
                --cnt[u];
                pq.push({cnt[u], u});
            }
        }
    }
    return R;
}

DensityResult peel_app_triangle(const Graph &G) {
    int n = G.n;
    vector<int> all(n); iota(all.begin(), all.end(), 0);
    LocalSubgraph whole = build_local(G, all);
    vector<Triangle> tris = enumerate_triangles(whole);

    vector<int> tdeg(n, 0);
    vector<vector<int>> inc(n);
    for (int tid = 0; tid < (int)tris.size(); ++tid) {
        auto &t = tris[tid];
        ++tdeg[t.a]; inc[t.a].push_back(tid);
        ++tdeg[t.b]; inc[t.b].push_back(tid);
        ++tdeg[t.c]; inc[t.c].push_back(tid);
    }

    ll tri_count = (ll)tris.size();
    ll rem_nodes = n;
    vector<char> alive(n, 1);
    vector<char> tri_alive(tris.size(), 1);

    double best_dens = n > 0 ? double(tri_count) / n : 0.0;
    vector<char> best_alive = alive;

    priority_queue<pair<int,int>,
                   vector<pair<int,int>>,
                   greater<pair<int,int>>> pq;
    for (int i = 0; i < n; ++i) pq.push({tdeg[i], i});

    while (!pq.empty()) {
        check_global_time();
        auto [d, v] = pq.top(); pq.pop();
        if (!alive[v] || d != tdeg[v]) continue;

        alive[v] = 0;
        --rem_nodes;

        for (int tid : inc[v]) {
            if (!tri_alive[tid]) continue;
            tri_alive[tid] = 0;
            --tri_count;
            auto &t = tris[tid];
            for (int u : {t.a, t.b, t.c}) {
                if (u == v || !alive[u]) continue;
                --tdeg[u];
                pq.push({tdeg[u], u});
            }
        }
        if (rem_nodes > 0) {
            double cur = double(tri_count) / double(rem_nodes);
            if (cur > best_dens + EPS) { best_dens = cur; best_alive = alive; }
        }
    }

    DensityResult res;
    res.density = best_dens;
    for (int i = 0; i < n; ++i) if (best_alive[i]) res.nodes.push_back(i);
    return res;
}


struct ComponentData {
    vector<int>            verts;       
    LocalSubgraph          sub;
    vector<int>            tri_deg;      
    vector<Triangle>       tris;
    vector<pair<int,int>>  used_edges;   
    vector<array<int,3>>   tri_eids;     
};

static ComponentData prepare_component(const Graph &G, const vector<int> &verts) {
    ComponentData D;
    D.verts = verts;
    D.sub   = build_local(G, verts);
    D.tris  = enumerate_triangles(D.sub);

    int n = (int)D.sub.adj.size();
    D.tri_deg.assign(n, 0);

    unordered_map<uint64_t, int> eid_map;
    eid_map.reserve(D.tris.size() * 4);
    D.tri_eids.resize(D.tris.size());

    auto ensure_edge = [&](int u, int v) -> int {
        uint64_t key = edge_key(u, v);
        auto it = eid_map.find(key);
        if (it != eid_map.end()) return it->second;
        int id = (int)D.used_edges.size();
        eid_map[key] = id;
        if (u > v) swap(u, v);
        D.used_edges.push_back({u, v});
        return id;
    };

    for (int tid = 0; tid < (int)D.tris.size(); ++tid) {
        auto &t = D.tris[tid];
        ++D.tri_deg[t.a]; ++D.tri_deg[t.b]; ++D.tri_deg[t.c];
        D.tri_eids[tid] = {
            ensure_edge(t.b, t.c),   
            ensure_edge(t.a, t.c),  
            ensure_edge(t.a, t.b)   
        };
    }
    return D;
}

static double tri_density_of_selected(const vector<Triangle> &tris,
                                       const vector<char>     &sel) {
    ll cnt = 0, nodes = 0;
    for (const auto &t : tris)
        if (sel[t.a] && sel[t.b] && sel[t.c]) ++cnt;
    for (char c : sel) if (c) ++nodes;
    return nodes ? double(cnt) / double(nodes) : 0.0;
}

static DensityResult core_exact_on_component(const ComponentData &D,
                                              double proven_lo,
                                              double time_budget_sec) {
    int n    = (int)D.sub.adj.size();
    int ecnt = (int)D.used_edges.size();
    DensityResult ans;
    if (n == 0) return ans;
    if (D.tris.empty()) { ans.nodes = D.verts; return ans; }

    int local_kmax = *max_element(D.tri_deg.begin(), D.tri_deg.end());
    double rho0 = double(D.tris.size()) / n;  // density of whole component

    double lo  = proven_lo;       
    double hi  = double(local_kmax);
    double eps = (n >= 2) ? 1.0 / (double(n) * (n - 1)) : 1e-12;

    vector<char> best_sel(n, 1);   
    double       best_dens = rho0;

    vector<char> last_lo_sel(n, 1);
    bool         found_feasible = false;

    double deadline = elapsed_sec() + time_budget_sec;

    int source     = 0;
    int vbase      = 1;
    int ebase      = vbase + n;
    int sink       = ebase + ecnt;
    int flow_nodes = sink + 1;

    while (hi - lo >= eps) {
        check_global_time();
        if (elapsed_sec() >= deadline) { ans.timed_out = true; break; }

        double alpha = (lo + hi) * 0.5;
        Dinic din(flow_nodes);

        for (int v = 0; v < n; ++v) {
            din.add_edge(source, vbase + v, double(D.tri_deg[v]));
            din.add_edge(vbase + v, sink,   3.0 * alpha);
        }

        for (int eid = 0; eid < ecnt; ++eid) {
            auto [u, v] = D.used_edges[eid];
            int en = ebase + eid;
            din.add_edge(en, vbase + u, INF_CAP);
            din.add_edge(en, vbase + v, INF_CAP);
        }

        for (int tid = 0; tid < (int)D.tris.size(); ++tid) {
            auto &t   = D.tris[tid];
            auto &ids = D.tri_eids[tid];
            din.add_edge(vbase + t.a, ebase + ids[0], 1.0);
            din.add_edge(vbase + t.b, ebase + ids[1], 1.0);
            din.add_edge(vbase + t.c, ebase + ids[2], 1.0);
        }

        din.max_flow(source, sink);
        auto reach = din.reachable_from(source);

        vector<char> sel(n, 0);
        int sel_cnt = 0;
        for (int v = 0; v < n; ++v)
            if (reach[vbase + v]) { sel[v] = 1; ++sel_cnt; }

        if (sel_cnt == 0) {
            hi = alpha;
        } else {
            lo = alpha;
            found_feasible = true;
            last_lo_sel = sel;   

            double cur = tri_density_of_selected(D.tris, sel);
            if (cur > best_dens + EPS) { best_dens = cur; best_sel = sel; }
        }
    }

    if (found_feasible) {
        double final_dens = tri_density_of_selected(D.tris, last_lo_sel);
        if (final_dens > best_dens + EPS) {
            best_dens = final_dens;
            best_sel  = last_lo_sel;
        }
        if (lo > best_dens + EPS) {
            best_dens = lo;
            best_sel  = last_lo_sel;
        }
    }

    ans.density = best_dens;
    for (int i = 0; i < n; ++i)
        if (best_sel[i]) ans.nodes.push_back(D.verts[i]);
    return ans;
}

DensityResult core_exact_triangle(const Graph &G, double peelapp_density,
                                   double time_budget_sec = 3600.0) {
    DensityResult best;
    best.density = peelapp_density;  

    int n = G.n;
    vector<int> all(n); iota(all.begin(), all.end(), 0);

    LocalSubgraph whole = build_local(G, all);
    CoreResult    core  = triangle_core_decomp(whole);
    int kmax = core.kmax;

    int threshold = (kmax + 2) / 3;    

    vector<char> eligible(n, 0);
    for (int i = 0; i < n; ++i)
        if (core.core[i] >= threshold) eligible[i] = 1;

    vector<char> vis(n, 0);
    for (int s = 0; s < n; ++s) {
        if (!eligible[s] || vis[s]) continue;

        queue<int> q;
        vector<int> comp;
        vis[s] = 1; q.push(s);
        while (!q.empty()) {
            check_global_time();
            int u = q.front(); q.pop();
            comp.push_back(u);
            for (int v : G.adj[u])
                if (eligible[v] && !vis[v]) { vis[v] = 1; q.push(v); }
        }
        if (comp.empty()) continue;

        ComponentData D = prepare_component(G, comp);

        DensityResult cur = core_exact_on_component(D, best.density, time_budget_sec);
        if (cur.timed_out) best.timed_out = true;
        if (cur.density > best.density + EPS) best = cur;
    }

    if (best.nodes.empty()) {
        ComponentData D = prepare_component(G, all);
        DensityResult cur = core_exact_on_component(D, 0.0, time_budget_sec);
        if (cur.timed_out) best.timed_out = true;
        if (cur.density > best.density + EPS) best = cur;
    }

    return best;
}

static Graph read_graph(const string &path) {
    ifstream in(path);
    if (!in) throw runtime_error("Cannot open: " + path);

    unordered_map<long long, int> id;
    id.reserve(1 << 20);
    vector<pair<int,int>> raw;
    raw.reserve(1 << 20);

    string line;
    while (getline(in, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '%' || line[0] == 'c')
            continue;
        long long u0, v0;
        if (sscanf(line.c_str(), "%lld %lld", &u0, &v0) != 2) continue;
        if (u0 == v0) continue;

        auto get = [&](long long x) -> int {
            auto it = id.find(x);
            if (it != id.end()) return it->second;
            int nxt = (int)id.size();
            id[x] = nxt;
            return nxt;
        };
        int u = get(u0), v = get(v0);
        if (u > v) swap(u, v);
        raw.push_back({u, v});
    }

    sort(raw.begin(), raw.end());
    raw.erase(unique(raw.begin(), raw.end()), raw.end());

    int n = (int)id.size();
    Graph G(n);
    for (auto [u, v] : raw) G.add_edge(u, v);
    G.finalize();
    return G;
}

static void print_result(const string &name, const DensityResult &R) {
    cout << "  [" << name << "]\n";
    cout << "    Triangle-density : " << fixed << setprecision(6) << R.density << "\n";
    cout << "    Subgraph size    : " << R.nodes.size() << " vertices\n";
    if (R.timed_out) cout << "    (WARNING: timed out during binary search)\n";
}

int main(int argc, char **argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <graph-file> [greedy_pp_passes=5]\n";
        return 1;
    }

    int passes = 5;
    if (argc >= 3) passes = max(1, atoi(argv[2]));

    PROG_START = chrono::steady_clock::now();
    size_t mem_before = peak_rss_bytes();

    try {
        auto t_load0 = chrono::steady_clock::now();
        Graph G = read_graph(argv[1]);
        auto t_load1 = chrono::steady_clock::now();
        double load_sec = chrono::duration<double>(t_load1 - t_load0).count();

        cout << "Dataset  : " << argv[1] << "\n";
        cout << "Nodes    : " << G.n << "\n";
        cout << "Edges    : " << G.edges.size() << "\n";
        cout << "Load time: " << fixed << setprecision(3) << load_sec << " s\n";

        cout << "Running Charikar's greedy peeling (edge-density) ...\n";
        auto t0 = chrono::steady_clock::now();
        DensityResult ch = charikar(G);
        auto t1 = chrono::steady_clock::now();
        double sec_ch = chrono::duration<double>(t1 - t0).count();

        cout << "Running Greedy++ (" << passes << " passes, edge-density) ...\n";
        auto t2 = chrono::steady_clock::now();
        DensityResult gpp = greedy_pp(G, passes);
        auto t3 = chrono::steady_clock::now();
        double sec_gpp = chrono::duration<double>(t3 - t2).count();

        cout << "Running PeelApp (triangle-density peeling approx) ...\n";
        auto t4 = chrono::steady_clock::now();
        DensityResult pa = peel_app_triangle(G);
        auto t5 = chrono::steady_clock::now();
        double sec_pa = chrono::duration<double>(t5 - t4).count();

        cout << "Running CoreExact (triangle-density, Dinic flow) ...\n";
        double ce_budget = max(10.0, TIME_LIMIT_SEC - elapsed_sec() - 60.0);
        auto t6 = chrono::steady_clock::now();
        DensityResult ce = core_exact_triangle(G, pa.density, ce_budget);
        auto t7 = chrono::steady_clock::now();
        double sec_ce = chrono::duration<double>(t7 - t6).count();

        size_t mem_after = peak_rss_bytes();

        cout << "  [Charikar greedy (edge-density, 1/2-approx)]\n";
        cout << "    Edge-density     : " << fixed << setprecision(6)
             << ch.density << "\n";
        cout << "    Subgraph size    : " << ch.nodes.size() << " vertices\n";

        cout << "  [Greedy++ (edge-density, " << passes << " passes)]\n";
        cout << "    Edge-density     : " << fixed << setprecision(6)
             << gpp.density << "\n";
        cout << "    Subgraph size    : " << gpp.nodes.size() << " vertices\n";

        print_result("PeelApp (triangle-density peeling, 1/3-approx)", pa);
        print_result("CoreExact (triangle-density, Dinic flow, exact)", ce);

        double total_sec = chrono::duration<double>(t7 - t_load0).count();
        cout << "  Graph load             : " << fixed << setprecision(3) << load_sec  << " s\n";
        cout << "  Charikar (Alg 1)       : " << fixed << setprecision(3) << sec_ch    << " s\n";
        cout << "  Greedy++ (Alg 2)       : " << fixed << setprecision(3) << sec_gpp   << " s\n";
        cout << "  PeelApp triangle       : " << fixed << setprecision(3) << sec_pa    << " s\n";
        cout << "  CoreExact (Alg 4)      : " << fixed << setprecision(3) << sec_ce    << " s";
        if (ce.timed_out) cout << "  [TIMED OUT – partial result]";
        cout << "\n";
        cout << "  TOTAL                  : " << fixed << setprecision(3) << total_sec  << " s\n";

        size_t used = (mem_after > mem_before) ? (mem_after - mem_before) : mem_after;
        cout << "  RSS before algorithms  : " << bytes_to_human(mem_before) << "\n";
        cout << "  RSS after  algorithms  : " << bytes_to_human(mem_after)  << "\n";
        cout << "  Approx algo memory use : " << bytes_to_human(used)       << "\n";

    } catch (const exception &e) {
        cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
