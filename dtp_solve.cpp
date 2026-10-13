

#include <bits/stdc++.h>
#include <unordered_map>
#include <fstream>
#include <chrono>
#include <iomanip>
using namespace std;
namespace fs = std::filesystem;
using TimePoint = std::chrono::time_point<std::chrono::high_resolution_clock>;
using ll = long long;
const double INF = 1e18;

vector<vector<double>> dist_mat;
vector<vector<int>> next_node;

mt19937 rng;
struct Edge {
    int u, v;
    double w;
};

struct DSU {
    int n;
    vector<int> p;
    DSU(int n_ = 0) { init(n_); }
    void init(int n_) {
        n = n_;
        p.resize(n);
        iota(p.begin(), p.end(), 0);
    }
    int find(int x) { return p[x] == x ? x : p[x] = find(p[x]); }
    bool unite(int a, int b) {
        a = find(a); b = find(b);
        if (a == b) return false;
        p[b] = a;
        return true;
    }
};


struct TrieNode {
    unordered_map<int, TrieNode*> children;
    TrieNode* parent;
    int vertex;
    bool is_end;

    TrieNode(TrieNode* p = nullptr, int v = -1)
        : parent(p), vertex(v), is_end(false) {}
};



class EliteTrie {
private:
    TrieNode* root;
    vector<TrieNode*> elite_leaves;
    int max_safe_capacity;


    vector<int> to_sequence(const vector<char>& inDS) const {
        vector<int> seq;
        for (int i = 0; i < inDS.size(); ++i) {
            if (inDS[i]) seq.push_back(i);
        }
        return seq;
    }


    void free_tree(TrieNode* node) {
        if (!node) return;
        for (auto& pair : node->children) {
            free_tree(pair.second);
        }
        delete node;
    }

public:
    EliteTrie(int safe_capacity = 10000) : max_safe_capacity(safe_capacity) {
        root = new TrieNode();
    }

    ~EliteTrie() {
        free_tree(root);
    }


    void clear() {
        free_tree(root);
        root = new TrieNode();
        elite_leaves.clear();
    }

    int size() const {
        return elite_leaves.size();
    }


    bool insert(const vector<char>& inDS) {

        if (elite_leaves.size() >= max_safe_capacity) return false;

        vector<int> seq = to_sequence(inDS);
        TrieNode* curr = root;

        for (int val : seq) {
            if (curr->children.find(val) == curr->children.end()) {
                curr->children[val] = new TrieNode(curr, val);
            }
            curr = curr->children[val];
        }

        if (!curr->is_end) {
            curr->is_end = true;
            elite_leaves.push_back(curr);
            return true;
        }
        return false;
    }

    vector<char> get_random_solution(int n) const {
        vector<char> res(n, 0);
        if (elite_leaves.empty()) return res;

        int random_idx = rand() % elite_leaves.size();
        TrieNode* curr = elite_leaves[random_idx];

        while (curr != root) {
            res[curr->vertex] = 1;
            curr = curr->parent;
        }

        return res;
    }

    bool contains(const vector<char>& inDS) const {
        if (elite_leaves.empty()) return false;

        vector<int> seq = to_sequence(inDS);
        TrieNode* curr = root;

        for (int val : seq) {
            auto it = curr->children.find(val);
            if (it == curr->children.end()) {
                return false;
            }
            curr = it->second;
        }

        return curr->is_end;
    }

};


struct Graph {
    int n;
    vector<vector<pair<int,double>>> adj;
    vector<Edge> edges;
    Graph(int n_ = 0) { init(n_); }
    void init(int n_) {
        n = n_;
        adj.assign(n, {});
        edges.clear();
    }
    void addEdge(int u, int v, double w) {
        if (u < 0 || u >= n || v < 0 || v >= n) return;
        adj[u].push_back({v,w});
        adj[v].push_back({u,w});
        edges.push_back({u,v,w});
    }
};


pair<double,int> buildMST_on_X(const Graph &G, const vector<char> &inDS){
    int n = G.n;
    int cntNodes = 0;
    for (int i = 0; i < n; ++i) if (inDS[i]) ++cntNodes;
    if (cntNodes == 0) return {INF, 0};

    vector<Edge> cand;
    cand.reserve(G.edges.size());
    for (const Edge &e : G.edges) {
        if (inDS[e.u] && inDS[e.v]) cand.push_back(e);
    }
    sort(cand.begin(), cand.end(), [](const Edge &a, const Edge &b){ return a.w < b.w; });

    DSU dsu(n);
    double total = 0;
    int used = 0;
    for (const Edge &e : cand) {
        if (dsu.unite(e.u, e.v)) {
            total += e.w;
            used++;
        }
    }

    int root = -1;
    for (int i = 0; i < n; ++i) if (inDS[i]) {
        if (root == -1) root = dsu.find(i);
        else if (dsu.find(i) != root) return {INF, used};
    }
    if (used != cntNodes - 1) return {INF, used};

    return {total, used};
}


void initial_solution_by_mst_prune(const Graph &G, vector<char> &inDS, vector<int> &DSlist, vector<int> &leafEdgeIdx){
    int n = G.n;
    inDS.assign(n, 0);
    leafEdgeIdx.assign(n, -1);

    vector<int> idx(G.edges.size());
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(), [&](int a, int b){ return G.edges[a].w < G.edges[b].w; });
    DSU dsu(n);
    vector<vector<pair<int,double>>> treeAdj(n);
    for (int id : idx) {
        const Edge &e = G.edges[id];
        if (dsu.unite(e.u, e.v)) {
            treeAdj[e.u].push_back({e.v, e.w});
            treeAdj[e.v].push_back({e.u, e.w});
        }
    }

    for (int i = 0; i < n; ++i) inDS[i] = 1;

    priority_queue<pair<double,int>> pq;
    vector<int> degree(n, 0);
    vector<double> maxIncidentEdge(n, -1);
    for (int u = 0; u < n; ++u) {
        degree[u] = (int)treeAdj[u].size();
        if (degree[u] == 1) {
            maxIncidentEdge[u] = treeAdj[u][0].second;
            pq.push({maxIncidentEdge[u], u});
        }
    }

    auto isDominating = [&](const vector<char> &curInDS)->bool{
        vector<char> covered(n, 0);
        for (int v = 0; v < n; ++v) if (curInDS[v]) {
            covered[v] = 1;
            for (auto &pr : G.adj[v]) covered[pr.first] = 1;
        }
        for (int v = 0; v < n; ++v) if (!covered[v]) return false;
        return true;
    };

    while (!pq.empty()) {
        auto pr = pq.top(); pq.pop();
        double w = pr.first;
        int u = pr.second;
        if (!inDS[u]) continue;
        inDS[u] = 0;
        if (!isDominating(inDS)) {
            inDS[u] = 1;
            continue;
        }
        for (auto &nbr : treeAdj[u]) {
            int v = nbr.first;
            if (degree[v] > 0) {
                degree[v]--;
                if (degree[v] == 1) {
                    double mw = -1;
                    for (auto &pr2 : treeAdj[v]) {
                        int z = pr2.first;
                        if (z == u && inDS[u] == 0) continue;
                        mw = max(mw, pr2.second);
                    }
                    if (mw < 0) mw = 0;
                    pq.push({mw, v});
                }
            }
        }
    }

    DSlist.clear();
    for (int i = 0; i < n; ++i) if (inDS[i]) DSlist.push_back(i);
}

void recompute_cover(const Graph &G, const vector<char> &inDS, vector<int> &cover){
    int n = G.n;
    cover.assign(n, 0);
    for (int u = 0; u < n; ++u) if (inDS[u]) {
        cover[u]++;
        for (auto &pr : G.adj[u]) cover[pr.first]++;
    }
}

int computeA1_delete(const Graph &G, int v, const vector<char> &inDS, const vector<int> &cover){
    if (!inDS[v]) return 0;
    int cnt = (cover[v] == 1 ? 1 : 0);
    for (auto &pr : G.adj[v]) {
        int nei = pr.first;
        if (!inDS[nei] && cover[nei] == 1) cnt++;
    }
    return cnt;
}

int computeA1_add(const Graph &G, int v, const vector<char> &inDS, const vector<int> &cover) {
    if (inDS[v]) return 0;
    int cnt = 0;
    if (cover[v] == 0) cnt--;
    for (auto &pr : G.adj[v]) {
        int nei = pr.first;
        if (cover[nei] == 0) cnt--;
    }
    return cnt;
}

double computeA2_by_rebuild(const Graph &G, const vector<char> &inDS_tmp){
    auto res = buildMST_on_X(G, inDS_tmp);
    return res.first;
}

int countUndominated(const Graph &G, const vector<char> &inDS) {
    int n = G.n;
    vector<char> covered(n, 0);

    for (int u = 0; u < n; ++u) if (inDS[u]) {
        covered[u] = 1;
        for (auto &pr : G.adj[u]) {
            covered[pr.first] = 1;
        }
    }

    int cnt = 0;
    for (int i = 0; i < n; ++i) {
        if (!covered[i]) cnt++;
    }
    return cnt;
}



Graph loadGraph(const string &filename) {
    ifstream fin(filename);
    if (!fin) {
        cerr << "Cannot open " << filename << "\n";
        exit(1);
    }
    int V, E;
    fin >> V >> E;
    Graph G(V);
    for (int i = 0; i < E; ++i) {
        int u, v; double w;
        fin >> u >> v >> w;
        G.addEdge(u, v, w);
    }
    return G;
}

double buildInitialSolution(const Graph &G, vector<char> &inDS, vector<int> &DSlist, vector<int> &tmp, vector<int> &cover) {
    int n = G.n;
    inDS.assign(n, 1);
    tmp.assign(n, -1);
    DSlist.clear();
    for (int i = 0; i < n; ++i) {
        DSlist.push_back(i);
    }
    recompute_cover(G, inDS, cover);
    auto res = buildMST_on_X(G, inDS);
    if (res.first == INF) {
        cerr << "Initial solution not feasible (G[X] not connected). Abort\n";
        exit(1);
    }
    return res.first;
}


void flip_to_1(int v, vector<char>& inDS, vector<int>& cover, int& undom_count, const Graph& G) {
    if (inDS[v] == 1) return;
    inDS[v] = 1;

    if (cover[v] == 0) undom_count--;
    cover[v]++;

    for (auto& edge : G.adj[v]) {
        int u = edge.first;
        if (cover[u] == 0) undom_count--;
        cover[u]++;
    }
}


void flip_to_0(int v, vector<char>& inDS, vector<int>& cover, int& undom_count, const Graph& G) {
    if (inDS[v] == 0) return;
    inDS[v] = 0;

    cover[v]--;
    if (cover[v] == 0) undom_count++;

    for (auto& edge : G.adj[v]) {
        int u = edge.first;
        cover[u]--;
        if (cover[u] == 0) undom_count++;
    }
}

void repair_connectivity(const Graph &G, vector<char> &inDS,vector<int> &cover, int &undom_count) {
    int n = G.n;
    int total_inDS = 0;
    int start_node = -1;

    for (int i = 0; i < n; ++i) {
        if (inDS[i]) {
            total_inDS++;
            if (start_node == -1) start_node = i;
        }
    }

    if (total_inDS <= 1) return;

    vector<char> in_C1(n, 0);
    vector<int> newly_added;
    int count_C1 = 0;


    queue<int> q;
    q.push(start_node);
    in_C1[start_node] = 1;
    count_C1++;
    newly_added.push_back(start_node);

    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (auto &edge : G.adj[u]) {
            int v = edge.first;
            if (inDS[v] && !in_C1[v]) {
                in_C1[v] = 1;
                count_C1++;
                q.push(v);
                newly_added.push_back(v);
            }
        }
    }

    vector<double> min_dist_to_C1(n, INF);
    vector<int> best_u_in_C1(n, -1);

    while (count_C1 < total_inDS) {

        for (int u : newly_added) {
            for (int v = 0; v < n; ++v) {
                if (inDS[v] && !in_C1[v]) {
                    if (dist_mat[u][v] < min_dist_to_C1[v]) {
                        min_dist_to_C1[v] = dist_mat[u][v];
                        best_u_in_C1[v] = u;
                    }
                }
            }
        }
        newly_added.clear();

        double best_d = INF;
        int target_v = -1;
        for (int v = 0; v < n; ++v) {
            if (inDS[v] && !in_C1[v]) {
                if (min_dist_to_C1[v] < best_d) {
                    best_d = min_dist_to_C1[v];
                    target_v = v;
                }
            }
        }

        if (target_v == -1) break;

        int curr = best_u_in_C1[target_v];
        while (curr != target_v) {
            curr = next_node[curr][target_v];

            if (!inDS[curr]) {

                flip_to_1(curr, inDS, cover, undom_count, G);
                total_inDS++;
            }

            if (!in_C1[curr]) {
                in_C1[curr] = 1;
                count_C1++;
                newly_added.push_back(curr);
            }
        }


        if (!in_C1[target_v]) {
            in_C1[target_v] = 1;
            count_C1++;
            newly_added.push_back(target_v);
        }

        q.push(target_v);
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (auto &edge : G.adj[u]) {
                int v = edge.first;
                if (inDS[v] && !in_C1[v]) {
                    in_C1[v] = 1;
                    count_C1++;
                    q.push(v);
                    newly_added.push_back(v);
                }
            }
        }
    }
}

vector<char> mixed_perturb(const Graph &G, const vector<char> &baseDS, int k1, int k2) {
    int n = G.n;
    vector<char> newDS = baseDS;
    cout<< "Perturbing: k1="<<k1<<", k2="<<k2<<endl;
    vector<pair<int, int>> ds_nodes;
    vector<int> non_ds_nodes;

    for (int i = 0; i < n; ++i) {
        if (newDS[i]) {
            ds_nodes.push_back({(int)G.adj[i].size(), i});
        } else {
            non_ds_nodes.push_back(i);
        }
    }


    sort(ds_nodes.rbegin(), ds_nodes.rend());

    int pool_size = min(10, (int)ds_nodes.size());
    int actual_k1 = min(k1, pool_size);

    vector<int> remove_indices;
    for (int i = 0; i < pool_size; ++i) remove_indices.push_back(i);
    std::shuffle(remove_indices.begin(), remove_indices.end(), rng);

    for (int i = 0; i < actual_k1; ++i) {
        int v = ds_nodes[remove_indices[i]].second;
        newDS[v] = 0;
    }


    std::shuffle(non_ds_nodes.begin(), non_ds_nodes.end(), rng);
    int actual_k2 = min(k2, (int)non_ds_nodes.size());


    for (int i = 0; i < actual_k2; ++i) {
        newDS[non_ds_nodes[i]] = 1;
    }

    return newDS;
}




double runDNS(const Graph &G, vector<char> &inDS, vector<int> &cover, double curWeight,
              int maxIter, int perturb_threshold,
              double min_k1_rate, double max_k1_rate,
              double min_k2_rate, double max_k2_rate,
              double &time_to_best, int &iter_to_best,
              double setup_time, TimePoint search_start_time,
              double max_time_seconds) {
    int n = G.n;
    int iter = 0;
    vector<int> add_tabu(n, 0), delete_tabu(n, 0);
    int best_DS_size = n;



    int dynamic_k1_bonus = 0;
    int dynamic_k2_bonus = 0;

    int no_improve_count = 0;
    double bignumber = 100000.0;


    int current_undom = 0;
    for (int i = 0; i < n; ++i) {
        if (cover[i] == 0) current_undom++;
    }

    int initA1 = current_undom;
    curWeight = initA1 * bignumber + curWeight;
    double best_weight = curWeight;
    vector<char> best_inDS = inDS;


    time_to_best = setup_time;
    iter_to_best = 0;


    EliteTrie elite_trie(10000);
    elite_trie.insert(inDS);


    while (iter < maxIter) {
        ++iter;


        if (iter % 10 == 0) {
            auto current_check_time = chrono::high_resolution_clock::now();
            double elapsed_search_time = chrono::duration<double>(current_check_time - search_start_time).count();


            if (setup_time + elapsed_search_time >= max_time_seconds) {

                break;
            }
        }
        // ==========================================

        vector<int> candidates;
        vector<vector<char>> candidate_states;
        int best_A1 = INT_MAX;
        double bestA2 = INF;


        if (no_improve_count >= perturb_threshold) {
           // cout << ">>> [Perturb] Stuck at " << best_weight << ". Reconstructing..." << endl;


            vector<char> baseDS = elite_trie.get_random_solution(n);


            int base_size = 0;
            for (char c : baseDS) {
                if (c) base_size++;
            }
            cout <<"base_size="<<base_size<<endl;

            int base_k1 = max(1, (int)(base_size * min_k1_rate));
            int base_k2 = max(1, (int)(base_size * min_k2_rate));


            int max_k1 = max(1, (int)(base_size * max_k1_rate));
            int max_k2 = max(1, (int)(base_size * max_k2_rate));


            int actual_k1 = min(base_k1 + dynamic_k1_bonus, max_k1);
            int actual_k2 = min(base_k2 + dynamic_k2_bonus, max_k2);
            cout << "dynamic_k1_bonus=" << dynamic_k1_bonus << ", dynamic_k2_bonus=" << dynamic_k2_bonus << endl;

            inDS = mixed_perturb(G, baseDS, actual_k1, actual_k2);


            recompute_cover(G, inDS, cover);
            current_undom = 0;
            for (int i = 0; i < n; ++i) {
                if (cover[i] == 0) current_undom++;
            }

            fill(add_tabu.begin(), add_tabu.end(), 0);
            fill(delete_tabu.begin(), delete_tabu.end(), 0);


            curWeight = INF;
            no_improve_count = 0;

            continue;
        }



        for (int v = 0; v < n; ++v) if (inDS[v]) {

            vector<char> temp = inDS;
            vector<int> temp_cover = cover;
            int temp_undom = current_undom;


            flip_to_0(v, temp, temp_cover, temp_undom, G);


            repair_connectivity(G, temp, temp_cover, temp_undom);

            int a1 = temp_undom;
            if (a1 > best_A1) continue;

            double newW = a1 * bignumber + computeA2_by_rebuild(G, temp);

            if (delete_tabu[v] > iter && newW >= best_weight) continue;

            if (a1 < best_A1) {
                best_A1 = a1;
                bestA2 = newW;
                candidates.clear();
                candidate_states.clear();
                candidates.push_back(v);
                candidate_states.push_back(temp);
            } else if (a1 == best_A1) {
                if (newW == bestA2) {
                    candidates.push_back(v);
                    candidate_states.push_back(temp);
                }else if (newW < bestA2) {
                    bestA2 = newW;
                    candidates.clear();
                    candidate_states.clear();
                    candidates.push_back(v);
                    candidate_states.push_back(temp);
                }
            }
        }


        for (int v = 0; v < n; ++v) if (!inDS[v]) {


            vector<char> temp = inDS;
            vector<int> temp_cover = cover;
            int temp_undom = current_undom;



            flip_to_1(v, temp, temp_cover, temp_undom, G);
            int a1 = temp_undom;
            if (a1 > best_A1) continue;

            double newW = a1 * bignumber + computeA2_by_rebuild(G, temp);

            if (add_tabu[v] > iter && newW >= best_weight) continue;

            if (a1 < best_A1) {
                best_A1 = a1;
                bestA2 = newW;
                candidates.clear();
                candidate_states.clear();
                candidates.push_back(v);
                candidate_states.push_back(temp);
            } else if (a1 == best_A1) {
                if (newW == bestA2) {
                    candidates.push_back(v);
                    candidate_states.push_back(temp);
                }else if (newW < bestA2) {
                    bestA2 = newW;
                    candidates.clear();
                    candidate_states.clear();
                    candidates.push_back(v);
                    candidate_states.push_back(temp);
                }
            }
        }


        if (candidates.empty()) {

            no_improve_count = perturb_threshold;
            continue;
        }


        int chosen_idx = rand() % candidates.size();
        int best_vertex = candidates[chosen_idx];
        bool is_add = !inDS[best_vertex];   // 根据当前状态判断操作类型




        inDS = candidate_states[chosen_idx];
        recompute_cover(G, inDS, cover);

        current_undom = 0;
        for (int i = 0; i < n; ++i) {
            if (cover[i] == 0) current_undom++;
        }
        curWeight = bestA2;

        // 计算当前支配集 (DS) 和非支配集的顶点个数
        int current_ds_size = 0;
        for (int i = 0; i < n; ++i) {
            if (inDS[i]) current_ds_size++;
        }
        int non_ds_size = n - current_ds_size;


        int dynamic_tabu_len = 0;
        if (inDS[best_vertex] == 1) {

            dynamic_tabu_len = current_ds_size;
        } else {

            dynamic_tabu_len = non_ds_size;
        }



        add_tabu[best_vertex] = iter + dynamic_tabu_len;
        delete_tabu[best_vertex] = iter + dynamic_tabu_len;



        if (curWeight < best_weight - 1e-7) {
            best_weight = curWeight;
            best_inDS = inDS;

            elite_trie.clear();
            elite_trie.insert(inDS);
            no_improve_count = 0;


            dynamic_k1_bonus = 0;
            dynamic_k2_bonus = 0;


           if (!best_inDS.empty()) best_DS_size = accumulate(best_inDS.begin(), best_inDS.end(), 0);

            auto current_time = chrono::high_resolution_clock::now();
            double search_duration = chrono::duration<double>(current_time - search_start_time).count();
            time_to_best = setup_time + search_duration;
            iter_to_best = iter;

        }
        else if (abs(curWeight - best_weight) <= 1e-7) {

            bool is_new_elite = elite_trie.insert(inDS);
            if (!is_new_elite) {

                dynamic_k1_bonus++;
                dynamic_k2_bonus++;
            }

            no_improve_count++;
        }
        else {
            no_improve_count++;
        }


        cout << "Iter " << iter << " " << (is_add ? "add" : "remove")
             << " " << best_vertex << " -> MST " << curWeight
             << "           best_result  " << best_weight << "   size  "  <<  best_DS_size << " best solution number : "<< elite_trie.size()<< "\n";


    }

    inDS = elite_trie.get_random_solution(n);

    return best_weight;
}


void precompute_apsp(const Graph &G) {
    int n = G.n;
    dist_mat.assign(n, vector<double>(n, INF));
    next_node.assign(n, vector<int>(n, -1));

    for (int i = 0; i < n; ++i) {
        dist_mat[i][i] = 0;
        next_node[i][i] = i;
        for (const auto& edge : G.adj[i]) {
            int v = edge.first;
            double w = edge.second;
            if (w < dist_mat[i][v]) {
                dist_mat[i][v] = w;
                next_node[i][v] = v;
            }
        }
    }


    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (dist_mat[i][k] + dist_mat[k][j] < dist_mat[i][j]) {
                    dist_mat[i][j] = dist_mat[i][k] + dist_mat[k][j];
                    next_node[i][j] = next_node[i][k];
                }
            }
        }
    }
}


void printResult(const vector<char> &inDS, double finalWeight) {

    int dsSize = 0;
    if (!inDS.empty()) dsSize = accumulate(inDS.begin(), inDS.end(), 0);
    cout << "Final MST weight: " << finalWeight << "\n";
    cout << "Final DS size: " << dsSize << "\n";
}

struct PerturbRateConfig {
    double min_k1;
    double max_k1;
    double min_k2;
    double max_k2;
};


int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);


    unsigned int global_seed = static_cast<unsigned>(time(nullptr));
    // global_seed = 123456789;

    srand(global_seed);
    rng.seed(global_seed);


    string filename = "C:\\work\\ins_500_2.txt";

    // 本地调试建议把时间调短（比如 30 秒或 60 秒），避免等待过久
    double time_limit_seconds = 1000.0;
    int fixed_maxIter = 2000000;
    // 直接定义你要死磕的这一组参数
    int th = 200;
    double min_k1 = 0.05, max_k1 = 0.15;
    double min_k2 = 0.30, max_k2 = 0.50;
    Graph G = loadGraph(filename);
    auto setup_start = chrono::high_resolution_clock::now();
    cout << "Graph loaded: V=" << G.n << " E=" << G.edges.size() << "\n";


    precompute_apsp(G);


    auto setup_end = chrono::high_resolution_clock::now();
    double setup_time = chrono::duration<double>(setup_end - setup_start).count();

    vector<char> inDS;
    vector<int> DSlist;
    vector<int> tmp;
    vector<int> cover;

    double init_Weight = buildInitialSolution(G, inDS, DSlist, tmp, cover);
    int initialDSsize = accumulate(inDS.begin(), inDS.end(), 0);
    cout << "Initial DS size: " << initialDSsize << ", initial MST weight: " << init_Weight << "\n";


    double time_to_best = 0.0;
    int iter_to_best = 0;
    auto search_start = chrono::high_resolution_clock::now();

    double finalWeight = runDNS(G, inDS, cover, init_Weight, fixed_maxIter, th,
                                min_k1, max_k1, min_k2, max_k2,
                                time_to_best, iter_to_best, setup_time, search_start, time_limit_seconds);;
    printResult(inDS, finalWeight);

    return 0;
}
