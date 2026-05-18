#pragma once
#include "graph.hpp"
#include "commons.hpp"
#include "random"
#include "SolutionChecking.hpp"
class UnionFind{
    private:
        std::vector<int> parent;
        std::vector<int> rank;
        int n;
    public:
        UnionFind(int n):parent(n,0),rank(n,0),n(n){
            for(int i = 0; i< n;i++){
                parent[i] = i;
            }
        }
        int find(int x){
            int rep = x;
            while(parent[rep] != rep){
                rep = parent[rep];
            }
            int temp = x;
            while(parent[temp] != temp){
                parent[temp] = rep;
                temp = parent[temp];
            }
            return rep;
        }
        //Returns true if a merge occured, i.e. x and y were previously in different sets
        bool merge(int x, int y){
            int parent_x = find(x);
            int parent_y = find(y);
            if(parent_x == parent_y)return false;
            if(rank[parent_x] < rank[parent_y]){
                parent[parent_x] = parent_y;
                return true;
            }
            parent[parent_y] = parent_x;
            if(rank[parent_x] == rank[parent_y]){
                rank[parent_x]++;
            }
            return true;
        }
        void print_parents(){
            for(int v = 0; v < n;v++){
                std::cout << v << ":" << parent[v] << ",";
            }
            std::cout << "\n";
        }
        void clear(){
            for(int i = 0; i< n;i++){
                parent[i] = i;
            }
            rank = std::vector<int>(n,0);
        }   

};
class SpanningTreeHeuristic{
    private:
        Graph G;
        int n,m,k;
        double lambda;
        double offset;
        std::vector<std::vector<double>> weights;
        std::vector<std::vector<Vertex>> partition;
        double value;
        UnionFind union_find;
        bool has_run;
        bool value_computed;
        std::vector<std::pair<Vertex,Vertex>> compute_spanning_tree(bool complete_graph,bool full_tree, bool randomized = false, int seed = 0){
            std::vector<std::tuple<Vertex,Vertex,double>> sorted_edges;
            if(complete_graph){
                for(int u = 0; u < n; u++){
                    for(int v = 0; v < u; v++){
                        sorted_edges.push_back({u,v,get_dissimilarity(u,v)});
                    }
                }

            }
            else{
                for(int u = 0; u < n; u++){
                    for(auto [vit,end] = G.getNeighbours(u);vit != end ;vit++){
                        Vertex v = *vit;
                        if(v < u){
                            sorted_edges.push_back({u,v,get_dissimilarity(u,v)});
                        }
                    }
                }
            }
            if(!randomized){
                std::sort(sorted_edges.begin(),sorted_edges.end(), [](const auto &a, const auto &b){ 
                        return std::get<2>(a) > std::get<2>(b); 
                });
            }
            else{
                //std::random_device rd;
                std::mt19937 g(seed);
                std::shuffle(sorted_edges.begin(),sorted_edges.end(),g);

            }
            union_find.clear();
            std::vector<std::pair<Vertex,Vertex>> spanning_tree;
            while(spanning_tree.size() < (full_tree ? n-1:n-k)){
                auto largest_edge = sorted_edges.back();
                sorted_edges.pop_back();
                int u = std::get<0>(largest_edge);
                int v = std::get<1>(largest_edge);
                bool no_cycle = union_find.merge(u,v);
                if(no_cycle){
                    spanning_tree.push_back(std::make_pair(u,v));
                }
            }
            return spanning_tree;
        }
        void print_list(const std::vector<std::tuple<Vertex,Vertex,double>>& list){
            for(auto [x,y,z]:list){
                std::cout << "Edge: " << x << " "<< y <<" " << z << "\n"; 
            }
        }

        std::pair<std::vector<Vertex>,std::vector<Vertex>> compute_cc(const Graph& Tree, Vertex a, Vertex b){
            std::vector<Vertex> component_a;
            std::vector<bool> visited(n,false);
            std::queue<Vertex> bfs;
            visited[a] = true;
            bfs.push(a);
            while(!bfs.empty()){
                Vertex u = bfs.front();
                bfs.pop();
                component_a.push_back(u);
                for(auto [wit,end] = Tree.getNeighbours(u); wit != end;wit++){
                    Vertex w = *wit;
                    if(visited[w] || w == b) continue;
                    bfs.push(w);
                    visited[w] = true;
                }
            }
            std::vector<Vertex> component_b;
            bfs = {};
            visited[b] = true;
            bfs.push(b);
            while(!bfs.empty()){
                Vertex u = bfs.front();
                bfs.pop();
                component_b.push_back(u);
                for(auto [wit,end] = Tree.getNeighbours(u); wit != end;wit++){
                    Vertex w = *wit;
                    if(visited[w]) continue;
                    bfs.push(w);
                    visited[w] = true;
                }
            }
            return std::make_pair(component_a,component_b);
        }
        double bicut(const std::vector<Vertex>& set_a, const std::vector<Vertex>& set_b){
            double cut = 0;
            for(Vertex a:set_a){
                for(Vertex b:set_b){
                    cut += get_dissimilarity(a,b);
                }
            }
            return cut;
        }
        void partition_from_cc(const Graph& Tree){
            assert(Tree.numVertices() == n);
            assert(Tree.numEdges() == n-k);
            int partition_counter = 0;
            std::vector<bool> visited(n,false);
            for(Vertex c = 0; c < n; c++){
                if(visited[c])continue;
                std::queue<Vertex> bfs;
                bfs.push(c);
                visited[c] = true;
                while(!bfs.empty()){
                    Vertex u = bfs.front();
                    bfs.pop();
                    partition[partition_counter].push_back(u); 
                    for(auto [wit,end] = Tree.getNeighbours(u); wit != end;wit++){
                        Vertex w = *wit;
                        if(visited[w]) continue;
                        bfs.push(w);
                        visited[w] = true;
                        assert(G.adjacent(u,w));
                    }
                }
                partition_counter++;
            }
        }
        double get_dissimilarity(Vertex u, Vertex v){
            if(u == v) return 0;
            return (u > v) ? weights[u][v]:weights[v][u]; 
        }        
    public:
        SpanningTreeHeuristic(const Graph& G, int k,double lambda):
        G(G),n(G.numVertices()),m(G.numEdges()),k(k),lambda(lambda),weights(n,std::vector<double>()),
        partition(k,std::vector<Vertex>()),value(0),union_find(n),has_run(false),value_computed(false){
            //Compute weights of weighted sum of vertex dissimilarity and border length
            for(Vertex u = 1; u < n;u++){
                weights[u] = std::vector<double>(u,0);
                for(Vertex v = 0; v < u;v++){
                    weights[u][v] = lambda*std::abs(G.getAttribute(u).feature - G.getAttribute(v).feature);
                }
            }
            for(Vertex u = 0; u < n;u++){
                for(auto[vit,end] = G.getNeighbours(u);vit != end;vit++){
                    Vertex v = *vit;
                    if(u < v) continue;
                    weights[u][v] += (lambda - 1) * G.borderLength(u,v);
                    offset += (1 - lambda) * G.borderLength(u,v);
                }
            }
        }
        void compute_solution(bool dp_cut = false, bool random = false, int seed = 0){
            std::vector<std::pair<Vertex,Vertex>> spanning_tree = compute_spanning_tree(false, dp_cut,random,seed);
            if(!dp_cut){
                std::map<int,int> partition_id;
                int id_counter = 0;
                for(Vertex v = 0; v < n; v++){
                    int rep = union_find.find(v);
                    if(partition_id.find(rep) == partition_id.end()){
                        partition_id[rep] = id_counter;
                        id_counter++;
                    }
                    partition[partition_id[rep]].push_back(v); 
                }
            }
            else{
                Graph Tree(n);
                for(auto [u,v]:spanning_tree){
                    Tree.addEdge(u,v);
                }
                assert(Tree.is_connected());
                while(Tree.numEdges() > n-k){
                    std::pair<Vertex,Vertex> best_pair;
                    double best_value = -1E10;
                    for(int u = 0; u < n; u++){
                        for(auto [vit,end] = Tree.getNeighbours(u);vit != end;vit++){
                            Vertex v = *vit;
                            if(u < v) continue;
                            auto [comp_u, comp_v] = compute_cc(Tree,u,v);
                            double cut_value = bicut(comp_u,comp_v);
                            if(cut_value > best_value){
                                cut_value = best_value;
                                best_pair = {u,v};
                            }
                        }
                    }
                    Tree.removeEdge(best_pair.first,best_pair.second);
                }
                partition_from_cc(Tree);
            }
            has_run = true;
            int total_vertices = 0;
            for(int i = 0; i < k;i++){
                total_vertices += partition[i].size();
            }
            assert(total_vertices == n);
            assert(checkConnectivity(G,partition_to_assignment(partition),k));
            assert(check_valid_partitioning(partition));
        }
        void compute_solution_noconnect(bool dp_cut = false){
            std::vector<std::pair<Vertex,Vertex>> spanning_tree = compute_spanning_tree(true, dp_cut);
            if(!dp_cut){
                std::map<int,int> partition_id;
                int id_counter = 0;
                for(Vertex v = 0; v < n; v++){
                    int rep = union_find.find(v);
                    if(partition_id.find(rep) == partition_id.end()){
                        partition_id[rep] = id_counter;
                        id_counter++;
                    }
                    partition[partition_id[v]].push_back(v); 
                }
            }
            has_run = true;
        }

        std::vector<std::vector<Vertex>> get_solution(){
            if(!has_run){
                std::cout << "Algorithm has not run yet \n";
            }
            return partition;
        }

        double get_value(){
            if(!has_run){
                std::cout << "Algorithm has not run yet \n";
            }
            if(!value_computed){
                value = 0;
                for(auto& region: partition){
                    for(Vertex u: region){
                        for(Vertex v: region){
                            if(u<v)value += get_dissimilarity(u,v);
                        }
                    }
                }
                value += offset;
                value_computed = true;
            }
            return value;
        }
        
};

class KernighanLin{
    private: 
        Graph G;
        int n,m,k;
        double lambda;
        double offset = 0;
        std::vector<std::vector<double>> weights;
        std::vector<std::vector<Vertex>> partition;
        std::vector<int> partition_map;
        std::vector<int> partition_sizes;
        std::vector<std::vector<double>> vxpt_cost; //Total weights between vertex and partition
        double value;
        bool has_run,value_computed;
        std::vector<bool> is_artic; //Maintain which vertices are articulation vertices of current partitions
        //Update articulation vertices of partition prt
        void update_articulation(int prt){
            std::vector<Vertex> stack;
            std::vector<int> discovery(n,-1);
            std::vector<int> parent(n,-1);
            int time = 1;
            std::vector<int> early_reach(n,-1); // Earliest vertex reachable by back edge
            Vertex r = 0;
            while(partition_map[r] != prt) r++; // Find arbitrary vertex in partition prt
            //Run DFS to obtain discovery times
            discovery[r] = 0;
            int root_children = 0;
            stack.push_back(r);
            while(! stack.empty()){
                Vertex v = stack.back();
                is_artic[v] = false;
                early_reach[v] = discovery[v];
                bool new_tree = false;
                for(auto [wit,end] = G.getNeighbours(v);wit != end;wit++){
                    Vertex w = *wit;
                    if(partition_map[w] != prt) continue;
                    if(parent[v] == w) continue; //Tree edge from other side
                    if(discovery[w] == -1){
                        //Newly discovered tree edge {v,w}
                        discovery[w] = time;
                        parent[w] = v;
                        stack.push_back(w);
                        time++;
                        if (v == r)root_children++;
                        new_tree = true;
                        break;
                    }
                    else if(parent[w] == v){
                        //Previously discovered tree edge {v,w}
                        if(early_reach[w] < early_reach[v]) early_reach[v] = early_reach[w];
                        if(discovery[v] <= early_reach[w]) {
                            is_artic[v] = true;
                            if(v == 6){
                                is_artic[v] = true;
                            }
                        }
                    }
                    else{
                        //Back edge {v,w}
                        if(discovery[w] < early_reach[v]) early_reach[v] = discovery[w];
                    }
                }
                if(new_tree) continue;
                stack.pop_back();
            }
            is_artic[r] = (root_children > 1);

        }

        void update_articulation2(int prt){
            std::vector<Vertex> stack;
            std::vector<int> discovery(n,-1);
            std::vector<int> parent(n,-1);
            std::vector<int> last_called(n,-1);
            int time = 0;
            //std::vector<int> early_reach(n,-1); // Earliest vertex reachable by back edge
            Vertex r = 0;
            while(partition_map[r] != prt) r++; // Find arbitrary vertex in partition prt
            //Run DFS to obtain discovery times
            stack.push_back(r);
            while(! stack.empty()){
                Vertex v = stack.back();
                stack.pop_back();
                parent[v] = last_called[v];
                discovery[v] = time;
                time++;
                for(auto [wit,end] = G.getNeighbours(v);wit != end;wit++){
                    Vertex w = *wit;
                    if(partition_map[w] != prt) continue;
                    if(discovery[w] != -1) continue;
                    if(last_called[w] == -1)stack.push_back(w);
                    last_called[w] = v;
                }
            }
        }
        double get_dissimilarity(Vertex u, Vertex v){
            if(u == v ) return 0;
            return (u > v) ? weights[u][v]:weights[v][u]; 
        }        
    public:
        KernighanLin(const Graph& G, int k,double lambda,const std::vector<std::vector<Vertex>>& partition):
        G(G),n(G.numVertices()),m(G.numEdges()),k(k),lambda(lambda),weights(n,std::vector<double>()),
        partition(partition),partition_map(n,0),partition_sizes(k,0),vxpt_cost(n,std::vector<double>(k,0)),
        value(0),has_run(false),value_computed(false),is_artic(n){
            //Compute weights of weighted sum of vertex dissimilarity and border length
            for(Vertex u = 1; u < n;u++){
                weights[u] = std::vector<double>(u,0);
                for(Vertex v = 0; v < u;v++){
                    weights[u][v] = lambda*std::abs(G.getAttribute(u).feature - G.getAttribute(v).feature);
                }
            }
            for(Vertex u = 0; u < n;u++){
                for(auto[vit,end] = G.getNeighbours(u);vit != end;vit++){
                    Vertex v = *vit;
                    if(u < v) continue;
                    weights[u][v] += (lambda - 1) * G.borderLength(u,v);
                    offset += (1 - lambda) * G.borderLength(u,v);
                }
            }
            
            for(int i = 0; i < k; i++){
                for(Vertex w:partition[i]){
                    partition_map[w] = i;
                    partition_sizes[i]++;
                    for(Vertex v = 0; v < n; v++){
                        vxpt_cost[v][i] += get_dissimilarity(v,w);
                    }
                }
            }
            //assert that initial partition is feasible
            assert(checkConnectivity(G,partition_map,k));
            //update_articulation2(0);
            for(int i = 0; i< k;i++){
                update_articulation(i);
            }

        }


        void improve(int iterations, bool verbose = false){

            Vertex min_vertex = 0;
            int origin = 0;
            int target = 0;
            double min_delta = 0;
            for(Vertex v = 0; v < n; v++){
                int part_v = partition_map[v];
                if(is_artic[v] || partition_sizes[part_v] == 1) continue; //Cannot reassign vertices that are artic. or isolated
                for(auto [wit,end] = G.getNeighbours(v); wit != end;wit++){
                    Vertex w = *wit;
                    int part_w = partition_map[w];
                    if(part_v == part_w) continue;
                    double delta = vxpt_cost[v][part_w]-vxpt_cost[v][part_v];
                    if(delta < min_delta){
                        min_delta = delta;
                        origin = part_v;
                        target = part_w;
                        min_vertex = v;
                    }
                }
            }
            for(int i = 0; i < iterations;i++){
                if(min_delta >= -1E-2){
                    if(verbose)std::cout << "No improvement possible, Iteration Count:" << i << "\n";
                    break;
                }
                if(verbose){
                    std::cout << "Reassigning " << min_vertex << " from " << origin << " to "<< target << ", Improvement: "<< -min_delta<<"\n"; 
                }
                partition_sizes[partition_map[min_vertex]]--;
                partition_sizes[target]++;
                partition_map[min_vertex] = target;
                //Update articulation vertices of changed partitions
                update_articulation(origin);
                update_articulation(target);
                for(Vertex v = 0; v < n; v++){
                    vxpt_cost[v][origin] -= get_dissimilarity(v,min_vertex);
                }

                for(Vertex v = 0; v < n; v++){
                    vxpt_cost[v][target] += get_dissimilarity(v,min_vertex);
                }

                //This can be optimized by not checking all reassignments again but storing them the first time and only checking recomputed
                min_delta = 0;
                for(Vertex v = 0; v < n; v++){
                    int part_v = partition_map[v];
                    if(is_artic[v] || partition_sizes[part_v] == 1) continue; //Cannot reassign articulation vertices as partition becomes disconnected
                    for(auto [wit,end] = G.getNeighbours(v); wit != end;wit++){
                        Vertex w = *wit;
                        int part_w = partition_map[w];
                        if(part_v == part_w) continue;
                        double delta = vxpt_cost[v][part_w]-vxpt_cost[v][part_v];
                        if(delta < min_delta){
                            min_delta = delta;
                            origin = part_v;
                            target = part_w;
                            min_vertex = v;
                        }
                    }
                }
                
            }
            partition = std::vector<std::vector<Vertex>>(k,std::vector<Vertex>());
            for(Vertex v = 0; v < n ;v++){
                partition[partition_map[v]].push_back(v);
            }
            assert(checkConnectivity(G,partition_map,k));
            assert(check_valid_partitioning(partition));
            has_run = true;
        }

        std::vector<std::vector<Vertex>> get_solution(){
            if(!has_run){
                std::cout << "Algorithm has not run yet \n";
            }
            return partition;
        }

        double get_value(){
            if(!has_run){
                std::cout << "Algorithm has not run yet \n";
            }
            if(!value_computed){
                value = 0;
                for(auto& region: partition){
                    for(Vertex u: region){
                        for(Vertex v: region){
                            if(u<v)value += get_dissimilarity(u,v);
                        }
                    }
                }
                value += offset;
                value_computed = true;
            }
            return value;
        }

};

double pair_cluster_cost(const Graph& G, Vertex a, Vertex b, Vertex c){
    std::queue<Vertex> bfs = {};
    const int n = G.numVertices();
    std::vector<bool> visited(n,false); 
    std::vector<Vertex> parent(n,n); 
    visited[a] = true;
    bfs.push(a);
    while(!bfs.empty()){
        Vertex u = bfs.front();
        bfs.pop();
        for(auto [wit,end] = G.getNeighbours(u); wit != end;wit++){
            Vertex w = *wit;
            if(visited[w]) continue;
            bfs.push(w);
            parent[w] = u;
            visited[w] = true;
        }
        if(visited[b] && visited[c])break;
    }
    double total_cost = 0;
    std::vector<Vertex> cluster = {b};
    Vertex on_path = parent[b];
    while(on_path != n){
        for(Vertex v:cluster){
            total_cost += std::abs(G.getAttribute(on_path).feature - G.getAttribute(v).feature);
        }
        cluster.push_back(on_path);
        on_path = parent[on_path];
    }
    on_path = c;
    while(std::find(cluster.begin(),cluster.end(),c) == cluster.end()){
        for(Vertex v:cluster){
            total_cost += std::abs(G.getAttribute(on_path).feature - G.getAttribute(v).feature);
        }
        cluster.push_back(on_path);
        on_path = parent[on_path];
    }
    return total_cost;

}

double pair_cluster_cost(const Graph& G, Vertex a, Vertex b){
    std::queue<Vertex> bfs = {};
    const int n = G.numVertices();
    std::vector<bool> visited(n,false); 
    std::vector<Vertex> parent(n,n); 
    visited[a] = true;
    bfs.push(a);
    while(!bfs.empty()){
        Vertex u = bfs.front();
        bfs.pop();
        for(auto [wit,end] = G.getNeighbours(u); wit != end;wit++){
            Vertex w = *wit;
            if(visited[w]) continue;
            bfs.push(w);
            parent[w] = u;
            visited[w] = true;
        }
        if(visited[b])break;
    }
    double total_cost = 0;
    std::vector<Vertex> cluster = {b};
    Vertex on_path = parent[b];
    while(on_path != n){
        for(Vertex v:cluster){
            total_cost += std::abs(G.getAttribute(on_path).feature - G.getAttribute(v).feature);
        }
        cluster.push_back(on_path);
        on_path = parent[on_path];
    }
    return total_cost;

}
class RandomKLHeuristic{
    private:
        Graph graph;
        int k;
        int iter;
        int seed;
        double best_obj = 1E9;
        std::vector<std::vector<Vertex>> best_partition;
    public:
        RandomKLHeuristic(const Graph& graph, int k, int iter = 100, int seed = 0):graph(graph),k(k),iter(iter),seed(seed){}
        void run(){
            std::mt19937 rng(seed);
            for(int i = 1; i <= 100;i++){
                SpanningTreeHeuristic random_mst(graph,k,1);
                random_mst.compute_solution(true,true,rng());
                KernighanLin local_search(graph,k,1,random_mst.get_solution());
                local_search.improve(100,false);
                auto solution = local_search.get_solution();
                if(local_search.get_value() < best_obj){
                    best_obj = local_search.get_value();
                    best_partition = local_search.get_solution();
                }
            }
        }
        double get_objective(){
            return best_obj;
        }
};
class ClusterDijkstra{
    private:
        Graph graph;
        int k;
        std::vector<double> distance;
        std::vector<Vertex> predecessor;
        double uvdist;
        double path_vertex_cut(Vertex a, Vertex b){
            Vertex pivot = a;
            double cut_value = 0;
            while(pivot != graph.numVertices()){
                cut_value += get_dissimilarity(pivot,b);
                pivot = predecessor[pivot];
            }
            return cut_value;
        }
        double get_dissimilarity(Vertex u, Vertex v){
            return std::abs(graph.getAttribute(u).feature - graph.getAttribute(v).feature); 
        }   
    public:
        ClusterDijkstra(const Graph& graph, int k):graph(graph),k(k),distance(graph.numVertices(),1E10),predecessor(graph.numVertices(),graph.numVertices()){}
        void run(Vertex u ,Vertex v){
            distance = std::vector<double>(graph.numVertices(),1E10);
            predecessor = std::vector<Vertex>(graph.numVertices(),graph.numVertices());
            std::priority_queue<std::pair<double,Vertex>> q;
            q.push({0,u});
            distance[u] = 0;
            while(!q.empty()){
                auto[prio,a] = q.top();
                q.pop();
                if(a == v){
                    uvdist = distance[v];
                    return;
                }
                if(prio > distance[a]) continue; // Duplicates in queue
                for(auto [it,end] = graph.getNeighbours(a); it != end;it++){
                    Vertex b = *it;
                    double pvc = distance[a] + path_vertex_cut(a,b);
                    if(pvc < distance[b]){
                        distance[b] = pvc;
                        predecessor[b] = a;
                        q.push({pvc,b});
                    }
                }
            }
        }
        double get_objective(){
            return uvdist;
        }
};