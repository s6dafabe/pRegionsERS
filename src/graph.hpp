#pragma once

#include <vector>
#include <map>
#include <limits>
#include "commons.hpp"
#include <utility>
#include "boost/graph/strong_components.hpp"
#include "boost/graph/push_relabel_max_flow.hpp"
#include "boost/graph/edmonds_karp_max_flow.hpp"
#include "boost/graph/boykov_kolmogorov_max_flow.hpp"
#include <queue>
#include <algorithm>
#include <random>
class Graph{
    private:
        BoostGraph G;
        enum VType{In = 0,Out = 1};
        int split_idx(Vertex v,VType t)const{
            return v*2+t;
        }
        int edge_idx(Vertex u, Vertex v,const std::map<std::pair<Vertex,Vertex>,int>& e_map)const{
            return numVertices()+e_map.at({u,v});
        }
    public:
        Graph(int n = 0):G(n){}

        void addVertex(const VertexProperties& p = {"",1,0,0,0}){
            boost::add_vertex(p,G);
        }
        void removeVertex(Vertex v){
            boost::clear_vertex(v,G);
            boost::remove_vertex(v,G);
        }
        void addEdge(Vertex u, Vertex v,double borderLength = 1){
            boost::add_edge(u,v,{borderLength},G);
        }
        void removeEdge(Vertex u, Vertex v){
            if(! adjacent(u,v)){
                std::cout << "Edge {"<<u<<","<<v<<"} does not exist \n";
            }
            boost::remove_edge(u,v,G);
        }
        std::pair<adjacency_iterator, adjacency_iterator> getNeighbours(Vertex u)const{
            return boost::adjacent_vertices(u,G);
        }
        std::pair<edge_iterator, edge_iterator> getEdges()const{
            return boost::edges(G);
        }
        std::pair<Vertex,Vertex> get_incident(const Edge& e)const{
            return {boost::source(e,G),boost::target(e,G)};
        }
        bool adjacent(Vertex u, Vertex v)const{
            for(auto [start,end] = boost::adjacent_vertices(u,G);start !=end;start++){
                if(*start == v)return true;
            }
            return false;
        }
        double borderLength(Vertex u, Vertex v) const{
            for(auto [start,end] = boost::out_edges(u,G);start !=end;start++){
                if(boost::target(*start,G) == v){
                    return G[*start].borderLength;
                }
            }
            return std::numeric_limits<double>::infinity();
        }
        int numVertices()const{
            return boost::num_vertices(G);
        }
        int numEdges()const{
            return boost::num_edges(G);
        }
        VertexProperties getAttribute(Vertex v) const{
            return G[v];
        }
        void setAttribute(Vertex v, VertexProperties f){
            G[v] = f;
        }
        int getDegree(Vertex v) const{
            return boost::degree(v,G);
        }
        bool is_connected(){
            if(numVertices() == 0) return true;
            std::vector<bool> visited(numVertices(),false);
            std::queue<Vertex> bfs;
            bfs.push(0);
            visited[0] = true;
            while(! bfs.empty()){
                Vertex v = bfs.front();
                bfs.pop();
                for(auto [wit,end] = getNeighbours(v);wit != end;wit++){
                    Vertex w = *wit;
                    if(visited[w]) continue;
                    visited[w] = true;
                    bfs.push(w);
                }
            }
            for(bool isvisit: visited){
                if(!isvisit)return false;
            }
            return true;

        }

        void permute_indices(int seed = 0){
            std::mt19937 g(seed);
            std::vector<int> shuffled_idx(numVertices());
            std::iota(shuffled_idx.begin(),shuffled_idx.end(),0);
            std::shuffle(shuffled_idx.begin(),shuffled_idx.end(),g);
            BoostGraph G_shuffled(numVertices());
            for(Vertex v = 0; v < numVertices();v++){
                G_shuffled[shuffled_idx[v]] = G[v];
            }
            for(auto [eit,end] = boost::edges(G); eit != end;eit++){
                Edge e = *eit;
                Vertex a = boost::source(e,G);
                Vertex b = boost::target(e,G);
                boost::add_edge(shuffled_idx[a],shuffled_idx[b],G[e],G_shuffled);
            }
            G = G_shuffled;
            
        }

        std::pair<double,std::vector<Vertex>> minimum_separator(Vertex u, Vertex v, const std::vector<double>& flow_values) const{
            FlowGraph flow_graph(2*numVertices());

            for(auto[eit,end] = boost::edges(G);eit!=end;eit++){
                Edge ex = *eit;
                Vertex a = boost::source(ex,G);
                Vertex b = boost::target(ex,G);
                if((u==a && v ==b)||(u==b&&v==a)){
                    std::vector<Vertex> empty = {};
                    return std::make_pair(numVertices(),empty);
                }
                FlowEdge e1 = boost::add_edge(split_idx(a,Out),split_idx(b,In),flow_graph).first;
                FlowEdge e2 = boost::add_edge(split_idx(b,Out),split_idx(a,In),flow_graph).first;
                FlowEdge e1r = boost::add_edge(split_idx(b,In),split_idx(a,Out),flow_graph).first;
                FlowEdge e2r = boost::add_edge(split_idx(a,In),split_idx(b,Out),flow_graph).first;
                flow_graph[e1].capacity = numVertices();
                flow_graph[e2].capacity = numVertices();
                flow_graph[e1r].capacity = 0;
                flow_graph[e2r].capacity = 0;
                flow_graph[e1].reverse_edge = e1r;
                flow_graph[e1r].reverse_edge = e1;
                flow_graph[e2].reverse_edge = e2r;
                flow_graph[e2r].reverse_edge = e2;
            }
            for(Vertex u = 0; u < numVertices(); u++){
                FlowEdge e = boost::add_edge(split_idx(u,In),split_idx(u,Out),flow_graph).first;
                FlowEdge er = boost::add_edge(split_idx(u,Out),split_idx(u,In),flow_graph).first;
                flow_graph[e].capacity = flow_values[u];
                flow_graph[er].capacity = 0;
                flow_graph[e].reverse_edge = er;
                flow_graph[er].reverse_edge = e;
            }
            /*for(auto[eit,end] = boost::edges(flow_graph);eit!=end;eit++){
                FlowEdge ex = *eit;
                std::cout << ex <<" "<<flow_graph[ex].capacity<< "\n";
            }
            std::cout <<"---------------------------\n";*/
            std::vector<FlowEdge> pred_storage(2*numVertices());
            std::vector<boost::default_color_type> color_storage(2*numVertices());
            std::vector<long> distance_storage(2*numVertices());
            auto index_map = boost::get(boost::vertex_index, flow_graph);
            auto predecessor_map = boost::make_iterator_property_map(pred_storage.begin(), index_map);
            auto color_map = boost::make_iterator_property_map(color_storage.begin(), index_map);
            auto distance_map = boost::make_iterator_property_map(distance_storage.begin(), index_map);
            double flow = boost::boykov_kolmogorov_max_flow(flow_graph,split_idx(u,Out),split_idx(v,In),
                capacity_map(boost::get(&FlowProperties::capacity,flow_graph)).
                residual_capacity_map(boost::get(&FlowProperties::residual_capacity,flow_graph)).
                reverse_edge_map(boost::get(&FlowProperties::reverse_edge,flow_graph)).
                predecessor_map(predecessor_map).
                color_map(color_map).
                distance_map(distance_map)
            );
            //Extracting the cut from The flow
            //First Step: BFS on the residual graph to determine the u-component
            std::queue<Vertex> q;
            std::vector<bool> visited(2*numVertices(),false);
            q.push(split_idx(u,Out));
            visited[split_idx(u,Out)] = true;
            while(!q.empty()){
                Vertex v = q.front();
                q.pop();
                for(auto [eit,end] = boost::out_edges(v,flow_graph);eit != end;eit++){
                    FlowEdge e = *eit;
                    Vertex w = boost::target(e,flow_graph);
                    if(visited[w] || flow_graph[e].residual_capacity < ERR_TOLERANCE) continue;
                    visited[w] = true;
                    q.push(w);
                }
            }
            //Second Step: The cut edges (going from S to T component) induce the seperator vertices
            std::unordered_set<Vertex> seperator;
            for(Vertex w = 0; w < numVertices();w++){
                if(w == u || w == v)continue;
                
                if(visited[split_idx(w,In)] != visited[split_idx(w,Out)]){
                    seperator.insert(w);
                }
                    
                
            }
            /*We want a minimal (minimal wrt. to vertex set) separator. 
            Not all edges with residual capacity of 0 are necessarily part of a minimal cut, as they can simply have 0 capacity
            We therefore perform a BFS from v on G/S and only include vertices from S that are adjacent to a vertex in 
            the CC of v as well as v
            */
            std::vector<Vertex> min_seperator;
            q = {};
            std::fill(visited.begin(), visited.end(), false);
            q.push(v);
            visited[v] = true;
            while(!q.empty()){
                Vertex v = q.front();
                q.pop();
                for(auto[vit,end] = getNeighbours(v);vit != end;vit++){
                    Vertex w = *vit;
                    if(visited[w])continue;
                    visited[w] = true;
                    if(seperator.find(w) != seperator.end()){
                        min_seperator.push_back(w);
                        continue;
                    }
                    q.push(w);
                }
            }
            return std::make_pair(flow,min_seperator);
        }
        std::pair<double,std::vector<std::pair<Vertex,Vertex>>> minimum_cut2(Vertex u, Vertex v, const std::vector<std::vector<double>>& flow_values, bool pr = true) const{
            FlowGraph flow_graph(numVertices());
            for(auto[eit,end] = boost::edges(G);eit!=end;eit++){
                Edge ex = *eit;
                Vertex a = boost::source(ex,G);
                Vertex b = boost::target(ex,G);
                if((u==a && v ==b)||(u==b&&v==a)){
                    std::vector<std::pair<Vertex,Vertex>> empty = {};
                    return std::make_pair(numVertices(),empty);
                }
                FlowEdge uv = boost::add_edge(a,b,flow_graph).first;
                FlowEdge vu = boost::add_edge(b,a,flow_graph).first;
                FlowEdge uv_rev = boost::add_edge(b,a,flow_graph).first;
                FlowEdge vu_rev = boost::add_edge(a,b,flow_graph).first;
                flow_graph[uv].capacity = (a > b) ? flow_values[a][b]:flow_values[b][a];
                flow_graph[uv].residual_capacity = flow_graph[uv].capacity;
                flow_graph[uv_rev].capacity = 0;
                flow_graph[uv_rev].residual_capacity = 0;
                flow_graph[vu].capacity = (a > b) ? flow_values[a][b]:flow_values[b][a];
                flow_graph[vu].residual_capacity = flow_graph[vu].capacity;
                flow_graph[vu_rev].capacity = 0;
                flow_graph[vu_rev].residual_capacity = 0;
                //flow_graph[uv].reverse_edge = vu;
                //flow_graph[vu].reverse_edge = uv;
                flow_graph[uv].reverse_edge = uv_rev;
                flow_graph[uv_rev].reverse_edge = uv;
                flow_graph[vu].reverse_edge = vu_rev;
                flow_graph[vu_rev].reverse_edge = vu;
            }
            assert(boost::num_edges(flow_graph) == 4*numEdges());
            double flow = 0;
            if(pr)flow = boost::push_relabel_max_flow(flow_graph,u,v,
                capacity_map(boost::get(&FlowProperties::capacity,flow_graph)).
                residual_capacity_map(boost::get(&FlowProperties::residual_capacity,flow_graph)).
                reverse_edge_map(boost::get(&FlowProperties::reverse_edge,flow_graph))
            );
            else flow = boost::edmonds_karp_max_flow(flow_graph,u,v,
                capacity_map(boost::get(&FlowProperties::capacity,flow_graph)).
                residual_capacity_map(boost::get(&FlowProperties::residual_capacity,flow_graph)).
                reverse_edge_map(boost::get(&FlowProperties::reverse_edge,flow_graph))
            );
            //Extracting the cut from The flow
            //First Step: BFS on the residual graph to determine the u-component
            std::queue<Vertex> q;
            std::vector<bool> visited(numVertices(),false);
            std::vector<Vertex> parent(numVertices(),numVertices());
            q.push(u);
            visited[u] = true;
            while(!q.empty()){
                Vertex a = q.front();
                q.pop();
                for(auto [eit,end] = boost::out_edges(a,flow_graph);eit != end;eit++){
                    FlowEdge e = *eit;
                    Vertex b = boost::target(e,flow_graph);
                    if(visited[b] || flow_graph[e].residual_capacity < ERR_TOLERANCE) continue;
                    visited[b] = true;
                    parent[b] = a;
                    q.push(b);
                }
            }
            //Second Step: The cut edges (going from S to T component) induce the seperator vertices
            std::vector<std::pair<Vertex,Vertex>> cut_edges;
            for(auto[eit,end] = boost::edges(G);eit!=end;eit++){
                Edge ex = *eit;
                Vertex a = boost::source(ex,G);
                Vertex b = boost::target(ex,G);
                if(visited[a] != visited[b]){
                    cut_edges.push_back(std::make_pair(a,b));
                }
            }
            
            return std::make_pair(flow,cut_edges);
        }


       
        std::pair<double,std::vector<std::pair<Vertex,Vertex>>> minimum_cut(Vertex u, Vertex v, const std::vector<std::vector<double>>& flow_values) const{
            BKGraph flow_graph(numVertices());
            for(auto[eit,end] = boost::edges(G);eit!=end;eit++){
                Edge ex = *eit;
                Vertex a = boost::source(ex,G);
                Vertex b = boost::target(ex,G);
                if((u==a && v ==b)||(u==b&&v==a)){
                    std::vector<std::pair<Vertex,Vertex>> empty = {};
                    return std::make_pair(numVertices(),empty);
                }
                FlowEdge uv = boost::add_edge(a,b,flow_graph).first;
                FlowEdge uv_rev = boost::add_edge(b,a,flow_graph).first;
                flow_graph[uv].capacity = (a > b) ? flow_values[a][b]:flow_values[b][a];
                flow_graph[uv].residual_capacity = flow_graph[uv].capacity;
                flow_graph[uv_rev].capacity = flow_graph[uv].capacity;
                flow_graph[uv_rev].residual_capacity = flow_graph[uv].capacity;
                flow_graph[uv].reverse_edge = uv_rev;
                flow_graph[uv_rev].reverse_edge = uv;
            }
            assert(boost::num_edges(flow_graph) == 2*numEdges());
            std::vector<FlowEdge> pred_storage(numVertices());
            std::vector<boost::default_color_type> color_storage(numVertices());
            std::vector<long> distance_storage(numVertices());
            auto index_map = boost::get(boost::vertex_index, flow_graph);
            auto predecessor_map = boost::make_iterator_property_map(pred_storage.begin(), index_map);
            auto color_map = boost::make_iterator_property_map(color_storage.begin(), index_map);
            auto distance_map = boost::make_iterator_property_map(distance_storage.begin(), index_map);
            double flow = boost::boykov_kolmogorov_max_flow(flow_graph,u,v,
                capacity_map(boost::get(&FlowProperties::capacity,flow_graph)).
                residual_capacity_map(boost::get(&FlowProperties::residual_capacity,flow_graph)).
                reverse_edge_map(boost::get(&FlowProperties::reverse_edge,flow_graph)).
                predecessor_map(predecessor_map).
                color_map(color_map).
                distance_map(distance_map)
            );

            //Extracting the cut from The flow
            //First Step: BFS on the residual graph to determine the u-component
            std::queue<Vertex> q;
            std::vector<bool> visited(numVertices(),false);
            std::vector<Vertex> parent(numVertices(),numVertices());
            q.push(u);
            visited[u] = true;
            while(!q.empty()){
                Vertex a = q.front();
                q.pop();
                for(auto [eit,end] = boost::out_edges(a,flow_graph);eit != end;eit++){
                    FlowEdge e = *eit;
                    Vertex b = boost::target(e,flow_graph);
                    if(visited[b] || flow_graph[e].residual_capacity < ERR_TOLERANCE) continue;
                    visited[b] = true;
                    parent[b] = a;
                    q.push(b);
                }
            }
            
            //Second Step: The cut edges (going from S to T component) induce the seperator vertices
            std::vector<std::pair<Vertex,Vertex>> cut_edges;
            for(auto[eit,end] = boost::edges(G);eit!=end;eit++){
                Edge ex = *eit;
                Vertex a = boost::source(ex,G);
                Vertex b = boost::target(ex,G);
                if(visited[a] != visited[b]){
                    cut_edges.push_back(std::make_pair(a,b));
                }
            }
            
            return std::make_pair(flow,cut_edges);
        }

    };

class StronglyConnectedComponents{
    private:
        DirectedGraph directed;
        int num = 0;
        std::vector<int> components;
    public:
        StronglyConnectedComponents(std::vector<std::pair<int,int>> activeEdges,int n):directed(n),components(n){
            for(auto [u,v]:activeEdges){
                boost::add_edge(u,v,directed);
            }
        }
        void run(){
            num = boost::strong_components(directed, &components[0]);
        }
        int getNumComponents() const{
            return num;
        }
        std::vector<int> getComponents() const{
            return components;
        }
        std::vector<int> getComponentSizes() const{
            std::vector<int> sizes(num);
            for(int c:components){
                sizes[c]++;
            }
            return sizes;
        }
};


class DenseSubgraph{
    private:
        FlowGraph flow_graph;
        int n;
        double previous_degree;
        Vertex previous_vertex;
        const Vertex source;
        const Vertex sink;
        double total_assign;
        double difference;
        std::vector<Vertex> source_component;
        void clear(){
            for(auto [eit,end] = boost::edges(flow_graph);eit != end;eit++){
                FlowEdge e = *eit;
                flow_graph[e].residual_capacity = flow_graph[e].capacity;
            }
            if(previous_vertex != sink+1){
                FlowEdge sv = boost::edge(source,previous_vertex,flow_graph).first;
                flow_graph[sv].capacity = previous_degree;
                flow_graph[sv].residual_capacity = previous_degree;
            }
        }
        double flow;
    public:
        DenseSubgraph(int n, const std::vector<std::tuple<Vertex,Vertex,double>>& edges):
        flow_graph(n+2),n(n),previous_vertex(n+2),source(n),sink(n+1),total_assign(0){
            std::vector<double> degree(n);
            for(auto [u,v,assign]:edges){
                FlowEdge uv = boost::add_edge(u,v,flow_graph).first;
                FlowEdge uv_rev = boost::add_edge(v,u,flow_graph).first;
                flow_graph[uv].capacity = assign/2;
                flow_graph[uv].residual_capacity = flow_graph[uv].capacity;
                flow_graph[uv_rev].capacity = flow_graph[uv].capacity;
                flow_graph[uv_rev].residual_capacity = flow_graph[uv].capacity;
                flow_graph[uv].reverse_edge = uv_rev;
                flow_graph[uv_rev].reverse_edge = uv;
                degree[u] += assign;
                degree[v] += assign;
                total_assign += assign;
                assert(flow_graph[uv].capacity > -ERR_TOLERANCE);
            }
            const int source = n;
            const int sink = n+1;
            for(Vertex v = 0; v < n;v++){
                FlowEdge sv = boost::add_edge(source,v,flow_graph).first;
                FlowEdge sv_rev = boost::add_edge(v,source,flow_graph).first;
                flow_graph[sv].capacity = degree[v]/2;
                flow_graph[sv].residual_capacity = flow_graph[sv].capacity;
                flow_graph[sv_rev].capacity = 0;
                flow_graph[sv_rev].residual_capacity = 0;
                flow_graph[sv].reverse_edge = sv_rev;
                flow_graph[sv_rev].reverse_edge = sv;

                FlowEdge vt = boost::add_edge(v,sink,flow_graph).first;
                FlowEdge vt_rev = boost::add_edge(sink,v,flow_graph).first;
                flow_graph[vt].capacity = 1;
                flow_graph[vt].residual_capacity = 1;
                flow_graph[vt_rev].capacity = 0;
                flow_graph[vt_rev].residual_capacity = 0;
                flow_graph[vt].reverse_edge = vt_rev;
                flow_graph[vt_rev].reverse_edge = vt;
            }

        }
        void run(Vertex v){
            clear();
            FlowEdge sv = boost::edge(source,v,flow_graph).first; 
            previous_degree = flow_graph[sv].capacity;
            previous_vertex = v;
            const int M = boost::num_edges(flow_graph)+n+1;
            flow_graph[sv].capacity = M;
            flow_graph[sv].residual_capacity = M;
            std::vector<FlowEdge> pred_storage(boost::num_vertices(flow_graph));
            std::vector<boost::default_color_type> color_storage(boost::num_vertices(flow_graph));
            std::vector<long> distance_storage(boost::num_vertices(flow_graph));
            auto index_map = boost::get(boost::vertex_index, flow_graph);
            auto predecessor_map = boost::make_iterator_property_map(pred_storage.begin(), index_map);
            auto color_map = boost::make_iterator_property_map(color_storage.begin(), index_map);
            auto distance_map = boost::make_iterator_property_map(distance_storage.begin(), index_map);
            
            flow = boost::boykov_kolmogorov_max_flow(flow_graph,source,sink,
                capacity_map(boost::get(&FlowProperties::capacity,flow_graph)).
                residual_capacity_map(boost::get(&FlowProperties::residual_capacity,flow_graph)).
                reverse_edge_map(boost::get(&FlowProperties::reverse_edge,flow_graph)).
                predecessor_map(predecessor_map).
                color_map(color_map).
                distance_map(distance_map)
            );

            //BFS on the residual graph to determine the source component
            std::queue<Vertex> q;
            std::vector<bool> visited(n+2,false);
            std::vector<Vertex> parent(n+2,n+2);
            std::vector<Vertex> s_component;
            q.push(source);
            visited[source] = true;
            while(!q.empty()){
                Vertex a = q.front();
                q.pop();
                for(auto [eit,end] = boost::out_edges(a,flow_graph);eit != end;eit++){
                    FlowEdge e = *eit;
                    Vertex b = boost::target(e,flow_graph);
                    if(visited[b] || flow_graph[e].residual_capacity < ERR_TOLERANCE) continue;
                    visited[b] = true;
                    s_component.push_back(b);
                    parent[b] = a;
                    q.push(b);
                }
            }
            difference = flow-total_assign; //Difference |S|-x(E(s))
            source_component = s_component;
        }
        std::vector<Vertex> get_source_component(){
            return source_component;
        }
        double get_difference(){
            return difference;
        }
};
