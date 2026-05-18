#pragma once
#include "gurobi_c++.h"
#include "graph.hpp"
#include <memory>
#include <cmath>
#include "commons.hpp"
#include <unordered_map>
#include <string>
#include "Models.hpp"
namespace filesystem = std::filesystem;
class CycleBreaking : public GRBCallback{
    private:
        GRBVar *x;
        GRBVar *t;
        const Graph &G;
        const int n;
        const int k;
        const std::unordered_map<int,int>& edgeMap;
        double getValue(int i,int j){
            return getSolution(x[edgeMap.at(i+j*n)]);
        }
        double getValueFractional(int i,int j){
            return getNodeRel(x[edgeMap.at(i+j*n)]);
        }
        GRBVar& getVar(int i,int j){
            return x[edgeMap.at(i+j*n)];
        }
        double getValueT(int i,int j){
            return getNodeRel(getVarT(i,j));
        }
        double getValueTInt(int i,int j){
            return getSolution(getVarT(i,j));
        }
        GRBVar& getVarT(int i,int j){
            if(i < j) return t[j*(j-1)/2+i];
            return t[i*(i-1)/2+j];
        }
        bool fractional_cut;
        bool strengthened;
    public:
        CycleBreaking(GRBVar *x, GRBVar *t, const Graph &G, int k, const std::unordered_map<int,int>& eMap, bool strngth = false) : 
        x(x),t(t), G(G),n(G.numVertices()),k(k),edgeMap(eMap),fractional_cut(false),strengthened(strngth){}
        void set_fractional(bool frac){
            fractional_cut = frac;
        }
    protected:
        void callback(){
            try{
                if (where == GRB_CB_MIPSOL){
                    const double eps = 0.001;
                    const int n = G.numVertices();
                    std::vector<std::pair<int, int>> activeEdges;
                    for (int u = 0; u < n; u++){
                        for (auto [vit, vend] = G.getNeighbours(u); vit != vend; vit++){
                            Vertex v = *vit;
                            if (getValue(u,v) >= 1 - eps){
                                activeEdges.push_back({u, v});
                            }
                        }
                    }
                    StronglyConnectedComponents scc(activeEdges, n);
                    scc.run();
                    std::vector<int> connectedComponents = scc.getComponents();
                    int num = scc.getNumComponents();
                    std::vector<int> componentSizes = scc.getComponentSizes();
                    //std::vector<int> areaSizes = scc.getComponentSizes();
                    std::vector<GRBLinExpr> componentSums(num, 0);
                    /*
                    //Compute component sets
                    std::vector<std::vector<Vertex>> component_sets(num,std::vector<Vertex>());
                    for(int v = 0; v < n;v++){
                        component_sets[connectedComponents[v]].push_back(v);
                    }
                    for(int c = 0; c < num;c++){
                        if(componentSizes[c] < 2) continue;
                        Vertex r = component_sets[c][0];
                        std::set<Vertex> area_set;
                        for(auto [eit,end] = G.getEdges(); eit !=end; eit++){
                            Edge e = *eit;
                            auto [a,b] = G.get_incident(e);
                            if((a == r || getValueTInt(a,r) > 0.5) && (b == r || getValueTInt(b,r) > 0.5)){
                                //componentSums[c] += getVar(a,b)+getVar(b,a);
                                area_set.insert(a);
                                area_set.insert(b);
                            }
                        }
                        areaSizes[c] = area_set.size();
                    }*/
                    std::vector<GRBLinExpr> coreSums(num, 0);
                    for (int u = 0; u < n; u++){
                        for (auto [vit, vend] = G.getNeighbours(u); vit != vend; vit++){
                            Vertex v = *vit;
                            if (connectedComponents[u] == connectedComponents[v]){
                                componentSums[connectedComponents[u]] += getVar(u,v);
                                if(u < v && 2*(getValue(u,v)+getValue(v,u)) >= getValueTInt(u,v)) 
                                coreSums[connectedComponents[u]] += 2*(getVar(u,v)+getVar(v,u))-getVarT(u,v);
                            }
                        }
                    }
                    for(int c = 0; c < num;c++){
                        if(componentSizes[c] < 2)continue;
                        if(! strengthened || componentSums[c].size() > componentSizes[c]){
                            //addLazy(componentSums[c] <= areaSizes[c]-1);  
                            addLazy(componentSums[c] <= componentSizes[c]-1);  
                        }
                        if(strengthened){
                            addLazy(coreSums[c] <= componentSizes[c]-2);  
                        }
                    }
                }
                else if(fractional_cut && where == GRB_CB_MIPNODE && getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL){
                    std::vector<std::tuple<Vertex,Vertex,double>> edge_list;
                    for(auto [eit,end] = G.getEdges(); eit !=end; eit++){
                        Edge e = *eit;
                        auto [a,b] = G.get_incident(e);
                        //We are interested in the total assigned edges in a subgraph, so we can consider the undirected graph with the sum of arcs
                        edge_list.push_back({a,b,getValueFractional(a,b)+getValueFractional(b,a)}); 
                    }
                    DenseSubgraph dns_sbgrph(n,edge_list);
                    for (Vertex v = 0; v < n; v++){
                        dns_sbgrph.run(v);
                        if (dns_sbgrph.get_difference() < 1*0.95){
                            std::vector<Vertex> s_comp = dns_sbgrph.get_source_component();
                            GRBLinExpr edge_sum = 0;
                            double edge_val = 0;
                            GRBLinExpr core_sum = 0;
                            double core_val = 0;
                            for (Vertex v : s_comp){
                                for (auto [wit, end] = G.getNeighbours(v); wit != end; wit++){
                                    Vertex w = *wit;
                                    // if(w < v) continue;
                                    if (std::find(s_comp.begin(), s_comp.end(), w) == s_comp.end())continue;
                                    edge_sum += getVar(v, w);
                                    edge_val += getValueFractional(v,w);
                                    if (v < w && 2*(getValue(v,w)+getValue(w,v)) >= getValueT(v,w)){
                                        core_sum += getVarT(v, w);
                                        core_val += getValueT(v,w);
                                    }
                                }
                            }
                            if(edge_sum.size() > s_comp.size()){
                                addCut(edge_sum - s_comp.size() <=  s_comp.size()-1);
                            }
                            if(core_val > (s_comp.size() - 2)*1.05){
                                addCut(core_sum <=  s_comp.size()-2);
                            }
                        }
                    }
                }
            }
            catch (GRBException e){
                std::cout << "Error number: " << e.getErrorCode() << std::endl;
                std::cout << e.getMessage() << std::endl;
            }
            catch (...) {
                std::cout << "Error during callback" << std::endl;
            }
        }
    };


class TreeModel_Duque:public BaseModel{
    private:
        GRBVar* t;
        GRBVar* x;
        GRBVar* o;
        CycleBreaking* callback;
        std::unordered_map<int,int> edgeMapDir;
        bool cycleCut;
        bool recursive;
        bool strengthened;
        //Bookkeeping for recursive strategy
        double total_runtime = 0;
        int total_nodes = 0;
        int total_iter = 0;
        GRBVar& getVar(int i,int j,Vars var = X)const{ 
            switch(var){
                case T:
                    return t[i*(i-1)/2+j];
                case X:
                    return x[edgeMapDir.at(idx(i,j))];
                default:
                    throw std::invalid_argument("Wrong Var idx \n");
            }
        }
        GRBVar& getVar(int i)const{
            return o[i];
        }
        int idx(int i, int j)const{
            return i+j*n;
        }
        double dissimilarity(Vertex u, Vertex v){
            return std::abs(G.getAttribute(u).feature-G.getAttribute(v).feature);
        }
    public:
        TreeModel_Duque(GRBEnv* env,const Graph& G, int k,double lambda = 1, bool strengthened = false):
        BaseModel(env,G,k,lambda),callback(nullptr),cycleCut(true),recursive(false),strengthened(strengthened){
            int counter = 0;
            for(int u = 0; u < n;u++){
                for(auto[vit,end] = G.getNeighbours(u);vit != end;vit++){
                    Vertex v = *vit;
                    edgeMapDir[idx(u,v)] = counter++;
                }
            }
            x = model.addVars(m*2,GRB_BINARY);
            t = model.addVars(n*(n-1)/2,GRB_BINARY);
            o = model.addVars(n,GRB_CONTINUOUS);

            const int M = n-k;
            GRBLinExpr edgeSum = 0;
            for(Vertex u = 0; u< n;u++){
                GRBLinExpr neighbourSum = 0;
                for(auto [vit,vend] = G.getNeighbours(u);vit != vend;vit++){
                    Vertex v = *vit;
                    edgeSum += getVar(u,v,X);
                    neighbourSum += getVar(u,v,X);
                    if(strengthened){
                        if(u > v){
                            model.addConstr(getVar(u,v,X)+getVar(v,u,X) <= getVar(u,v,T),"Pair-Edge");
                        }
                    }
                    else{
                        model.addConstr(getVar(u,v,X) <= (u > v ?getVar(u,v,T):getVar(v,u,T)),"Pair-Edge");
                    }
                    if(!cycleCut){
                        model.addConstr(getVar(u)-getVar(v) + M*getVar(u,v,X)+(M-2)*getVar(v,u,X) <= M-1);
                    }
                }
                model.addConstr(neighbourSum <= 1,"Outgoing-Arcs");
                for(Vertex v = 0; v< u;v++){
                    getVar(u,v,T).set(GRB_DoubleAttr_Obj,weights[u][v]);
                    //getVar(u,v).set(GRB_IntAttr_BranchPriority,5);
                    for(Vertex w = 0; w < v; w++){
                            model.addConstr(getVar(u,v,T) >= getVar(u,w,T)+getVar(v,w,T)-1,"Triangle");
                            model.addConstr(getVar(u,w,T) >= getVar(u,v,T)+getVar(v,w,T)-1,"Triangle");
                            model.addConstr(getVar(v,w,T) >= getVar(u,v,T)+getVar(u,w,T)-1,"Triangle");
                    }
                }
                if(! cycleCut){
                    getVar(u).set(GRB_DoubleAttr_LB,1);
                    getVar(u).set(GRB_DoubleAttr_LB,M);
                }
            }
            model.addConstr(edgeSum == M,"Edge-Sum");
            if(cycleCut && !recursive){
                model.set(GRB_IntParam_LazyConstraints,1);
                callback = new CycleBreaking(x,t,G,k,edgeMapDir,strengthened);
                model.setCallback(callback);
            }
        }
        double getVarVal(int i, int j,Vars v = X)const{
            return getVar(i,j,v).get(GRB_DoubleAttr_X);
        }
        void optimize(const std::vector<std::pair<std::string,std::string>>& parameters = {
            std::make_pair("TimeLimit","1000.0"),std::make_pair("Threads","1")
        }) override{
            double timelimit = 1E10;
            for(auto p:parameters){
               model.set(p.first,p.second);
               if(p.first == "TimeLimit") timelimit = std::stod(p.second);
            }
            total_runtime = total_nodes = total_iter = 0;
            double callback_time = 0;
            int callback_count = 0;
            if(cycleCut && recursive){
                filesystem::path log = model.get(GRB_StringParam_LogFile);
                std::ofstream log_file;
                log_file.open(log,std::ios_base::app);
                while(true){
                    model.set(GRB_DoubleParam_TimeLimit, timelimit-total_runtime);
                    model.optimize();
                    total_runtime += model.get(GRB_DoubleAttr_Runtime);
                    total_nodes += model.get(GRB_DoubleAttr_NodeCount);
                    total_iter += model.get(GRB_DoubleAttr_IterCount);
                    const double eps = 0.001;
                    //Time separation of constraints
                    auto t1 = std::chrono::high_resolution_clock::now();
                    std::vector<std::pair<int, int>> activeEdges;
                    for (int u = 0; u < n; u++){
                        for (auto [vit, vend] = G.getNeighbours(u); vit != vend; vit++){
                            Vertex v = *vit;
                            if (getVarVal(u,v,X) >= 1 - eps){
                                activeEdges.push_back({u, v});
                            }
                        }
                    }
                    StronglyConnectedComponents scc(activeEdges, n);
                    scc.run();
                    std::vector<int> connectedComponents = scc.getComponents();
                    int num = scc.getNumComponents();
                    std::vector<int> componentSizes = scc.getComponentSizes();
                    std::vector<GRBLinExpr> componentSums(num, 0);
                    for (int u = 0; u < n; u++){
                        for (auto [vit, vend] = G.getNeighbours(u); vit != vend; vit++){
                            Vertex v = *vit;
                            if (connectedComponents[u] == connectedComponents[v]){
                                componentSums[connectedComponents[u]] += getVar(u,v);
                            }
                        }
                    }
                    int cycle_count = 0;
                    for(int c = 0; c < num;c++){
                        if(componentSizes[c] < 2)continue;
                        model.addConstr(componentSums[c] <= componentSizes[c]-1);  
                        cycle_count ++;
                    }
                    auto t2 = std::chrono::high_resolution_clock::now();
                    const double single_callback = std::chrono::duration_cast<std::chrono::seconds>(t2-t1).count();
                    callback_time += single_callback;
                    total_runtime += single_callback;
                    callback_count ++;
                    if(cycle_count == 0)break;
                    std::cout << "Found " << cycle_count << " cycles, resolving \n";
                    log_file << "Found " << cycle_count << " cycles, resolving" << std::endl;
                }
                std::cout << std::fixed << std::setprecision(2) << 
                "User-callback calls "<< callback_count <<", time in user-callback " << std::setprecision(2)<< callback_time <<" sec\n";
                log_file << std::fixed << std::setprecision(2) << 
                "User-callback calls "<< callback_count <<", time in user-callback " << std::setprecision(2)<< callback_time <<" sec" << std::endl;
                log_file.close();
            }
            else{
                model.optimize();
            }
        }

        //Output solution as assigment to regions
        std::vector<int> return_solution()const override{
            std::vector<int> solution(n,-1);
            std::vector<int> represented_by(n,-1);
            int counter = 0;
            for(Vertex u = 0;u< n;u++){
                for(int v = u-1;v >= 0;v--){
                    if(getVarVal(u,v,T) == 1){
                        represented_by[u] = v;
                    }
                
                }
                if(represented_by[u] == -1){
                    solution[u] = counter++;
                }
            }
            for(Vertex u = 0; u< n;u++){
                if(represented_by[u] == -1) continue;
                solution[u] = solution[represented_by[u]];
            }
            return solution;
        }
        void set_fractional(bool frac){
            callback->set_fractional(frac);
        }
        /*void optimize(std::vector<GRBparameter> parameters = {
            std::make_pair(GRB_DoubleParam_TimeLimit,1000.0),
            std::make_pair(GRB_IntParam_Threads,1)
        }){
            for(auto p:parameters){
                if(std::holds_alternative<std::pair<GRB_IntParam,int>>(p)){
                    auto [param,value] = std::get<std::pair<GRB_IntParam,int>>(p);
                    model.set(param,value);
                }
                if(std::holds_alternative<std::pair<GRB_DoubleParam,double>>(p)){
                    auto [param,value] = std::get<std::pair<GRB_DoubleParam,double>>(p);
                    model.set(param,value);
                }
                if(std::holds_alternative<std::pair<GRB_StringParam,std::string>>(p)){
                    auto [param,value] = std::get<std::pair<GRB_StringParam,std::string>>(p);
                    model.set(param,value);
                }
            }
            model.optimize();
        }*/
        GRBattribute get_attributes(GRBattributeType attr){
            if(recursive && std::holds_alternative<GRB_DoubleAttr>(attr)){
                //Special case for attributes that are aggregated over all resolvings
                GRB_DoubleAttr a = std::get<GRB_DoubleAttr>(attr);
                if(a == GRB_DoubleAttr_Runtime) return total_runtime;
                if(a == GRB_DoubleAttr_NodeCount) return static_cast<double>(total_nodes);
                if(a == GRB_DoubleAttr_IterCount) return static_cast<double>(total_iter);
            }
            if(std::holds_alternative<GRB_IntAttr>(attr)){
                GRB_IntAttr a = std::get<GRB_IntAttr>(attr);
                return model.get(a);
            }
            if(std::holds_alternative<GRB_DoubleAttr>(attr)){
                GRB_DoubleAttr a = std::get<GRB_DoubleAttr>(attr);
                return model.get(a);
            }
            if(std::holds_alternative<GRB_StringAttr>(attr)){
                GRB_StringAttr a = std::get<GRB_StringAttr>(attr);
                return model.get(a);
            }
            return 0;
        }
        bool check_validity(){
            for(Vertex u = 0; u < n ;u++){
                for(Vertex v = 0; v < u; v++){
                    for(Vertex w = 0; w < v; w++){
                        if(getVarVal(u,v,T) < getVarVal(u,w,T)+getVarVal(v,w,T)-1-ERR_TOLERANCE)return false;
                        if(getVarVal(u,w,T) < getVarVal(u,v,T)+getVarVal(v,w,T)-1-ERR_TOLERANCE)return false;
                        if(getVarVal(v,w,T) < getVarVal(u,v,T)+getVarVal(u,w,T)-1-ERR_TOLERANCE)return false;
                        //TODO: missing checks for uv=0,uw=1 => vw = 0 etc
                    }

                } 
            }
            
            return true;
        }
        void setWarmStart(int u, int v, int value,Vars var){
            getVar(u,v,var).set(GRB_DoubleAttr_Start,value);
        }
        void set_init_solution([[maybe_unused]] const std::vector<std::vector<Vertex>>& init)override{
            int num_edges = 0;
            for(auto p:init){
                std::vector<bool> visited(n,true);
                for(auto u:p){
                    visited[u] = false;
                    for(auto v:p)
                    if(u > v) setWarmStart(u,v,1,T);
                }
                std::queue<Vertex> bfs;
                bfs.push(p[0]);
                visited[p[0]] = true;
                while(!bfs.empty()){
                    Vertex v = bfs.front();
                    bfs.pop();
                    for(auto [wit,end] = G.getNeighbours(v);wit != end;wit++){
                        Vertex w = *wit;
                        if(visited[w]) continue;
                        visited[w] = true;
                        //model.addConstr(getVar(v,w,X) == 1);
                        setWarmStart(w,v,1,X);
                        num_edges++;
                        bfs.push(w);
                    }
                }
            }
            assert(num_edges == n-k);


        }
        ~TreeModel_Duque(){
            delete[] x;
            delete[] t;
            delete[] o;
            delete callback;
        }

};



class OrderModel_Duque:public BaseModel{
    private:
        GRBVar* t;
        GRBVar* x;
        const int q;
        std::unordered_map<int,int> edgeMapDir;
        GRBVar& getVar(int i,int j)const{
                return t[i*(i-1)/2+j];
        }
        GRBVar& getVar(int u,int i,int o)const{
            return x[u*(q+1)*k+i*(q+1)+o];
        }
        double dissimilarity(Vertex u, Vertex v){
            return std::abs(G.getAttribute(u).feature-G.getAttribute(v).feature);
        }
        public:
            OrderModel_Duque(GRBEnv* env,const Graph& G, int k, double lambda = 1):BaseModel(env,G,k,lambda),q(G.numVertices()-k+1){
                t = model.addVars(n*(n-1)/2,GRB_BINARY);
                x = model.addVars(n*k*(q+1),GRB_BINARY);

                //Seed:
                model.addConstr(getVar(0,0,0) == 1);

                for(int i = 0; i < k; i++){
                    GRBLinExpr sum = 0;
                    for(int u = 0; u< n;u++){
                        sum += getVar(u,i,0);
                    }
                    model.addConstr(sum == 1);
                }

                for(int u = 0; u< n;u++){
                    GRBLinExpr sum = 0;
                    for(int i = 0; i < k; i++){
                        for(int o = 0; o <= q; o++){
                            sum += getVar(u,i,o);
                            if(o == 0) continue;
                            GRBLinExpr neighbourSum = 0;
                            for(auto [vit,vend] = G.getNeighbours(u);vit != vend;vit++){
                                Vertex v = *vit;
                                neighbourSum += getVar(v,i,o-1);
                            }
                            model.addConstr(getVar(u,i,o) <= neighbourSum);
                        }
                    }
                    model.addConstr(sum == 1);
                    for(int v = 0; v< u;v++){
                        for(int i = 0; i < k; i++){
                            GRBLinExpr orderSum = 0;
                            for(int o = 0; o <= q; o++){
                                orderSum += getVar(u,i,o)+ getVar(v,i,o);
                            }
                            model.addConstr(getVar(u,v) >= orderSum - 1);
                        }
                        getVar(u,v).set(GRB_DoubleAttr_Obj,weights[u][v]);
                    }
                }
            }
            double getVarVal(int u,int i,int o)const{
                return getVar(u,i,o).get(GRB_DoubleAttr_X);
            }
            double getVarVal(int u,int i)const{
                return getVar(u,i).get(GRB_DoubleAttr_X);
            }
            std::vector<int> return_solution()const override{
                std::vector<int> solution(n,-1);

                for(Vertex v = 0;v< n;v++){
                    for(int i = 0;i< k;i++){
                        for(int o = 0;o< q;o++){
                            if(getVarVal(v,i,o) > 0.5){
                                solution[v] = i;
                                break;
                            }
                        }
                        if(solution[v] > -1)break;
                    }
                }
                return solution;
                
            }
            /*
            void optimize(std::vector<GRBparameter> parameters = {
                std::make_pair(GRB_DoubleParam_TimeLimit,1000.0),
                std::make_pair(GRB_IntParam_Threads,1)
            }){
                for(auto p:parameters){
                    if(std::holds_alternative<std::pair<GRB_IntParam,int>>(p)){
                        auto [param,value] = std::get<std::pair<GRB_IntParam,int>>(p);
                        model.set(param,value);
                    }
                    if(std::holds_alternative<std::pair<GRB_DoubleParam,double>>(p)){
                        auto [param,value] = std::get<std::pair<GRB_DoubleParam,double>>(p);
                        model.set(param,value);
                    }
                    if(std::holds_alternative<std::pair<GRB_StringParam,std::string>>(p)){
                        auto [param,value] = std::get<std::pair<GRB_StringParam,std::string>>(p);
                        model.set(param,value);
                    }
                }
                model.optimize();
            }*/
            GRBattribute get_attributes(GRBattributeType attr){
                if(std::holds_alternative<GRB_IntAttr>(attr)){
                    GRB_IntAttr a = std::get<GRB_IntAttr>(attr);
                    return model.get(a);
                }
                if(std::holds_alternative<GRB_DoubleAttr>(attr)){
                    GRB_DoubleAttr a = std::get<GRB_DoubleAttr>(attr);
                    return model.get(a);
                }
                if(std::holds_alternative<GRB_StringAttr>(attr)){
                    GRB_StringAttr a = std::get<GRB_StringAttr>(attr);
                    return model.get(a);
                }
                return 0;
            }

            bool check_validity(){
                for(Vertex u = 0; u < n ;u++){
                    double sum = 0;
                    for(Vertex i = 0; i < k; i++){
                        for(int o = 0; o < q; o++){
                            sum += getVarVal(u,i,o);
                        }
                    } 
                    if(!DOUBLE_EQUALITY(sum,1))return false;
                }
                /*for(Vertex i = 0; i < k; i++){
                    std::vector<int> orders;
                    for(int o = 0; o < q; o++){
                        for(Vertex u = 0; u < n ;u++){
                            if(DOUBLE_EQUALITY(getVarVal(u,i,o),1))orders.push_back(o);
                        }
                    }
                    std::sort(orders.begin(),orders.end());
                    for(int i = 0; i< orders.size();i++){
                        std::cout << orders[i] <<",";
                        if(orders[i] != i){std::cout <<"aaaaaaaaaaaaaaaaaaaaaaa\n";return false;}
                    }
                }*/
                return true;
            }
            
            ~OrderModel_Duque(){
                delete[] t;
                delete[] x;
            }
     
};