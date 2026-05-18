#pragma once
#include "Models.hpp"
#include "VariableManager.hpp"
#include <random>
#include "IO.hpp"
struct CutInformation{
    GRBLinExpr lhs;
    GRBLinExpr rhs;
    double violation;
    double efficiacy;
    int num_tight;
};

class BaseCallback:public GRBCallback{
    protected:
        std::vector<bool> cutsActive;
        VariableManager& modelVars;
        const Graph& G;
        int n,m,k;
        const double tolerance = 0.001;
        ModelType m_type;
        std::vector<std::vector<CutInformation>> added_cuts;
        double getValueInt(int i, int j, Vars var = X){
            return getSolution(modelVars.getVar(i,j,var));
        }
        double getValueFrac(int i, int j, Vars var = X){
            return getNodeRel(modelVars.getVar(i,j,var));
        }
        friend int clique_heuristic(BaseCallback& cb, Vars type);
        friend int general_clique_heuristic_locsrch(BaseCallback& cb, Vars type, int limit, bool filter,bool locsrch);
        int iter_count = 0;
    public:
        BaseCallback(VariableManager& modelVars, const Graph& G,int ncuts,ModelType type):
        cutsActive(ncuts,false),
        modelVars(modelVars),
        G(G),n(G.numVertices()),m(G.numEdges()),k(modelVars.getK()),m_type(type),
        added_cuts(ncuts,std::vector<CutInformation>()){}
        std::vector<CutInformation> get_added_cuts(Cuts c){
            return added_cuts[c];
        }
};
/*
Cuts: 
Clique: x(S) >= 1 |S| = k+1
General Clique: x(S) >= (q+1)qr/2+q*(q-1)*(k-r)/2 |S| = kq+r
2-Partition: x(S:T)-x(S)-x(T) <= |S|, for |S=1|: x(S:T)-x(T) <= 1 
Wheel: x(c:W) - x(Wheel(W)) <= floor(|W|/2)
*/
/*Greedy separation heuristic for the general clique inequality (see computational study of valid inequalities for the maximum k-cut problem, rodrigues de sousa et al)
greedily adds vertices and adds  best cuts from k+1,...,limit if inequality violated -> at most 1 inequalty*/
int clique_heuristic(BaseCallback& cb, Vars type){
    int count = 0;
    const double cutoff = 1;
    auto get_value = [&cb,&type](Vertex u, Vertex v){return cb.getValueFrac(u,v,type);};
    
    std::unordered_set<Vertex> k1Set;
    double currentMin = cutoff;
    std::pair<Vertex,Vertex> bestEdge;
    for(Vertex u = 0; u < cb.n;u++){
        for(Vertex v = 0; v < u ;v++){
            if(get_value(u,v) < currentMin){
                currentMin = get_value(u,v);
                bestEdge = {u,v};
            }
            if(currentMin == 0)break;
        }
        if(currentMin == 0)break;
    }
    if(currentMin >= cutoff) return count;
    k1Set.insert(bestEdge.first);
    k1Set.insert(bestEdge.second);
    double totalVal = currentMin;
    GRBLinExpr sum = cb.modelVars.getVar(bestEdge.first,bestEdge.second,type);
    while(k1Set.size() < cb.k+1){
        double minVal = 1;  
        Vertex bestVertex;
        GRBLinExpr best_partial = 0;
        for(Vertex u = 0; u < cb.n;u++){
            double currentVal = 0;
            if(k1Set.find(u) != k1Set.end())continue;
            GRBLinExpr partial_sum = 0;
            for(Vertex v:k1Set){
                currentVal += get_value(u,v);
                partial_sum += cb.modelVars.getVar(u,v,type);
            }
            if(currentVal < minVal){
                minVal = currentVal;
                bestVertex = u;
                best_partial = partial_sum;
            }
        }
        totalVal += minVal;
        k1Set.insert(bestVertex);
        sum += best_partial;
        if(totalVal >= cutoff) return count;
    }
    if(totalVal < 1 && (1-totalVal)/std::sqrt(cb.k+1) > 0.05*0){
        cb.addCut(sum >= 1);
        cb.added_cuts[Clique].push_back({sum,GRBLinExpr(1),1-totalVal,(1-totalVal)/std::sqrt(sum.size()),0});
        count++;
    }
    return count;
}

int general_clique_heuristic_locsrch(BaseCallback& cb, Vars type, int limit,bool filter, bool locsrch){
    auto get_value = [&cb,&type](Vertex u, Vertex v){return cb.getValueFrac(u,v,type);};
    int count = 0;
    std::vector<CutInformation> candidate_cuts;
    const int p_limit = std::min(limit,cb.n-1);
    const int q_limit = static_cast<int>(std::floor(p_limit/cb.k));
    const int r_limit = p_limit % cb.k;
    const double cutoff = (q_limit+1)*q_limit*r_limit/2+q_limit*(q_limit-1)*(cb.k-r_limit)/2;
    std::unordered_set<Vertex> k1Set;
    double currentMin = cutoff;
    const int n = cb.n;
    const int k = cb.k;
    std::vector<double> clique_cut(n,0);
    //Find the edge with lowest assignment
    std::pair<Vertex,Vertex> bestEdge;
    for(Vertex u = 0; u < cb.n;u++){
        for(Vertex v = 0; v < u ;v++){
            if(get_value(u,v) < currentMin){
                currentMin = get_value(u,v);
                bestEdge = {u,v};
            }
            if(currentMin == 0)break;
        }
        if(currentMin == 0)break;
    }
    k1Set.insert(bestEdge.first);
    k1Set.insert(bestEdge.second);
    double totalVal = currentMin;
    for(Vertex u = 0; u < n;u++){
        if(u != bestEdge.first)
        clique_cut[u] += get_value(u,bestEdge.first);
        if(u != bestEdge.second)
        clique_cut[u] += get_value(u,bestEdge.second);
    }
    GRBLinExpr sum = cb.modelVars.getVar(bestEdge.first,bestEdge.second,type);
    while(k1Set.size() < p_limit){
        //search for vertex with best cut
        Vertex bestVertex;
        double best_cut = n;
        for(Vertex u = 0; u < cb.n;u++){
            if(k1Set.find(u) != k1Set.end()) continue;
            if(clique_cut[u] < best_cut){
                best_cut = clique_cut[u];
                bestVertex = u;
            }
        }
        //update the clique term
        for(Vertex u:k1Set){
            sum += cb.modelVars.getVar(bestVertex,u,type);

        }
        //update the clique set and cut values
        totalVal += best_cut;
        k1Set.insert(bestVertex);
        for(Vertex u = 0; u < n;u++){
            if(u == bestVertex) continue;
            clique_cut[u] += get_value(bestVertex,u);
        }
        //Local search swap
        if(locsrch){
            std::pair<Vertex,Vertex> swap;
            double best_delta = 1E-4;
            for(Vertex u = 0; u < n; u++){
                if(k1Set.find(u) != k1Set.end()) continue;
                for(Vertex v:k1Set){
                    const double delta = clique_cut[u]-clique_cut[v] - get_value(u,v);
                    if(delta < best_delta){
                        best_delta = delta;
                        swap = {v,u};
                    }
                }
            }
            
            if(best_delta <= 0){
                k1Set.erase(swap.first);
                totalVal += best_delta;
                for(Vertex u = 0; u < n; u++){
                    if(u != swap.second)
                    clique_cut[u] += get_value(u,swap.second);
                    if(u != swap.first)
                    clique_cut[u] -= get_value(u,swap.first);
                }
                for(Vertex u:k1Set){
                    sum += cb.modelVars.getVar(u,swap.second,type) - cb.modelVars.getVar(u,swap.first,type);
                }
                k1Set.insert(swap.second);
            }
        }
        if(k1Set.size() < k+1) continue;
        //Check if the cut is violated
        const int q = static_cast<int>(std::floor(k1Set.size()/k));
        const int r = k1Set.size() % cb.k;
        const int bound = (q+1)*q*r/2+q*(q-1)*(k-r)/2;
        //We want the violation to be atleast by 5 %
        if(totalVal < 0.95 * bound){
            candidate_cuts.push_back({sum,
                GRBLinExpr(bound),
                bound-totalVal,
                (bound-totalVal)/std::sqrt(sum.size()),
                0
            });
        }
    }

    if(filter){
        std::sort(candidate_cuts.begin(),candidate_cuts.end(),[](const CutInformation& a, const CutInformation& b){
            return a.efficiacy > b.efficiacy;
        });
        for(int i = 0; i < std::ceil(0.3*candidate_cuts.size());i++){
            cb.addCut(candidate_cuts[i].lhs >= candidate_cuts[i].rhs);
            cb.added_cuts[GeneralClique].push_back(candidate_cuts[i]);
            count++;
        }
    }
    else{
        for(const auto& cut:candidate_cuts){
            cb.addCut(cut.lhs >= cut.rhs);
            cb.added_cuts[GeneralClique].push_back(cut);
            count++;
        }
    }
    return count;
}


class CallbackER:public BaseCallback{
    private:
        void clique();
        void general_clique();
        void connectivityCut();
        void connectivityCut2();
        void connectivityCutFractional(bool filter = true);
        void cycle_cut();

        double getValueInt(int i,Vars var){
            return getSolution(modelVars.getVar(i,var));
        }
        double getValueInt(int i, int j, Vars var){
            if(var == X && i < j)return getSolution(modelVars.getVar(j,i,var));
            return getSolution(modelVars.getVar(i,j,var));
        }
        double getValueFrac(int i, int j, Vars var){
            if(var == X && i < j)return getNodeRel(modelVars.getVar(j,i,var));
            return getNodeRel(modelVars.getVar(i,j,var));
        }
        double getValueFrac(int i, Vars var){
            return getNodeRel(modelVars.getVar(i,var));
        }
        std::vector<int> cut_count;

        double get_expr_value(GRBLinExpr expr){
            double value = 0;
            for(int i = 0; i < expr.size();i++){
                value += expr.getCoeff(i)*getNodeRel(expr.getVar(i));
            }
            value += expr.getConstant();
            return value;
        }
        bool pre_filter;
        bool locsrch;
        int a,b;
        bool nested_cuts;
    public:
        CallbackER(VariableManager& modelVars, const Graph& G,bool pre_filter = true, bool locsrch = true,int a = 3, int b = 3):
        BaseCallback(modelVars,G,NUMBER_CUTS,EdgeRepresentative),cut_count(NUMBER_CUTS,0),
        pre_filter(pre_filter),locsrch(locsrch),a(a),b(b),nested_cuts(true){}
        void activateCut(Cuts cut){
            cutsActive[cut] = true;
        }
        int get_cut_count(Cuts cut){
            return cut_count[cut];
        }
        void set_filter(bool filter){
            pre_filter = filter;
        }
        void set_locsrch(bool ls){
            locsrch = ls;
        }
        void set_nested(bool nested){
            nested_cuts = nested;
        }
        void set_clqlimit(const std::string& clq_limit){
            auto [as,bs] = IO::seperate_comma(clq_limit);
            a = std::stoi(as);
            b = std::stoi(bs);
        }
    protected:
        void callback(){
            if(where == GRB_CB_MIPNODE && getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL)iter_count++;
            if(cutsActive[ConnectivityFractional])connectivityCutFractional(pre_filter);
            if(cutsActive[GeneralClique])general_clique();
            if(cutsActive[Clique])clique();
            if(cutsActive[Cycles])cycle_cut();
        }

};

void CallbackER::clique(){
    try {
        if(where == GRB_CB_MIPNODE && getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL){
            int ccount = clique_heuristic(*this,X);
            cut_count[Clique] += ccount;

        }
    }   
    catch (GRBException e) {
        std::cout << "Error number: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } 
    catch (...) {
        std::cout << "Error during callback" << std::endl;
    }
}





//Limiting Size of Clique Set to 2k-1, as the inequalities describe facets for k+1,...,2k-1
void CallbackER::general_clique(){
    try {
        if(where == GRB_CB_MIPNODE && getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL){
            int ccount = general_clique_heuristic_locsrch(*this,X,a*k+b,pre_filter,locsrch);
            cut_count[GeneralClique] += ccount;
        }
    }   
    catch (GRBException e) {
        std::cout << "Error number: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } 
    catch (...) {
        std::cout << "Error during callback" << std::endl;
    }
}


void CallbackER::cycle_cut(){
    if(where == GRB_CB_MIPNODE && getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL){
        std::vector<std::tuple<Vertex,Vertex,double>> edge_list;
        for(auto [eit,end] = G.getEdges(); eit !=end; eit++){
            Edge e = *eit;
            auto [a,b] = G.get_incident(e);
            //We are interested in the total assigned edges in a subgraph, so we can consider the undirected graph with the sum of arcs
            edge_list.push_back({a,b,getValueFrac(a,b,Y)+getValueFrac(b,a,Y)}); 
        }
        DenseSubgraph dns_sbgrph(n,edge_list);
        std::vector<CutInformation> candidates;
        for(Vertex v = 0; v < n; v++){
            dns_sbgrph.run(v);
            std::vector<Vertex> s_comp = dns_sbgrph.get_source_component();
            if(dns_sbgrph.get_difference() < 1*0.95 && s_comp.size() > 1){
                GRBLinExpr edge_sum = 0;
                GRBLinExpr pairing_sum = 0;
                double edge_value = 0;
                double pairing_value = 0;
                GRBLinExpr core = 0;
                double core_val = 0;
                std::vector<bool> in_core(s_comp.size());
                int core_size = 0;
                for(Vertex v:s_comp){
                    for(auto [wit,end] = G.getNeighbours(v);wit != end;wit++){
                        Vertex w = *wit;
                        if(std::find(s_comp.begin(),s_comp.end(),w) == s_comp.end())continue;
                        edge_sum += modelVars.getVar(v,w,Y);
                        edge_value += getValueFrac(v,w,Y);
                        if(v < w){
                            pairing_sum += modelVars.getVar(v,w,X);
                            pairing_value += getValueFrac(v,w,X);
                            if(2*(getValueFrac(v,w,Y)+getValueFrac(w,v,Y)) >= getValueFrac(v,w,X)){
                                core += 2*(modelVars.getVar(v,w,Y)+modelVars.getVar(w,v,Y)) - modelVars.getVar(v,w,X);
                                core_val += 2*(getValueFrac(v,w,Y)+getValueFrac(w,v,Y)) - getValueFrac(v,w,X);
                                core_size += (1-in_core[v])+(1-in_core[w]);
                                in_core[v] = true;
                                in_core[w] = true;
                            }
                        }
                    }
                }
                double violation = 1-dns_sbgrph.get_difference();
                //If the subset forms a cycle, the new cuts are strictly stronger than classical cycle elimination
                if(true && pairing_sum.size() > s_comp.size()){
                    candidates.push_back({edge_sum,GRBLinExpr(s_comp.size()-1),violation,violation/std::sqrt(edge_sum.size()),0});
                }
                /*
                violation = 2 - dns_sbgrph.get_difference();
                if(2*edge_value - pairing_value > 1.05*(s_comp.size()-2)){
                    candidates.push_back({2*edge_sum-pairing_sum,GRBLinExpr(s_comp.size()-2),violation,violation/std::sqrt(edge_sum.size()),0});
                }*/
                if(s_comp.size() < 2)continue;
                violation = 2 + core_val-core_size;
                if(core_val > core_size-2*1.05){
                    candidates.push_back({core,GRBLinExpr(core_size-2),violation,violation/std::sqrt(core.size()),0});
                }
            }

        }
        if(!pre_filter){
                for(auto cut:candidates){
                    addCut(cut.lhs <= cut.rhs);
                    added_cuts[Cycles].push_back(cut);
                    cut_count[Cycles]+=1;
                }
                return;
            }
        std::sort(candidates.begin(),candidates.end(),[](const CutInformation& a, const CutInformation& b){
                return a.efficiacy > b.efficiacy;}
            );
            for(int i = 0; i< std::ceil(0.3*candidates.size());i++){
                addCut(candidates[i].lhs <= candidates[i].rhs);
                added_cuts[Cycles].push_back(candidates[i]);
                cut_count[Cycles]+=1;
                
            }
    }
}
