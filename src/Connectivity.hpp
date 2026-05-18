#pragma once
#include "Models.hpp"
#include <queue>
#include "Cuts.hpp"
#include <random>

void FormulationNC::connectivity_scf(bool biDirConst){
    //Add region root variables
    if(!modelVars.variableExists(R)){
        modelVars.addVariableType(model.addVars(n*k,GRB_BINARY),R,UnorderedPairwise);
        for(Vertex i = 0; i < k;i++){
            GRBLinExpr sum = 0;
            for(Vertex v = 0; v < n;v++){
                sum += get_var(v,i,R);
                model.addConstr(get_var(v,i,R) <= get_var(v,i,X));
            }
            model.addConstr(sum == 1);
        }
    }

    if(!modelVars.variableExists(Y)){
        modelVars.addVariableType(model.addVars(m,GRB_BINARY),Y,UndirectedEdges);
            for(Vertex u = 0; u < n;u++){
                for(auto[vit,end] = G.getNeighbours(u);vit != end;vit++){
                    Vertex v = *vit;
                    if(u > v) continue;
                    for(int w = 0; w < k;w++){
                        model.addConstr(get_var(u,w,X)-get_var(v,w,X) <= get_var(u,v,Y));
                        if(biDirConst){
                            model.addConstr(get_var(v,w,X)-get_var(u,w,X) <= get_var(u,v,Y));
                        }
                    }
                }
            }
    }

    if(!modelVars.variableExists(Fscf)){
        modelVars.addVariableType(model.addVars(m*2,GRB_CONTINUOUS),Fscf,DirectedEdges);
    }
    int M = n-k+1;
    for(Vertex u = 0; u < n;u++){
        GRBLinExpr insum = 0;
        GRBLinExpr outsum = 0;
        for(auto [vit,vend] = G.getNeighbours(u);vit != vend;vit++){
            Vertex v = *vit;
            outsum += get_var(u,v,Fscf);
            insum += get_var(v,u,Fscf);
            if(u < v){
                model.addConstr(get_var(u,v,Fscf) + get_var(v,u,Fscf) <= M*(1-get_var(u,v,Y)));
            }
        }
        GRBLinExpr rsum = 0;
        for(int i = 0; i < k ;i++){
            rsum += get_var(u,i,R);
        }
        model.addConstr(insum-outsum >= 1 - M*rsum);
    }
    model.update();
}


void FormulationNC::connectivity_shir(){
    /*NOTE: This uses the shir model from the original p-regions paper by Duque. et al.
    Validi et al. present a slightly different model with additional variables that model the 
    amount of generated flow
    */
    //Add region root variables
    if(!modelVars.variableExists(R)){
        modelVars.addVariableType(model.addVars(n*k,GRB_BINARY),R,UnorderedPairwise);
        for(Vertex i = 0; i < k;i++){
            GRBLinExpr sum = 0;
            for(Vertex v = 0; v < n;v++){
                sum += get_var(v,i,R);
                model.addConstr(get_var(v,i,R) <= get_var(v,i,X));
            }
            model.addConstr(sum == 1);
        }
    }

    if(!modelVars.variableExists(Fshir)){
        modelVars.addVariableType(model.addVars(k*m*2,GRB_CONTINUOUS),Fshir,EdgesRegionPair);
    }
    int M = n-k;
    for(int i = 0; i < k ;i++){
        for(Vertex u = 0; u < n;u++){
            GRBLinExpr insum = 0;
            GRBLinExpr outsum = 0;
            for(auto [vit,vend] = G.getNeighbours(u);vit != vend;vit++){
                Vertex v = *vit;
                model.addConstr(get_var(i,u,v,Fshir) <= M*get_var(u,i,X));
                model.addConstr(get_var(i,u,v,Fshir) <= M*get_var(v,i,X));
                outsum += get_var(i,u,v,Fshir);
                insum += get_var(i,v,u,Fshir);    
            }
            model.addConstr(outsum-insum >= get_var(u,i,X)-M*get_var(u,i,R));
        }
    }
    model.update();
}

void CallbackER::connectivityCut(){
    try {
        if (where == GRB_CB_MIPSOL) {
            const int n = modelVars.getN();
            for(Vertex b = 0; b < n; b++){
                std::vector<bool> visited(n,false);
                //BFS loop
                std::queue<Vertex> bfs;
                std::vector<bool> border(n,false);
                //Compute C_b using BFS
                bfs.push(b);
                visited[b] = true;
                while(!bfs.empty()){
                    Vertex v = bfs.front();
                    bfs.pop();
                    for(auto [wit,wend] = G.getNeighbours(v);wit != wend;wit++){
                        Vertex w = *wit;
                        if(w == b)continue;
                        if((b>w ? getValueInt(b,w,X):getValueInt(w,b,X)) >= 1-tolerance){
                            if(!visited[w]){
                                bfs.push(w);
                                visited[w] = true;
                            }
                        }
                        else{
                            border[w] = true;
                        }
                        
                    }
                }
                for(Vertex a = 0; a < n; a++){
                    if(visited[a]||border[a]||(b>a ?getValueInt(b,a,X):getValueInt(a,b,X)) < 1-tolerance)continue;
                    bfs = {};
                    GRBLinExpr sepVertices = 0;
                    std::vector<Vertex> assignedVertices = {a};
                    bfs.push(a);
                    visited[a] = true;
                    while(!bfs.empty()){
                        Vertex v = bfs.front();
                        bfs.pop();
                        for(auto [wit,wend] = G.getNeighbours(v);wit != wend;wit++){
                            Vertex w = *wit;
                            if(border[w]){
                                sepVertices += (b>w ? modelVars.getVar(b,w,X):modelVars.getVar(w,b,X));
                                continue;
                            }
                            if(!visited[w]){
                                if((w>v ? getValueInt(w,v,X):getValueInt(v,w,X))>=1-tolerance){
                                    assignedVertices.push_back(w);
                                }
                                bfs.push(w);
                                visited[w] = true;
                            }
                            
                        }
                    }
                    for(Vertex v:assignedVertices){
                        if(v > b){
                            addLazy(modelVars.getVar(v,b,X) <= sepVertices);
                        }
                        else{
                            addLazy(modelVars.getVar(b,v,X) <= sepVertices);
                        }
                        cut_count[Connectivity]+=1;
                    }
                }
        }
    
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

void CallbackER::connectivityCut2(){
    try {
        if (where == GRB_CB_MIPSOL) {
            const int n = modelVars.getN();
            std::vector<bool> main_component(n,false);
            for(Vertex b = 0; b < n; b++){
                if(main_component[b]) continue;
                std::vector<bool> visited(n,false);
                //BFS loop
                std::queue<Vertex> bfs;
                std::vector<bool> border(n,false);
                //Compute C_b using BFS
                bfs.push(b);
                visited[b] = true;
                main_component[b] = true;
                while(!bfs.empty()){
                    Vertex v = bfs.front();
                    bfs.pop();
                    for(auto [wit,wend] = G.getNeighbours(v);wit != wend;wit++){
                        Vertex w = *wit;
                        if(w == b)continue;
                        if((b>w ? getValueInt(b,w,X):getValueInt(w,b,X)) >= 1-tolerance){
                            if(!visited[w]){
                                bfs.push(w);
                                visited[w] = true;
                                main_component[w] = true;
                            }
                        }
                        else{
                            border[w] = true;
                        }
                        
                    }
                }
                //While considering Connected components disconnected from b, we add a constraint between every vertex in that component and b
                //As the separator, we use the front of the connected component of a 
                for(Vertex a = 0; a < n; a++){
                    if(visited[a]||border[a]||(b>a ?getValueInt(b,a,X):getValueInt(a,b,X)) < 1-tolerance)continue;
                    bfs = {};
                    GRBLinExpr sepVertices = 0;
                    std::vector<Vertex> assignedVertices = {a};
                    bfs.push(a);
                    visited[a] = true;
                    while(!bfs.empty()){
                        Vertex v = bfs.front();
                        bfs.pop();
                        for(auto [wit,wend] = G.getNeighbours(v);wit != wend;wit++){
                            Vertex w = *wit;
                            if(border[w]){
                                sepVertices += (b>w ? modelVars.getVar(b,w,X):modelVars.getVar(w,b,X));
                                continue;
                            }
                            if(!visited[w]){
                                if((w>v ? getValueInt(w,v,X):getValueInt(v,w,X))>=1-tolerance){
                                    assignedVertices.push_back(w);
                                }
                                bfs.push(w);
                                visited[w] = true;
                            }
                            
                        }
                    }
                    
                    for(Vertex v:assignedVertices){
                        if(v > b){
                            addLazy(modelVars.getVar(v,b,X) <= sepVertices);
                        }
                        else{
                            addLazy(modelVars.getVar(b,v,X) <= sepVertices);
                        }
                        cut_count[Connectivity]+=1;
                    }
                    //break;
                }
        }
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

class DebugCallback:public GRBCallback{
    protected:
        void callback(){std::cout << "test\n";}
};

void FormulationER::connectivity_scf(){


    if(!modelVars.variableExists(Fscf)){
        modelVars.addVariableType(model.addVars(m*2,GRB_CONTINUOUS),Fscf,DirectedEdges);
    }
    int M = n-k+1;
    for(Vertex u = 0; u < n;u++){
        GRBLinExpr insum = 0;
        GRBLinExpr outsum = 0;
        for(auto [vit,vend] = G.getNeighbours(u);vit != vend;vit++){
            Vertex v = *vit;
            outsum += get_var(u,v,Fscf);
            insum += get_var(v,u,Fscf);
            if(u < v){
                model.addConstr(get_var(u,v,Fscf) + get_var(v,u,Fscf) <= M*get_var(u,v,X));
            }
        }
        model.addConstr(insum-outsum >= 1 - M*get_var(u,R));
    }
    model.update();
}

void CallbackER::connectivityCutFractional(bool filter){
    const bool creep_flows = false;
    try {
        if((where == GRB_CB_MIPNODE && getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL) || (where == GRB_CB_MIPSOL)){
            if(where == GRB_CB_MIPSOL) filter = false;
            const int n = modelVars.getN();
            std::vector<CutInformation> candidates;
            //We search for violated u,v separators
            std::vector<double> frac_values(n);
            for(Vertex u = 0; u< n; u++){
                for(Vertex v = 0; v < n ;v++){
                    //Collecting assignment values of vertices to u 
                    for(Vertex w = 0; w < n; w++){
                        if(u == w)continue;
                        frac_values[w] = ((where == GRB_CB_MIPNODE) ? getValueFrac(u,w,X):getValueInt(u,w,X))+creep_flows*1E-4; 
                    }
                    if(u == v) continue;
                    double value_uv = (where == GRB_CB_MIPNODE) ?getValueFrac(u,v,X):getValueInt(u,v,X);
                    if(value_uv < 1.0/k || G.adjacent(u,v)) continue;
                    while(true){
                        auto [flow_value,separator] = G.minimum_separator(u,v,frac_values);
                        //We want the violation to be atleast by 5 %
                        if(flow_value >= value_uv-0.05) break;
                        GRBLinExpr sum = 0;
                        for(Vertex z:separator){
                            sum+= modelVars.getVar(u,z,X);
                            if(nested_cuts)frac_values[z] = 1; //Nested cuts
                        }
                        //Either one of the seperator vertices must be assigned to u
                        candidates.push_back({
                            sum-modelVars.getVar(u,v,X),
                            GRBLinExpr(0),
                            value_uv-flow_value,
                            (value_uv-flow_value)/std::sqrt(separator.size()),
                            0
                        });
                        if(!nested_cuts)break;
                    }
                }
            }
            if(!filter){
                for(auto cut:candidates){
                    if(where == GRB_CB_MIPNODE){
                        addCut(cut.lhs >= cut.rhs);
                        added_cuts[ConnectivityFractional].push_back(cut);
                        cut_count[ConnectivityFractional]+=1;
                    }
                    else{
                        addLazy(cut.lhs >= cut.rhs);
                        added_cuts[Connectivity].push_back(cut);
                        cut_count[Connectivity]+=1;
                    }
                }
                return;
            }
            std::sort(candidates.begin(),candidates.end(),[](const CutInformation& a, const CutInformation& b){
                return a.efficiacy > b.efficiacy;}
            );
            for(int i = 0; i< std::ceil(0.3*candidates.size());i++){
                addCut(candidates[i].lhs >= candidates[i].rhs);
                added_cuts[ConnectivityFractional].push_back(candidates[i]);
                cut_count[ConnectivityFractional]+=1;
                
            }
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