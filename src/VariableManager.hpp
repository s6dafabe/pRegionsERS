#pragma once
#include "gurobi_c++.h"
#include "commons.hpp"
#include <unordered_map>
#include "graph.hpp"

/*
Class for managing model variables
*/
class VariableManager{
    private:
        std::vector<GRBVar*> modelVars;
        //Store indexing type of variables, i.e. vertex pairs, edges, etc.
        std::vector<IndexingTypes> indexing;
        int n,m,k;
        //Maps mapping from idx(u,v)->idx(e), where idx(e) in {0,...,m-1} respectively {0,...,2m-1} for directed edges
        std::unordered_map<int,int> edgeMapDir;
        std::unordered_map<int,int> edgeMapUn;
        /*
        get correct index using get_index, using the correct associated indexing type and (if applicable) edges indices stored for each variable type 
        */
        int getVarIdxInternal(int i, Vars var)const{
            return getIndex(i,indexing[var]);
        }
        int getVarIdxInternal(int i, int j,Vars var)const{
            return getIndex(i,j,indexing[var],n,edgeMapDir,edgeMapUn);
        }
        int getVarIdxInternal(int i, int j,int l,Vars var)const{
            return getIndex(i,j,l,indexing[var],n,k,edgeMapDir,edgeMapUn);
        }
        int pairHashInteral(int i, int j)const{
            return pairHash(i,j,n);
        }
    public:
        VariableManager(const Graph& G, unsigned k);
        void addVariableType(GRBVar* vars, Vars varType, IndexingTypes idxType);
        GRBVar& getVar(int i, Vars var)const;
        GRBVar& getVar(int i,int j,Vars var)const;
        GRBVar& getVar(int i,int j,int l,Vars var)const;
        double getValue(int i, int j,Vars v)const;
        bool variableExists(Vars var)const{
            return modelVars[var] != nullptr;
        }
        unsigned getN()const;
        unsigned getK()const;
        unsigned getM()const;
        ~VariableManager();
};
//Constructor: Set up edge indices for the graph
VariableManager::VariableManager(const Graph& G, unsigned k):modelVars(NUMBER_MODEL_VARS,nullptr),indexing(NUMBER_MODEL_VARS,UnorderedPairwise),n(G.numVertices()),m(G.numEdges()),k(k){
    int counter = 0;
    int counterUn = 0;
    for(int u = 0; u < n;u++){
        for(auto[vit,end] = G.getNeighbours(u);vit != end;vit++){
            Vertex v = *vit;
            edgeMapDir[pairHashInteral(u,v)] = counter++;
            if(u < v){
                edgeMapUn[pairHashInteral(u,v)] = counterUn++;
            }
        }
    }
}
//Adds model variables <vars> as variables <varType> with indexing type <idxType>
void VariableManager::addVariableType(GRBVar* vars, Vars varType, IndexingTypes idxType){
    modelVars[varType] = vars;
    indexing[varType] = idxType;
}
//Returns a variable with single index
GRBVar& VariableManager::getVar(int i,Vars var)const{
    return modelVars[var][getVarIdxInternal(i,var)];
}
//Returns a variable with double index
GRBVar& VariableManager::getVar(int i,int j,Vars var)const{
    return modelVars[var][getVarIdxInternal(i,j,var)];
}
//Returns a variable with triple index
GRBVar& VariableManager::getVar(int i,int j,int l, Vars var)const{
    return modelVars[var][getVarIdxInternal(i,j,l,var)];
}
//Returns the value of a variable with double index
double VariableManager::getValue(int i, int j,Vars v)const{
    return getVar(i,j,v).get(GRB_DoubleAttr_X);
}

unsigned VariableManager::getN()const{
    return n;
}
unsigned VariableManager::getK()const{
    return k;
}
unsigned VariableManager::getM()const{
    return m;
}
//Destructor: delete all model variables
VariableManager::~VariableManager(){
    for(int var = 0; var < NUMBER_MODEL_VARS;var++){
        if(modelVars[var] != nullptr){
            delete[] modelVars[var];
        }
    }
}


class PairingVariables:public VariableManager{};