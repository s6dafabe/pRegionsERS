#pragma once
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/graph_traits.hpp"
#include <string>
#include <variant>
#include <utility>
#include "gurobi_c++.h"

/*
===================================== Math operations==========================================================
*/
#define square(x) ((x)*(x))
#define ERR_TOLERANCE 1E-6
#define DOUBLE_EQUALITY(x,y) (std::abs(x-y) < ERR_TOLERANCE)
/*
===================================== Graph definitions==========================================================
*/
struct VertexProperties{
    std::string NUTId;
    int population;
    double xCoord;
    double yCoord;
    double feature;

};
struct EdgeProperties{
    double borderLength;
};
struct VertexValue{
    double val;
};
struct EdgeValue{
    double val;
};
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,VertexProperties,EdgeProperties> BoostGraph;
typedef boost::graph_traits<BoostGraph>::vertex_descriptor Vertex;
typedef boost::graph_traits<BoostGraph>::edge_descriptor Edge;
typedef boost::graph_traits<BoostGraph>::adjacency_iterator adjacency_iterator;
typedef boost::graph_traits<BoostGraph>::edge_iterator edge_iterator;
typedef boost::graph_traits<BoostGraph>::out_edge_iterator out_edge_iterator;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,boost::no_property,boost::no_property> DirectedGraph;
typedef boost::graph_traits<DirectedGraph>::vertex_descriptor DVertex;
typedef boost::graph_traits<DirectedGraph>::edge_descriptor DEdge;
typedef boost::graph_traits<DirectedGraph>::adjacency_iterator Dadjacency_iterator;
typedef boost::graph_traits<DirectedGraph>::out_edge_iterator Dedge_iterator;


typedef boost::adjacency_list_traits<boost::vecS, boost::vecS, boost::directedS> Traits;
typedef Traits::edge_descriptor FlowEdge;
struct FlowProperties{
    double capacity;
    double residual_capacity;
    FlowEdge reverse_edge;
};
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,boost::no_property,FlowProperties> FlowGraph;


typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,boost::no_property,FlowProperties> BKGraph;

/*
===================================== VARIABLE INDEXING==========================================================
*/
/*Model variable definitions
X: Centroid/Region/Ordering Assignment variables
Y: Edge cut variables
R: Centroid variables
T: Pairwise variables
A: Objective variables
Fshir: Flow variables for SHIR model
Fscf: Flow variables for SCF model
*/
enum Vars{X,Y,R,T,A,Aprime,Fshir,Fscf,NUMBER_MODEL_VARS};
enum IndexingTypes{OrderedPairwise,UnorderedPairwise,DirectedEdges,UndirectedEdges,EdgesRegionPair,Vertices,NUMBER_INDEXING_TYPES};
inline int getIndex(int i,IndexingTypes idxType){
    switch(idxType){
        case Vertices:
            return i;
        default:
            std::cout << "Wrong number of indices" << "\n";
            return -1;

    }
}
inline int getIndex(int i, int j,IndexingTypes idxType, int n, const std::unordered_map<int,int>& eMapDir, const std::unordered_map<int,int>& eMapUn){
    const int idx = i+j*n;
    switch(idxType){
        case UnorderedPairwise:
            return idx;
        case OrderedPairwise:
            if(i < j) return j*(j-1)/2+i; 
            return i*(i-1)/2+j;
        case DirectedEdges:
            return eMapDir.at(idx);
        case UndirectedEdges:
            return eMapUn.at(idx);
        default:
            std::cout << "Wrong number of indices" << "\n";
            return -1;

    }
}
inline int getIndex(int i, int j,int l, IndexingTypes idxType, int n, int k, [[maybe_unused]]const std::unordered_map<int,int>& eMapDir, [[maybe_unused]]const std::unordered_map<int,int>& eMapUn){
    switch(idxType){
        case EdgesRegionPair:
            return i+eMapDir.at(j+l*n)*k;
        default:
            std::cout << "Wrong number of indices" << "\n";
            return -1;

    }
}
//Variable indexing for variables with two indices
inline int getVarIdx(int i, int j,Vars vars, int n, const std::unordered_map<int,int>& eMapDir, const std::unordered_map<int,int>& eMapUn){
    const int idx = i+j*n;
    switch(vars){
        case X:
            return idx;
        case Y:
            return eMapUn.at(idx);
        case R:
            return idx;
        case A:
            return idx;
        case T:
            if(i < j) return j*(j-1)/2+i;
            return i*(i-1)/2+j;
        case Aprime:
            return idx;
        case Fscf:
            return eMapDir.at(idx);
        default:
            throw std::invalid_argument("Wrong number of indices\n");
            

    }
}
//Variable indexing for variables with three indices
inline int getVarIdx(int i, int j,int l,Vars vars, int n,const std::unordered_map<int,int>& eMapDir,[[maybe_unused]]std::unordered_map<int,int>& eMapUn){
    switch(vars){
        case Fshir:
            return i+eMapDir.at(j+l*n)*n;
        default:
            throw std::invalid_argument("Wrong number of indices\n");
    }
}

//Hash for vertex pair

inline int pairHash(unsigned i, unsigned j, unsigned n){
    return i+j*n;
}

/*
===================================== CUT DEFINITIONS ==========================================================
*/


enum Cuts{Clique,GeneralClique,Wheel,TwoPartition,InDegree,Connectivity,ConnectivityFractional,Cycles,NUMBER_CUTS};


/*
===================================== OTHER DEFINITIONS ==========================================================
*/

//typedef std::variant<std::pair<GRB_IntParam,int>,std::pair<GRB_DoubleParam,double>,std::pair<GRB_StringParam,std::string>> GRBparameter;
typedef std::pair<std::string,std::string> GRBparameter;
typedef std::variant<int,double,std::string> GRBattribute;
typedef std::variant<GRB_IntAttr,GRB_DoubleAttr,GRB_StringAttr> GRBattributeType;
typedef std::variant<GRB_IntParam,GRB_DoubleParam,GRB_StringParam> GRBparamType;
typedef std::variant<int,double,std::string> GRBparamReturn;
enum SymmetryBreaking{None,Seed,Weak,Strong};
enum kPartition{Standard,Extended};
enum ModelType{Base,NodeCluster,EdgeRepresentative};
enum Transitivity{Full,Lazy};

/*
===================================== STRING DEFINITIONS FOR OUTPUT ==========================================================
*/
std::map<std::string,SymmetryBreaking> symmetrymap = {
    {"Seed",Seed},
    {"Weak",Weak},
    {"Strong",Strong}
};


std::map<std::string,Cuts> cutmap = {
    {"Clique",Clique},
    {"GeneralClique",GeneralClique},
    {"Wheel",Wheel},
    {"2-Partition",TwoPartition},
    {"Indegree",InDegree},
    {"Connectivity",Connectivity},
    {"ConnectivityFractional",ConnectivityFractional},
    {"Cycles",Cycles}
};

std::vector<std::string> cut_strings = {
    "Clique",
    "GeneralClique",
    "Wheel",
    "2-Partition",
    "Indegree",
    "Connectivity",
    "ConnectivityFractional",
    "Cycles"
};
