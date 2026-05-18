#pragma once
#include "commons.hpp"
#include "graph.hpp"
#include <vector>
#include <queue>
std::vector<int> hessToAssignment(const std::vector<std::vector<int>> & x){
    int n = x.size();
    std::vector<int> compIdx(n,-1);
    int idxCounter = 0;
    std::vector<int> assignments(n);
    for(int i = 0;i< n;i++){
        for(int j = 0;j< n;j++){
            if(x[i][j] == 1){
                if(compIdx[j] == -1){
                    compIdx[j] = idxCounter++;
                }
                assignments[i] = compIdx[j];
            }
        }
    }
    return assignments;
}
std::vector<int> labelingtoAssignment(const std::vector<std::vector<int>> & x){
    int k = x.size();
    int n = x[0].size();
    std::vector<int> assignments(n);
    for(int i = 0;i< k;i++){
        for(int v = 0;v< n;v++){
            if(x[i][v] == 1){
                assignments[v] = i;
            }
        }
    }
    return assignments;
}

std::vector<int> partition_to_assignment(const std::vector<std::vector<Vertex>> & partitioning){
    int k = partitioning.size();
    int n = 0;
    for(int i = 0; i< k;i++){
        n += partitioning[i].size();
    }
    std::vector<int> assignments(n);
    for(int i = 0; i< k;i++){
        for(Vertex v: partitioning[i]){
            assignments[v] = i;
        }
    }
    return assignments;
}
std::vector<std::vector<int>> get_regs_from_ass(const std::vector<int>& assignment){
    const int max_region = *std::max_element(assignment.begin(),assignment.end());
    std::vector<std::vector<int>> regions(max_region+1,std::vector<int>());
    for(Vertex u = 0; u < assignment.size();u++){
        regions[assignment[u]].push_back(u);
    }
    return regions;

}
bool objectives_equal(double checked_obj, double true_obj){
    const double rel_tolerance = 0.0001;
    const double abs_tolerance = rel_tolerance*(std::abs(true_obj)+1);
    return (std::abs(checked_obj-true_obj) < abs_tolerance);
}
//Checks if the reported objective value corresponds to the value of the solution
bool check_pregions_objective(const std::vector<std::vector<int>>& regions, double obj, const Graph& G){
    double true_obj = 0;
    for(const std::vector<int>& reg:regions){
        for(int i = 0; i< reg.size();i++){
            for(int j = 0; j< i;j++){
                true_obj += std::abs(G.getAttribute(reg[i]).feature-G.getAttribute(reg[j]).feature);
            }
        }
    }
    return objectives_equal(true_obj,obj);
}
//Sort partitions according to smallest vertex index
void sort_partition(std::vector<std::vector<int>>& partition){
    std::sort(partition.begin(),partition.end(),[](const std::vector<int>& a, const std::vector<int>& b){
        return *std::min_element(a.begin(),a.end()) < *std::min_element(b.begin(),b.end());
    });
}
//Checks if two provided solutions are equal, under the assumption the regions of both partitions are sorted
bool check_equality_sorted(const std::vector<std::vector<int>>& partition1, const std::vector<std::vector<int>>& partition2){
    if(partition1.size() != partition2.size()){
        return false;
    }
    for(int i = 0; i< partition1.size();i++){
        const std::set<int> set1(partition1[i].begin(),partition1[i].end());
        const std::set<int> set2(partition2[i].begin(),partition2[i].end());
        if(set1 != set2) return false;
    }
    return true;
}
//Checks if provided solutions are equal, under the assumption the regions of both partitions are sorted
bool check_equality_sorted(const std::vector<std::vector<std::vector<int>>>& partitions){
    for(int reg_idx = 0; reg_idx< partitions.size()-1;reg_idx++){
        if(!check_equality_sorted(partitions[reg_idx],partitions[reg_idx+1])){
            return false;
        }
    }
    return true;
}

bool checkConnectivity(const Graph& G, const std::vector<int>& assignments,int k){
    int n = assignments.size();
    for(int i = 0; i< k;i++){
        std::queue<Vertex> bfs;
        std::vector<bool> visited(n,false);
        bool compFound = false;
        for(Vertex a = 0; a < n; a++){
            if(visited[a]||assignments[a] != i)continue;
            if(compFound)
            return false;
            compFound = true;
            bfs.push(a);
            visited[a] = true;
            while(!bfs.empty()){
                Vertex v = bfs.front();
                bfs.pop();
                for(auto [wit,wend] = G.getNeighbours(v);wit != wend;wit++){
                    Vertex w = *wit;
                    if(!visited[w] && assignments[w] == i){
                        visited[w] = true;
                        bfs.push(w);
                    }
                    
                }
            }
        }
    }
   return true;
}

void check_computed_solutions(
    const std::vector<std::vector<int>>& model_assignments,
    const std::vector<double>& objective_values, 
    const Graph& G, int k,
    bool connectivity = true,[[maybe_unused]]SymmetryBreaking symm = None){
    std::vector<std::vector<std::vector<int>>> sorted_partitions;
    for(const std::vector<int>& m: model_assignments){
        //Push partitions to vector
        sorted_partitions.push_back(get_regs_from_ass(m));
        //Check if each partition contains exactly k regions
        assert(sorted_partitions.back().size() == k);
        //Check for connectivity if it is enforced
        if(connectivity){
            assert(checkConnectivity(G,m,k));
        }
    }
    for(int i = 0; i< model_assignments.size();i++){
        //Check if objectives are equal
        if(i < model_assignments.size()-1){
            bool are_equal = objectives_equal(objective_values[i],objective_values[i+1]);
            if(! are_equal){
                std::cout << "Not equal: " << objective_values[i] << " vs " << objective_values[i+1] << "\n";
            }
            assert(objectives_equal(objective_values[i],objective_values[i+1]));
        }
        //Check if reported objectives correspond to actual objectives
        assert(check_pregions_objective(sorted_partitions[i],objective_values[i],G));
        sort_partition(sorted_partitions[i]);
    }
    if(sorted_partitions.size() > 1){
        //Notify if solutions differ
        bool all_solutions_equal = check_equality_sorted(sorted_partitions);
        if(!all_solutions_equal){
            std::cout << "Notification: Solutions are different \n";
        }
    }
    //TODO: check if symmetry breaking works
}

void check_obj_bound_consistency(const std::vector<std::pair<double,double>>& bounds){
    if(bounds.size()<1)return;
    double tightest_lb = bounds[0].first;
    double tightest_ub = bounds[0].second;
    for(auto [lb,ub]:bounds){
        assert(ub >= tightest_lb - ERR_TOLERANCE);
        assert(lb <= tightest_ub + ERR_TOLERANCE);
        if(ub < tightest_ub){
            tightest_ub = ub;
        }
        if(lb > tightest_lb){
            tightest_lb = lb;
        }
    }
}

bool check_valid_partitioning(const std::vector<std::vector<Vertex>>& partition){
    int total_vertices = 0;
    for(const auto& p:partition){
        if(p.size() < 1) return false;
        total_vertices += p.size();
    }
    std::vector<bool> assigned(total_vertices,false);
    for(const auto& p:partition){
        for(Vertex v:p){
            if(v > total_vertices-1) return false;
            if(assigned[v]) return false;
            assigned[v] = true;
        }
    }
    return true;
}