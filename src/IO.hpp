#pragma once

#include "graph.hpp"
#include <filesystem>
#include <fstream>
#include <string>
namespace filesystem = std::filesystem;

namespace IO{
    Graph readDimacs(filesystem::path inputFile,bool unit_weights = false){
        std::ifstream inFile(inputFile);
        std::string line;
        Graph G;
        while (std::getline(inFile, line)){
                if(line[0] == 'c'){
                    continue;
                }
                std::istringstream iss(line);
                if(line[0] == 'p'){
                    std::string c;
                    int nodes, edges;
                    iss >> c >> nodes >> edges;  
                    G = Graph(nodes);
                    continue;
                }
                if(line[0] == 'e'){
                    std::string c;
                    int u, v;
                    double border;
                    iss >> c >> u >> v >> border;
                    G.addEdge(u,v,unit_weights ? 1 : border);
                    continue;
                }
                if(line[0] == 'n'){
                    std::string c;
                    std::string idx;
                    Vertex v;
                    int pop;
                    double xCoord;
                    double yCoord;
                    double feature;
                    iss >> c >> v >> idx >> pop >> xCoord >> yCoord >> feature;
                    G.setAttribute(v,VertexProperties{idx,pop,xCoord,yCoord,feature});
                }
            }
            inFile.close();
            return G;
    }

    void write_solution(const std::vector<int> assign,const Graph& G, const filesystem::path& output_file){
        int max = 0;
        for(int a:assign) max = std::max(a,max);
        std::vector<std::vector<std::string>> partitioning(max+1,std::vector<std::string>());
        for(int i = 0; i< assign.size();i++){
            partitioning[assign[i]].push_back(G.getAttribute(i).NUTId);
        }
        
        std::ofstream out_file(output_file);
        for(const auto& partition:partitioning){
            for(const std::string& id:partition) out_file << id << " ";
            out_file << "\n";
        }
        out_file.close();

    }

    void matrix_to_file(const std::vector<std::vector<double>>& matrix, filesystem::path path){
        std::ofstream out_file(path);
        for(const auto& row:matrix){
            for(const auto& value:row){
                out_file << value << " ";
            }
            out_file << "\n";
        }
        out_file.close();
    }
    std::pair<std::string,std::string> seperate_comma(const std::string& s){
        int delimIdx = s.find(",");
        if(delimIdx == std::string::npos){
            return std::make_pair(s,"");
        }
        std::string param = s.substr(0,s.find(","));
        std::string value = s.substr(s.find(",")+1);
        return std::make_pair(param,value);
    }

    Graph readEdgeList(filesystem::path inputFile){
        std::ifstream inFile(inputFile);
        std::string line;
        Graph G = Graph();
        while (std::getline(inFile, line)){
                if(line[0] == '#'){
                    continue;
                }
                std::istringstream iss(line);

                std::string c;
                int u, v;
                double weight;
                iss >> u >> v >> weight;
                G.addEdge(u-1,v-1,weight);
            }
            inFile.close();
            return G;
    }
}