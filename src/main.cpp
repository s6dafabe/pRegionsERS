#include <iostream>
#include <string>
#include <filesystem> // Requires C++17 or later
#include "SolutionChecking.hpp"
#include "Models.hpp"
#include "DuqueModels.hpp"
#include "Cuts.hpp"
#include "Connectivity.hpp"
int main(int argc, char *argv[]) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <path/to/instance> <k> <Flow:Tree:ERS:ERSTree> <timelimit>" << std::endl;
        return 1; // Indicate an error
    }

    std::string directoryPath = argv[1];
    int k = 0;
    int timelimit = 0;
    try{
        k = std::stoi(argv[2]);
        timelimit = std::stoi(argv[4]);
    }
    catch(std::exception e){
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    std::string model_type = argv[3];
    std::filesystem::path path(directoryPath);

    if (!std::filesystem::is_regular_file(path)) {
        std::cerr << "Error: '" << directoryPath << "' is not a valid file." << std::endl;
        return 1; // Indicate an error
    }
    std::cout << "Reading instance: " << directoryPath << std::endl;
    try{
        Graph graph = IO::readDimacs(directoryPath);
        if(!graph.is_connected()){
            std::cerr << "Input graph must be connected \n";
            return 1;
        }
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "experiments.log");
        env.start();
        BaseModel* model;
        if(model_type == "Flow"){
            model = new FormulationNC(&env,graph,k);
            static_cast<FormulationNC*>(model)->dissimilarity_standard();
            static_cast<FormulationNC*>(model)->connectivity_shir();
            static_cast<FormulationNC*>(model)->add_symmetry_breaking(Seed);
        }
        else if(model_type == "Tree"){
            model = new TreeModel_Duque(&env,graph,k);
        }
        else if(model_type == "ERS"){
            model = new FormulationER(&env,graph,k,1,Standard);
            static_cast<FormulationER*>(model)->addCut(ConnectivityFractional);
            static_cast<FormulationER*>(model)->addCut(GeneralClique);
            static_cast<FormulationER*>(model)->full_clique_constraint();
            static_cast<FormulationER*>(model)->total_edge_constraint();
            static_cast<FormulationER*>(model)->cut_filter(true);
            static_cast<FormulationER*>(model)->cut_locsrch(true);
            static_cast<FormulationER*>(model)->nested_cut(true);
        }
        else if (model_type == "ERSTree"){
            model = new FormulationER(&env,graph,k,1,Standard);
            static_cast<FormulationER*>(model)->addCut(ConnectivityFractional);
            static_cast<FormulationER*>(model)->addCut(GeneralClique);
            static_cast<FormulationER*>(model)->full_clique_constraint();
            static_cast<FormulationER*>(model)->addCut(Cycles);
            static_cast<FormulationER*>(model)->cut_filter(true);
            static_cast<FormulationER*>(model)->cut_locsrch(true);
            static_cast<FormulationER*>(model)->nested_cut(true);
        }
        else{
            std::cerr << "Unknown model type \n";
            return 1;
        }
        const std::vector<GRBparameter> parameters = {
            {"Timelimit",std::to_string(timelimit)},
            {"Threads","1"},
            {"PreCrush","1"}
        };
        model->optimize(parameters);
        int status = std::get<int>(model->get_attributes(GRB_IntAttr_Status));
        std::map<int,std::string> status_name;
        status_name[2] = "Optimal";
        status_name[3] = "Infeasible";
        status_name[9] = "Time limit reached";
        std::cout << "Model status: " << status_name[status] << std::endl;
        std::cout << "Primal bound: " << std::get<double>(model->get_attributes(GRB_DoubleAttr_ObjVal)) << std::endl;
        std::cout << "Dual bound: " << std::get<double>(model->get_attributes(GRB_DoubleAttr_ObjBound))<< std::endl;
        std::cout << "Runtime: " << std::get<double>(model->get_attributes(GRB_DoubleAttr_Runtime))<< std::endl;
        std::cout << "B&B Nodes: " << std::get<double>(model->get_attributes(GRB_DoubleAttr_NodeCount))<< std::endl;
        //If at least one solution found, check if valid:
        if(std::get<int>(model->get_attributes(GRB_IntAttr_SolCount)) > 0){
            assert(model->check_validity());
            auto solution = get_regs_from_ass(model->return_solution());
            std::cout << "========== Best found solution ============ \n";
            int i = 0;
            for(auto s:solution){
                std::cout << "Partition " << i++ << ": ";
                for(auto v:s){
                std::cout << v << " ";
                }
                std::cout << "\n";
            }

        }
    }
    catch (GRBException e) {
        std::cout << "Error number: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;

        return 1;
    }
    return 0;
}