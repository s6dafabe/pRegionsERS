#pragma once
#include "gurobi_c++.h"
#include "graph.hpp"
#include <memory>
#include <cmath>
#include "commons.hpp"
#include <unordered_map>
#include <string>
#include "VariableManager.hpp"
#include "Cuts.hpp"
#include <algorithm>
#include <fstream>
//Placeholder callback
class EmptyCallback:public GRBCallback{
    protected:
        void callback(){}
};

class BaseModel{
    protected:
        GRBModel model;
        VariableManager modelVars;
        Graph G;
        int n,m,k;
        ModelType m_type;
        double lambda;
        std::vector<std::vector<double>> weights;
        double dissimilarity(Vertex u, Vertex v){
            return std::abs(G.getAttribute(u).feature - G.getAttribute(v).feature);
        }
    public:
        GRBVar& get_var(int i,Vars var)const{
            return modelVars.getVar(i,var);
        }
        GRBVar& get_var(int i,int j,Vars var)const{
            return modelVars.getVar(i,j,var);
        }
        GRBVar& get_var(int i,int j,int l,Vars var)const{
            return modelVars.getVar(i,j,l,var);
        }
        BaseModel(GRBEnv* env,const Graph& G,int k,double lambda):
        model(env),
        modelVars(G,k),
        G(G),n(G.numVertices()),m(G.numEdges()),k(k),lambda(lambda),weights(n,std::vector<double>()){
            m_type=Base;
            //Compute weights of weighted sum of vertex dissimilarity and border length
            for(Vertex u = 1; u < n;u++){
                weights[u] = std::vector<double>(u,0);
                for(Vertex v = 0; v < u;v++){
                    weights[u][v] = lambda*dissimilarity(u,v);
                }
            }
            double edgeweight_sum = 0;
            for(Vertex u = 0; u < n;u++){
                for(auto[vit,end] = G.getNeighbours(u);vit != end;vit++){
                    Vertex v = *vit;
                    if(u < v) continue;
                    weights[u][v] += (lambda - 1) * G.borderLength(u,v);
                    edgeweight_sum += (1 - lambda) * G.borderLength(u,v);
                }
            }
            model.set(GRB_DoubleAttr_ObjCon,edgeweight_sum);
        }
        virtual void optimize(const std::vector<std::pair<std::string,std::string>>& parameters = {
            std::make_pair("TimeLimit","1000.0"),
            std::make_pair("Threads","1")
        }){
            for(auto p:parameters){
               model.set(p.first,p.second);
            }
            model.optimize();
        }
        void set_attribute(const std::pair<std::string,std::string>& p){
            model.set(p.first,p.second);
        }
        double get_value(int i, int j,Vars v = X)const {
            return get_var(i,j,v).get(GRB_DoubleAttr_X);
        }
        int getN(){
            return n;
        }
        int getK(){
            return k;
        }
        std::vector<std::vector<double>> get_matrix(){
            int a = model.get(GRB_IntAttr_NumVars);
            int b = model.get(GRB_IntAttr_NumConstrs);
            std::vector<std::vector<double>> matrix;
            GRBVar* all_vars = model.getVars();
            GRBConstr* all_constr = model.getConstrs();
            for(int i = 0; i< b ;i++){
                char sense = all_constr[i].get(GRB_CharAttr_Sense);
                if(sense == '='|| sense == '<'){
                    matrix.push_back(std::vector<double>(a,0));
                    for(int j = 0; j< a ;j++){
                        matrix.back()[j] = model.getCoeff(all_constr[i],all_vars[j]);
                    }
                }
                if(sense == '='|| sense == '>'){
                    matrix.push_back(std::vector<double>(a,0));
                    for(int j = 0; j< a ;j++){
                        matrix.back()[j] = -model.getCoeff(all_constr[i],all_vars[j]);
                    }
                }
                    
            }
            free(all_vars);
            free(all_constr);
            return matrix;
        }
        GRBattribute get_attributes(GRBattributeType attr)const{
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
        void fix_variable(double val, int i, int j,Vars type){
            model.addConstr(get_var(i,j,type) == val);
        }
        GRBparamReturn get_attributes(GRBparamType attr)const{
            if(std::holds_alternative<GRB_IntParam>(attr)){
                GRB_IntParam a = std::get<GRB_IntParam>(attr);
                return model.get(a);
            }
            if(std::holds_alternative<GRB_DoubleParam>(attr)){
                GRB_DoubleParam a = std::get<GRB_DoubleParam>(attr);
                return model.get(a);
            }
            if(std::holds_alternative<GRB_StringParam>(attr)){
                GRB_StringParam a = std::get<GRB_StringParam>(attr);
                return model.get(a);
            }
            return 0;
        }
        const Graph& get_graph()const{
            return G;
        }
        virtual int get_cut_count([[maybe_unused]] int cut)const{
            return 0;
        }
        virtual int number_tight_cuts([[maybe_unused]] Cuts cut)const{
            return 0;
        }
        virtual void set_init_solution([[maybe_unused]] const std::vector<std::vector<Vertex>>& init){}
        virtual bool check_validity()const{return true;}
        virtual std::vector<int> return_solution()const{return {};}
        virtual ~BaseModel() = default;
        ModelType get_type()const{
            return m_type;
        }
};


class FormulationNC:public BaseModel{
    protected:
        int rg;
    public:
        FormulationNC(GRBEnv* env,Graph& G,int k,double lambda = 1):BaseModel(env,G,k,lambda){
            m_type = NodeCluster;
            rg = k;
            modelVars.addVariableType(model.addVars(n*k,GRB_BINARY),X,UnorderedPairwise);
            //Each vertex is assigned to exactly 1 district
            for(int u = 0; u < n;u++){
                GRBLinExpr sum;
                for(int i = 0; i<k;i++){
                    sum += get_var(u,i,X);
                }
                model.addConstr(sum == 1);
            }
            //Each district must contain at least one vertex
            for(int i = 0;i < k;i++){
                GRBLinExpr sum = 0;
                for(Vertex u = 0; u < n;u++){
                    sum += get_var(u,i,X);
                }
                model.addConstr(sum >= 1);
            }
            model.update();

        }
        void dissimilarity_standard(){
            modelVars.addVariableType(model.addVars(n*(n-1)/2,GRB_BINARY),T,OrderedPairwise);
            for(int u = 0; u < n;u++){
                for(int v = 0; v < u;v++){
                    get_var(u,v,T).set(GRB_DoubleAttr_Obj,weights[u][v]);
                    get_var(u,v,T).set(GRB_IntAttr_BranchPriority,0);
                    for(int w = 0; w < rg;w++){
                        model.addConstr(get_var(u,v,T) >= get_var(u,w,X)+get_var(v,w,X)-1);
                    }
                }
            }
            model.update();
        }
        void dissimilarity_quadratic_objective(){
            GRBQuadExpr sum;
            for(int i = 0; i < rg;i++){
                for(int u = 0; u < n;u++){
                    for(int v = 0; v < u;v++){
                        sum += weights[u][v]*get_var(u,i,X)*get_var(v,i,X);
                    }
                }
            }
            model.setObjective(model.getObjective()+sum,GRB_MINIMIZE);
            model.update();
        }
        void dissimilarity_quadratic_variables(){
            modelVars.addVariableType(model.addVars(n*(n-1)/2,GRB_BINARY),T,OrderedPairwise);
            for(int u = 0; u < n;u++){
                for(int v = 0; v < u;v++){
                    get_var(u,v,T).set(GRB_DoubleAttr_Obj,weights[u][v]);
                    get_var(u,v,T).set(GRB_IntAttr_BranchPriority,0);
                    GRBQuadExpr sum;
                    for(int w = 0; w < rg;w++){
                        sum += get_var(u,w,X)*get_var(v,w,X);
                    }
                    model.addQConstr(get_var(u,v,T) == sum);
                }
            }
            model.update();
        }
        void add_symmetry_breaking(SymmetryBreaking symm){
            if(symm == None)return;

            model.addConstr(get_var(0,0,X) == 1);
            if(modelVars.variableExists(R)){
                model.addConstr(get_var(0,0,R) == 1);
            }
            
            if(symm==Seed) return;
            for(int v = 0; v < n; v++){
                for(int i = v+1; i < k;i++){
                    model.addConstr(get_var(v,i,X) == 0);

                }
                if(symm == Strong){
                    for(int i = 1; i < std::min(k,v+1);i++){
                        GRBLinExpr sum = 0;
                        for(int w = i-1;w < v;w++){
                            sum += get_var(w,i-1,X);
                        }
                        model.addConstr(get_var(v,i,X) <= sum);
                    }
                    if(modelVars.variableExists(R)){
                        for(int i = 0; i < k;i++){
                            if(v == 0)continue;
                            model.addConstr(get_var(v,i,R) <= 1-get_var(v-1,i,R));
                        
                        }
                    }
                }
            }
        }
        void add_fvsgc_constraint(){
            if(!modelVars.variableExists(T)){
                std::cout << "Does not contain pairing variables \n";
            }
            GRBLinExpr sum = 0;
            for(Vertex u = 0; u < n;u++){
                for(Vertex v = 0; v < u;v++){
                    sum += get_var(u,v,T);
                }
            }
            const int q = static_cast<int>(n/k);
            const int r = n % k;
            const int limit = (q+1)*q*r/2+q*(q-1)*(k-r)/2;
            model.addConstr(sum >= limit);
        }
        void connectivity_scf(bool biDirConst = false);
        void connectivity_shir();
        //Output solution as assigment to regions
        std::vector<int> return_solution()const override{
            std::vector<int> solution(n,-1);
            for(Vertex u = 0;u< n;u++){
                for(int i = 0;i< k;i++){
                    if(get_value(u,i,X) > 0.5)solution[u] = i;
                }
                assert(solution[u] > -1);
            }
            return solution;
            
        }
        bool check_validity()const override{
            for(int i = 0; i< k;i++){
                int sum = 0;
                for(Vertex v = 0; v < n;v++){
                    sum += get_value(v,i); 
                }
                if(sum < 1-ERR_TOLERANCE || sum > n+ERR_TOLERANCE){
                    return false;
                }
            }
            for(Vertex v = 0; v < n;v++){
                int sum = 0;
                for(int i = 0; i< k;i++){
                    sum += get_value(v,i);
                } 
                if(!DOUBLE_EQUALITY(sum,1)){
                    return false;
                }
            }
            return true;
        }

        void setWarmStart(int u, int v, int value){
            get_var(u,v,X).set(GRB_DoubleAttr_Start,value);
        }
        void set_init_solution([[maybe_unused]] const std::vector<std::vector<Vertex>>& init)override{
            std::vector<std::pair<Vertex,int>> lowest_vertices;
            for(int i = 0; i < init.size(); i++){
                Vertex lowest = n;
                for(auto u:init[i]){
                    if(u < lowest)lowest = u;
                }
                lowest_vertices.push_back({lowest,i});
            }
            
            std::sort(lowest_vertices.begin(),lowest_vertices.end());
            for(int i = 0; i < init.size(); i++){
                for(auto u:init[i]){
                    setWarmStart(u,lowest_vertices[i].second,1);
                }
            }
        }
};

class FormulationER:public BaseModel{
    private:
        CallbackER* callback;
        kPartition partition_constrs;
        GRBEnv* env;
        std::map<std::tuple<Vertex,Vertex,Vertex>,GRBConstr> stored_constraints;
    public:
        int idx(int i, int j){
            return i*(i-1)/2+j;
        }
        double get_varVal(int i, int j, Vars type) const {
            return modelVars.getValue(i,j,type);
        }
        FormulationER(GRBEnv* env,const Graph& G,int k, double lambda = 1, kPartition part_const = Extended,Transitivity trans = Full):
        BaseModel(env,G,k,lambda),partition_constrs(part_const),env(env){
            m_type = EdgeRepresentative;

            modelVars.addVariableType(model.addVars(n*(n-1)/2,GRB_BINARY),X,OrderedPairwise);
            modelVars.addVariableType(model.addVars(n,GRB_BINARY),R,Vertices);

            //Transitivity constraints
            for(Vertex u = 0; u < n ;u++){
                for(Vertex v = 0; v < u; v++){

                    get_var(u,v,X).set(GRB_DoubleAttr_Obj,weights[u][v]);
                    for(Vertex w = 0; w < v; w++){
                        GRBConstr cstr1 = model.addConstr(get_var(u,v,X) >= get_var(u,w,X)+get_var(v,w,X)-1);
                        GRBConstr cstr2 = model.addConstr(get_var(u,w,X) >= get_var(u,v,X)+get_var(v,w,X)-1);
                        GRBConstr cstr3 = model.addConstr(get_var(v,w,X) >= get_var(u,v,X)+get_var(u,w,X)-1);
                        if(trans == Full)continue;
                        else{
                            cstr1.set(GRB_IntAttr_Lazy,3);
                            cstr2.set(GRB_IntAttr_Lazy,3);
                            cstr3.set(GRB_IntAttr_Lazy,3);
                        
                        }
                    } 
                }
            }
            GRBLinExpr csum = 0;
            switch(partition_constrs){
                case Standard:{
                    /*
                    Models rep by two inequalities: 
                    - u cannot be rep if paired with lower vertex v
                    - u must be rep if paired with no lower vertex v
                    */
                    for(Vertex u = 0; u < n ;u++){
                        GRBLinExpr sum = 0;
                        for(Vertex v = 0; v < u; v++){
                            model.addConstr(get_var(u,R) <= 1-get_var(u,v,X));
                            sum += get_var(u,v,X);
                        }
                        model.addConstr(get_var(u,R) >= 1-sum);
                        csum += get_var(u,R);
                    }
                    break;
                }
                case Extended:{
                    /*
                    Models rep by one quad inequality: 
                    - u is rep iff no lower vertex is rep and assigned to u
                    */
                    modelVars.addVariableType(model.addVars(n*(n-1)/2,GRB_BINARY),T,OrderedPairwise);
                    for(Vertex u = 0; u < n ;u++){
                        GRBLinExpr sum = 0;
                        for(Vertex v = 0; v < u; v++){
                            model.addConstr(get_var(u,v,T) >= get_var(v,R)+get_var(u,v,X)-1);
                            model.addConstr(get_var(u,v,T) <= get_var(v,R));
                            model.addConstr(get_var(u,v,T) <= get_var(u,v,X));
                            model.addConstr(get_var(u,R) <= 1-get_var(u,v,X));//Necessary for theoretical dominance over Standard
                            sum += get_var(u,v,T);
                        }
                        model.addConstr(get_var(u,R) == 1-sum);
                        csum += get_var(u,R);
                    }
                    break;
                }
            }
            model.addConstr(csum == k);
            callback = new CallbackER(modelVars,G);
            model.setCallback(callback);
            model.update();
            
        }
        void full_clique_constraint(){
            GRBLinExpr sum = 0;
            for(Vertex u = 0; u < n;u++){
                for(Vertex v = 0; v < u;v++){
                    sum += get_var(u,v,X);
                }
            }
            const int q = static_cast<int>(n/k);
            const int r = n % k;
            const int limit = (q+1)*q*r/2+q*(q-1)*(k-r)/2;
            model.addConstr(sum >= limit);
        }
        void total_edge_constraint(){
            GRBLinExpr sum = 0;
            for(Vertex u = 0; u < n;u++){
                for(auto [vit,end] = G.getNeighbours(u); vit != end;vit++){
                    Vertex v = *vit;
                    if(u < v)sum += get_var(u,v,X);
                }
            }
            model.addConstr(sum >= n-k);
        }
        void neighboorhood_constraint(double upper_bound){
            std::vector<Vertex> non_isolated;
            for(Vertex u = 0; u < n;u++){
                FormulationER relaxation_isolated(env,G,k,lambda);
                relaxation_isolated.set_attribute({"OutputFlag","0"});
                relaxation_isolated.set_attribute({"LogFile",""});
                relaxation_isolated.set_attribute({"LogToConsole","0"});
                for(Vertex v = 0; v < n;v++){
                    if(u == v) continue;
                    relaxation_isolated.fix_variable(0,u,v,X);
                }
                relaxation_isolated.solve_relaxation();

                if(std::get<double>(relaxation_isolated.get_attributes(GRB_DoubleAttr_ObjVal)) > upper_bound)
                non_isolated.push_back(u);
            }
            for(Vertex u:non_isolated){
                GRBLinExpr sum = 0;
                for(auto [vit,end] = G.getNeighbours(u); vit != end;vit++){
                    Vertex v = *vit;
                    sum += get_var(u,v,X);
                }
                model.addConstr(sum >= 1);
            }
        }
        void add_edge_vars(){
            modelVars.addVariableType(model.addVars(m*2,GRB_BINARY),Y,DirectedEdges);
            GRBLinExpr total_sum = 0;
            for(auto [eit,end] = G.getEdges(); eit != end; eit++){
                Edge e = *eit;
                auto [u,v] = G.get_incident(e);
                model.addConstr(get_var(u,v,Y) + get_var(v,u,Y) <= get_var(u,v,X));
                total_sum += get_var(u,v,Y)+ get_var(v,u,Y);
            }
            model.addConstr(total_sum == n-k);
            for(Vertex u = 0; u < n; u++){
                GRBLinExpr out_sum = 0;
                for(auto [vit,end] = G.getNeighbours(u); vit != end; vit++){
                    Vertex v = *vit;
                    out_sum += get_var(u,v,Y);
                }
                model.addConstr(out_sum == 1-get_var(u,R));
            }
        }
        void setWarmStart(int u, int v, int value){
            get_var(u,v,X).set(GRB_DoubleAttr_Start,value);
        }
        void relax_model(bool relax){
            if(relax){
                for(Vertex u = 0; u < n;u++){
                    get_var(u,R).set(GRB_CharAttr_VType,GRB_CONTINUOUS);
                    for(Vertex v = 0; v < u;v++){
                        get_var(u,v,X).set(GRB_CharAttr_VType,GRB_CONTINUOUS);
                        get_var(u,v,T).set(GRB_CharAttr_VType,GRB_CONTINUOUS);
                    }
                }
                model.addVar(0,1,0,GRB_BINARY);

                model.reset();
            }
            else{
                for(Vertex u = 0; u < n;u++){
                    get_var(u,R).set(GRB_CharAttr_VType,GRB_BINARY);
                    for(Vertex v = 0; v < u;v++){
                        get_var(u,v,X).set(GRB_CharAttr_VType,GRB_BINARY);
                        get_var(u,v,T).set(GRB_CharAttr_VType,GRB_BINARY);
                    }
                }
                model.reset();
            }
        }
        void lp_lazy();
        void solve_relaxation();
        void relaxation_experiment();
        void bound_separator_preprocessing();
        //Cutting Plane section
        void addCut(Cuts cut){
            const std::unordered_set<Cuts> containLazyConstraints = {Connectivity,ConnectivityFractional};
            if(containLazyConstraints.find(cut) != containLazyConstraints.end()){
                model.set(GRB_IntParam_LazyConstraints,1);
            }
            if(cut == Cycles && ! modelVars.variableExists(Y)){
                std::cout << "Adding edge variables to model \n";
                add_edge_vars();
            }
            callback->activateCut(cut);
        }
        void triangle_inequalities(){
            if(!modelVars.variableExists(Y)){
                std::cout << "Adding edge variables to model \n";
                add_edge_vars();
            }
            for(auto [eit,end] = G.getEdges();eit != end;eit++){
                auto [a,b] = G.get_incident(*eit);
                for(auto [vit,vend] = G.getNeighbours(a);vit != vend;vit++){
                    Vertex v = *vit;
                    if(G.adjacent(b,v)){
                        GRBLinExpr arc_sum = get_var(a,b,Y)+get_var(b,v,Y)+get_var(v,a,Y);
                        arc_sum += get_var(a,v,Y)+get_var(v,b,Y)+get_var(b,a,Y);
                        model.addConstr(2*arc_sum <= get_var(a,b,X)+get_var(b,v,X)+get_var(v,a,X)+1);
                    }

                }
            }
        }
        //Output solution as assigment to regions
        std::vector<int> return_solution()const override{
            std::vector<int> solution(n,-1);
            std::vector<int> represented_by(n,-1);
            int counter = 0;
            for(Vertex u = 0;u< n;u++){
                for(int v = u-1;v >= 0;v--){
                    if(get_varVal(u,v,X) > 0.5){
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
        bool check_validity()const override{
            for(Vertex u = 0; u < n ;u++){
                for(Vertex v = 0; v < u; v++){
                    for(Vertex w = 0; w < v; w++){
                        if(get_varVal(u,v,X) < get_varVal(u,w,X)+get_varVal(v,w,X)-1-ERR_TOLERANCE)return false;
                        if(get_varVal(u,w,X) < get_varVal(u,v,X)+get_varVal(v,w,X)-1-ERR_TOLERANCE)return false;
                        if(get_varVal(v,w,X) < get_varVal(u,v,X)+get_varVal(u,w,X)-1-ERR_TOLERANCE)return false;
                        //TODO: missing checks for uv=0,uw=1 => vw = 0 etc
                    }

                } 
            }
            
            return true;
        }
        void tune_model(){
            model.set("TuneResults","2");
            model.set("TuneTimeLimit","300");
            model.set("TimeLimit","60");
            model.tune();
            model.getTuneResult(1);
            model.write("tune.prm");

        }
        void connectivity_scf();
        void connectivity_neighseparator();
        void add_uv_neighseparator(Vertex u,Vertex v);
        ~FormulationER(){
            delete callback;
        }
        int get_cut_count(Cuts cut){
            return callback->get_cut_count(cut);
        }
        int get_cut_count(int cut)const override{
            return callback->get_cut_count(static_cast<Cuts>(cut));
        }
        void callback_info(Cuts cut){
            std::ofstream log_file;
            log_file.open("Log_"+cut_strings[cut]+".txt",std::ios_base::app);
            log_file << "Left-hand " 
                        << "Right-hand " 
                        <<"Size "
                        <<"Violation "
                        <<"Efficiacy "
                        <<"Nodes-Tight "
                        <<"Coeff-Sum \n";
            for(auto [lhs,rhs,violation,efficiacy,tght]:callback->get_added_cuts(cut)){
                int coeff_sum = 0;
                for(int i = 0; i< lhs.size();i++){
                    coeff_sum += lhs.getVar(i).get(GRB_DoubleAttr_Obj);
                }
                for(int i = 0; i< rhs.size();i++){
                    coeff_sum += rhs.getVar(i).get(GRB_DoubleAttr_Obj);
                }
                log_file << lhs.getValue() << " "
                        << rhs.getValue() << " "
                        << lhs.size() << " "
                        << violation<< " "
                        << efficiacy<< " "
                        << tght << " "
                        << coeff_sum << "\n";

            }
            log_file.close();   
        }
        void cut_filter(bool filter){
            callback->set_filter(filter);
        }
        void cut_locsrch(bool locsrch){
            callback->set_locsrch(locsrch);
        }
        void cut_clqlimit(const std::string& clqlimit){
            callback->set_clqlimit(clqlimit);
        }
        void nested_cut(bool nested_cut){
            callback->set_nested(nested_cut);
        }
        int number_tight_cuts(Cuts cut)const override{
            int tight = 0;
            for(auto [lhs,rhs,violation,efficiacy,tght]:callback->get_added_cuts(cut)){
                if(std::abs(lhs.getValue()-rhs.getValue()) < 1E-4)tight++;
            }   
            return tight;
        }

        void extended_callback_info(Cuts cut){
            int tight = 0;
            double av_violation = 0;
            double av_efficiacy= 0;
            double av_coeffs= 0;
            double av_violation_tight= 0;
            double av_efficiacy_tight= 0;
            double av_coeffs_tight= 0;
            for(auto [lhs,rhs,violation,efficiacy,tght]:callback->get_added_cuts(cut)){
                if(std::abs(lhs.getValue()-rhs.getValue()) < 1E-4){
                    tight++;
                    av_violation_tight += violation;
                    av_efficiacy_tight += efficiacy;
                    av_coeffs_tight += lhs.size();
                }
                av_violation += violation;
                av_efficiacy += efficiacy;
                av_coeffs += lhs.size();
            }
            const int num_cuts = callback->get_added_cuts(cut).size();
            std::cout << "Avg Violation: " << av_violation/num_cuts << " Avg Efficiacy: " << av_efficiacy/num_cuts <<" Avg Coeffs: "<<av_coeffs/num_cuts << "\n";
            std::cout << "Number Tight: " << tight 
                    << " Avg Violation Tgt: " << av_violation_tight/tight 
                    << " Avg Efficiacy Tgt: " << av_efficiacy_tight/tight 
                    <<" Avg Coeffs Tgt: "<<av_coeffs_tight/tight << "\n";
        }

        void set_init_solution([[maybe_unused]] const std::vector<std::vector<Vertex>>& init)override{
            for(auto p:init){
                for(auto u:p){
                    for(auto v:p)
                    if(u < v) setWarmStart(u,v,1);
                }
            }
        }

        void tofile(const std::string& path){
            model.write(path);
        }
};