#pragma once
#include "type.h"
#include "graph.h"
#include "separation.h"
#include <ilcplex/ilocplex.h>
#include <map>
#include <memory>
#include <unordered_map>
#include <string>

using namespace std;

typedef unordered_map<NODE, IloNumVar> NODE_VARS;
typedef map<NODE, vector<NODE>> ADJ_LIST;

class SmpSolver
{
public:
	SmpSolver() {};
	SmpSolver(IloEnv env,
		std::shared_ptr<Graph>g_ptr,
		SmpForm formulation,
		double epsilon,
		int time_limit,
		int max_cuts);

	void update_problem(const unordered_map<NODE, double> &obj_coeff);

	void add_constraint();

	void solve();
	void solveLP_Steiner();

private:
	void build_problem_scf();
	void build_problem_mcf();
	void build_problem_steiner();
	void build_problem_ns();

	void printInfo(IloNumVar var);

	/* All used variable */
	unordered_map<NODE, IloNumVar> primal_node_vars;					  // x_i
	unordered_map<pair<NODE, INDEX>, IloNumVar> partition_node_vars;	  // x_i^k
	IloNumVarArray x_vararray;

	/* SCF */
	unordered_map<INDEX, IloNumVar> source_node_vars;					  // x_s_k
	unordered_map<pair<NODE_PAIR, INDEX>, IloNumVar> partition_flow_vars; // y_ij^k

	/* MCF */
	unordered_map<pair<NODE_PAIR, NODE_PAIR>, IloNumVar> multi_flow_vars; //y_ij_km

	/* STRINER */
	unordered_map<pair<NODE_PAIR, INDEX>, IloNumVar> edge_vars;			  //y_ij_k
	unordered_map <pair< NODE_PAIR, INDEX>, int> x_varindex;
	unordered_map<INDEX, NODE> Steiner_root;

	IloModel model;
	IloCplex cplex;
	IloObjective objective;

	int time_limit;
	int ncuts;
	int max_cuts;
	double tol;
	SmpForm formulation;
	std::shared_ptr<Graph>	G;
	double elapsed_time;
	double elapsed_ticks;
};
