#include "smp.h"
#include "graph.h"
#include "separation.h"
//#include "type.h"
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <numeric>
#define LOG if(false) cerr

using namespace std;
using namespace lemon;

// Lazy constraint: called only when integer feasible incumbent is found

class StrongComponentLazyCallbackI : public IloCplex::LazyConstraintCallbackI
{   
	std::shared_ptr<Graph> G;
	//unordered_map<pair<NODE_PAIR, INDEX>, IloNumVar> edge_vars;
	unordered_map<INDEX, NODE> root;
	unordered_map<NODE, IloNumVar> primal_node_vars;

	IloNumVarArray x_vararray;
	//unordered_map<pair<NODE_PAIR, INDEX>, int> x_varindex_Steiner;

	const double tol;
	const int max_cuts;
	const SmpForm form;

public:
	/*ILOCOMMONCALLBACKSTUFF(StrongComponentLazyCallback)
		StrongComponentLazyCallbackI(IloEnv env, std::shared_ptr<Graph>graph, unordered_map<pair<NODE_PAIR, INDEX>, IloNumVar>edge_vars_,
			IloNumVarArray x_vararray_, unordered_map<pair<NODE_PAIR, INDEX>, int>x_varindex_Steiner_, double tol_, int max_cuts_, SmpForm form_,
			unordered_map<INDEX, NODE>root_, unordered_map<NODE, IloNumVar>primal_node_vars_)
		: IloCplex::LazyConstraintCallbackI(env), G(graph), edge_vars(edge_vars_), x_vararray(x_vararray_), x_varindex_Steiner(x_varindex_Steiner_),
		tol(tol_), max_cuts(max_cuts_), form(form_), root(root_), primal_node_vars(primal_node_vars_) {}*/

	ILOCOMMONCALLBACKSTUFF(StrongComponentLazyCallback)
		StrongComponentLazyCallbackI(IloEnv env, std::shared_ptr<Graph>graph, IloNumVarArray x_vararray_, 
			double tol_, int max_cuts_, SmpForm form_)
		: IloCplex::LazyConstraintCallbackI(env), G(graph), x_vararray(x_vararray_), tol(tol_), max_cuts(max_cuts_), form(form_) {}

	void main();
};

/*
IloCplex::Callback StrongComponentLazyCallback(IloEnv env, std::shared_ptr<Graph>graph, unordered_map<pair<NODE_PAIR, INDEX>, IloNumVar>edge_vars,
	IloNumVarArray x_vararray, unordered_map<pair<NODE_PAIR, INDEX>, int>x_varindex_Steiner, double tol, int max_cuts, SmpForm form,
	unordered_map<INDEX, NODE>root, unordered_map<NODE, IloNumVar> primal_node_vars)
{
	return (IloCplex::Callback(new (env) StrongComponentLazyCallbackI(env, graph, edge_vars, x_vararray, x_varindex_Steiner,
		tol, max_cuts, form, root, primal_node_vars)));
}
*/

IloCplex::Callback StrongComponentLazyCallback(IloEnv env, std::shared_ptr<Graph>graph, IloNumVarArray x_vararray, double tol ,int max_cuts,SmpForm form)
{
	return (IloCplex::Callback(new (env) StrongComponentLazyCallbackI(env, graph, x_vararray, tol, max_cuts, form)));
}


void StrongComponentLazyCallbackI::main()
{
	/*IloEnv masterEnv = getEnv();
	LOG << "--STRONGCOMPONENT-LAZYCONSTR--" << endl;

	IloNumArray val = IloNumArray(masterEnv, edge_vars.size());
	getValues(val, x_vararray);

	pair<NODE_PAIR, INDEX> pair_ij_k;
	SUB_Graph subG;
	unordered_map<pair<NODE_PAIR, INDEX>, double>xSol;
	for (auto k : G->p_set())
	{
		subG = G->get_subgraph()[k];
		pair_ij_k.second = k;
		for (auto arc : subG.arcs())
		{
			pair_ij_k.first.first = arc.first;
			pair_ij_k.first.second = arc.second;
			xSol[pair_ij_k] = val[x_varindex_Steiner[pair_ij_k]];
		}
	}

	vector<IloExpr> cutLhs, cutRhs;
	vector<double> violation;
	vector<IloRange> cons;
	switch (form)
	{
	case STEINER:
	{
		separate_sc_Steiner(masterEnv, xSol, G, edge_vars, cutLhs, cutRhs, violation);

		// Only need to get the max_cuts maximally-violated inequalities
		vector<int> p(violation.size()); / * vector with indices * /
		iota(p.begin(), p.end(), 0);     / * increasing * /
		bool sorted = false;

		int attempts = 0;
		if (max_cuts < 0)
			attempts = violation.size();
		else
		{
			attempts = min(max_cuts, int(violation.size()));
			partial_sort(p.begin(), p.begin() + attempts, p.end(), [&](int i, int j)
			{ return violation[i] > violation[j]; });/ * sort indices according to violation * /
			sorted = true;
		}

		for (unsigned int i = 0; i < attempts; ++i)
		{
			LOG << violation[p[i]] << endl;
			if (violation[p[i]] >= tol)
			{
				LOG << "Adding user cut for the " << i + 1 << "-th maximally violated constraint. Violation: "
					<< violation[p[i]] << endl;
				try
				{
					LOG << (cutLhs[p[i]] >= cutRhs[p[i]]) << endl;
					add(cutLhs[p[i]] >= cutRhs[p[i]]);
				}
				catch (IloException e)
				{
					cerr << "Cannot add cut" << endl;
				}
			}
			else / * sorted, so no further violated ineq exist * /
				if (sorted)
					break;
		}
		for (unsigned int i = 0; i < cutLhs.size(); ++i)
		{
			cutLhs[i].end();
			cutRhs[i].end();
		}
	}
	break;
	}

	LOG << "--END ILOLAZYCONSTRAINTCALLBACK--" << endl;
	return;*/
}

/*

// User cut: called also on fractional solutions
class SmpCutCallbackI : public IloCplex::UserCutCallbackI
{
	std::shared_ptr<Graph> G;
	unordered_map<pair<NODE_PAIR, INDEX>, IloNumVar> edge_vars;
	unordered_map<INDEX, NODE> root;
	unordered_map<NODE, IloNumVar> primal_node_vars;

	IloNumVarArray x_vararray;
	unordered_map<pair<NODE_PAIR, INDEX>, int> x_varindex_Steiner;

	const double tol;
	const int max_cuts;
	const SmpForm form;

public:
	ILOCOMMONCALLBACKSTUFF(SmpCutCallback)
		SmpCutCallbackI(IloEnv env, std::shared_ptr<Graph>graph, unordered_map<pair<NODE_PAIR, INDEX>, IloNumVar>edge_vars_,
			IloNumVarArray x_vararray_, unordered_map<pair<NODE_PAIR, INDEX>, int>x_varindex_Steiner_, double tol_, int max_cuts_, SmpForm form_,
			unordered_map<INDEX, NODE>root_, unordered_map<NODE, IloNumVar>primal_node_vars_)
		: IloCplex::UserCutCallbackI(env), G(graph), edge_vars(edge_vars_), x_vararray(x_vararray_), x_varindex_Steiner(x_varindex_Steiner_),
		tol(tol_), max_cuts(max_cuts_), form(form_), root(root_), primal_node_vars(primal_node_vars_) {}

	void main();
};

IloCplex::Callback SmpCutCallback(IloEnv env, std::shared_ptr<Graph>graph, unordered_map<pair<NODE_PAIR, INDEX>, IloNumVar>edge_vars,
	IloNumVarArray x_vararray, unordered_map<pair<NODE_PAIR, INDEX>, int>x_varindex_Steiner, double tol, int max_cuts, SmpForm form,
	unordered_map<INDEX, NODE>root, unordered_map<NODE, IloNumVar> primal_node_vars)
{
	return (IloCplex::Callback(new (env) SmpCutCallbackI(env, graph, edge_vars, x_vararray, x_varindex_Steiner,
		tol, max_cuts, form, root, primal_node_vars)));
}

void SmpCutCallbackI::main()
{
	// Skip the separation if not at the end of the cut loop
	if (!isAfterCutLoop())
		return;

	LOG << "--SMP USERCUT--" << endl;

	IloEnv masterEnv = getEnv();
	IloNumArray val = IloNumArray(masterEnv, edge_vars.size());
	getValues(val, x_vararray);

	pair<NODE_PAIR, INDEX> pair_ij_k;
	SUB_Graph subG;
	unordered_map<pair<NODE_PAIR, INDEX>, double>xSol;
	for (auto k : G->p_set())
	{
		subG = G->get_subgraph()[k];
		pair_ij_k.second = k;
		for (auto arc : subG.arcs())
		{
			pair_ij_k.first.first = arc.first;
			pair_ij_k.first.second = arc.second;
			xSol[pair_ij_k] = val[x_varindex_Steiner[pair_ij_k]];
		}
	}

	vector<IloExpr> cutLhs, cutRhs;
	vector<double> violation;
	vector<IloRange> cons;
	switch (form)
	{
	case STEINER:
	{
		if (!separate_sc_Steiner(masterEnv, xSol, G, edge_vars, cutLhs, cutRhs, violation))
			seperate_min_cut_Steiner(masterEnv, xSol, G, edge_vars, cutLhs, cutRhs, violation, root, primal_node_vars);

		// Only need to get the max_cuts maximally-violated inequalities
		vector<int> p(violation.size()); / * vector with indices * /
		iota(p.begin(), p.end(), 0);     / * increasing * /
		bool sorted = false;

		int attempts = 0;
		if (max_cuts < 0)
			attempts = violation.size();
		else
		{
			attempts = min(max_cuts, int(violation.size()));
			partial_sort(p.begin(), p.begin() + attempts, p.end(), [&](int i, int j)
			{ return violation[i] > violation[j]; });/ * sort indices according to violation * /
			sorted = true;
		}

		for (unsigned int i = 0; i < attempts; ++i)
		{
			LOG << violation[p[i]] << endl;
			if (violation[p[i]] >= tol)
			{
				LOG << "Adding user cut for the " << i + 1 << "-th maximally violated constraint. Violation: "
					<< violation[p[i]] << endl;
				try
				{
					LOG << (cutLhs[p[i]] >= cutRhs[p[i]]) << endl;
					add(cutLhs[p[i]] >= cutRhs[p[i]]);
				}
				catch (IloException e)
				{
					cerr << "Cannot add cut" << endl;
				}
			}
			else / * sorted, so no further violated ineq exist * /
				if (sorted)
					break;
		}
		for (unsigned int i = 0; i < cutLhs.size(); ++i)
		{
			cutLhs[i].end();
			cutRhs[i].end();
		}
	}
	break;
	}

	LOG << "---END ILOUSERCUTCALLBACK---" << endl;
	return;
}*/