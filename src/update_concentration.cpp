#include "lung_model_discrete_branching.hpp"
//Checked and commented 16/05/2018

int update_conc(Tree &stree, vector<Tree> &ltree, double time, Options &o) //update concentration on all trees
{
	int error = 1;
	
	//----Assign old values of concentation-----//
	stree.masstotold = stree.masstot;
	for (unsigned i = 0; i < stree.st.size(); i++)
	{
		stree.st[i].totIGvolold = stree.st[i].totIGvol;
		for (unsigned j = stree.st[i].StartGen; j <= stree.st[i].EndGen; j++)
		{
			for (unsigned k = 0; k < stree.st[i].gn[j].p.size(); k++)
			{
				stree.st[i].gn[j].p[k].cold = stree.st[i].gn[j].p[k].c;
			}
		}
	}
	error = update_c_stree(stree, time, o);   //update symmetric case
	
	for (unsigned n = 0; n < o.Ntrees; n++)
	{
		ltree[n].masstotold = ltree[n].masstot;
		for (unsigned i = 0; i < ltree[n].st.size(); i++)
		{
			ltree[n].st[i].totIGvolold = ltree[n].st[i].totIGvol;
			for (unsigned j = ltree[n].st[i].StartGen; j <= ltree[n].st[i].EndGen; j++)
			{
				for (unsigned k = 0; k < ltree[n].st[i].gn[j].p.size(); k++)
				{
					ltree[n].st[i].gn[j].p[k].cold = ltree[n].st[i].gn[j].p[k].c;
				}
			}
		}
		error += update_c_ltree(ltree[n], stree, time, o);   //update symmetric case
	}
	
	return error;
}


int update_c_stree(Tree &stree, double time, Options &o) //update concentration on stree
{
	unsigned long Ntot = o.kmtot;  //total number of points
	
	Eigen::SparseMatrix<double, Eigen::RowMajor> ACrow(Ntot,Ntot);
	Eigen::VectorXd BC(Ntot), XC(Ntot);
	
	fill_Ab(stree,ACrow,XC,BC,Ntot,time,o);
	if(o.SolverOp == DIRECT)
	{
		Eigen::SparseMatrix<double, Eigen::ColMajor> AC = ACrow;
		AC.makeCompressed();
		if(time == 0) stree.solver_dir->analyzePattern(AC);
		if (LU_solver(stree.solver_dir, AC, XC, BC)) return 1;
	}
	else
	{
		ACrow.makeCompressed();
		if(time == 0) stree.solver_iter->analyzePattern(ACrow);
		if (iterative_solver(stree.solver_iter, ACrow, XC, BC)) return 1;
	}

	calc_gas_concs(stree,XC,o);
	
	return 0;
}

int update_c_ltree(Tree &ltree, Tree &stree, double time, Options &o) //update concentration on ltree
{
	//----Concentration update in linear pert calc----//
	unsigned long Ntot = ltree.kmtot;  //total number of points
	
	if (o.InputOp == VOLUME || o.VDM == 0)  //as long as dead-space fills up as before (or isn't there)
	{
		//matrices for update
		Eigen::SparseMatrix<double, Eigen::RowMajor> ACrow(Ntot,Ntot);
		Eigen::VectorXd BC(Ntot), XC(Ntot);
		
		fill_Ab_lp(ltree,stree,ACrow,XC,BC,Ntot,o);
		if(o.SolverOp == DIRECT)
		{
			Eigen::SparseMatrix<double, Eigen::ColMajor> AC = ACrow;
			AC.makeCompressed();
			if(time == 0) ltree.solver_dir->analyzePattern(AC);
			if (LU_solver(ltree.solver_dir, AC, XC, BC)) return 1;
		}
		else
		{
			ACrow.makeCompressed();
			if(time == 0) ltree.solver_iter->analyzePattern(ACrow);
			if (iterative_solver(ltree.solver_iter, ACrow, XC, BC)) return 1;
		}
		
		calc_gas_concs_lp(ltree,stree,XC,o);
	}
	else
	{
		cout << "Can't currently do linear perturbations with pressure input and mouth dead-space, discontinuous boundary conditions cause problems\n";
		cerr << "Pressure input and Mouth DS currently not compatible\n";
		return 1;
	}
	
	return 0;
}


