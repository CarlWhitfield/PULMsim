#include "lung_model_discrete_branching.hpp"

void subtree_resistance(Tree &tree, unsigned i, Options &o)
{
	//calculate resistance of subtree
	double aj, lj;
	unsigned j, jn;
	unsigned imp = tree.st[i].imeanpath;
	tree.st[i].Rtree = 0;

	if (tree.st[i].EndGen > tree.st[imp].Ncond + tree.st[imp].StartGen) jn = tree.st[imp].Ncond + tree.st[imp].StartGen;
	else jn = tree.st[i].EndGen;
	if (tree.st[i].blocked == false)   //only calculate if branch isn't 100% blocked
	{
		double qmag = fabs((tree.st[i].Vnew - tree.st[i].Vold) / o.dt);
		for (j = tree.st[i].StartGen; j <= jn; j++)
		{
			aj = tree.st[i].gn[j].p[0].al;
			if(j>0) lj = tree.st[i].gn[j].p.size()*tree.st[i].gn[j].dx;
			else lj = (tree.st[i].gn[j].p.size()-1)*tree.st[i].gn[j].dx;
			double Rpois = 8.0*lj*o.Viscosity*M_PI / (tree.st[i].gn[j].Nb*aj*aj);   //poiseuille resistance of branch
			if (o.FlowType == PEDLEY)    //pedley resistance if used
			{
				double Z = (0.5*PEDLEY_C)*sqrt(qmag*o.Density / (2.0*tree.st[i].gn[j].Nb*M_PI*o.Viscosity*lj));
				if (Z > 1) tree.st[i].Rtree += Z*Rpois;
				else tree.st[i].Rtree += Rpois;
			}
			else tree.st[i].Rtree += Rpois;  //if just poiseuille
		}
	}
}

void subtree_lp_symm_resistance(subtree &st, subtree &symm, Options &o)
{
	double aj, lj;
	unsigned j, jn;
	//calculate 0th order resistance of this branch in ltree (not perturbed resistance)
	st.Rtree = 0;

	if (st.EndGen > st.Ncond + st.StartGen) jn = st.Ncond + st.StartGen;
	else jn = st.EndGen;
	if (st.blocked == false)   //don't use Rtree otherwise
	{
		double qmag = pow(2.0, symm.StartGen - st.StartGen)*fabs((symm.Vnew - symm.Vold) / o.dt);  //equivalent qmag here
		for (j = st.StartGen; j <= jn; j++)
		{
			//get properties from symm tree
			aj = symm.gn[j].p[0].al;
			if (j == 0) lj = (symm.gn[j].p.size() - 1)*symm.gn[j].dx;
			else  lj = symm.gn[j].p.size()*symm.gn[j].dx;

			double Rpois = 8.0*lj*o.Viscosity*M_PI / (st.gn[j].Nb*aj*aj);
			if (o.FlowType == PEDLEY)
			{
				double Z = (0.5*PEDLEY_C)*sqrt(qmag*o.Density / (2.0*st.gn[j].Nb*M_PI*o.Viscosity*lj));
				if (Z > 1) st.Rtree += 1.5*Z*Rpois;     //in linear case resistance goes up by factor 1.5 for pedley branches
				else st.Rtree += Rpois;
			}
			else st.Rtree += Rpois;  //if it is first gen, add viscosity pert too
		}
	}
}

int full_tree_fluxcalc(Tree &tree, double time, Options &o)
{
	if (o.RespFunc == LINEAR_RESP)
	{
		unsigned NT = ((unsigned) tree.EndSubtrees.size());

		int *it = new int[tree.st.size()];     //matrix for re-numbering trees
		Eigen::VectorXd Vprev(NT+1);     //stores volume from previous iteration
		Eigen::VectorXd Vold(NT+1);     //stores volume from previous iteration
		Eigen::VectorXd Bvec = Eigen::VectorXd::Zero(NT+1);    //RHS of equation

		for (unsigned m = 0; m < NT; m++)  //loop over subtrees -- initialise these
		{
			unsigned im = tree.EndSubtrees[m];    //current terminating subtree
			//--Make Vold and S vector--//
			Vold(m) = tree.st[im].Vold;   //stores vold for outlets
		}
		if(o.InputOp == VOLUME)
		{
			Vold(NT) = tree.Ppl*o.dt;   //last entry is - pressure (old) * dt
		}
		else
		{
			Vold(NT) = o.dt*pressure_func(time + 0.5*o.dt, o);   //last entry is - dt * pressure at t + dt/2
		}

		double residual = 1;
		while (residual > PTOL)  //iterate linear problem to solve non-linear problem
		{
			if (o.FlowType != POISEUILLE || o.RespFunc != LINEAR_RESP)   //iterative process if not linear
			{
				for (unsigned i = 0; i < tree.st.size(); i++)   //loop over all subtrees
				{
					subtree_resistance(tree, i, o);  //calculate resistance of subtrees from Vnew
				}
				if (build_resistance_matrix(tree, o)) return 1;    //calculate resistance and invert A matrix
			}
			
			//build b for Ax = b
			Bvec = tree.Bmat*Vold + o.dt*tree.KVs;
			if(o.InputOp == VOLUME) Bvec(NT) += 0.5*o.dt*pressure_func(time + 0.5*o.dt, o);  //pressure func returns flow rate

			//solve Ax = b
			Eigen::VectorXd Vnew = tree.AmatLU->solve(Bvec);
			
			//if not using poiseuille check how close this is to actual solution
			residual = 0;
			for (unsigned m = 0; m < NT; m++)  //loop over subtrees
			{
				unsigned im = tree.EndSubtrees[m];    //current terminating subtree
				Vprev(m) = tree.st[im].Vnew;
				tree.st[im].Vnew = Vnew(m);    //final expression for Vnew
				if (o.FlowType != POISEUILLE) residual += fabs((tree.st[im].Vnew - Vprev(m)) / tree.Vs(m));   //difference between solution and
			}
			tree.Ppl = Vnew(NT)/o.dt;  //pleural pressure

			for (unsigned i = 0; i < tree.st.size(); i++)  //loop over subtrees
			{
				double Vtot = 0;
				for (unsigned m = 0; m < tree.st[i].EndSubtrees.size(); m++)
				{
					Vtot += tree.st[tree.st[i].EndSubtrees[m]].Vnew;   //parent branches are sum of subtree volumes
				}
				tree.st[i].Vnew = Vtot;
			}
		}

		delete[] it;

		return 0;
	}
	else
	{
		cout << "Not coded non-linear resp-function yet.\n";
		return 1;
	}
}

int full_lp_tree_fluxcalc(Tree &ltree, Tree &stree, Options &o)
{
	//linearly perturbed version of full_tree_fluxcalc
	if (o.RespFunc == LINEAR_RESP)
	{
		unsigned NT = ((unsigned) ltree.EndSubtrees.size());     //number of terminating subtrees
		double sVnew=0, sVold=0, sV0=0;
		
		Eigen::VectorXd DVold = Eigen::VectorXd::Zero(NT+1);   //stores old volume values
		Eigen::VectorXd Bvec = Eigen::VectorXd::Zero(NT+1);
		
		unsigned ip = ltree.pert_ij[0];    //ltree subtree containing pert
		unsigned isp = ltree.st[ip].i_stree;     //stree subtree containing pert
		unsigned mp = 0;   //matrix index of subtree containing pert
		
		for (unsigned m = 0; m < NT; m++)  //loop over subtrees -- initialise these
		{
			unsigned im = ltree.EndSubtrees[m];    //current terminating subtree
			if (im == ip)
			{
				mp = m;
				//unperturbed quantities
				sVnew = stree.st[ltree.st[im].i_stree].Vnew / ((double) stree.st[ltree.st[im].i_stree].gn[ltree.st[im].StartGen].Nb);
				sVold = stree.st[ltree.st[im].i_stree].Vold / ((double) stree.st[ltree.st[im].i_stree].gn[ltree.st[im].StartGen].Nb);
				sV0 = stree.st[ltree.st[im].i_stree].V0 / ((double) stree.st[ltree.st[im].i_stree].gn[ltree.st[im].StartGen].Nb);
			}
			//--Make Vold and S vector--//
			DVold(m) = ltree.st[im].Vold;   //stores vold for outlets
		}
		DVold(NT) = -ltree.Ppl*o.dt;
		
		if (o.FlowType != POISEUILLE || o.RespFunc != LINEAR_RESP)   //iterative process if not linear
		{
			for (unsigned i = 0; i < ltree.st.size(); i++)   //loop over all subtrees
			{
				if(ltree.st[i].Ntot >= o.Ngen2 - o.Ngen) subtree_lp_symm_resistance(ltree.st[i], stree.st[ltree.st[i].i_stree], o);  //calculate resistance of symm subtrees in perturbed network
			}
			if (build_resistance_matrix(ltree, o)) return 1;    //calculate resistance and invert A matrix
		}
		
		double dr = 0, dK = 0;  //resistance and K perts
		unsigned jp = ltree.pert_ij[1];             //generation perturbed
		if(ltree.st[ip].Ntot + o.Ngen >= o.Ngen2)   //only if pert is in conducting tree
		{
			if (ltree.pert_type == 'K')
			{
				double KNk = ltree.st[ip].E;
				dK = KNk*ltree.st[ip].Ep;       //actual change in stiffness
			}
			else
			{
				double aj = stree.st[isp].gn[jp].p[0].al;     //original area
				double lj = stree.st[isp].gn[jp].dx*stree.st[isp].gn[jp].p.size();    //original length
				double rjk = 8.0*lj*o.Viscosity*M_PI / (ltree.st[ip].gn[jp].Nb*aj*aj);
				double qmag = (sVnew - sVold) / o.dt;    //velocity magnitude in perturbed tree
				double Z = (0.5*PEDLEY_C)*sqrt(qmag*o.Density / (2.0*M_PI*o.Viscosity*lj));
				if (o.FlowType == PEDLEY && Z > 1) dr = Z*rjk*(0.5*ltree.st[ip].gn[jp].Lp - 2.0*ltree.st[ip].gn[jp].Ap);  //pedley pert
				else dr = rjk*(ltree.st[ip].gn[jp].Lp - 2.0*ltree.st[ip].gn[jp].Ap);    //poiseuille pert
			}
		}
		
		//build b vector for Ax = b
		Bvec = ltree.Bmat*DVold;
		Bvec(mp) += -dr*(sVnew - sVold) - 0.5*o.dt*dK*(sVnew + sVold - 2.0*sV0);

		//solve Ax = b
		Eigen::VectorXd DVnew = ltree.AmatLU->solve(Bvec);

		for (unsigned m = 0; m < NT; m++)  //loop over subtrees
		{
			unsigned im = ltree.EndSubtrees[m];    //current terminating subtree
			ltree.st[im].Vnew = DVnew(m);
		}
		ltree.Ppl = DVnew(NT)/o.dt;   //otherwise pressure is the input
		
		for (unsigned i = 0; i < ltree.st.size(); i++)  //loop over subtrees
		{
			double Vtot = 0;
			for (unsigned m = 0; m < ltree.st[i].EndSubtrees.size(); m++)
			{
				Vtot += ltree.st[ltree.st[i].EndSubtrees[m]].Vnew;   //parent branches are sum of subtree volumes
			}
			ltree.st[i].Vnew = Vtot;
		}
		
		return 0;
	}
	else
	{
		cout << "Not coded non-linear resp-function yet.\n";
		return 1;
	}
}

int build_resistance_matrix(Tree &tree, Options &o)
{
	//build resistance matrix for given tree structure
	unsigned NT = ((unsigned)tree.EndSubtrees.size());
	Eigen::MatrixXd RM = Eigen::MatrixXd::Zero(NT,NT);
	Eigen::VectorXd One_vec = Eigen::VectorXd::Ones(NT);

	//Re-zero the matrices
	tree.Amat = Eigen::MatrixXd::Zero(NT+1,NT+1);   //extra entry for pressure
	tree.Bmat = Eigen::MatrixXd::Zero(NT+1,NT+1);
	tree.KVs = Eigen::VectorXd::Zero(NT+1);    //last entry zero
	tree.Kmat = Eigen::MatrixXd::Zero(NT,NT);   //contains K values
	tree.Vs = Eigen::VectorXd::Zero(NT);  //this just contains volumes
	
	for (unsigned m = 0; m < NT; m++)  //loop over end subtrees
	{
		unsigned im = tree.EndSubtrees[m];    //current terminating subtree
		tree.Kmat(m,m) = tree.st[im].E;     //make vector of stiffnesses
		tree.Vs(m) = tree.st[im].V0;      //make vector of resting volumes
	}
	tree.KVs.block(0,0,NT,1) = tree.Kmat*tree.Vs;

	for (unsigned i = 0; i < tree.st.size(); i++)   //loop over all subtrees
	{
		if(tree.st[i].Ntot >= o.Ngen2- o.Ngen)  //not in acinus
		{
			for (unsigned m = 0; m < tree.st[i].EndSubtrees.size(); m++)   //subtree i contributes to entries for all terminating subtrees descended from it
			{
				unsigned im = tree.st[i].EndSubtrees[m];
				unsigned ESTno_x = tree.STno_to_ESTno[im];
				for (unsigned n = 0; n < tree.st[i].EndSubtrees.size(); n++)   
				{
					unsigned in = tree.st[i].EndSubtrees[n];
					unsigned ESTno_y = tree.STno_to_ESTno[in];
					RM(ESTno_x, ESTno_y) +=  tree.st[i].Rtree; //add resistance from parent branch to related outflows (normalised by portion of tree)
					if (i == 0) RM(ESTno_x, ESTno_y) += o.Rmouth;  //add mouth resistance
				}
				if(i==im) RM(ESTno_x, ESTno_x) += tree.st[im].Rb;  //Add diagonal resistance (bag resistance) for terminal trees (when i is end sub tree im)
			}
		}
	}
	
	Eigen::VectorXd TestSol = RM.partialPivLu().solve(One_vec);
	tree.Amat.block(0,0,NT,NT) = RM + 0.5*o.dt*tree.Kmat;
	tree.Bmat.block(0,0,NT,NT) = RM - 0.5*o.dt*tree.Kmat;
	if(o.InputOp == PRESSURE)   //last row is just P = b[NT]  -- pressure fixed
	{
		tree.Amat(NT,NT) = 1;
		tree.Bmat(NT,NT) = 1;
	}
	else   //last row is 1^T . Vnew = b[NT]  -- flux fixed
	{
		tree.Amat.block(0,NT,NT,1) = 0.5*One_vec;
		tree.Amat.block(NT,0,1,NT) = 0.5*One_vec.transpose();   //keeps matrix symmetric 1/2 dv/dt = fixed
		tree.Bmat.block(0,NT,NT,1) = -0.5*One_vec;
		tree.Bmat.block(NT,0,1,NT) = 0.5*One_vec.transpose();
	}
	
	//	cout << "K\n" << tree.Kmat << "\nKVs\n" << tree.KVs << "\nRes mat\n" << RM << "\nA mat\n" << tree.Amat << "\nB mat\n" << tree.Bmat << '\n';
	//tree.AmatQR->compute(tree.Amat);
	tree.AmatLU->compute(tree.Amat);

	tree.R0 = 1.0/(One_vec.transpose()*TestSol);     //effective resistance

	return 0;
}

bool isparent(Tree &tree, unsigned ip, unsigned id)  //returns true if subtree id is daughter of subtree ip
{
    bool parent = false;
    int ih = id;
    while(parent == false && ih >= 0)
    {
        if(ih == ip) parent = true;
        ih = tree.st[ih].treein;
    }
    
    return parent;
}
