#include "lung_model_discrete_branching.hpp"
//--Checked and commented 16/05/18---//

//----Discretisation functions----//

void upwind_coeffs(vector<double> ucfl[3], vector<double> ucfr[3], vector<double> dxh[3], double u0, Options &o)
{
	//---Works out upwind interpolation from nearest neighbours---//
	//stores upwind coefficients of neighbours at given node
	unsigned mmax = ((unsigned) dxh[2].size());
	switch(o.UpwindOp)
	{
		case FIRST_ORDER_UPWIND:
		{
			for(unsigned n=0; n<3; n++)   //to be calculated
			{
				ucfl[n].clear();
				ucfr[n].clear();
				for(unsigned m=0; m<mmax; m++)
				{
					if(n==0)
					{
						if(u0>0) ucfl[n].push_back(1.0);
						else ucfl[n].push_back(0.0);
						ucfr[n].push_back(0.0);
					}
					if(n==1)
					{
						if(u0<=0)
						{
							ucfl[n].push_back(1.0);
							ucfr[n].push_back(0.0);
						}
						else
						{
							ucfl[n].push_back(0.0);
							ucfr[n].push_back(1.0);
						}
					}
					if(n==2)
					{
						ucfl[n].push_back(0.0);
						if(u0<=0) ucfr[n].push_back(1.0);
						else ucfr[n].push_back(0.0);
					}
				}
			}
		} break;
		
		case CENTRAL_UPWIND:
		{
			for(unsigned n=0; n<3; n++)   //to be calculated
			{
				ucfl[n].clear();
				ucfr[n].clear();
				for(unsigned m=0; m<mmax; m++)
				{
					if(n==0)
					{
						ucfl[n].push_back(dxh[0][0]/(dxh[1][0]+dxh[0][0]));
						ucfr[n].push_back(0.0);
					}
					if(n==1)
					{
						ucfl[n].push_back(dxh[1][0]/(dxh[1][0]+dxh[0][0]));
						ucfr[n].push_back(dxh[2][m]/(dxh[1][0]+dxh[2][m]));
					}
					if(n==2)
					{
						ucfl[n].push_back(0.0);
						ucfr[n].push_back(dxh[1][0]/(dxh[1][0]+dxh[2][m]));
					}
				}
			}
		} break;
		
		default: cout << "Something went wrong with Upwind Option.\n";
	}
}

void lp_upwind_coeffs(vector<double> ducfl[3], vector<double> ducfr[3], vector<double> dx0[3], vector<double> dxh[3], double u0, Options &o)
{
	unsigned n,m;
	
	//-----Linear perturbations of upwind interpolation due to length change -----//
	
	//ducfl contains linear pert to ucfl etc.
	//dx0 contains unperturbed node length, dxh contains linear perturbation
	
	switch(o.UpwindOp)
	{
		case CENTRAL_UPWIND:
		{
			for(n=0; n<3; n++)
			{
				ducfl[n].clear();    //initialise
				ducfr[n].clear();
				for(m=0; m<dxh[2].size(); m++)   //there can be more than one node at k+1, so we loop over these
				{
					if(n==0)   //central interpolation only changes at boundary between generation
					{
						ducfl[n].push_back((dxh[1][0]*dx0[0][0] - dx0[1][0]*dxh[0][0])/((dx0[0][0]+dx0[1][0])*(dx0[0][0]+dx0[1][0])));
						ducfr[n].push_back(0.0);
					}
					if(n==1)
					{
						ducfl[n].push_back((dx0[1][0]*dxh[0][0] - dxh[1][0]*dx0[0][0])/((dx0[0][0]+dx0[1][0])*(dx0[0][0]+dx0[1][0])));
						ducfr[n].push_back((dxh[2][m]*dx0[1][0] - dx0[2][0]*dxh[1][0])/((dx0[1][0]+dx0[2][0])*(dx0[1][0]+dx0[2][0])));
					}
					if(n==2)
					{
						ducfr[n].push_back((dx0[2][0]*dxh[1][0] - dxh[2][m]*dx0[1][0])/((dx0[1][0]+dx0[2][0])*(dx0[1][0]+dx0[2][0])));
						ducfl[n].push_back(0.0);
					}
				}
			}
		} break;
		
		case FIRST_ORDER_UPWIND:
		{
			for(n=0; n<3; n++)
			{
				ducfl[n].clear();     //initialise
				ducfr[n].clear();
				for(m=0; m<dxh[2].size(); m++)    //upwind is unchanged to first order
				{
					ducfl[n].push_back(0.0);
					ducfr[n].push_back(0.0);
				}
			}
		} break;
	}
}

void fd_coeffs(vector<double> dcfl[3], vector<double> dcfr[3], vector<double> dxh[3])
{
	unsigned mmax = ((unsigned) dxh[2].size());
	//---Works out central gradient from nearest neighbours---//
	//stores central difference coefficients of neighbours at given node
	for(unsigned n=0; n<3; n++)
	{
		dcfl[n].clear();
		dcfr[n].clear();
		for(unsigned m=0; m<mmax; m++)
		{
			if(n==0)
			{
				dcfl[n].push_back(-2.0/(dxh[0][0]+dxh[1][0]));
				dcfr[n].push_back(0.0);
			}
			if(n==1)
			{
				dcfl[n].push_back(2.0/(dxh[0][0]+dxh[1][0]));
				dcfr[n].push_back(-2.0/(dxh[1][0]+dxh[2][m]));
			}
			if(n==2)
			{
				dcfl[n].push_back(0.0);
				dcfr[n].push_back(2.0/(dxh[1][0]+dxh[2][m]));
			}
		}
	}
}

void lp_fd_coeffs(vector<double> ddcfl[3], vector<double> ddcfr[3], vector<double> dx0[3], vector<double> dxh[3])
{
	//linear perturbation to finite difference coeffs
	unsigned n,m;
	
	for(n=0; n<3; n++)
	{
		ddcfl[n].clear();
		ddcfr[n].clear();
		for(m=0; m<dxh[2].size(); m++)
		{
			if(n==0)
			{
				ddcfl[n].push_back(2.0*(dxh[1][0] + dxh[0][0])/((dx0[0][0]+dx0[1][0])*(dx0[0][0]+dx0[1][0])));
				ddcfr[n].push_back(0.0);
			}
			if(n==1)
			{
				ddcfl[n].push_back(-2.0*(dxh[0][0] + dxh[1][0])/((dx0[0][0]+dx0[1][0])*(dx0[0][0]+dx0[1][0])));
				ddcfr[n].push_back(2.0*(dxh[2][m] + dxh[1][0])/((dx0[1][0]+dx0[2][0])*(dx0[1][0]+dx0[2][0])));
			}
			if(n==2)
			{
				ddcfl[n].push_back(0.0);
				ddcfr[n].push_back(-2.0*(dxh[1][0] + dxh[2][m])/((dx0[1][0]+dx0[2][0])*(dx0[1][0]+dx0[2][0])));
			}
			
		}
	}
}

int iterative_solver(Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::IncompleteLUT<double>> *solver, Eigen::SparseMatrix<double, Eigen::RowMajor> &AC, Eigen::VectorXd &x, Eigen::VectorXd &b)
{
	//uses Eigen libraries to iteratively solve Ax = b for sparse matrix
	// Compute the numerical factorization
	solver->factorize(AC);
	if(solver->info()!=Eigen::Success)
	{
		cout << "Error in factorise.\n";
		return 1;
	}
	//Use the factors to solve the linear system
	x = solver->solveWithGuess(b, x);
	if(solver->info() != Eigen::Success)
	{
		cout << "Error in solve.\n";
		return 1;
	}
	return 0;
}

int LU_solver(Eigen::SparseLU<Eigen::SparseMatrix<double>> *solver, Eigen::SparseMatrix<double> &AC, Eigen::VectorXd &x, Eigen::VectorXd &b)
{
	//uses Eigen libraries to iteratively solve Ax = b for sparse matrix
	// Compute the numerical factorization
	solver->factorize(AC);
	if(solver->info()!=Eigen::Success)
	{
		cout << "Error in factorise.\n";
		return 1;
	}
	//Use the factors to solve the linear system
	x = solver->solve(b);
	if(solver->info()!=Eigen::Success)
	{
		cout << "Error in solve.\n";
		return 1;
	}
	return 0;
}

void tree_point_pos(Tree &tree, unsigned i, unsigned j, unsigned k, int n, vector<unsigned> &in, vector<unsigned> &jn, vector<unsigned> &kn)
{
	//--increments to k+n th element in symmetric grid--//
	unsigned Nj;
	vector<int> nn;
	in.clear();    //remove anything stored in these vectors
	jn.clear();
	kn.clear();

	in.push_back(i);
	jn.push_back(j);
	kn.push_back(k);        //initialise starting point
	nn.push_back(n);
	for (unsigned nh = 0; nh < in.size(); nh++)
	{
		Nj = ((unsigned) tree.st[in[nh]].gn[jn[nh]].p.size());
		while (((int) kn[nh]) + nn[nh] > ((int) Nj) - 1 || nn[nh] + ((int)kn[nh]) < 0)   //iterate until at correct point
		{
			if (nn[nh] + ((int)kn[nh]) > ((int)Nj) - 1)            //if kn + n is below this generation
			{
				if (jn[nh] + 1 > tree.st[in[nh]].EndGen)    //if it goes outside tree range, move to next tree
				{
					if (tree.st[in[nh]].treeout[0] < 0)  //if it is outside all trees, start to subtract points back off
					{
						nn[nh] = -(nn[nh] - (((int)Nj) - ((int)kn[nh])));   //n becomes negative
						kn[nh] = Nj - 1;
					}
					else
					{
						nn[nh] += ((int)kn[nh]) - ((int)Nj);   //n reduces by size of points left in this generation
						nn.push_back(nn[nh]);    //add new counter to carry on here
						in.push_back(tree.st[in[nh]].treeout[1]);  //mve to these trees
						in[nh] = tree.st[in[nh]].treeout[0];
						jn[nh] += 1;
						jn.push_back(jn[nh]);
						kn[nh] = 0;
						kn.push_back(0);
					}
				}
				else                   //it is in the tree, move down a gen
				{
					nn[nh] += ((int)kn[nh]) - ((int)Nj);    //n reduces by size of points left in this generation
					jn[nh] += 1;
					kn[nh] = 0;
				}
			}
			else                       //kn+n above this generation
			{
				if (jn[nh] == tree.st[in[nh]].StartGen)           //if point is not in this tree
				{
					if (tree.st[in[nh]].treein >= 0) //if tree exists
					{
						nn[nh] += ((int)kn[nh]) + 1;
						in[nh] = tree.st[in[nh]].treein;
						jn[nh] = tree.st[in[nh]].EndGen;
						kn[nh] = ((unsigned) tree.st[in[nh]].gn[jn[nh]].p.size())-1;  //move to last point in tree in
					}
					else   //just go to zero, must be base tree
					{
						nn[nh] = 0;
						kn[nh] = 0;
					}
				}
				else                   //move up a generation
				{
					nn[nh] += ((int)kn[nh]) + 1;
					jn[nh] -= 1;
					kn[nh] = ((unsigned) tree.st[in[nh]].gn[jn[nh]].p.size()) - 1; //move to last point in gen up
				}
			}
			Nj = ((unsigned) tree.st[in[nh]].gn[jn[nh]].p.size());
		}
		kn[nh] = kn[nh] + nn[nh];                     //final value of k
	}
}