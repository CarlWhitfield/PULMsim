//
//  initialisation.cpp
//
//
//  Created by Carl Whitfield on 09/03/2017.
//
//

//--Checked and commented 16/05/18---//

//-----Initialisation functions-----//

#include "lung_model_discrete_branching.hpp"
int initialise_system(Tree &stree, Conversions &cons, Options &o)
{
	//-- Output file of simulation summary --//
	ofstream summary_file;
	stringstream fnss;

	//build file name from options

	fnss << o.TreeOp % 100 << o.BcOp % 100;
	fnss << o.UpwindOp % 100 << o.FlowType % 100 << o.RespFunc % 100 << o.PressOp % 100 << o.ShapeOp % 100 << o.SolverOp % 100 << o.InitOp % 100;
	o.SimID = fnss.str();  //simulation ID made up of sim options	
	fnss.str("");
	fnss << "summary_" << o.SimID << ".txt";
	summary_file.open(fnss.str().c_str());

	summary_file << "Simulation Options\n";
	if (print_simops(summary_file, o)) return 1;
	if(o.BcOp == BAG) o.Ngen2 = o.Ngen;  //zero gens in acinus
	double LL;
	
	//choose which lung structure to build
	switch (o.ShapeOp)
	{
	case ALT_LOBE_GEOMETRIC: LL = build_alternative_lobe_branches(stree, o);    //returns path length from mouth to end of RL lobe
		break;

	default: LL = build_symmetric_branches(stree, o);    //returns path length from mouth to end
		break;
	}
	
	//print infor to simulation summary file
	summary_file << "\nSimulation Parameters (Physiological units)\n";
	if (print_params(summary_file, o)) return 1;
	
    if(convert_units(stree, cons, LL, o)) return 1;     //normalise units -- lengthscale LL

	summary_file << "\nSimulation Parameters (Simulation units)\n";
	if (print_params(summary_file, o)) return 1;
    
	if(o.def.size() > 0) initialise_defects(stree, o);

	//---Calc resistance and Approximate Peclet number in trachea (not accounting for Taylor disp)-----/
	if (o.InputOp == VOLUME)
	{
		o.Peclet = o.VT*stree.st[0].L0 / (stree.st[0].A0*o.Diffusion);  //assuming flow rate in trachea is ~o.VT/A0
	}
	else
	{
		cout << "Need to work out Peclet for this case\n";
        if(o.PressOp == LINEAR) o.Peclet = o.P0/(o.Diffusion*o.E*stree.st[0].A0);
        else o.Peclet = o.P0/(stree.st[0].A0*o.Diffusion*stree.R0);
	}
	summary_file << "\nPeclet number: " << o.Peclet << '\n';
	
	//---- Initialise the system ----//
	if (initialise_tree(stree, o)) return 1;    //error if returns 1
	//------------------------------//

	//Check total volume is correct (for debugging)
//    double Vcheck = 0, Vcheck2 = 0, Vcheck3 = 0;
//    for (unsigned i = 0; i < stree.st.size(); i++)
//    {
//        for (unsigned j = stree.st[i].StartGen; j <= stree.st[i].EndGen; j++)
//        {
//            unsigned k0;
//            if (j == 0) k0 = 1;
//            else k0 = 0;
//            for (unsigned k = k0; k < stree.st[i].gn[j].p.size(); k++)
//            {
//                if (j<stree.st[i].StartGen + stree.st[i].Ncond) Vcheck += stree.st[i].gn[j].p[k].al*stree.st[i].gn[j].Nb*stree.st[i].gn[j].dx;
//                else Vcheck2 += stree.st[i].gn[j].p[k].al*stree.st[i].gn[j].Nb*stree.st[i].gn[j].dx;
//                Vcheck3 += stree.st[i].gn[j].p[k].Anew*stree.st[i].gn[j].Nb*stree.st[i].gn[j].dx;
//            }
//        }
//    }

	summary_file.close();
	
	return 0;
}

int initialise_defects(Tree &stree, Options &o)
{
	vector<unsigned> stnextgen;
	subtree sth;
	unsigned osts = ((unsigned) stree.st.size());   //original number of subtrees
	
	//--Build base tree--//
	for(unsigned im = 0; im < osts; im++) //loop over mean paths
	{
		stnextgen.push_back(im);
		while (stnextgen.size() > 0)    //while branches left unterminated
		{
			vector<unsigned> stthisgen = stnextgen;  //branches to look at in this gen
			stnextgen.clear();   //reset
			for (unsigned q = 0; q < stthisgen.size(); q++)  //all unterminated subtrees in next gen
			{
				bool symmcheck = true;
				unsigned ic = stthisgen[q];     //current subtree
				//use unordered_map to store defects (vector of defects for each map point)
				//loop over elements of unordered map
				//if defect is in this subtree -> check if all equiv branches have same defect -> otherwise symmcheck = false
                unsigned n = 0;       //counts defects
				while (symmcheck == true && n < o.def[im].size())   //check whether this subtree is symmetric
				{
                    unsigned jc = o.def[im][n].jn + stree.st[im].StartGen - stree.st[ic].StartGen;   //defect gen relative to tree ic startgen
					if (o.def[im][n].jn + stree.st[im].StartGen <= stree.st[ic].StartGen || o.def[im][n].kstart < stree.st[ic].StartBranch*pow(2, jc))
					{
						n++;   //before this subtree
					}
					else
					{
						unsigned jh = o.def[im][n].jn + stree.st[im].StartGen;   //gen of first defect in whole tree
						unsigned long khs = o.def[im][n].kstart;
						if (jh > stree.st[ic].StartGen && khs == stree.st[ic].StartBranch*((unsigned long) pow(2, jc)))
						{
							//defect at left of subtree
							unsigned long k = o.def[im][n].kend + 1;  //kstart until kend
                            if(o.def[im][n].kend >= (stree.st[ic].StartBranch+1)*((unsigned long) pow(2, jc)))
                            {
                                //this stretch of defects ends outside this tree
                                //insert next stretch of defects into defect list
                                Defect overflow;
                                overflow.jn = o.def[im][n].jn;
                                overflow.mag = o.def[im][n].mag;
                                overflow.type = o.def[im][n].type;
                                overflow.kstart = (stree.st[ic].StartBranch+1)*((unsigned long) pow(2, jc));
                                overflow.kend = o.def[im][n].kend;
                                o.def[im][n].kend = overflow.kstart - 1;
                                o.def[im].insert(o.def[im].begin()+n+1,overflow);
                                o.def_map[ijk_index(im,o.def[im][n+1].jn,o.def[im][n+1].kstart,o)].push_back(overflow);
                                k = overflow.kstart;
                            }
							unsigned long nup; //index of next branch
							//check whole row of defects is the same
                            if(k==(stree.st[ic].StartBranch + 1)*pow(2, jc)) n++;   //if end branch is at end of this tree
                            else
                            {
                                while (symmcheck == true && k < (stree.st[ic].StartBranch + 1)*pow(2, jc))
                                {
                                    nup = ijk_index(im,o.def[im][n].jn,k,o);
                                    n += o.def_map[nup - 1].size();   //can now count past all defects in last branch
                                    auto inext = o.def_map.find(nup);  //find if point exists
                                    if (inext != o.def_map.end())
                                    {
                                        if (o.def_map[nup].size() == o.def_map[nup - 1].size()) //are there the same number of defects here?
                                        {
                                            for (unsigned p = 0; p < o.def_map[nup].size(); p++)
                                            {
                                                if (o.def_map[nup][p] == o.def_map[nup - 1][p]);   //are the defects equal
                                                else symmcheck = false;
                                            }
                                            if (symmcheck == true) k=o.def_map[nup][0].kend;   //move to next branch
                                        }
                                        else symmcheck = false;
                                    }
                                    else symmcheck = false;
                                }
                            }
						}
						else
						{
							if (khs >= (stree.st[ic].StartBranch + 1)*pow(2, jc))   //defect is outside this tree
							{
								n++;
							}
							else symmcheck = false;
						}
					}    //end of inner while
					if (symmcheck == false)
					{
						//tree is not symmetric
						stree.st[ic].EndGen = stree.st[ic].StartGen;    //end this sub tree
						stree.st[ic].treeout[0] = ((unsigned) stree.st.size());
						stree.st[ic].treeout[1] = ((unsigned) stree.st.size())+1;
						stnextgen.push_back(stree.st[ic].treeout[0]);
						stnextgen.push_back(stree.st[ic].treeout[1]);
						
						subtree st1;               //create new subtrees with start at next gen
						st1.StartGen = stree.st[ic].StartGen + 1;
						st1.treein = ic;           //tree in is ic
						st1.imeanpath = im;         //it is in mean path im
                        if(stree.st[ic].Ncond > 0) st1.Ncond = stree.st[ic].Ncond - 1;
                        else st1.Ncond = 0;
                        if(stree.st[ic].Ntot > 0) st1.Ntot = stree.st[ic].Ntot - 1;
                        else
                        {
                            stree.st[ic].Ntot = 0;
                            cout << "Error: Should not be triggered\n";
                            cerr << "Less than 0 gens left to build.\n";
                            return 1;
                        }
						stree.st.push_back(st1);    //add these to tree
						stree.mpsubtrees[im].push_back(((unsigned) stree.st.size())-1);
						stree.st.push_back(st1);
						stree.mpsubtrees[im].push_back(((unsigned) stree.st.size())-1);
						stree.st[stree.st[ic].treeout[0]].allocate_gn(stree.st[im].StartGen + stree.st[im].Ntot + 1);
						stree.st[stree.st[ic].treeout[1]].allocate_gn(stree.st[im].StartGen + stree.st[im].Ntot + 1);
						stree.st[stree.st[ic].treeout[0]].StartBranch = 2 * stree.st[ic].StartBranch;   //left and right branch numbers
						stree.st[stree.st[ic].treeout[1]].StartBranch = 2 * stree.st[ic].StartBranch + 1;
						stree.st[stree.st[ic].treeout[0]].z0 = stree.st[ic].z0 - 0.5 / pow(2,stree.st[stree.st[ic].treeout[0]].StartGen);   //left and right branch numbers
						stree.st[stree.st[ic].treeout[1]].z0 = stree.st[ic].z0 + 0.5 / pow(2,stree.st[stree.st[ic].treeout[1]].StartGen);
					}
					else
					{
						//subtree is symetric
						stree.st[ic].EndGen = stree.st[im].StartGen + stree.st[im].Ntot;    //subtree terminates at tree end
						stree.st[ic].treeout[0] = -1;     //no tree out
						stree.st[ic].treeout[1] = -1;
					}
				}
			}       //end of for
		}   //end of while
	}   //end of outer for
		//----Tree structure built-----//
		
	return 0;
}

int pre_run_breaths(Tree &stree, vector<Tree> &ltree, double &time, Options &o)
{
	//run ventilation solver with transport solver to equilibrate
	unsigned long t;
	
	for (t = 0; t < o.PrerunTime / o.dt; t++)
	{
		time = t*o.dt;
		if (update_flux(stree, ltree, time, o))
		{
			cout << "Error in flux update\n.";
			return 1;
		}                            //work out new volume (symmetric case)
		
		for (unsigned long n = 0; n < o.Ntrees; n++)   
		{   //this checks if perturbed tree flow changes sign at different time point
			if ((stree.st[0].Vnew - stree.st[0].Vold + ltree[n].st[0].Vnew - ltree[n].st[0].Vold)*(stree.st[0].Vnew - stree.st[0].Vold) < 0)  cout << "Warning: boundary conditions not perfectly aligned " << n << '\n';
		}
	}
	
	return 0;
}

int initialise_lung_volumes(Tree &stree, vector<Tree> &ltree, Options &o)
{
	//calculate initial values to the total inert gas in acinar units totIGvol

	for (unsigned m = 0; m < stree.EndSubtrees.size(); m++)  //loop trees that reach acinus
	{
		unsigned i = stree.EndSubtrees[m];   //current subtree
		stree.st[i].totIGvol = 0;
		for (unsigned p = 0; p < stree.st[i].isub.size(); p++)  //loop over subtrees
		{
			if (i != stree.st[i].isub[p])stree.st[stree.st[i].isub[p]].totIGvol = 0;
			for (unsigned j = stree.st[stree.st[i].isub[p]].StartGen + stree.st[stree.st[i].isub[p]].Ncond; j <= stree.st[stree.st[i].isub[p]].EndGen; j++)
			{
				for (unsigned k = 0; k < stree.st[stree.st[i].isub[p]].gn[j].p.size(); k++)
				{
					stree.st[stree.st[i].isub[p]].totIGvol += stree.st[stree.st[i].isub[p]].gn[j].p[k].c*stree.st[stree.st[i].isub[p]].gn[j].p[k].Anew*stree.st[stree.st[i].isub[p]].gn[j].dx*stree.st[stree.st[i].isub[p]].gn[j].Nb;
				}
			}
			if (i != stree.st[i].isub[p]) stree.st[i].totIGvol += stree.st[stree.st[i].isub[p]].totIGvol;
		}
	}

	for (unsigned long n = 0; n < o.Ntrees; n++)
	{
		for (unsigned m = 0; m < ltree[n].EndSubtrees.size(); m++)  //loop trees that reach acinus
		{
			unsigned i = ltree[n].EndSubtrees[m];   //current subtree
			ltree[n].st[i].totIGvol = 0;
			for (unsigned p = 0; p < ltree[n].st[i].isub.size(); p++)  //loop over subtrees
			{
				if (i != ltree[n].st[i].isub[p]) ltree[n].st[ltree[n].st[i].isub[p]].totIGvol = 0;
				unsigned isst = ltree[n].st[ltree[n].st[i].isub[p]].i_stree;
				for (unsigned j = ltree[n].st[ltree[n].st[i].isub[p]].StartGen + ltree[n].st[ltree[n].st[i].isub[p]].Ncond; j <= ltree[n].st[ltree[n].st[i].isub[p]].EndGen; j++)
				{
					for (unsigned k = 0; k < ltree[n].st[ltree[n].st[i].isub[p]].gn[j].p.size(); k++)
					{
						ltree[n].st[ltree[n].st[i].isub[p]].totIGvol += ltree[n].st[ltree[n].st[i].isub[p]].gn[j].Nb
							*(ltree[n].st[ltree[n].st[i].isub[p]].gn[j].p[k].c*stree.st[isst].gn[j].p[k].Anew*stree.st[isst].gn[j].dx
							+ stree.st[isst].gn[j].p[k].c*ltree[n].st[ltree[n].st[i].isub[p]].gn[j].p[k].Anew*stree.st[isst].gn[j].dx
							+ stree.st[isst].gn[j].p[k].c*stree.st[isst].gn[j].p[k].Anew*ltree[n].st[ltree[n].st[i].isub[p]].gn[j].dx);
					}
				}
				if (i != ltree[n].st[i].isub[p]) ltree[n].st[i].totIGvol += ltree[n].st[ltree[n].st[i].isub[p]].totIGvol;
			}
		}
		
		ltree[n].masstot = 0;
		vector<double> dmtotold;
		for (unsigned i = 0; i < ltree[n].st.size(); i++)
		{
			unsigned istree = ltree[n].st[i].i_stree;
			dmtotold.push_back(ltree[n].st[i].masstot);
			ltree[n].st[i].masstot = 0;
			for (unsigned j = ltree[n].st[i].StartGen; j <= ltree[n].st[i].EndGen; j++)
			{
				unsigned Nj = ((unsigned)ltree[n].st[i].gn[j].p.size());
				for (unsigned k = 0; k < Nj; k++)
				{
					if (j>0 || k > 0) ltree[n].st[i].masstot += ((double)ltree[n].st[i].gn[j].Nb)*(stree.st[istree].gn[j].p[k].Anew*stree.st[istree].gn[j].p[k].c*ltree[n].st[i].gn[j].dx + stree.st[istree].gn[j].p[k].Anew*ltree[n].st[i].gn[j].p[k].c*stree.st[istree].gn[j].dx + ltree[n].st[i].gn[j].p[k].Anew*stree.st[istree].gn[j].p[k].c*stree.st[istree].gn[j].dx);
				}
			}
			ltree[n].masstot += ltree[n].st[i].masstot;
		}
	}
	return 0;
}
