#include "lung_model_discrete_branching.hpp"
//Checked and commented 16/05/18

int update_flux(Tree &stree, vector<Tree> &ltree, double time, Options &o)   //updates volume of all trees
{
	if (tree_flux(stree, time, o)) return 1;
	for (unsigned n = 0; n < o.Ntrees; n++)
	{
		if (tree_lp_flux(ltree[n], stree, o)) return 1;
	}

	return 0;
}

int tree_flux(Tree &stree, double time, Options &o)    //update volumes on stree
{
	double Vtot;

	stree.Vlungold = stree.Vlung;
	stree.Pplold = stree.Ppl;
	for (unsigned i = 0; i < stree.st.size(); i++)  //loop over subtrees
	{
		stree.st[i].Vold = stree.st[i].Vnew;    //move volume to old value
		stree.st[i].qend = 0;   //set to 0 for all trees
	}

	if (full_tree_fluxcalc(stree, time, o)) return 1;  //do full calculation of tree flux

	for (unsigned m = 0; m < stree.EndSubtrees.size(); m++)  //loop trees that reach acinus
	{
		unsigned i = stree.EndSubtrees[m];   //current subtree
		Vtot = stree.st[i].Vnew / stree.st[i].Valv0;     //total volume to be spread over tree isub[n]

        double Vcheck = 0, Vcheck2 = 0 ;
		for (unsigned n = 0; n < stree.st[i].isub.size(); n++)  //loop over subtrees
		{
			stree.st[stree.st[i].isub[n]].tot_area_function(Vtot);
            for(unsigned j=stree.st[stree.st[i].isub[n]].StartGen; j<=stree.st[stree.st[i].isub[n]].EndGen; j++)
            {
                for(unsigned k=0; k<stree.st[stree.st[i].isub[n]].gn[j].p.size(); k++)
                {
                    Vcheck += (stree.st[stree.st[i].isub[n]].gn[j].p[k].Anew - stree.st[stree.st[i].isub[n]].gn[j].p[k].al)*stree.st[stree.st[i].isub[n]].gn[j].dx*stree.st[stree.st[i].isub[n]].gn[j].Nb;
					Vcheck2 += (stree.st[stree.st[i].isub[n]].gn[j].p[k].Aold - stree.st[stree.st[i].isub[n]].gn[j].p[k].al)*stree.st[stree.st[i].isub[n]].gn[j].dx*stree.st[stree.st[i].isub[n]].gn[j].Nb;
				}
            }
		}
		Vcheck = Vcheck;
	}
	stree.Vlung = stree.Vairways + stree.st[0].Vnew;
	//Trees are ordered by starting generation -- loop backwards
    unsigned ic,i;
	for (ic = 0; ic < stree.st.size(); ic++)
	{
		i = ((unsigned) (stree.st.size() - 1 - ic));     //list of subtrees at this level
		stree.st[i].velocity_calc(stree.st[i].qend, o);
		if (stree.st[i].StartGen > 0) stree.st[stree.st[i].treein].qend += stree.st[i].gn[stree.st[i].StartGen].p[0].ul*stree.st[i].gn[stree.st[i].StartGen].p[0].al*stree.st[i].gn[stree.st[i].StartGen].Nb;   //flux in/out subtree
	}
	
	return 0;
}

int tree_lp_flux(Tree &ltree, Tree &stree, Options &o)    //update volumes on ltree
{
	double DVtot;

	ltree.Vlungold = ltree.Vlung;
	ltree.Pplold = ltree.Ppl;
	for (unsigned i = 0; i < ltree.st.size(); i++)  //loop over subtrees
	{
		ltree.st[i].Vold = ltree.st[i].Vnew;    //move volume to old value
		ltree.st[i].qend = 0;   //set to 0 for all trees
	}

	if (ltree.pert_ij[1] <= ltree.st[ltree.pert_ij[0]].StartGen + ltree.st[ltree.pert_ij[0]].Ncond)   //fluxes only need updating if resistance or stiffness is changed
	{
		if (full_lp_tree_fluxcalc(ltree, stree, o)) return 1;  //do full calculation of tree flux
	}

	for (unsigned m = 0; m < ltree.EndSubtrees.size(); m++)  //loop trees that reach acinus
	{
		unsigned i = ltree.EndSubtrees[m];   //current subtree
        unsigned is = ltree.st[i].i_stree;   //corresponding subtree in stree
		DVtot = (ltree.st[i].Vnew - ltree.st[i].Valv0*(stree.st[is].Vnew / stree.st[is].Valv0))
			/ ((((double)ltree.st[i].gn[ltree.st[i].StartGen].Nb) / ((double)stree.st[is].gn[ltree.st[i].StartGen].Nb))*stree.st[is].Valv0);
		//double Vcheck = 0;
		for (unsigned n = 0; n < ltree.st[i].isub.size(); n++)  //loop over subtrees
		{
			unsigned isst = ltree.st[ltree.st[i].isub[n]].i_stree;  //number if symm st
			ltree.st[ltree.st[i].isub[n]].tot_lp_area_function(stree.st[isst], DVtot);
		//	for (unsigned j = ltree.st[ltree.st[i].isub[n]].StartGen; j <= ltree.st[ltree.st[i].isub[n]].EndGen; j++)
		//	{
		//		for (unsigned k = 0; k<ltree.st[ltree.st[i].isub[n]].gn[j].p.size(); k++)
		//		{
		//			Vcheck += (ltree.st[ltree.st[i].isub[n]].gn[j].p[k].Anew - 0.5*(ltree.st[ltree.st[i].isub[n]].gn[j].p[k].al + ltree.st[ltree.st[i].isub[n]].gn[j].p[k].ar))*ltree.st[ltree.st[i].isub[n]].gn[j].Nb*stree.st[isst].gn[j].dx
		//				+ (stree.st[isst].gn[j].p[k].Anew - 0.5*(stree.st[isst].gn[j].p[k].al + stree.st[isst].gn[j].p[k].ar))*ltree.st[ltree.st[i].isub[n]].gn[j].Nb*ltree.st[ltree.st[i].isub[n]].gn[j].dx;
		//		}
		//	}
		}
		//cout << Vcheck << '\n';
	}
	ltree.Vlung = ltree.Vairways + ltree.st[0].Vnew;
	
	//Trees are ordered by starting generation -- loop backwards
	for (unsigned ic = 0; ic < ltree.st.size(); ic++)
	{
		unsigned i = ((unsigned) ltree.st.size()) - 1 - ic;     //list of subtrees at this level
		unsigned istree = ltree.st[i].i_stree;
		ltree.st[i].velocity_lp_calc(stree.st[ltree.st[i].i_stree], ltree.st[i].qend, o);
		if (ltree.st[i].StartGen > 0) ltree.st[ltree.st[i].treein].qend += ltree.st[i].gn[ltree.st[i].StartGen].Nb*
			(ltree.st[i].gn[ltree.st[i].StartGen].p[0].ul*stree.st[istree].gn[ltree.st[i].StartGen].p[0].al
			 + stree.st[istree].gn[ltree.st[i].StartGen].p[0].ul*ltree.st[i].gn[ltree.st[i].StartGen].p[0].al);   //flux in/out subtree
	}

	return 0;
}



