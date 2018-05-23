#include "lung_model_discrete_branching.hpp"

//--Checked and commented 16/05/18---//
int initialise_tree(Tree &stree, Options &o)   //setup the symmetric tree
{
    //-----Fill tree-----//
    unsigned p, Nxj;
    node ph1;
    double x0 = 0;
    
    p = 0;
    
    stree.Vtot.push_back(o.V0);  //store total lung sac volume
	stree.Ppl = 0;
	stree.Pplold = 0;
    for (unsigned i = 0; i < stree.st.size(); i++) //loop over subtrees
    {
        unsigned imp = stree.st[i].imeanpath;   //stores mean path number
        stree.st[i].fluxin = 0;   //initialise - to be calculated
        stree.st[i].Valv0 = 0;  
        stree.st[i].Vold = stree.st[imp].V0 / pow(2.0, ((double)stree.st[i].StartGen) - ((double)stree.st[imp].StartGen));
        stree.st[i].Vnew = stree.st[i].Vold;     //volume in associated parenchymal sub-unit - (if all at end of tree)
        stree.st[i].V0 = stree.st[i].Vold;   //all sac volumes equal at time = 0
        if (i == 0)      //initialise
		{
			stree.st[i].Ep = 1.0;   //size of elasticity defect
			stree.st[i].Rbp = 1.0;   //size of resistance defect
		}
        else
		{
			stree.st[i].Ep = stree.st[stree.st[i].treein].Ep;   //copy from tree in, so defects of elasticity are multiplicative
			stree.st[i].Rbp = stree.st[stree.st[i].treein].Rbp;
		}
		
        if (i>0)
        {
            x0 = stree.st[i].gn[stree.st[i].StartGen].x0;    //shift x0 to start of this subtree
        }
        for (unsigned j = stree.st[i].StartGen; j <= stree.st[i].EndGen; j++)   //loop over generations
        {
            unsigned jmp = j - stree.st[imp].StartGen;       // j relative to subtree starting gen
			if (j<stree.st[i].StartGen + stree.st[i].Ncond)     //if in conducting tree
			{
				double BranchPeclet = o.Peclet;   
				switch(o.TaylorDisp)    //calculate max branch peclet number based on dispersion type
				{
					case NONE:
					{
						BranchPeclet = o.Peclet*((stree.st[i].V0*stree.st[0].A0)/(o.V0*pow(2.0,jmp)*stree.st[i].area_function(0,jmp,o)));
					} break;
						
					case TAYLOR:
					{
						BranchPeclet = stree.st[i].length_function(jmp,o)/(2.0*sqrt(stree.st[i].area_function(0,jmp,o)/(192 * M_PI)));  //max possible peclet (at u = d/sqrt(a/192pi))
					} break;
						
					case SCHERER:
					{
						BranchPeclet = stree.st[i].length_function(jmp,o)/(0.37*sqrt(stree.st[i].area_function(0,jmp,o)*4/M_PI)); //max possible peclet -- on exhalation -- limit u->infinity
					} break;
				}
				Nxj = max(((unsigned)(BranchPeclet/o.MaxPeclet)), ((unsigned)(stree.st[imp].length_function(jmp, o) / o.dxmax)));   //min spacing (max N points) from either MaxPeclet or dx settings
				Nxj = max(Nxj,o.MinGenSize);   //min spacing could be superceded by min gen size parameter
			}
			else    //in acinus
			{
				Nxj = max(((unsigned)(stree.st[imp].length_function(jmp, o) / o.dxmax)),o.MinAcinGenSize);   //min spacing could be superceded by min acin gen size parameter
			}
            //--work out effect of defects--//
            double Ap = 1.0;                        //factors to multiply a and l by from symm case
            double Lp = 1.0;
            unsigned long jk_map = ijk_index(imp, j - stree.st[imp].StartGen, ((unsigned long)stree.st[i].StartBranch)*((unsigned long)pow(2, j - stree.st[i].StartGen)), o);   //index for first branch in gen
            auto it = o.def_map.find(jk_map);  //find if point exists
            if (it != o.def_map.end())   //apply user inputted defects
            {
                for (unsigned long m = 0; m < o.def_map[jk_map].size(); m++)
                { //loop over defects at this generation
                    switch (o.def_map[jk_map][m].type)
                    {
                        case BLOCKAGE:
                        {
                            Ap = 0;     //zero cross section for blockage
                            stree.st[i].blocked = true;
                        } break;
                        case AREA:
                        {
                            Ap *= (1 + o.def_map[jk_map][m].mag);  //(1 + mag)*a = area of branch
                        } break;
                        case LENGTH:
                        {
                            Lp *= (1 + o.def_map[jk_map][m].mag);  //(1 + mag)*l = length of branch
                        } break;
                        case ELASTICITY:
                        {
                            stree.st[i].Ep *= (1 + o.def_map[jk_map][m].mag);   //bag elastiticity
                        } break;
						case BAG_RESISTANCE:
						{
							stree.st[i].Rbp *= (1 + o.def_map[jk_map][m].mag);    //bag resistance
						}
                    }
                }
            }
            
            stree.st[i].gn[j].dx = Lp*stree.st[imp].length_function(jmp, o) / Nxj;      //voxel length
            stree.st[i].gn[j].Nb = ((unsigned long)(pow(2, j - stree.st[i].StartGen)));   //assign number of branches
            if (j == 0)
            {
                Nxj++;
                x0 = 0;   //this node is external to the system
                stree.st[i].gn[0].x0 = 0;
            }
            else stree.st[i].gn[j].x0 = x0;                  //stores starting position of generation
            for (unsigned k = 0; k < Nxj; k++)              //each node contains c at x[k] and velocity value at x[k-1/2]
            {
				//initialise concentration field
                if (o.InitOp == EMPTY)
                {
                    ph1.c = 0.0;                                            //0 inside
                    ph1.cold = 0.0; //old concentration
                }
                if (o.InitOp == FULL)
                {
                    ph1.c = 1.0;
                    ph1.cold = 1.0;
                }
                if (k == 0 && j == 0)
                {
                    if (o.VDM < 1E-09) ph1.c = c_stim(0, o);                       //concentration initially 1 outside, if no extra ds volume
                    stree.cmouth = c_stim(0, o);
					stree.cmouthold = stree.cmouth;
                    stree.ctop.push_back(ph1.c);
                }
                ph1.ul = 0.0;                            //flow speed (to be calculated)
                ph1.ur = 0.0;                            //flow speed central (to be calculated)
                ph1.Dl = o.Diffusion;                     //Diffusion left as constant
                ph1.Dr = o.Diffusion;
                
                ph1.km = p;                                                    //position in A matrix
                
                ph1.al = Ap*stree.st[imp].area_function(x0 - stree.st[i].gn[j].x0, jmp, o);            //cross sectional areadefined at x[k-1/2]
                ph1.ar = Ap*stree.st[imp].area_function(x0 + stree.st[i].gn[j].dx - stree.st[i].gn[j].x0, jmp, o);     //c-s area defined at x[k+1/2]
                
                ph1.Anew = 0.5*(ph1.al + ph1.ar);                                            //total cs area
                ph1.Aold = ph1.Anew;                                              //to be calculated
                ph1.DA = stree.st[imp].alveolar_density(jmp, k, o)*sqrt(0.5*(ph1.al + ph1.ar));
                
                if (j > 0 || k > 0) x0 += 0.5*stree.st[i].gn[j].dx;                  //position of node centre
                ph1.x = x0;
                ph1.y0 = j + (ph1.x - stree.st[i].gn[j].x0) / (Nxj*stree.st[i].gn[j].dx);   //normalised generation length
                
                stree.st[i].gn[j].p.push_back(ph1);             //add point to gen
                
                if (j > 0 || k > 0) x0 += 0.5*stree.st[i].gn[j].dx;                 //start of next node position
                
                p++;
            }
        }
        stree.st[i].E = o.E*stree.st[i].Ep*(o.V0 / stree.st[i].Vnew);  //elasticity scales with 1/resting volume
        stree.st[i].Rb = o.Rb*stree.st[i].Rbp*(o.V0 / stree.st[i].Vnew);   
        subtree_resistance(stree, i, o);  //calculate resistance of each subtree
        
        if (stree.st[i].EndGen < stree.st[imp].StartGen + stree.st[imp].Ntot)   //if subtree terminates before end of tree
        {
            stree.st[stree.st[i].treeout[0]].gn[stree.st[stree.st[i].treeout[0]].StartGen].x0 = x0;
            stree.st[stree.st[i].treeout[1]].gn[stree.st[stree.st[i].treeout[1]].StartGen].x0 = x0;
        }
        if (stree.st[i].StartGen <= stree.st[imp].StartGen + stree.st[imp].Ncond && stree.st[i].EndGen >= stree.st[imp].StartGen + stree.st[imp].Ncond)  //tree terminates at acinus -- used in flow calc
        {
            int iup = i;   //starting at current tree
            while (iup >= 0)   //while subtree exists
            {
                stree.st[iup].EndSubtrees.push_back(i);  //add to end sub trees (whether blocked or not)
                if (stree.st[iup].blocked) stree.st[i].blocked = true;   //if parent tree is blocked, so is this one
                iup = stree.st[iup].treein;    //move up
            }
            vector<unsigned> idwn;
            if (stree.st[i].treeout[0] > -1) idwn.push_back(stree.st[i].treeout[0]);   //keep track of all trees distal to the terminal bronchiole
            if (stree.st[i].treeout[1] > -1) idwn.push_back(stree.st[i].treeout[1]);
            for (unsigned n = 0; n < ((unsigned)idwn.size()); n++)
            {
                if (stree.st[idwn[n]].treeout[0] > -1) idwn.push_back(stree.st[idwn[n]].treeout[0]);
                if (stree.st[idwn[n]].treeout[1] > -1) idwn.push_back(stree.st[idwn[n]].treeout[1]);
                stree.st[idwn[n]].EndSubtrees.push_back(i);     //all sub branches also have tree i as corresponsing endtree
            }
            if (stree.st[i].blocked == false) stree.EndSubtrees.push_back(i);    //stores end trees for resistance calc
        }
    }
    o.kmtot = p;   //total number of finite volume elements
	stree.kmtot = o.kmtot;

    cout << "Number of points on base tree: " << p << '\n';
    
	stree.Vairways = 0;
    for (unsigned i = 0; i < stree.st.size(); i++)     //fill in properties of neighbours that do not change
    {
        for (unsigned j = stree.st[i].StartGen; j <= stree.st[i].EndGen; j++)
        {
            for (unsigned k = 0; k < stree.st[i].gn[j].p.size(); k++)
            {
                unsigned ih, jh, kh;
				if(j>0 || k>0) stree.Vairways += stree.st[i].gn[j].Nb*stree.st[i].gn[j].dx*stree.st[i].gn[j].p[k].Anew;
                for (unsigned n = 0; n < 3; n++)                                //loop over k+n-1 (k-1, k and k+1)
                {
                    tree_point_pos(stree, i, j, k, ((int)n) - 1, stree.st[i].gn[j].p[k].iup[n], stree.st[i].gn[j].p[k].jup[n], stree.st[i].gn[j].p[k].kup[n]);
                    for (unsigned m = 0; m < ((unsigned)stree.st[i].gn[j].p[k].iup[n].size()); m++) //fill in spacing
                    {
                        ih = stree.st[i].gn[j].p[k].iup[n][m];
                        jh = stree.st[i].gn[j].p[k].jup[n][m];  //store indices of neighbours
                        kh = stree.st[i].gn[j].p[k].kup[n][m];
						if(j==0 && k+n<2)
						{
							stree.st[i].gn[j].p[k].dxup[n].push_back(0);    //stores neighbouring dxs
						}
						else
						{
							stree.st[i].gn[j].p[k].dxup[n].push_back(stree.st[ih].gn[jh].dx);    //stores neighbouring dxs
						}
                    }
                }
				unsigned nodes_right= ((unsigned) stree.st[i].gn[j].p[k].iup[2].size());
                stree.st[i].gn[j].p[k].sl = stree.st[i].gn[j].p[k].al*stree.st[i].gn[j].Nb;    //set sl
                ih = stree.st[i].gn[j].p[k].iup[0][0];
                jh = stree.st[i].gn[j].p[k].jup[0][0];
                kh = stree.st[i].gn[j].p[k].kup[0][0];
                if (ih != i)  //if tree down is different to tree here then left boundary is at a bifurcation
                {
                    unsigned iother;
                    if (stree.st[ih].treeout[0] == ((int)i)) iother = ((unsigned)stree.st[ih].treeout[1]); //find tree number of other tree at bifurcation
                    else iother = ((unsigned)stree.st[ih].treeout[0]);
                    if (stree.st[i].gn[j].p[k].sl == 0. && stree.st[iother].gn[j].p[k].al == 0.) stree.st[i].gn[j].p[k].sdr = 0.5*stree.st[ih].gn[jh].p[kh].ar*stree.st[ih].gn[jh].Nb;
                    else stree.st[i].gn[j].p[k].sdr = stree.st[ih].gn[jh].p[kh].ar*stree.st[ih].gn[jh].Nb*stree.st[i].gn[j].p[k].sl / (stree.st[i].gn[j].p[k].sl + stree.st[iother].gn[j].p[k].al*stree.st[iother].gn[j].Nb); //flux splits asymmetrically at bifurc.
                }
                else
                {
                    stree.st[i].gn[j].p[k].sdr = stree.st[ih].gn[jh].p[kh].ar*stree.st[ih].gn[jh].Nb;   //right edge of element down
                }
                for (unsigned m = 0; m < nodes_right; m++) //flux(es) at right boundary
                {
                    ih = stree.st[i].gn[j].p[k].iup[2][m];
                    jh = stree.st[i].gn[j].p[k].jup[2][m];
                    kh = stree.st[i].gn[j].p[k].kup[2][m];
                    stree.st[i].gn[j].p[k].sul.push_back(stree.st[ih].gn[jh].p[kh].al*stree.st[ih].gn[jh].Nb);
                }

                for (unsigned m = 0; m < nodes_right; m++) //flux(es) at right boundary
                {
                    if (nodes_right == 1) stree.st[i].gn[j].p[k].sr.push_back(stree.st[i].gn[j].p[k].ar*stree.st[i].gn[j].Nb);  //c-s area of edge of cell
                    else stree.st[i].gn[j].p[k].sr.push_back(stree.st[i].gn[j].p[k].ar*stree.st[i].gn[j].Nb*stree.st[i].gn[j].p[k].sul[m] / (stree.st[i].gn[j].p[k].sul[m] + stree.st[i].gn[j].p[k].sul[bitflip(m)]));
                }

				
				upwind_coeffs(stree.st[i].gn[j].p[k].ucflpos, stree.st[i].gn[j].p[k].ucfrpos, stree.st[i].gn[j].p[k].dxup, 1.0, o);    //coefficients for upwind scheme based on spacing
				upwind_coeffs(stree.st[i].gn[j].p[k].ucflneg, stree.st[i].gn[j].p[k].ucfrneg, stree.st[i].gn[j].p[k].dxup, -1.0, o);    //coefficients for upwind scheme based on spacing
				fd_coeffs(stree.st[i].gn[j].p[k].dcfl, stree.st[i].gn[j].p[k].dcfr, stree.st[i].gn[j].p[k].dxup);            //coefficients for fd scheme based on spacing
            }
        }
    }
	stree.Vlungold = stree.Vairways + stree.st[0].Vold;
	stree.Vlung = stree.Vairways + stree.st[0].Vnew;
    //first voxel has zero size
    stree.st[0].gn[0].p[0].dxup[1][0] = 0.0;
    stree.st[0].gn[0].p[1].dxup[0][0] = 0.0;

    //-------Calculate volume in alveolar region--------//
    double Vtot;
    //    double V0check = 0;
    //    double V0check2 = 0;
    for (unsigned m = 0; m < stree.EndSubtrees.size(); m++) //loop over subtrees
    {
        unsigned i = stree.EndSubtrees[m];   //terminating subtree number
		stree.STno_to_ESTno[i] = m;              //inverse map from number of ST in network, to numbered terminal ST in res matrix
        stree.st[i].Valv0 = 0;
		stree.st[i].Vacinduct = 0;
        stree.st[i].isub.clear();    //vector for storing sub-trees of i
        stree.st[i].isub.push_back(i);
        for (unsigned n = 0; n < stree.st[i].isub.size(); n++)  //make list of subtrees
        {
            if (stree.st[stree.st[i].isub[n]].treeout[0] > -1) 		stree.st[i].isub.push_back(stree.st[stree.st[i].isub[n]].treeout[0]);
            if (stree.st[stree.st[i].isub[n]].treeout[1] > -1) 		stree.st[i].isub.push_back(stree.st[stree.st[i].isub[n]].treeout[1]);
        }

        for (unsigned n = 0; n < stree.st[i].isub.size(); n++)  //loop over subtrees in acinus
        {
            unsigned ih = stree.st[i].isub[n];   //number of subtree
            stree.st[ih].Valv0 = 0;
			stree.st[ih].Vacinduct = 0;
			for (unsigned j = stree.st[ih].StartGen + stree.st[ih].Ncond; j <= stree.st[ih].EndGen; j++)
            {
                for (unsigned k = 0; k < stree.st[ih].gn[j].p.size(); k++)
                {
					stree.st[ih].Valv0 += stree.st[ih].gn[j].p[k].DA*stree.st[ih].gn[j].dx*stree.st[ih].gn[j].Nb;    //tot up surface area available for alveoli
					stree.st[ih].Vacinduct += 0.5*(stree.st[ih].gn[j].p[k].al + stree.st[ih].gn[j].p[k].ar)*stree.st[ih].gn[j].dx*stree.st[ih].gn[j].Nb;  //tot up acin duct volume
                }
            }
			if (ih != i)
			{
				stree.st[i].Valv0 += stree.st[ih].Valv0;   //add all subtree valv0 to valv0
				stree.st[i].Vacinduct += stree.st[ih].Vacinduct;
			}
        }

        Vtot = stree.st[i].Vnew / stree.st[i].Valv0;     //total volume to be spread over tree isub[n]
        //        V0check2 += stree.st[i].Vnew;

        //double Vcheck = 0;
        for (unsigned n = 0; n < stree.st[i].isub.size(); n++)  //loop over subtrees
        {
            stree.st[stree.st[i].isub[n]].tot_area_function(Vtot);
            for (unsigned j = stree.st[stree.st[i].isub[n]].StartGen; j <= stree.st[stree.st[i].isub[n]].EndGen; j++)
            {
                for (unsigned k = 0; k < stree.st[stree.st[i].isub[n]].gn[j].p.size(); k++)
                {
                    stree.st[stree.st[i].isub[n]].gn[j].p[k].Aold = stree.st[stree.st[i].isub[n]].gn[j].p[k].Anew;
                    //                    V0check += (stree.st[stree.st[i].isub[n]].gn[j].p[k].Anew - 0.5*(stree.st[stree.st[i].isub[n]].gn[j].p[k].al + stree.st[stree.st[i].isub[n]].gn[j].p[k].ar))*stree.st[stree.st[i].isub[n]].gn[j].Nb*stree.st[stree.st[i].isub[n]].gn[j].dx;
                }
            }
        }
    }

    stree.masstot = 0;
    for (unsigned i = 0; i < stree.st.size(); i++)     //fill in properties of neighbours that do not change
    {
        stree.st[i].masstot = 0;    //initialisation
        for (unsigned j = stree.st[i].StartGen; j <= stree.st[i].EndGen; j++)
        {
            for (unsigned k = 0; k < stree.st[i].gn[j].p.size(); k++)
            { 
                if (j>0 || k>0) stree.st[i].masstot += stree.st[i].gn[j].Nb*stree.st[i].gn[j].dx*stree.st[i].gn[j].p[k].Anew*stree.st[i].gn[j].p[k].c;  //add up total gas volume
            }
        }
        stree.masstot += stree.st[i].masstot;
    }
	stree.masstotold = stree.masstot;

    return 0;
}

void setup_pert_trees(vector<Tree> &ltree, Tree &stree, Options &o)
{
    o.Ntrees = 0;    //count how many ltrees we need
    for (unsigned im = 0; im < stree.meanpath.size(); im++)
    {
        for(unsigned m = 0; m < stree.mpsubtrees[im].size(); m++)
        {
            unsigned i = stree.mpsubtrees[im][m];   //loop over every subtree of each mean path
            if (stree.st[i].StartGen <= stree.st[im].StartGen + stree.st[im].Ncond && stree.st[i].EndGen >= stree.st[im].StartGen + stree.st[im].Ncond) o.Ntrees += 1;                    //1 for each elasticity pert (only for terminal sub-trees)
            o.Ntrees += 2 * (stree.st[i].EndGen + 1 - stree.st[i].StartGen);      //2 for each area and length pert
        }
    }
    ltree.resize(o.Ntrees);  //declared perturbed trees
    
    calc_perts(ltree, stree);
	
    //-----Setup perturbed tree-----//
    if (o.TreeOp == LINEARPERT)
    {
        for (unsigned n = 0; n < o.Ntrees; n++)
        {
            for (unsigned i = 0; i < stree.st.size(); i++)
            {
                ltree[n].st.push_back(stree.st[i]);    //copy subtree from base tree
            }
        }
        setup_lp_tree(ltree, stree, o);
    }
}


void calc_perts(vector<Tree> &ltree, Tree &stree)  //assigns perturbation values to trees
{
    unsigned n = 0;
    for (unsigned im = 0; im < stree.meanpath.size(); im++)
    {
        for(unsigned m = 0; m < stree.mpsubtrees[im].size(); m++)
        {
            unsigned i = stree.mpsubtrees[im][m];   //loop over every subtree of each mean path
            if (stree.st[i].StartGen <= stree.st[im].StartGen + stree.st[im].Ncond && stree.st[i].EndGen >= stree.st[im].StartGen + stree.st[im].Ncond)
            {
                ltree[n].pert_ij[0] = i;           //Add elasticity pert
                ltree[n].pert_ij[1] = stree.st[im].StartGen + stree.st[im].Ncond;     //elasticity perturbs final generation of res tree
                ltree[n].pert_type = 'K';
                n++;          //this perturbation added
            }
            for (unsigned j = stree.st[i].StartGen; j <= stree.st[i].EndGen; j++)
            {
                ltree[n].pert_ij[0] = i;           //Add area pert
                ltree[n].pert_ij[1] = j;
                ltree[n].pert_type = 'A';
                n++;          //this perturbation added
                ltree[n].pert_ij[0] = i;           //Add area pert
                ltree[n].pert_ij[1] = j;
                ltree[n].pert_type = 'L';
                n++;          //this perturbation added
            }
        }
    }  //perturbations associated with ltrees
}

void setup_lp_tree(vector<Tree> &ltree, Tree &stree, Options &o)    //set up perturbed tree
{	
    int i_prev = ((int) stree.st.size()) - 1;  //max i value in stree
    for (unsigned n = 0; n < o.Ntrees; n++)    //looping over ltrees
    {
        ltree[n].EndSubtrees.clear();
        ltree[n].meanpath = stree.meanpath;
        ltree[n].mpsubtrees = stree.mpsubtrees;
		ltree[n].cmouth = 0;
		ltree[n].cmouthold = 0;
		ltree[n].Ppl = 0;
		ltree[n].Pplold = 0;
        for (unsigned i = 0; i < stree.st.size(); i++)
        {
            unsigned im = stree.st[i].imeanpath;
            ltree[n].st[i].set_quantities_zero();    //set all perturbed quantities to zero
            ltree[n].st[i].EndSubtrees.clear();   //needs to be reset
            if (ltree[n].pert_ij[0] == i)   //perturbation is in this sub tree
            {
                if (stree.st[i].StartGen == ltree[n].pert_ij[1])  //perturbation is in first gen of tree
                {
                    //do pert on single branch
                    ltree[n].st[i] = stree.st[i];  //set equal to stree
                    ltree[n].st[i].EndSubtrees.clear();   //needs to be reset
                    ltree[n].st[i].i_stree = i;
                    ltree[n].st[i].Ncond = stree.st[i].Ncond;  //copy ncond and ntot values
                    ltree[n].st[i].Ntot = stree.st[i].Ntot;
                    ltree[n].st[i].set_quantities_zero();    //set all perturbed quantities to zero
                    for (unsigned j = ltree[n].st[i].StartGen; j <= ltree[n].st[i].EndGen; j++)
                    {
                        ltree[n].st[i].gn[j].set_quantities_zero();
                        for (unsigned k = 0; k < ltree[n].st[i].gn[j].p.size(); k++)   //count through generations up to perturbation
                        {
                            ltree[n].st[i].gn[j].p[k].set_quantities_zero();
                        }
                    }
                    //Rtree and E should automatically be copied (not zeroed)
                    //add pertubations
                    if (ltree[n].pert_type == 'A') ltree[n].st[i].gn[ltree[n].pert_ij[1]].Ap = DEFPERTSIZE;
                    if (ltree[n].pert_type == 'L') ltree[n].st[i].gn[ltree[n].pert_ij[1]].Lp = DEFPERTSIZE;
                    if (ltree[n].pert_type == 'K') ltree[n].st[i].Ep = DEFPERTSIZE;
                }
                else  //tree needs resizing
                {
                    ltree[n].st[i].EndGen = ltree[n].st[i].StartGen;      //branch terminates in same generation now
                    ltree[n].st[i].treeout[0] = i_prev + 1;     //connected to trees iprev and iprev+1
                    ltree[n].st[i].treeout[1] = i_prev + 2;
                    ltree[n].st[i].i_stree = i;
                    ltree[n].st[i].Ncond = stree.st[i].Ncond;  //copy ncond and ntot values
                    ltree[n].st[i].Ntot = stree.st[i].Ntot;
                    ltree[n].st[i].gn[ltree[n].st[i].StartGen].set_quantities_zero();
                    for (unsigned k = 0; k < ltree[n].st[i].gn[ltree[n].st[i].StartGen].p.size(); k++)   //zero all parts
                    {
                        ltree[n].st[i].gn[ltree[n].st[i].StartGen].p[k].set_quantities_zero();
                    }
                    //-----Change iup for endpoint of prev branch-------//
                    ltree[n].st[i].gn[ltree[n].st[i].StartGen].p[ltree[n].st[i].gn[ltree[n].st[i].StartGen].p.size() - 1].iup[2][0] = i_prev + 1;
                    if (ltree[n].st[i].gn[ltree[n].st[i].StartGen].p[ltree[n].st[i].gn[ltree[n].st[i].StartGen].p.size() - 1].iup[2].size() > 1)
                    {
                        ltree[n].st[i].gn[ltree[n].st[i].StartGen].p[ltree[n].st[i].gn[ltree[n].st[i].StartGen].p.size() - 1].iup[2][1] = i_prev + 2;
                    }
					else
					{
						ltree[n].st[i].gn[ltree[n].st[i].StartGen].p[ltree[n].st[i].gn[ltree[n].st[i].StartGen].p.size() - 1].iup[2].push_back(i_prev + 2);
						ltree[n].st[i].gn[ltree[n].st[i].StartGen].p[ltree[n].st[i].gn[ltree[n].st[i].StartGen].p.size() - 1].jup[2].push_back(ltree[n].st[i].gn[ltree[n].st[i].StartGen].p[ltree[n].st[i].gn[ltree[n].st[i].StartGen].p.size() - 1].jup[2][0]); //need resizing: values same as 0 entry
						ltree[n].st[i].gn[ltree[n].st[i].StartGen].p[ltree[n].st[i].gn[ltree[n].st[i].StartGen].p.size() - 1].kup[2].push_back(ltree[n].st[i].gn[ltree[n].st[i].StartGen].p[ltree[n].st[i].gn[ltree[n].st[i].StartGen].p.size() - 1].kup[2][0]);
						ltree[n].st[i].gn[ltree[n].st[i].StartGen].p[ltree[n].st[i].gn[ltree[n].st[i].StartGen].p.size() - 1].sr.push_back(0);   //these need to be right size
						ltree[n].st[i].gn[ltree[n].st[i].StartGen].p[ltree[n].st[i].gn[ltree[n].st[i].StartGen].p.size() - 1].sul.push_back(0);
						ltree[n].st[i].gn[ltree[n].st[i].StartGen].p[ltree[n].st[i].gn[ltree[n].st[i].StartGen].p.size() - 1].dxup[2].push_back(0);
					}
                    //---------------------------------------------------//
                    
                    subtree_lp_symm_resistance(ltree[n].st[i], stree.st[i], o);        //calculate resistance of each subtree
					if(ltree[n].st[i].Ntot < o.Ngen2-o.Ngen) ltree[n].st[i].Rtree = 0;  //st starts in acinus: 0 effective resistance
                    unsigned Pt = i;                     //Parent branch number for next iteration
                    for (unsigned j = ltree[n].st[i].StartGen + 1; j <= ltree[n].pert_ij[1]; j++)   //count through generations up to perturbation
                    {
                        subtree sth;                                    //left (even) branches are on perturbed route, right branches are not
                        unsigned jr = j - ltree[n].st[i].StartGen;
                        unsigned ir = (i_prev + 2 * jr);            //right trees (unperturbed)
                        unsigned il = (i_prev + 2 * jr - 1);
                        sth.treein = Pt;
                        sth.StartGen = j;
                        sth.i_stree = i;    //stores subtree number in stree
                        sth.imeanpath = im;    //stores mean path number
                        sth.allocate_gn(stree.st[i].EndGen + 1);
                        sth.gn[j] = stree.st[i].gn[j];                  //set equal to previous
                        sth.gn[j].Nb = 1;                               //always starts with single branch
                        sth.set_quantities_zero();
                        sth.gn[j].set_quantities_zero();                //zero all perts
						if((j-stree.st[i].StartGen) > stree.st[i].Ncond) sth.Ncond = 0;
                        else sth.Ncond = stree.st[i].Ncond - (j-stree.st[i].StartGen);  //update ncond and ntot values
                        sth.Ntot = stree.st[i].Ntot  - (j-stree.st[i].StartGen);
                        ltree[n].st[i].gn[j].p.clear();                 //erase points that used to be in tree i
                        ltree[n].st.push_back(sth);                      //add new trees
                        ltree[n].mpsubtrees[im].push_back(((unsigned) ltree[n].st.size())-1);
                        ltree[n].st.push_back(sth);
                        ltree[n].mpsubtrees[im].push_back(((unsigned) ltree[n].st.size())-1);
                        for (unsigned k = 0; k < sth.gn[j].p.size(); k++)
                        {
                            ltree[n].st[il].gn[j].p[k].set_quantities_zero();
                            ltree[n].st[ir].gn[j].p[k].set_quantities_zero();
                            //----correct iup (jup and kup unchanged)---//
                            for (unsigned mm = 0; mm < 3; mm++)
                            {
                                ltree[n].st[il].gn[j].p[k].iup[mm][0] = il;
                                ltree[n].st[ir].gn[j].p[k].iup[mm][0] = ir;
                            }
                            if (k == 0)
                            {
                                ltree[n].st[il].gn[j].p[k].iup[0][0] = ltree[n].st[il].treein;
                                ltree[n].st[ir].gn[j].p[k].iup[0][0] = ltree[n].st[ir].treein;
                            }
                        }
                        ltree[n].st[il].StartBranch = stree.st[i].StartBranch*((unsigned long)pow(2, j - stree.st[i].StartGen));
                        ltree[n].st[ir].StartBranch = stree.st[i].StartBranch*((unsigned long)pow(2, j - stree.st[i].StartGen)) + 1;
						ltree[n].st[il].z0 = ltree[n].st[ltree[n].st[il].treein].z0 - 0.5/pow(2.0,ltree[n].st[il].StartGen);
						ltree[n].st[ir].z0 = ltree[n].st[ltree[n].st[ir].treein].z0 + 0.5/pow(2.0,ltree[n].st[il].StartGen);
                        if (j == ltree[n].pert_ij[1])                   //final trees
                        {
                            ltree[n].st[il].EndGen = stree.st[i].EndGen;
                            ltree[n].st[il].treeout[0] = stree.st[i].treeout[0];         //-1 indicates no tree out
                            ltree[n].st[il].treeout[1] = stree.st[i].treeout[1];
                            for (unsigned jp = j + 1; jp <= ltree[n].st[il].EndGen; jp++)   //assign rest of generations
                            {
                                ltree[n].st[il].gn[jp] = stree.st[i].gn[jp];
                                ltree[n].st[il].gn[jp].Nb = ((unsigned long)pow(2, jp - j));
                                ltree[n].st[il].gn[jp].set_quantities_zero();
                                for (unsigned k = 0; k < ltree[n].st[il].gn[jp].p.size(); k++)   //count through generations up to perturbation
                                {
                                    ltree[n].st[il].gn[jp].p[k].set_quantities_zero();
                                    //----------Set iup for all points in this new tree------//
                                    for (unsigned mm = 0; mm < 3; mm++) ltree[n].st[il].gn[jp].p[k].iup[mm][0] = il;
                                    //-------------------------------------------------------//
                                }
                            }
                            if (ltree[n].pert_type == 'A') ltree[n].st[il].gn[j].Ap = DEFPERTSIZE;
                            if (ltree[n].pert_type == 'L') ltree[n].st[il].gn[j].Lp = DEFPERTSIZE;
                            if (ltree[n].pert_type == 'K') ltree[n].st[il].Ep = DEFPERTSIZE;
                            ltree[n].pert_ij[0] = il;    //update location of perturbation
                        }
                        else
                        {
                            ltree[n].st[il].EndGen = j;                  //left trees terminate in same gen
                            ltree[n].st[il].treeout[0] = il + 2;
                            ltree[n].st[il].treeout[1] = il + 3;
                            //----------Set iup for end point in this new tree------//
                            ltree[n].st[il].gn[j].p[ltree[n].st[il].gn[j].p.size() - 1].iup[2][0] = il + 2;
                            if (ltree[n].st[il].gn[j].p[ltree[n].st[il].gn[j].p.size() - 1].iup[2].size() > 1)
                            {
                                ltree[n].st[il].gn[j].p[ltree[n].st[il].gn[j].p.size() - 1].iup[2][1] = il + 3;
                            }
							else
							{
								ltree[n].st[il].gn[j].p[ltree[n].st[il].gn[j].p.size() - 1].iup[2].push_back(il + 3);
								ltree[n].st[il].gn[j].p[ltree[n].st[il].gn[j].p.size() - 1].jup[2].push_back(ltree[n].st[il].gn[j].p[ltree[n].st[il].gn[j].p.size() - 1].jup[2][0]);
								ltree[n].st[il].gn[j].p[ltree[n].st[il].gn[j].p.size() - 1].kup[2].push_back(ltree[n].st[il].gn[j].p[ltree[n].st[il].gn[j].p.size() - 1].kup[2][0]);
								ltree[n].st[il].gn[j].p[ltree[n].st[il].gn[j].p.size() - 1].sr.push_back(0);
								ltree[n].st[il].gn[j].p[ltree[n].st[il].gn[j].p.size() - 1].sul.push_back(0);
								ltree[n].st[il].gn[j].p[ltree[n].st[il].gn[j].p.size() - 1].dxup[2].push_back(0);
							}
                            //-------------------------------------------------------//
                        }
                        ltree[n].st[ir].treeout[0] = stree.st[i].treeout[0];        //right trees connected to base
                        ltree[n].st[ir].treeout[1] = stree.st[i].treeout[1];
                        ltree[n].st[ir].EndGen = stree.st[i].EndGen;
                        for (unsigned jp = j + 1; jp <= ltree[n].st[ir].EndGen; jp++)   //assign rest of generations
                        {
                            ltree[n].st[ir].gn[jp] = stree.st[i].gn[jp];
                            ltree[n].st[ir].gn[jp].Nb = ((unsigned long)pow(2, jp - j));
                            ltree[n].st[ir].gn[jp].set_quantities_zero();
                            for (unsigned k = 0; k < ((unsigned) ltree[n].st[ir].gn[jp].p.size()); k++)   //count through generations up to perturbation
                            {
                                ltree[n].st[ir].gn[jp].p[k].set_quantities_zero();
                                //----------Set iup for all points in this new tree------//
                                for (unsigned mm = 0; mm < 3; mm++) ltree[n].st[ir].gn[jp].p[k].iup[mm][0] = ir;
                                //-------------------------------------------------------//
                            }
                        }
                        subtree_lp_symm_resistance(ltree[n].st[il], stree.st[ltree[n].st[il].i_stree], o);        //calculate resistance of each subtree
                        subtree_lp_symm_resistance(ltree[n].st[ir], stree.st[ltree[n].st[ir].i_stree], o);        //calculate resistance of each subtree
						if(ltree[n].st[il].Ntot < o.Ngen2-o.Ngen) ltree[n].st[il].Rtree = 0;
						if(ltree[n].st[ir].Ntot < o.Ngen2-o.Ngen) ltree[n].st[ir].Rtree = 0;
                        ltree[n].st[il].E = pow(2.0, (j - stree.st[ltree[n].st[il].i_stree].StartGen))*stree.st[ltree[n].st[il].i_stree].E;  //stiffness calculation
                        ltree[n].st[ir].E = pow(2.0, (j - stree.st[ltree[n].st[ir].i_stree].StartGen))*stree.st[ltree[n].st[ir].i_stree].E;  //stiffness calculation
                        ltree[n].st[il].Rb = pow(2.0, (j - stree.st[ltree[n].st[il].i_stree].StartGen))*stree.st[ltree[n].st[il].i_stree].Rb;  //stiffness calculation
                        ltree[n].st[ir].Rb = pow(2.0, (j - stree.st[ltree[n].st[ir].i_stree].StartGen))*stree.st[ltree[n].st[ir].i_stree].Rb;  //stiffness calculation
                        Pt = il;                                       //Parent branch number for next iteration
                    }
                    for (unsigned j = ltree[n].pert_ij[1]; j <= stree.st[i].EndGen; j++) ltree[n].st[i].gn[j].p.clear();   //delete all excess points from ltree
                }
            }
            else //not perturbed subtree
            {
                ltree[n].st[i] = stree.st[i];  //set equal to stree
                ltree[n].st[i].i_stree = i;
                ltree[n].st[i].set_quantities_zero();    //set all perturbed quantities to zero
                for (unsigned j = ltree[n].st[i].StartGen; j <= ltree[n].st[i].EndGen; j++)
                {
                    ltree[n].st[i].gn[j].set_quantities_zero();
                    for (unsigned k = 0; k < ltree[n].st[i].gn[j].p.size(); k++)   //count through generations up to perturbation
                    {
                        ltree[n].st[i].gn[j].p[k].set_quantities_zero();
                    }
                }
            }
        }
		
		ltree[n].ctop.push_back(0);
        apply_pert(ltree[n], stree);
        //ltree[n] is now built
        
        //re-build endsubtree structure
		unsigned long kmcount = 0;
        for (unsigned i = 0; i < ltree[n].st.size(); i++)
        {
            unsigned im = ltree[n].st[i].imeanpath;

            if (ltree[n].st[i].Ntot + o.Ngen >= o.Ngen2 && ltree[n].st[i].EndGen >= stree.st[im].StartGen + stree.st[im].Ncond)  //tree terminates at acinus -- used in flow calc
            {
                int iup = i;   //starting at current tree
                while (iup >= 0)   //while subtree exists
                {
                    ltree[n].st[iup].EndSubtrees.push_back(i);  //add to end sub trees (whether blocked or not)
                    //if (stree.st[iup].blocked) ltree[n].st[i].blocked = true;   //if parent tree is blocked, so is this one
                    iup = ltree[n].st[iup].treein;    //move up
                }
                vector<unsigned long> idwn;
                if (ltree[n].st[i].treeout[0] > -1) idwn.push_back(ltree[n].st[i].treeout[0]);
                if (ltree[n].st[i].treeout[1] > -1) idwn.push_back(ltree[n].st[i].treeout[1]);
                for (unsigned z = 0; z < idwn.size(); z++)
                {
                    if (ltree[n].st[idwn[z]].treeout[0] > -1) idwn.push_back(ltree[n].st[idwn[z]].treeout[0]);
                    if (ltree[n].st[idwn[z]].treeout[1] > -1) idwn.push_back(ltree[n].st[idwn[z]].treeout[1]);
                    ltree[n].st[idwn[z]].EndSubtrees.push_back(i);     //all sub branches also have tree i as corresponsing endtree
                }
                if (stree.st[ltree[n].st[i].i_stree].blocked == false) ltree[n].EndSubtrees.push_back(i);    //stores end trees for resistance calc
            }
			for(unsigned j = ltree[n].st[i].StartGen; j <= ltree[n].st[i].EndGen; j++)
			{
				for(unsigned k = 0; k < ltree[n].st[i].gn[j].p.size(); k++)   //count number of points and build upwind coeffs map
				{
					unsigned istree = ltree[n].st[i].i_stree;
					lp_upwind_coeffs(ltree[n].st[i].gn[j].p[k].ucflpos, ltree[n].st[i].gn[j].p[k].ucfrpos, stree.st[istree].gn[j].p[k].dxup, ltree[n].st[i].gn[j].p[k].dxup, 1.0, o);    //coefficients for upwind scheme based on spacing
					lp_upwind_coeffs(ltree[n].st[i].gn[j].p[k].ucflneg, ltree[n].st[i].gn[j].p[k].ucfrneg, stree.st[istree].gn[j].p[k].dxup, ltree[n].st[i].gn[j].p[k].dxup, -1.0, o);    //coefficients for upwind scheme based on spacing
					lp_fd_coeffs(ltree[n].st[i].gn[j].p[k].dcfl, ltree[n].st[i].gn[j].p[k].dcfr, stree.st[istree].gn[j].p[k].dxup, ltree[n].st[i].gn[j].p[k].dxup);            //coefficients for fd scheme based on spacing
					ltree[n].st[i].gn[j].p[k].km = kmcount;
					kmcount++;
				}
			}
        }
		ltree[n].kmtot = kmcount;
        
        //-------Calculate volume in alveolar region--------//
        //double Vtot;
        //    double V0check = 0;
        //    double V0check2 = 0;
        for (unsigned m = 0; m < ltree[n].EndSubtrees.size(); m++) //loop over subtrees
        {
            unsigned i = ltree[n].EndSubtrees[m];   //terminating subtree number
			ltree[n].STno_to_ESTno[i] = m;              //inverse map from number of ST in network, to numbered terminal ST in res matrix
            unsigned is = ltree[n].st[i].i_stree;
            ltree[n].st[i].isub.clear();    //vector for storing sub-trees of i
            ltree[n].st[i].isub.push_back(i);
            for (unsigned p = 0; p < ltree[n].st[i].isub.size(); p++)  //make list of subtrees
            {
                if (ltree[n].st[ltree[n].st[i].isub[p]].treeout[0] > -1) 		ltree[n].st[i].isub.push_back(ltree[n].st[ltree[n].st[i].isub[p]].treeout[0]);
                if (ltree[n].st[ltree[n].st[i].isub[p]].treeout[1] > -1) 		ltree[n].st[i].isub.push_back(ltree[n].st[ltree[n].st[i].isub[p]].treeout[1]);
				if (i != ltree[n].st[i].isub[p]) ltree[n].st[i].Valv0 += ltree[n].st[ltree[n].st[i].isub[p]].Valv0;  //total Valv0 is sum of all sub Valv0s
			}
            
            double sVtot = stree.st[is].Vnew / stree.st[is].Valv0;     //total volume to be spread over tree 
            double DVtot = (ltree[n].st[i].Vnew - ltree[n].st[i].Valv0*sVtot)
				/ ((((double)ltree[n].st[i].gn[ltree[n].st[i].StartGen].Nb) / ((double)stree.st[is].gn[ltree[n].st[i].StartGen].Nb))*stree.st[is].Valv0); 

			for (unsigned p = 0; p < ltree[n].st[i].isub.size(); p++)  //make list of subtrees
			{
				unsigned isp = ltree[n].st[ltree[n].st[i].isub[p]].i_stree;
				ltree[n].st[ltree[n].st[i].isub[p]].tot_lp_area_function(stree.st[isp], DVtot);
				for (unsigned j = ltree[n].st[ltree[n].st[i].isub[p]].StartGen; j <= ltree[n].st[ltree[n].st[i].isub[p]].EndGen; j++)
				{
					for (unsigned k = 0; k < ltree[n].st[ltree[n].st[i].isub[p]].gn[j].p.size(); k++)
					{
						ltree[n].st[ltree[n].st[i].isub[p]].gn[j].p[k].Aold = ltree[n].st[ltree[n].st[i].isub[p]].gn[j].p[k].Anew;
					}
				}
			}

            //recalc sValv0
            //        V0check2 += stree.st[i].Vnew;
        }
        //if (ltree[n].pert_ij[1] >= o.Ngen && ltree[n].pert_type != 'K')   //if branch pert is in acinus, need to re-calculate sValv0
    }
}

void apply_pert(Tree &ltree, Tree &stree)
{
    unsigned ih = ltree.pert_ij[0];
    unsigned jh = ltree.pert_ij[1];
    unsigned is = ltree.st[ih].i_stree;
    unsigned il,jl,kl;
    vector<unsigned> ir,jr,kr;
    ///Add changes to Vairways, Vlung, IG vols and mass tot to this
    switch (ltree.pert_type)
    {
        case 'A':
        {
			ltree.Vairways = 0;  //change in DS volume
            for (unsigned k = 0; k < ltree.st[ih].gn[jh].p.size(); k++)
            {
                il = ltree.st[ih].gn[jh].p[k].iup[0][0];  //tree no. left
                ir = ltree.st[ih].gn[jh].p[k].iup[2];     //tree no. right (vector)
                jl = ltree.st[ih].gn[jh].p[k].jup[0][0];  //gen left
                jr = ltree.st[ih].gn[jh].p[k].jup[2];     //gen right (vector)
                kl = ltree.st[ih].gn[jh].p[k].kup[0][0];  //branch no. left
                kr = ltree.st[ih].gn[jh].p[k].kup[2];     //branch no. right (vector)
                ltree.st[ih].gn[jh].p[k].al =  ltree.st[ih].gn[jh].Ap*stree.st[is].gn[jh].p[k].al;   //perturbation to area
				ltree.st[ih].gn[jh].p[k].ar = ltree.st[ih].gn[jh].Ap*stree.st[is].gn[jh].p[k].ar;   //perturbation to area
				ltree.st[ih].gn[jh].p[k].Aold = 0.5*(ltree.st[ih].gn[jh].p[k].ar + ltree.st[ih].gn[jh].p[k].al);
				ltree.st[ih].gn[jh].p[k].Anew = 0.5*(ltree.st[ih].gn[jh].p[k].ar + ltree.st[ih].gn[jh].p[k].al);   //perturbation to area
                ltree.st[ih].gn[jh].p[k].sl = ltree.st[ih].gn[jh].p[k].al*ltree.st[ih].gn[jh].Nb;    //perturbation to total duct cross section
                if(ltree.st[il].gn[jl].p[kl].iup[2].size() > 1)    //deal with bifurcation upstream -- relative sr's change
                {
                    unsigned ils = ltree.st[il].i_stree;  //tree no. left on stree
                    unsigned iother, isother;
                    for(unsigned m = 0; m < ltree.st[il].gn[jl].p[kl].iup[2].size(); m++)
                    {
                        if(ltree.st[il].gn[jl].p[kl].iup[2][m] == ih)   //if branch m is this branch
                        {
                            iother = ltree.st[il].gn[jl].p[kl].iup[2][bitflip(m)];   //other treeout of point to left
                            isother = ltree.st[ltree.st[il].gn[jl].p[kl].iup[2][bitflip(m)]].i_stree; //other treeout of point to left in stree
							ltree.st[il].gn[jl].p[kl].sr[m] = stree.st[ils].gn[jl].p[kl].ar*(ltree.st[ih].gn[jh].p[k].al / (stree.st[is].gn[jh].p[k].al + stree.st[isother].gn[jh].p[k].al))
								*(1.0 - stree.st[is].gn[jh].p[k].al / (stree.st[is].gn[jh].p[k].al + stree.st[isother].gn[jh].p[k].al)); //sr for left point changed
                            ltree.st[il].gn[jl].p[kl].sul[m] = ltree.st[ih].gn[jh].p[k].sl;
							ltree.st[ih].gn[jh].p[k].sdr = ltree.st[il].gn[jl].p[kl].sr[m];  //value of sdr for this point also changes
                        }
                        else   //sr for other branch also changed, since total cross section is same, only splitting has changed
                        {
                            iother = ltree.st[il].gn[jl].p[kl].iup[2][m];   //other treeout of point to left
                            isother = ltree.st[ltree.st[il].gn[jl].p[kl].iup[2][m]].i_stree; //other treeout of point to left in stree
							ltree.st[il].gn[jl].p[kl].sr[m] = -stree.st[ils].gn[jl].p[kl].ar*(ltree.st[ih].gn[jh].p[k].al / (stree.st[is].gn[jh].p[k].al + stree.st[isother].gn[jh].p[k].al))
								*(stree.st[is].gn[jh].p[k].al / (stree.st[is].gn[jh].p[k].al + stree.st[isother].gn[jh].p[k].al)); //sr for left point changed
							ltree.st[iother].gn[jh].p[k].sdr = ltree.st[il].gn[jl].p[kl].sr[m];  //value of sdr for this point also changes
                        }
                    }
                }
				else ltree.st[il].gn[jl].p[kl].sul[0] = ltree.st[ih].gn[jh].p[k].sl;
                if(ir.size() > 1)  //splits at end point
                {
                    for(unsigned m = 0; m < ir.size(); m++)
                    {
                        if(stree.st[is].gn[jh].p[k].iup[2].size() > 1) //used to be a splitting here
                        {
							ltree.st[ih].gn[jh].p[k].sr[m] = (stree.st[is].gn[jh].p[k].sr[m] / (stree.st[is].gn[jh].p[k].sr[m] + stree.st[is].gn[jh].p[k].sr[bitflip(m)]))*ltree.st[ih].gn[jh].p[k].ar; //same splitting ratio as before
                        }
                        else ltree.st[ih].gn[jh].p[k].sr[m] = (1.0/2.0)*ltree.st[ih].gn[jh].p[k].ar; //otherwise even splitting
                        ltree.st[ir[m]].gn[jr[m]].p[kr[m]].sdr = ltree.st[ih].gn[jh].p[k].sr[m]; //sr reference for cell to right
                    }
                }
                else
                {
                    ltree.st[ih].gn[jh].p[k].sr[0] = ltree.st[ih].gn[jh].p[k].ar*ltree.st[ih].gn[jh].Nb;
                    ltree.st[ir[0]].gn[jr[0]].p[kr[0]].sdr = ltree.st[ih].gn[jh].p[k].sr[0];   //sr reference for cell to right
                    //sdr for next point
                }
                if(jh>=ltree.st[ih].StartGen + ltree.st[ih].Ncond)  //if in acinus
                {
                    ltree.st[ih].gn[jh].p[k].DA =  0.5*ltree.st[ih].gn[jh].Ap*stree.st[is].gn[jh].p[k].DA;   //DA scales as sqrt area pert
                    ltree.st[ih].Valv0 += ltree.st[ih].gn[jh].p[k].DA*stree.st[is].gn[jh].dx*ltree.st[ih].gn[jh].Nb;
					ltree.st[ih].Vacinduct += 0.5*(ltree.st[ih].gn[jh].p[k].al + ltree.st[ih].gn[jh].p[k].ar)*stree.st[is].gn[jh].dx*ltree.st[ih].gn[jh].Nb;
                }
				if (jh > 0 || k>0) ltree.Vairways += 0.5*(ltree.st[ih].gn[jh].p[k].al + ltree.st[ih].gn[jh].p[k].ar)*stree.st[is].gn[jh].dx*ltree.st[ih].gn[jh].Nb;
            }
        } break;
            
        case 'L':
        {
            ltree.st[ih].gn[jh].dx = stree.st[is].gn[jh].dx*ltree.st[ih].gn[jh].Lp;   //perturbed length of cells
            double xshift = (ltree.st[ih].gn[jh].p.size() - 1)*(1 + ltree.st[ih].gn[jh].Lp); //shift in x position following this branch
            for (unsigned j = jh; j <= ltree.st[ih].EndGen; j++) //shift all up
            {
                ltree.st[ih].gn[j].x0 = xshift;      //shift all subsequent generations along
            }
            for (unsigned k = 0; k < ltree.st[ih].gn[jh].p.size(); k++)  //change values in dxup
            {
                il = ltree.st[ih].gn[jh].p[k].iup[0][0];  //tree no. left
                ir = ltree.st[ih].gn[jh].p[k].iup[2];     //tree no. right (vector)
                jl = ltree.st[ih].gn[jh].p[k].jup[0][0];  //gen left
                jr = ltree.st[ih].gn[jh].p[k].jup[2];     //gen right (vector)
                kl = ltree.st[ih].gn[jh].p[k].kup[0][0];  //branch no. left
                kr = ltree.st[ih].gn[jh].p[k].kup[2];     //branch no. right (vector)
                ltree.st[ih].gn[jh].p[k].dxup[1][0] =  ltree.st[ih].gn[jh].dx;  //dx of this point is changed
                if(k<ltree.st[ih].gn[jh].p.size()-1) ltree.st[ih].gn[jh].p[k].dxup[2][0] =  ltree.st[ih].gn[jh].dx;  //if not final point, then dx of point to right is changed
                else   //if final point, then dxup[0] for point to right is changed
                {
                    for(unsigned m = 0; m < ltree.st[ih].gn[jh].p[k].iup[2].size(); m++)
                    {
                        ltree.st[ir[m]].gn[jr[m]].p[kr[m]].dxup[0][0] = ltree.st[ih].gn[jh].dx;
                    }
                }
                if(k>0) ltree.st[ih].gn[jh].p[k].dxup[0][0] =  ltree.st[ih].gn[jh].dx;  //if not first point, dxup[0] is changed
                else   //if first point, dxup[2] of point to left is changed
                {
                    for(unsigned m = 0; m < ltree.st[il].gn[jl].p[kl].iup[2].size(); m++)
                    {
                        if(ltree.st[il].gn[jl].p[kl].iup[2][m] == ih) ltree.st[il].gn[jl].p[kl].dxup[2][m] = ltree.st[ih].gn[jh].p[k].dxup[1][0];
                    }
                }
                if(jh>=ltree.st[ih].StartGen + ltree.st[ih].Ncond)  //if in acinus
                {
                    ltree.st[ih].Valv0 += stree.st[is].gn[jh].p[k].DA*ltree.st[ih].gn[jh].dx*ltree.st[ih].gn[jh].Nb;
					ltree.st[ih].Vacinduct += 0.5*(stree.st[is].gn[jh].p[k].al + stree.st[is].gn[jh].p[k].ar)*ltree.st[ih].gn[jh].dx*ltree.st[ih].gn[jh].Nb;
                }
				if (jh > 0 || k>0) ltree.Vairways += 0.5*(stree.st[is].gn[jh].p[k].al + stree.st[is].gn[jh].p[k].ar)*ltree.st[ih].gn[jh].dx*ltree.st[ih].gn[jh].Nb;
            }
        } break;
            
        default:
            break;    //Ep stores change in K
    }
}
