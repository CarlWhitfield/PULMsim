//
//  build_lobe_models.cpp
//
//
//  Created by Carl Whitfield on 03/01/2018.
//
//  Commented and checked 16/05/18

#include "lung_model_discrete_branching.hpp"

double build_symmetric_branches(Tree &stree, Options &o)  //build basic symmetric tree model
{
	stree.st.resize(1);  //just one subtree
	stree.st[0].treein = -1;
	stree.meanpath.push_back(0);
	stree.st[0].Ncond = o.Ngen;
	stree.st[0].Ntot = o.Ngen2;
	stree.st[0].StartBranch = 0;
	stree.st[0].z0 = 0;   //position for output
	stree.st[0].imeanpath = 0;    //number of meanpath

	if (o.ShapeOp == WEIBEL)
	{
		//TODO fix DS and things in weibel
		stree.st[0].L0 = 0.12;   //in Weibel, length of trachea is 12cm
		stree.st[0].A0 = 0.000254;   //CS area in 2.54cm^2
		o.VD = 0.06;   //need to update
		o.V0 = o.VFRC - 0.5;  //need to update
	}
	else
	{
		o.V0 = o.VFRC - o.VD - o.VDM;   //volume of acini
		stree.st[0].L0 = pow((4.0*o.VD*o.LDratio*o.LDratio*(1 - 2 * o.lambda*o.lambda*o.lambda) / (M_PI*(1 - pow(2 * o.lambda*o.lambda*o.lambda, o.Ngen)))), 1.0 / 3.0);
		stree.st[0].A0 = 0.25*M_PI*stree.st[0].L0*stree.st[0].L0/(o.LDratio*o.LDratio);
        stree.st[0].L1 = pow((4.0*o.V0*o.Vductvol*o.LDratio2*o.LDratio2*(1 - 2 * o.lambda2*o.lambda2*o.lambda2) / (pow(2.0,o.Ngen)*M_PI*(1 - pow(2 * o.lambda2*o.lambda2*o.lambda2, o.Ngen2 - o.Ngen + 1)))), 1.0 / 3.0);
        stree.st[0].A1 = 0.25*M_PI*stree.st[0].L1*stree.st[0].L1/(o.LDratio2*o.LDratio2);
		
	}
//    double Vlungcheck = o.VDM;
//    for(unsigned j=0; j<o.Ngen; j++)
//    {
//        Vlungcheck += pow(2.0*o.lambda*o.lambda*o.lambda,j)*stree.st[0].L0*stree.st[0].A0;
//    }
//    for(unsigned j=o.Ngen; j<=o.Ngen2; j++)
//    {
//        Vlungcheck += pow(2.0,j)*pow(o.lambda2,3.0*(j-o.Ngen))*stree.st[0].L1*stree.st[0].A1;
//    }
	stree.st[0].Ncond = o.Ngen;   //final gen of conducting tree
	stree.st[0].Ntot = o.Ngen2;   //final gen of acini
	stree.st[0].StartGen = 0;
	stree.st[0].EndGen = stree.st[0].Ntot;
	stree.st[0].allocate_gn(stree.st[0].Ntot + 1);
    stree.st[0].V0 = (1-o.Vductvol)*o.V0;
    o.V0 = stree.st[0].V0;
//    Vlungcheck += o.V0;
	double LL = 0;
	for (unsigned j = 0; j <= stree.st[0].Ntot; j++) LL += stree.st[0].length_function(j, o);

	stree.mpsubtrees.resize(1);
	stree.mpsubtrees[0].push_back(0);

	return LL;
}

double build_alternative_lobe_branches(Tree &stree, Options &o) //build lobar regional model 
{
	double VDtop, LL;
//    double Lscale = pow(o.VD / 0.00012, 1.0/3.0);   //branch size scaled by 120ml DS
	stree.st.resize(13);             //baseline number of strees
	stree.mpsubtrees.resize(13);

	//i=0 trachea
	stree.st[0].treein = -1;        //add connections to trees -1 if no subtree
	stree.st[0].treeout[0] = BR;
	stree.st[0].treeout[1] = BL;
	stree.st[0].z0 = 0;   //position for output
	stree.st[BR].z0 = stree.st[0].z0 - 0.5;   //left and right branch positions
	stree.st[BL].z0 = stree.st[0].z0 + 0.5;
	stree.st[0].StartGen = 0;
	stree.st[0].EndGen = 0;
	stree.st[0].allocate_gn(1);
	stree.st[0].L0 = 10.0*cm_m;// *Lscale;
	stree.st[0].A0 = 0.25*M_PI*1.6*1.6*cm_m*cm_m;// *Lscale*Lscale;

	//right primary bronchiole
	stree.st[BR].treein = 0;
	stree.st[BR].treeout[0] = BRU;
	stree.st[BR].treeout[1] = BRML;
	stree.st[BRU].z0 = stree.st[BR].z0 - 0.25;   //left and right branch numbers
	stree.st[BRML].z0 = stree.st[BR].z0 + 0.25;
	stree.st[BR].StartGen = 1;
	stree.st[BR].EndGen = 1;
	stree.st[BR].L0 = 2.2*cm_m;// *Lscale;
	stree.st[BR].A0 = 0.25*M_PI*1.11*1.11*cm_m*cm_m;// *Lscale*Lscale;
	stree.st[BR].allocate_gn(2);

	//left primary bronchiole
	stree.st[BL].treein = 0;
	stree.st[BL].treeout[0] = BLU;
	stree.st[BL].treeout[1] = BLL;
	stree.st[BLL].z0 = stree.st[BL].z0 - 0.25;
	stree.st[BLU].z0 = stree.st[BL].z0 + 0.25;   //left and right branch numbers
	stree.st[BL].StartGen = 1;
	stree.st[BL].EndGen = 1;
	stree.st[BL].allocate_gn(2);
	stree.st[BL].L0 = 5.0*cm_m;// *Lscale;
	stree.st[BL].A0 = 0.25*M_PI*1.2*1.2*cm_m*cm_m;// *Lscale*Lscale;

	//right intermediate bronchiole
	stree.st[BRML].treein = BR;
	stree.st[BRML].treeout[0] = BRM;
	stree.st[BRML].treeout[1] = BRL;
	stree.st[BRM].z0 = stree.st[BRML].z0 - 0.125;
	stree.st[BRL].z0 = stree.st[BRML].z0 + 0.125;   //left and right branch numbers
	stree.st[BRML].StartGen = 2;
	stree.st[BRML].EndGen = 2;
	stree.st[BRML].allocate_gn(3);
	stree.st[BRML].L0 = 2.6*cm_m;// *Lscale;
	stree.st[BRML].A0 = 0.25*M_PI*0.89*0.89*cm_m*cm_m;// *Lscale*Lscale;

	//right lower lobar bronchiole
	stree.st[BRL].treein = BRML;
	stree.st[BRL].treeout[0] = BRL1;
	stree.st[BRL].treeout[1] = BRL2;
	stree.st[BRL1].z0 = stree.st[BRL].z0 - 0.0625;
	stree.st[BRL2].z0 = stree.st[BRL].z0 + 0.0625;   //left and right branch numbers
	stree.st[BRL].StartGen = 3;
	stree.st[BRL].EndGen = 3;
	stree.st[BRL].allocate_gn(4);
	stree.st[BRL].L0 = 0.8*cm_m;// *Lscale;
	stree.st[BRL].A0 = 0.25*M_PI*0.64*0.64*cm_m*cm_m;// *Lscale*Lscale;

	//left lower lobar bronchiole
	stree.st[BLL].treein = BL;
	stree.st[BLL].treeout[0] = BLL1;
	stree.st[BLL].treeout[1] = BLL2;
	stree.st[BLL1].z0 = stree.st[BLL].z0 - 0.125;
	stree.st[BLL2].z0 = stree.st[BLL].z0 + 0.125;   //left and right branch numbers
	stree.st[BLL].StartGen = 2;
	stree.st[BLL].EndGen = 2;
	stree.st[BLL].allocate_gn(3);
	stree.st[BLL].L0 = 1.1*cm_m;// *Lscale;
	stree.st[BLL].A0 = 0.25*M_PI*0.8*0.8*cm_m*cm_m;// *Lscale*Lscale;

	unsigned long Nacintot=0;
	for (unsigned i = 0; i < stree.st.size(); i++)
	{
		stree.st[i].imeanpath = i;    //store mean paths
		stree.st[i].StartBranch = 0;   //all start branches in mean paths 0
		stree.st[i].Ntot = lobe_gens(i, o);
		stree.st[i].Ncond = stree.st[i].Ntot - (o.Ngen2 - o.Ngen);
	}

	VDtop = stree.st[BRL].L0*stree.st[BRL].A0 + stree.st[BLL].L0*stree.st[BLL].A0 
		  + stree.st[BRML].L0*stree.st[BRML].A0 + stree.st[BL].L0*stree.st[BL].A0 
		  + stree.st[BR].L0*stree.st[BR].A0 + stree.st[0].L0*stree.st[0].A0;   //proximal airway vol in m^3
	LL = stree.st[BRL].L0 + stree.st[BRML].L0 + stree.st[BR].L0 + stree.st[0].L0;     //stores length of BRL major path 

	//create list of mean path models
	stree.meanpath.push_back(BRU);
	stree.st[BRU].treein = BR;
	stree.meanpath.push_back(BRM);
	stree.st[BRM].treein = BRML;
	stree.meanpath.push_back(BRL1);
	stree.st[BRL1].treein = BRL;
	stree.meanpath.push_back(BRL2);
	stree.st[BRL2].treein = BRL;
	stree.meanpath.push_back(BLU);
	stree.st[BLU].treein = BL;
	stree.meanpath.push_back(BLL1);
	stree.st[BLL1].treein = BLL;
	stree.meanpath.push_back(BLL2);
	stree.st[BLL2].treein = BLL;
	Nacintot = 0;

	double VDtot = 0;
	//calculate fit for VD and total number of acini
	for (unsigned m = 0; m < stree.meanpath.size(); m++)
	{
		unsigned im = stree.meanpath[m];
		stree.st[im].StartGen = stree.st[stree.st[im].treein].EndGen + 1;
		stree.st[im].EndGen = stree.st[im].StartGen + stree.st[im].Ntot;
		stree.st[im].StartBranch = 0;
		VDtot += pow(o.lambda*o.lambda*o.lambda, o.Ngen - 2 - stree.st[im].Ncond)
			     *(1.0 - pow((2 * o.lambda*o.lambda*o.lambda), stree.st[im].Ncond)) 
				 /(1.0 - 2 * o.lambda*o.lambda*o.lambda);
		stree.st[im].allocate_gn(stree.st[im].EndGen + 1);
		Nacintot += ((unsigned long)pow(2, stree.st[im].Ncond));
	}
	o.V0 = o.VFRC - o.VD - o.VDM;  //V0 is currently volume of acini
	for (unsigned m = 0; m < stree.meanpath.size(); m++)
	{
		unsigned im = stree.meanpath[m];
		stree.mpsubtrees[im].push_back(im);
		//length and areas fitted to VD and Vacinduct volumes
        stree.st[im].L0 = pow(o.lambda, o.Ngen - 2 - stree.st[im].Ncond)*pow((4.0*o.LDratio*o.LDratio*(o.VD-VDtop) / (M_PI*VDtot)), 1.0 / 3.0);
		stree.st[im].A0 = 0.25*M_PI*stree.st[im].L0*stree.st[im].L0 / (o.LDratio*o.LDratio);      
		stree.st[im].L1 = pow((4.0 / M_PI)*o.LDratio2*o.LDratio2*(o.V0*o.Vductvol/ Nacintot)*((1 - 2 * o.lambda2*o.lambda2*o.lambda2) / (1 - pow(2 * o.lambda2*o.lambda2*o.lambda2, stree.st[im].Ntot - stree.st[im].Ncond + 1))), 1.0 / 3.0);
		stree.st[im].A1 = 0.25*M_PI*stree.st[im].L1*stree.st[im].L1 / (o.LDratio2*o.LDratio2);      
		if (im == BRL1)   //calc lung length scale LL
		{
			LL += stree.st[im].L0*(1 - pow(o.lambda, stree.st[im].Ncond)) / (1 - o.lambda);
			LL += stree.st[im].L1*(1 - pow(o.lambda2, stree.st[im].Ntot - stree.st[im].Ncond+1)) / (1 - o.lambda2);
		}
		stree.st[im].imeanpath = im;
	}
	o.V0 = (1.0 - o.Vductvol)*(o.VFRC - o.VD - o.VDM);  //V0 updated to volume of acinar sacs
	//initial sac volume distribution
	stree.st[BRL1].V0 = 0.2*o.V0;
	stree.st[BRL2].V0 = 0.05*o.V0;
	stree.st[BLL1].V0 = 0.2*o.V0;
	stree.st[BLL2].V0 = 0.05*o.V0;
	stree.st[BRM].V0 = 0.1*o.V0;
	stree.st[BRU].V0 = 0.2*o.V0;
	stree.st[BLU].V0 = 0.2*o.V0;

	//All of these trees count as mean paths
	stree.meanpath.insert(stree.meanpath.begin(), BRL);
	stree.mpsubtrees[BRL].push_back(BRL);
	stree.st[BRL].V0 = stree.st[BRL1].V0 + stree.st[BRL2].V0;

	stree.meanpath.insert(stree.meanpath.begin(), BLL);
	stree.mpsubtrees[BLL].push_back(BLL);
	stree.st[BLL].V0 = stree.st[BLL1].V0 + stree.st[BLL2].V0;

	stree.meanpath.insert(stree.meanpath.begin(), BRML);
	stree.mpsubtrees[BRML].push_back(BRML);
	stree.st[BRML].V0 = stree.st[BRM].V0 + stree.st[BRL].V0;

	stree.meanpath.insert(stree.meanpath.begin(), BR);
	stree.mpsubtrees[BR].push_back(BR);
	stree.st[BR].V0 = stree.st[BRML].V0 + stree.st[BRU].V0;

	stree.meanpath.insert(stree.meanpath.begin(), BL);
	stree.mpsubtrees[BL].push_back(BL);
	stree.st[BL].V0 = stree.st[BLL].V0 + stree.st[BLU].V0;

	stree.meanpath.insert(stree.meanpath.begin(), 0);
	stree.mpsubtrees[0].push_back(0);
	stree.st[0].V0 = stree.st[BL].V0 + stree.st[BR].V0;

	return LL;
}
