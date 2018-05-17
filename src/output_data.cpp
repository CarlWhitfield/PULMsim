//
//  output_data.cpp
//  
//
//  Created by Carl Whitfield on 10/03/2017.
//
//  Commented and checked 16/05/2018

#include "lung_model_discrete_branching.hpp"

int printfunc(Tree &stree, vector<Tree> &ltree, Options &o, Conversions cons, double t)
{
	int error = 0;
	{
		stringstream ss;
		string filename;
		
		//options for conc output format
		if (o.OutputOp == VTK || o.OutputOp == VTK_CSV)
		{
			ss.str("");
			ss << "conc_" << o.SimID << "_stree_" << t << ".vtk";
			filename=ss.str();
			error += stree.printfunc_tree(filename, cons);
		}
		if (o.OutputOp == CSV || o.OutputOp == VTK_CSV)
		{
			ss.str("");
			ss << "conc_" << o.SimID << "_stree_" << t << ".csv";
			filename=ss.str();
			error += stree.printfunc_tree_csv(filename, cons);
		}
	}


	if (o.output_perts)   //output perturbation conc files
	{
		for (unsigned n=0; n < ltree.size(); n++)
		{
			stringstream ss;
			string filename;
			if (o.OutputOp == VTK || o.OutputOp == VTK_CSV)
			{
				ss.str("");
				ss << "conc_" << o.SimID << "_ltree_" << ltree[n].pert_type << "_"
				<< read_lobe_no(ltree[n].st[ltree[n].pert_ij[0]].i_stree,o) << "_j" << (ltree[n].pert_ij[1] - ltree[n].st[ltree[n].st[ltree[n].pert_ij[0]].i_stree].StartGen) << '_' << t << ".vtk";
				filename=ss.str();
				error += ltree[n].printfunc_lp_tree(filename, stree, cons);
			}
			if (o.OutputOp == CSV || o.OutputOp == VTK_CSV)
			{
				ss.str("");
				ss << "conc_" << o.SimID << "_ltree_" << ltree[n].pert_type << "_"
				<< read_lobe_no(ltree[n].st[ltree[n].pert_ij[0]].i_stree,o) << "_j" << (ltree[n].pert_ij[1] - ltree[n].st[ltree[n].st[ltree[n].pert_ij[0]].i_stree].StartGen) << '_' << t << ".csv";
				//cout << ss.str().c_str() << '\n';
				filename=ss.str();
				error += ltree[n].printfunc_lp_tree_csv(filename, stree, cons);
				
			}
		}
	}
	return error;
}

int Tree::printfunc_tree(string filename, Conversions cons) //print stree snapshot (VTK)
{
	unsigned long Nx = kmtot;  //total number of cells in this tree
	ofstream output;
	
	output.open(filename.c_str(),std::ofstream::out);
	
	if(!output.good())    //check file
	{
		cout << "Could not open output file.\n";
		return 1;
	}
	
	output << fixed << setprecision(OUTPUT_PRECISION_CONC);
	output <<  "# vtk DataFile Version 3.0\nConc_data\nASCII\nDATASET UNSTRUCTURED_GRID\n";
	output <<  "POINTS " << Nx + 1 << " float\n";
	
	output << (st[0].gn[0].p[0].x - 0.5*st[0].gn[0].dx)*cons.LL_to_cm << " 0 0\n";   //0 point
	for (unsigned i = 0; i < st.size(); i++)
	{
		for (unsigned j = st[i].StartGen; j <= st[i].EndGen; j++)
		{  
			for (unsigned k = 0; k < st[i].gn[j].p.size(); k++)
			{
				output << (st[i].gn[j].p[k].x + 0.5*st[i].gn[j].dx)*cons.LL_to_cm << ' ' << st[i].gn[j].p[k].y0 << ' ' << st[i].z0*cons.LL_to_cm << '\n';
			}
		}
	}

	output <<  "\n\nCELLS " << Nx << ' ' << 3*Nx << '\n';
	
	output <<  "2 0 1\n";
	for (unsigned i = 0; i < st.size(); i++)
	{
		for (unsigned j = st[i].StartGen; j <= st[i].EndGen; j++)
		{
			unsigned k0 = 0;
			if(j==0) k0 = 1;
			for (unsigned k = k0; k < st[i].gn[j].p.size(); k++)
			{
				unsigned il = st[i].gn[j].p[k].iup[0][0];
				unsigned jl = st[i].gn[j].p[k].jup[0][0];
				unsigned kl = st[i].gn[j].p[k].kup[0][0];
				output << "2 " << st[il].gn[jl].p[kl].km + 1 << ' ' << st[i].gn[j].p[k].km + 1 << '\n';
			}
		}
	}

	output <<  "\n\nCELL_TYPES " << Nx << '\n';
	
	for (unsigned k = 0; k<Nx; k++) output << "3\n";
	
	output <<  "\n\nPOINT_DATA " << Nx + 1 << "\n";
	
	output <<  "\nSCALARS Velocity_cms-1 float\nLOOKUP_TABLE default\n";

	output << st[0].gn[0].p[0].ul*cons.LL_to_cm / cons.t_to_s << '\n';
	for (unsigned i = 0; i < st.size(); i++)
	{
		for (unsigned j = st[i].StartGen; j <= st[i].EndGen; j++)
		{
			for (unsigned k = 0; k < st[i].gn[j].p.size(); k++)
			{
				output << st[i].gn[j].p[k].ur*cons.LL_to_cm/cons.t_to_s << '\n';
			}
		}
	}

	output <<  "\nSCALARS Diffusion_Constant_cm2s-1 float\nLOOKUP_TABLE default\n";
	
	output << st[0].gn[0].p[0].Dl*cons.LL_to_cm*cons.LL_to_cm / cons.t_to_s << '\n';
	for (unsigned i = 0; i < st.size(); i++)
	{
		for (unsigned j = st[i].StartGen; j <= st[i].EndGen; j++)
		{
			for (unsigned k = 0; k < st[i].gn[j].p.size(); k++)
			{
				output << st[i].gn[j].p[k].Dr*cons.LL_to_cm*cons.LL_to_cm/cons.t_to_s << '\n';
			}
		}
	}
	
	output <<  "\n\nCELL_DATA " << Nx << "\n\n";

	output <<  "\nSCALARS Concentration float\nLOOKUP_TABLE default\n";
	
	for (unsigned i = 0; i < st.size(); i++)
	{
		for (unsigned j = st[i].StartGen; j <= st[i].EndGen; j++)
		{
			for (unsigned k = 0; k < st[i].gn[j].p.size(); k++)
			{
				output << 0.5*(st[i].gn[j].p[k].c + st[i].gn[j].p[k].cold) << '\n';
			}
		}
	}
		
	output <<  "\nSCALARS Total_Outer_CS_Area_cm2 float\nLOOKUP_TABLE default\n";
	
	for (unsigned i = 0; i < st.size(); i++)
	{
		for (unsigned j = st[i].StartGen; j <= st[i].EndGen; j++)
		{
			for (unsigned k = 0; k < st[i].gn[j].p.size(); k++)
			{
				output << 0.5*st[i].gn[j].Nb*(st[i].gn[j].p[k].Anew + st[i].gn[j].p[k].Aold)*cons.LL_to_cm*cons.LL_to_cm << '\n';
			}
		}
	}
	
	output <<  "\nSCALARS Branch_Outer_CS_Area_cm2 float\nLOOKUP_TABLE default\n";
	
	for (unsigned i = 0; i < st.size(); i++)
	{
		for (unsigned j = st[i].StartGen; j <= st[i].EndGen; j++)
		{
			for (unsigned k = 0; k < st[i].gn[j].p.size(); k++)
			{
				output << 0.5*(st[i].gn[j].p[k].Anew + st[i].gn[j].p[k].Aold)*cons.LL_to_cm*cons.LL_to_cm << '\n';
			}
		}
	}

	output << "\nSCALARS Branch_Inner_CS_Area_cm2 float\nLOOKUP_TABLE default\n";

	for (unsigned i = 0; i < st.size(); i++)
	{
		for (unsigned j = st[i].StartGen; j <= st[i].EndGen; j++)
		{
			for (unsigned k = 0; k < st[i].gn[j].p.size(); k++)
			{
				output << 0.5*(st[i].gn[j].p[k].al + st[i].gn[j].p[k].ar)*cons.LL_to_cm*cons.LL_to_cm << '\n';
			}
		}
	}
	output.close();
	return 0;
}

int Tree::printfunc_lp_tree(string filename, Tree &stree, Conversions cons)   //print ltree snapshot (VTK)
{
	unsigned long Nx = kmtot;  //total number of cells in this tree
	ofstream output;
	
	output.open(filename.c_str(),std::ofstream::out);
	
	if(!output.good())    //check file
	{
		cout << "Could not open output file.\n";
		return 1;
	}
	
	output << fixed << setprecision(OUTPUT_PRECISION_CONC);
	output <<  "# vtk DataFile Version 3.0\nConc_data\nASCII\nDATASET UNSTRUCTURED_GRID\n";
	output <<  "POINTS " << Nx + 1 << " float\n";
	
	output <<  stree.st[0].gn[0].p[0].x - 0.5*stree.st[0].gn[0].dx << " 0 0\n";   //0 point
	for (unsigned i = 0; i < st.size(); i++)
	{
		unsigned istree = st[i].i_stree;
		for (unsigned j = st[i].StartGen; j <= st[i].EndGen; j++)
		{
			for (unsigned k = 0; k < st[i].gn[j].p.size(); k++)
			{
				output << (stree.st[istree].gn[j].p[k].x + 0.5*stree.st[istree].gn[j].dx)*cons.LL_to_cm << ' ' << stree.st[istree].gn[j].p[k].y0 << ' ' << st[i].z0*cons.LL_to_cm << '\n';
			}
		}
	}
	
	output <<  "\n\nCELLS " << Nx << ' ' << 3*Nx << '\n';
	
	output <<  "2 0 1\n";
	for (unsigned i = 0; i < st.size(); i++)
	{
		for (unsigned j = st[i].StartGen; j <= st[i].EndGen; j++)
		{
			unsigned k0 = 0;
			if(j==0) k0 = 1;
			for (unsigned k = k0; k < st[i].gn[j].p.size(); k++)
			{
				unsigned il = st[i].gn[j].p[k].iup[0][0];
				unsigned jl = st[i].gn[j].p[k].jup[0][0];
				unsigned kl = st[i].gn[j].p[k].kup[0][0];
				output << "2 " << st[il].gn[jl].p[kl].km + 1 << ' ' << st[i].gn[j].p[k].km + 1 << '\n';
			}
		}
	}
	
	output <<  "\n\nCELL_TYPES " << Nx << '\n';
	
	for (unsigned k = 0; k<Nx; k++) output << "3\n";
	
	output <<  "\n\nPOINT_DATA " << Nx + 1 << "\n";
	
	output <<  "\nSCALARS Velocity0_cms-1 float\nLOOKUP_TABLE default\n";
	
	output <<  st[0].gn[0].p[0].ul << '\n';
	for (unsigned i = 0; i < st.size(); i++)
	{
		unsigned istree = st[i].i_stree;
		for (unsigned j = st[i].StartGen; j <= st[i].EndGen; j++)
		{
			for (unsigned k = 0; k < st[i].gn[j].p.size(); k++)
			{
				output << stree.st[istree].gn[j].p[k].ur*cons.LL_to_cm/cons.t_to_s << '\n';
			}
		}
	}
	
	output <<  "\nSCALARS Delta_Velocity_cms-1 float\nLOOKUP_TABLE default\n";
	
	output <<  st[0].gn[0].p[0].ul << '\n';
	for (unsigned i = 0; i < st.size(); i++)
	{
		for (unsigned j = st[i].StartGen; j <= st[i].EndGen; j++)
		{
			for (unsigned k = 0; k < st[i].gn[j].p.size(); k++)
			{
				output << st[i].gn[j].p[k].ur*cons.LL_to_cm/cons.t_to_s << '\n';
			}
		}
	}
	
	output <<  "\nSCALARS Diffusion_constant0_cm2s-1 float\nLOOKUP_TABLE default\n";
	
	output << st[0].gn[0].p[0].Dl << '\n';
	for (unsigned i = 0; i < st.size(); i++)
	{
		unsigned istree = st[i].i_stree;
		for (unsigned j = st[i].StartGen; j <= st[i].EndGen; j++)
		{
			for (unsigned k = 0; k < st[i].gn[j].p.size(); k++)
			{
				output << stree.st[istree].gn[j].p[k].Dr*cons.LL_to_cm*cons.LL_to_cm/cons.t_to_s << '\n';
			}
		}
	}
	
	output <<  "\nSCALARS Delta_Diffusion_Constant_cm2s-1 float\nLOOKUP_TABLE default\n";
	
	output << st[0].gn[0].p[0].Dl << '\n';
	for (unsigned i = 0; i < st.size(); i++)
	{
		for (unsigned j = st[i].StartGen; j <= st[i].EndGen; j++)
		{
			for (unsigned k = 0; k < st[i].gn[j].p.size(); k++)
			{
				output << st[i].gn[j].p[k].Dr*cons.LL_to_cm*cons.LL_to_cm/cons.t_to_s << '\n';
			}
		}
	}
	
	output <<  "\n\nCELL_DATA " << Nx << "\n\n";
	
	output <<  "\nSCALARS Concentration0 float\nLOOKUP_TABLE default\n";
	
	for (unsigned i = 0; i < st.size(); i++)
	{
		unsigned istree = st[i].i_stree;
		for (unsigned j = st[i].StartGen; j <= st[i].EndGen; j++)
		{
			for (unsigned k = 0; k < st[i].gn[j].p.size(); k++)
			{
				output << 0.5*(stree.st[istree].gn[j].p[k].c + stree.st[istree].gn[j].p[k].cold) << '\n';
			}
		}
	}
	
	output <<  "\nSCALARS Delta_Concentration float\nLOOKUP_TABLE default\n";
	
	for (unsigned i = 0; i < st.size(); i++)
	{
		for (unsigned j = st[i].StartGen; j <= st[i].EndGen; j++)
		{
			for (unsigned k = 0; k < st[i].gn[j].p.size(); k++)
			{
				output << 0.5*(st[i].gn[j].p[k].c + st[i].gn[j].p[k].cold) << '\n';
			}
		}
	}
	
	output <<  "\nSCALARS Tot_Outer_CS_Area0_cm2 float\nLOOKUP_TABLE default\n";
	
	for (unsigned i = 0; i < st.size(); i++)
	{
		unsigned istree = st[i].i_stree;
		for (unsigned j = st[i].StartGen; j <= st[i].EndGen; j++)
		{
			for (unsigned k = 0; k < st[i].gn[j].p.size(); k++)
			{
				output << 0.5*st[i].gn[j].Nb*(stree.st[istree].gn[j].p[k].Anew + stree.st[istree].gn[j].p[k].Aold)*cons.LL_to_cm*cons.LL_to_cm << '\n';
			}
		}
	}
	
	output <<  "\nSCALARS Delta_Tot_Outer_CS_Area_cm2 float\nLOOKUP_TABLE default\n";
	
	for (unsigned i = 0; i < st.size(); i++)
	{
		for (unsigned j = st[i].StartGen; j <= st[i].EndGen; j++)
		{
			for (unsigned k = 0; k < st[i].gn[j].p.size(); k++)
			{
				output << 0.5*st[i].gn[j].Nb*(st[i].gn[j].p[k].Anew + st[i].gn[j].p[k].Aold)*cons.LL_to_cm*cons.LL_to_cm << '\n';
			}
		}
	}
	
	output <<  "\nSCALARS Branch_Outer_CS_Area0_cm2 float\nLOOKUP_TABLE default\n";
	
	for (unsigned i = 0; i < st.size(); i++)
	{
		unsigned istree = st[i].i_stree;
		for (unsigned j = st[i].StartGen; j <= st[i].EndGen; j++)
		{
			for (unsigned k = 0; k < st[i].gn[j].p.size(); k++)
			{
				output << 0.5*(stree.st[istree].gn[j].p[k].Anew + stree.st[istree].gn[j].p[k].Aold)*cons.LL_to_cm*cons.LL_to_cm << '\n';
			}
		}
	}
	
	output <<  "\nSCALARS Delta_Branch_Outer_CS_Area_cm2 float\nLOOKUP_TABLE default\n";
	
	for (unsigned i = 0; i < st.size(); i++)
	{
		for (unsigned j = st[i].StartGen; j <= st[i].EndGen; j++)
		{
			for (unsigned k = 0; k < st[i].gn[j].p.size(); k++)
			{
				output << 0.5*(st[i].gn[j].p[k].Anew + st[i].gn[j].p[k].Aold)*cons.LL_to_cm*cons.LL_to_cm << '\n';
			}
		}
	}
	
	output << "\nSCALARS Branch_Inner_CS_Area0_cm2 float\nLOOKUP_TABLE default\n";
	
	for (unsigned i = 0; i < st.size(); i++)
	{
		unsigned istree = st[i].i_stree;
		for (unsigned j = st[i].StartGen; j <= st[i].EndGen; j++)
		{
			for (unsigned k = 0; k < st[i].gn[j].p.size(); k++)
			{
				output << 0.5*(stree.st[istree].gn[j].p[k].al + stree.st[istree].gn[j].p[k].ar)*cons.LL_to_cm*cons.LL_to_cm << '\n';
			}
		}
	}
	
	output << "\nSCALARS Delta_Branch_Inner_CS_Area_cm2 float\nLOOKUP_TABLE default\n";
	
	for (unsigned i = 0; i < st.size(); i++)
	{
		for (unsigned j = st[i].StartGen; j <= st[i].EndGen; j++)
		{
			for (unsigned k = 0; k < st[i].gn[j].p.size(); k++)
			{
				output << 0.5*(st[i].gn[j].p[k].al + st[i].gn[j].p[k].ar)*cons.LL_to_cm*cons.LL_to_cm << '\n';
			}
		}
	}
	output.close();
	return 0;
}

int Tree::printfunc_tree_csv(string filename, Conversions cons)  //print stree content (csv)
{
	ofstream output;
	output.open(filename.c_str(),std::ofstream::out);
	
	if(!output.good())    //check file
	{
		cout << "Could not open output file.\n";
		return 1;
	}
	
	output << fixed << setprecision(OUTPUT_PRECISION_CONC);
	output << "X_cm , Y_gen , Z_na , Velocity_cms-1 , Diffusion_Constant_cm2s-1, Concentration, "
		   << "Total_Outer_CS_Area_cm2 , Branch_Outer_CS_Area_cm2, Branch_Inner_CS_Area_cm2\n";

	for (unsigned i = 0; i < st.size(); i++)
	{
		for (unsigned j = st[i].StartGen; j <= st[i].EndGen; j++)
		{
			for (unsigned k = 0; k < st[i].gn[j].p.size(); k++)
			{
				output << (st[i].gn[j].p[k].x)*cons.LL_to_cm << " , " << st[i].gn[j].p[k].y0 << " , " << st[i].z0*cons.LL_to_cm << " , "
					<< 0.5*(st[i].gn[j].p[k].ul + st[i].gn[j].p[k].ur)*cons.LL_to_cm / cons.t_to_s << " , "
					<< 0.5*(st[i].gn[j].p[k].Dl + st[i].gn[j].p[k].Dr)*cons.LL_to_cm*cons.LL_to_cm / cons.t_to_s << " , " 
					<< 0.5*(st[i].gn[j].p[k].c + st[i].gn[j].p[k].cold) << " , " 
					<< 0.5*st[i].gn[j].Nb*(st[i].gn[j].p[k].Anew + st[i].gn[j].p[k].Aold)*cons.LL_to_cm*cons.LL_to_cm << " , "
					<< 0.5*(st[i].gn[j].p[k].Anew + st[i].gn[j].p[k].Aold)*cons.LL_to_cm*cons.LL_to_cm << " , "
					<< 0.5*(st[i].gn[j].p[k].al + st[i].gn[j].p[k].ar)*cons.LL_to_cm*cons.LL_to_cm << '\n';
			}
		}
	}
	output.close();
	return 0;
}

int Tree::printfunc_lp_tree_csv(string filename, Tree &stree, Conversions cons)   //print ltree content (csv)
{
	ofstream output;
	output.open(filename.c_str(),std::ofstream::out);
	
	if(!output.good())    //check file
	{
		cout << "Could not open output file.\n";
		return 1;
	}
	
	output << fixed << setprecision(OUTPUT_PRECISION_CONC);

	output << "X_cm , Y_gen , Z_na , Delta_X_cm, Velocity0_cms-1 , Delta_Velocity_cms-1 , Diffusion_Constant0_cm2s-1 , "
		<< "Delta_Diffusion_Constant_cm2s-1 , Concentration0 , Delta_Concentration , Total_Outer_CS_Area0_cm2 , "
		<< "Delta_Total_Outer_CS_Area_cm2 , Branch_Outer_CS_Area0_cm2 , Delta_Branch_Outer_CS_Area_cm2 , "
		<< "Branch_Inner_CS_Area0_cm2 , Delta_Branch_Inner_CS_Area_cm2\n";

	for (unsigned i = 0; i < st.size(); i++)
	{
		unsigned istree = st[i].i_stree;
		for (unsigned j = st[i].StartGen; j <= st[i].EndGen; j++)
		{
			for (unsigned k = 0; k < st[i].gn[j].p.size(); k++)
			{
				output << (stree.st[istree].gn[j].p[k].x)*cons.LL_to_cm << " , " << stree.st[istree].gn[j].p[k].y0 << " , "
					<< st[i].z0*cons.LL_to_cm << " , " << st[i].gn[j].p[k].x*cons.LL_to_cm << " , "
					<< 0.5*(stree.st[istree].gn[j].p[k].ul + stree.st[istree].gn[j].p[k].ur)*cons.LL_to_cm / cons.t_to_s << " , "
					<< 0.5*(st[i].gn[j].p[k].ul + st[i].gn[j].p[k].ur)*cons.LL_to_cm / cons.t_to_s << " , "
					<< 0.5*(stree.st[istree].gn[j].p[k].Dl + stree.st[istree].gn[j].p[k].Dr)*cons.LL_to_cm*cons.LL_to_cm / cons.t_to_s << " , "
					<< 0.5*(st[i].gn[j].p[k].Dl + st[i].gn[j].p[k].Dr)*cons.LL_to_cm*cons.LL_to_cm / cons.t_to_s << " , "
					<< 0.5*(stree.st[istree].gn[j].p[k].c + stree.st[istree].gn[j].p[k].cold) << " , "
					<< 0.5*(st[i].gn[j].p[k].c + st[i].gn[j].p[k].cold) << " , "
					<< 0.5*st[i].gn[j].Nb*(stree.st[istree].gn[j].p[k].Anew + stree.st[istree].gn[j].p[k].Aold)*cons.LL_to_cm*cons.LL_to_cm << " , "
					<< 0.5*st[i].gn[j].Nb*(st[i].gn[j].p[k].Anew + st[i].gn[j].p[k].Aold)*cons.LL_to_cm*cons.LL_to_cm << " , "
					<< 0.5*(stree.st[istree].gn[j].p[k].Anew + stree.st[istree].gn[j].p[k].Aold)*cons.LL_to_cm*cons.LL_to_cm << " , "
					<< 0.5*(st[i].gn[j].p[k].Anew + st[i].gn[j].p[k].Aold)*cons.LL_to_cm*cons.LL_to_cm << " , "
					<< 0.5*(stree.st[istree].gn[j].p[k].al + stree.st[istree].gn[j].p[k].ar)*cons.LL_to_cm*cons.LL_to_cm << " , "
					<< 0.5*(st[i].gn[j].p[k].al + st[i].gn[j].p[k].ar)*cons.LL_to_cm*cons.LL_to_cm << '\n';
			}
		}
	}
	output.close();
	return 0;
}

int append_masterout(string filename, double time, Tree &tree, Options &o, Conversions cons)   //update "flux" file
{	
	ofstream mf;
	mf.open(filename.c_str(),std::ofstream::app);
	//cout << mf.good() << ' ' << mf.eof() << ' ' << mf.fail() <<  ' ' << mf.bad() << '\n';
//	if(!mf.good()) return 1;
	
	mf << setprecision(OUTPUT_PRECISION_FLUX);
	mf << (time + o.dt)*cons.t_to_s << " , " << tree.Vlung*cons.V_to_Litres << " , " << tree.st[0].Vnew*cons.V_to_Litres << " , "
		<< tree.masstot*cons.V_to_Litres << " , " << tree.cmouth << " , " << tree.Ppl*cons.P_to_cmH20 << " , "
		<< (time + 0.5*o.dt)*cons.t_to_s << " , " << (tree.st[0].Vnew - tree.st[0].Vold)*cons.V_to_Litres / (o.dt*cons.t_to_s) << " , "
		<< tree.fluxin*cons.V_to_Litres;
	//old terms used for volumes since volume has been updated
	if (o.N2leak.exists)   //if there is a leak
	{
		if (time > o.StimTime && time - o.StimTime >= o.N2leak.start && time - o.StimTime <= o.N2leak.end)
		{
			if (o.N2leak.type != EXPIRATORY_LEAK && tree.st[0].Vnew - tree.st[0].Vold >= 0)
			{
				mf << " , " << 0 << " , " << (1-o.N2leak.size)*(tree.st[0].Vnew - tree.st[0].Vold)*cons.V_to_Litres / (o.dt*cons.t_to_s);   //if there is a leak, it isn't seen in gas conc
			}
			else
			{
				if (o.N2leak.type != INSPIRATORY_LEAK && tree.st[0].Vnew - tree.st[0].Vold <= 0)
				{
					mf << " , " << tree.cmouth << " , " << (1 - o.N2leak.size)*(tree.st[0].Vnew - tree.st[0].Vold)*cons.V_to_Litres / (o.dt*cons.t_to_s);   //if there is a leak, it isn't seen in gas conc
				}
				else mf << " , " << tree.cmouth << " , " << (tree.st[0].Vnew - tree.st[0].Vold)*cons.V_to_Litres / (o.dt*cons.t_to_s);
			}

		}
		else
		{
			mf << " , " << tree.cmouth << " , " << (tree.st[0].Vnew - tree.st[0].Vold)*cons.V_to_Litres / (o.dt*cons.t_to_s);
		}
	}

	if(o.output_lung_volumes)
	{
		for (unsigned m = 0; m < tree.EndSubtrees.size(); m++)
		{
			mf << " , " << (tree.st[tree.EndSubtrees[m]].Vnew + tree.st[tree.EndSubtrees[m]].Vacinduct)*cons.V_to_Litres;
			mf << " , " << tree.st[tree.EndSubtrees[m]].totIGvol*cons.V_to_Litres;
		}
	}
	mf << '\n';
	mf.close();
	
	return 0;
}

int print_simops(ofstream &summary_file, Options &o)
{
	if(!summary_file.good()) return 1;
	
	summary_file << "Simualtion ID: " << o.SimID << '\n';
	summary_file << "Tree Option: " << o.read_tree_option() << '\n';
	summary_file << "Boundary Condition Option: " << o.read_bc_option() << '\n';
	summary_file << "Upwind Interpolation Option: " << o.read_upwind_option() << '\n';
	summary_file << "Flow Type Option: " << o.read_flow_option() << '\n';
	summary_file << "Elastic Response Function Option: " << o.read_resp_option() << '\n';
	summary_file << "Pressure Option: " << o.read_pressure_option() << '\n';
	summary_file << "Tree Size and Shape Option: " << o.read_shape_option() << '\n';
	summary_file << "Solver Option: " << o.read_solver_option() << '\n';
	summary_file << "Tree Initial Condition: " << o.read_init_option() << '\n';
	summary_file << "Taylor Dispersion: " << o.read_taylor_option() << '\n';
	summary_file << "Output Option: " << o.read_output_option() << '\n';
	summary_file << "Number of breaths simulated: " << o.RunTime << '\n';
	//for(i=0; i<o.Nperts; i++)
	//{
	//	summary_file << "Perturbation generation " << i << ": " << o.GenPert[i] << '\n';
	//	if(o.TreeOp != AMPUTATED)
	//	{
	//		summary_file << "Perturbation size " << i << ": " << o.PertSize[i] << '\n';
	//		summary_file << "Perturbation type " << i << ": " << o.PertType << '\n';
	//	}
	//}
	
	return 0;
}

int print_params(ofstream &summary_file, Options &o)
{
    if(!summary_file.good()) return 1;
    
    summary_file << "Air Viscosity = " << o.Viscosity << '\n';
    summary_file << "Air Density = " << o.Density << '\n';
    summary_file << "Elastance = " << o.E << '\n';

    if(o.ShapeOp != WEIBEL)
    {
        summary_file << "Geometric ratio = " << o.lambda << '\n';
    }
    summary_file << "First breath duration: " << o.Tin << '\n';
    summary_file << "Diffusion constant =  " << o.Diffusion << '\n';
    
	return 0;
}
