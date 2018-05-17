//
//  Lung_model_discrete_branching.cpp
//  
//
//  Created by Carl Whitfield on 09/03/2017.
//
//
#include "lung_model_discrete_branching.hpp"

int main (int argc, char** argv)   //main function
{
	Options o;           //Contains simulation options
	Tree stree;             //Contains unperturbed tree information
	vector<Tree> ltree;   //Contains perturbed tree info
	Conversions cons;    //stores conversion factors
	//--Start simulation timer--//
	o.tstart = chrono::system_clock::now();	
	//------------------------------//
	
	
	//--Direct input file to function//
	if(argc>1)
	{
		string opfile = (string) argv[1];
		if(parse_options(argv[1], o)) return 1;      //error if returns 1
	}
	//------------------------------//
	
	//--set up tree shape and convert physiological units into simulation units--//
    check_option_consistency(o);
	if(initialise_system(stree, cons, o)) return 1;
	//-------------------------------------//

	unsigned long t;
	double time=0;

	//----Setup ltrees----//
	if(o.TreeOp != NOPERT) setup_pert_trees(ltree, stree, o);
	else o.Ntrees = 0;
	
	//-------------------//
	
	//Compute resistance matrices if linear -- only has to be done once//

	stree.AmatQR = new Eigen::HouseholderQR<Eigen::MatrixXd>;
	stree.AmatLU = new Eigen::PartialPivLU<Eigen::MatrixXd>;
	if(o.SolverOp == DIRECT)
	{
		stree.solver_dir = new Eigen::SparseLU<Eigen::SparseMatrix<double>>;
	}
	else
	{
		stree.solver_iter = new Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::IncompleteLUT<double>>;
		stree.solver_iter->setMaxIterations(o.kmtot*o.kmtot);
		stree.solver_iter->setTolerance(TOL);
	}
	if (o.FlowType == POISEUILLE && o.RespFunc == LINEAR_RESP)   //if linear, can build resistance matrix immediately
	{
		if (build_resistance_matrix(stree, o)) return 1;
	}
	for (unsigned n = 0; n < o.Ntrees; n++)
	{
		ltree[n].AmatQR = new Eigen::HouseholderQR<Eigen::MatrixXd>;
		ltree[n].AmatLU = new Eigen::PartialPivLU<Eigen::MatrixXd>;
		if(o.SolverOp == DIRECT)
		{
			ltree[n].solver_dir = new Eigen::SparseLU<Eigen::SparseMatrix<double>>;
		}
		else
		{
			ltree[n].solver_iter = new Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::IncompleteLUT<double>>;
			ltree[n].solver_iter->setMaxIterations(ltree[n].kmtot*ltree[n].kmtot);
			ltree[n].solver_iter->setTolerance(PTOL);
		}
		if (o.FlowType == POISEUILLE && o.RespFunc == LINEAR_RESP)   //if linear, can invert flow equations immediately
		{
			if (build_resistance_matrix(ltree[n], o)) return 1;
		}
	}

	//output to simulation summary file
    ofstream summary_file;
    stringstream fnss;
    fnss.str("");
    fnss << "summary_" << o.SimID << ".txt";
    summary_file.open(fnss.str().c_str(),ofstream::app);
    //--Output Resistances and Volumes--//
    summary_file << "\nCalculated Parameters (Physiological Units)\n";
    summary_file << "Linear Resistance R0 = " << stree.R0*(cons.P_to_cmH20)*(cons.t_to_s)/(cons.LL_to_cm*cons.LL_to_cm*cons.LL_to_cm/1000) << "cmH20sL^(-1)\n";
    summary_file << "Initial alveolar volume: " << o.V0*(cons.LL_to_cm*cons.LL_to_cm*cons.LL_to_cm/1000) << "L\n";

    summary_file << "\nCalculated Parameters (Simulation Units)\n";
    summary_file << "Linear Resistance R0 = " << stree.R0 << '\n';
    summary_file << "Initial alveolar volume: " << o.V0 << '\n';

    //--Print out timestep--//
    unsigned long Nt = ((unsigned long)(o.RunTime / o.dt));
    cout << "Time steps: " << Nt << ' ' << " Step size: " << o.dt << '\n';
    summary_file << "\nTime steps: " << Nt << ' ' << " Step size: " << o.dt << '\n';
    time_t time_start = chrono::system_clock::to_time_t(o.tstart);
    summary_file << "\nSimulation initialised: " << ctime(&time_start) << '\n';
    summary_file.close();

	//---- Setting up master output file for each tree----//
	string *master_filenames;
	master_filenames = new string[o.Ntrees + 1];

	for (unsigned n = 0; n <= o.Ntrees; n++)   
	{
		ofstream masterout;
		stringstream ss;
		ss.str("");

		//build "flux" output files for each tree

		if(n==0) ss << "flux_" << o.SimID << "_stree.csv";
		else ss << "flux_" << o.SimID << "_ltree_" << ltree[n - 1].pert_type << '_' << read_lobe_no(ltree[n-1].st[ltree[n - 1].pert_ij[0]].i_stree,o) << "_j" <<  (ltree[n-1].pert_ij[1] - ltree[n-1].st[ltree[n-1].st[ltree[n-1].pert_ij[0]].i_stree].StartGen) << ".csv";
		master_filenames[n]=ss.str();
		masterout.open(ss.str().c_str(),std::ofstream::out);
		//cout << masterout.good() << ' ' << masterout.eof() << ' ' << masterout.fail() << ' ' << masterout.bad() << '\n';
		masterout << setprecision(OUTPUT_PRECISION_FLUX);

		ofstream output_lungunits;
		ss.str("");
		if (n == 0) ss << "lung_unit_map_" << o.SimID << "_stree.csv";
		else ss << "lung_unit_map_" << o.SimID << "_ltree_" << ltree[n - 1].pert_type << '_' << read_lobe_no(ltree[n - 1].st[ltree[n - 1].pert_ij[0]].i_stree, o) << "_j" << (ltree[n - 1].pert_ij[1] - ltree[n - 1].st[ltree[n - 1].st[ltree[n - 1].pert_ij[0]].i_stree].StartGen) << ".csv";
		output_lungunits.open(ss.str().c_str(), std::ofstream::out);
		output_lungunits << "Lung_unit_number , Stree_unit_number , No_of_acini\n";
		if (n == 0)
		{
			masterout << "Time_Secs , Tot_Lung_Vol_Litres, Tot_Lung_Acin_Vol_Litres, Tot_Lung_Gas_Vol_Litres , Inert_Gas_Conc_Mouth_Normalised , Pleural_Pressure_cmH20, Flow_Time_Secs, Flow_Rate_Mouth_Litres_per_Sec , Trachea_Flux_Sum_Litres";
			if (o.N2leak.exists) masterout << " , Inert_Gas_Conc_Mouth_Measured , Flow_Rate_Mouth_Litres_per_Sec_Measured";
			for (unsigned m = 0; m < stree.EndSubtrees.size(); m++)
            {
				unsigned im = stree.EndSubtrees[m];
				unsigned jend = stree.st[im].StartGen + stree.st[im].Ncond;
				if(o.output_lung_volumes)
				{
					masterout << " , Volume_Litres_Unit" << m;
					masterout << " , IG_Acin_Volume_Litres_Unit" << m;
				}
				output_lungunits << m << " , " << m << " , " << stree.st[im].gn[jend].Nb << '\n';
            }
		}
		else
		{
			masterout << "Time_Secs , Delta_Tot_Lung_Vol_Litres, Delta_Tot_Lung_Acin_Vol_Litres, Delta_Tot_Lung_Gas_Vol_Litres , Delta_Inert_Gas_Conc_Mouth_Normalised , Delta_Pleural_Pressure_cmH20 , Flow_Time_Secs, Delta_Flow_Rate_Mouth_Litres_per_Sec ,  Delta_Trachea_Flux_Sum_Litres";
			if (o.N2leak.exists) masterout << " , Delta_Inert_Gas_Conc_Mouth_Measured , Delta_Flow_Rate_Mouth_Litres_per_Sec_Measured";
			for (unsigned m = 0; m < ltree[n-1].EndSubtrees.size(); m++)
            {
				unsigned im = ltree[n-1].EndSubtrees[m];
				unsigned jend = ltree[n - 1].st[im].StartGen + ltree[n - 1].st[im].Ncond;
				if(o.output_lung_volumes)
				{
					masterout << " , Delta_Volume_Litres_Unit" << m;
					masterout << " , Delta_IG_Acin_Volume_Litres_Unit" << m;
				}
				output_lungunits << m << " , " << stree.STno_to_ESTno[ltree[n - 1].st[im].i_stree] << " , " << ltree[n - 1].st[im].gn[jend].Nb << '\n';
            }
		}
		masterout << '\n';
		masterout.close();
		output_lungunits.close();
	}
	//-------------------------------------//

	//-------------- Pre-run breath dynamics-------------//
	if(pre_run_breaths(stree,ltree,time,o)) return 1;
	if(initialise_lung_volumes(stree, ltree, o)) return 1;
	//---------------------------------------------------//

	//-------------- Iterate over time -----------------//

	double at = 0;          //counts time between full outputs
	//print initial conditions 
	for (unsigned n = 0; n <= o.Ntrees; n++)
	{
		if (n == 0)
		{
			if (append_masterout(master_filenames[n], -o.dt, stree, o, cons)) cout << "Error writing to master output file.\n";
		}//master output (mass and fluxes)
		else
		{
			if (append_masterout(master_filenames[n], -o.dt, ltree[n - 1], o, cons)) cout << "Error writing to master output file.\n";    //master output (mass and fluxes)
		}
	}
	if (time >= at)                                                  //print out full concentration profile
	{
		printfunc(stree, ltree, o, cons, at / o.printerval);
		at = at + o.printerval;
	}
	//start updates
	for(t=0; t<Nt; t++)
	{
		time = t*o.dt;
		//--Update flow--//
		if (update_flux(stree, ltree, time, o))
		{
			cout << "Problem updating flux. Aborting.\n";
			return 1;
		}                            //work out new volume (symmetric case)
		//--Update concentration--//
		if (update_conc(stree, ltree, time, o))
		{
			cout << "Problem updating conc. Aborting.\n";
			return 1;        //update concentration profile
		}
		//--Print out after velocity and conc is calculated--//
		for (unsigned n = 0; n <= o.Ntrees; n++)
		{
			if (n == 0)
			{
				if (append_masterout(master_filenames[n], time, stree, o, cons)) cout << "Error writing to master output file.\n";
			}//master output (mass and fluxes)
			else
			{
				if (append_masterout(master_filenames[n], time, ltree[n - 1], o, cons)) cout << "Error writing to master output file.\n";    //master output (mass and fluxes)
			}
		}
		if (time + o.dt >= at)                                                  //print out full concentration profile
		{
			printfunc(stree, ltree, o, cons, at / o.printerval);
			at = at + o.printerval;
		}

	}

	//--Time simulation end--//
    summary_file.open(fnss.str().c_str(),ofstream::app);
	o.tend = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = o.tend-o.tstart;
	summary_file << "\nTime taken: " << elapsed_seconds.count() << "s\n";
	//-----------------------//

	//--Clear memory--//
	delete[] master_filenames;
	summary_file.close();

	for (unsigned i = 0; i < stree.st.size(); i++)
	{
		stree.st[i].deallocate_gn();
		for (unsigned n = 0; n < o.Ntrees; n++)
		{
			ltree[n].st[i].deallocate_gn();
		}
	}
	ltree.clear();

	return 0;
}
