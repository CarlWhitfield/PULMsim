//
//  lung_model_discrete_branching.hpp
//  
//
//  Created by Carl Whitfield on 09/03/2017.
//
//
// Checked and commented 16/05/2018
#if defined(_WIN32) || defined(_WIN64)
#define _USE_MATH_DEFINES
#include <windows.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <string.h>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <chrono>
#include <ctime>
#include <Eigen/Eigen/Dense>
#include <Eigen/Eigen/LU>
#include <Eigen/Eigen/Sparse>
#include <Eigen/Eigen/SparseLU>
#include <Eigen/Eigen/IterativeLinearSolvers>

//Options for tree config
#define LINEARPERT 1       //option TREE l
#define NOPERT 2    //option TREE n

//options for BCs
#define BAG 101            //option BC b
#define NOFLUX 102         //option BC n
#define SINK 103

//options for Defects
#define BLOCKAGE 201
#define AREA 202
#define LENGTH 203
#define ELASTICITY 204
#define BAG_RESISTANCE 205

//options for upwind scheme
#define CENTRAL_UPWIND 301
#define FIRST_ORDER_UPWIND 302      //First order upwind interpolation

//options for flow type
#define POISEUILLE 401   //--Poiseuille Flow Assumption--//
#define PEDLEY 402       //--See Pedley 1970--//

//option for elastic response function
#define LINEAR_RESP 501
#define NONLINEAR_RESP 502

//options for pressure profile
#define STEP_FUNCTION 601
#define SIGMOIDAL 602
#define SINUSOIDAL 603
#define LINEAR 604

//options for parameters to use
#define GEOMETRIC    701
#define HOMOGENISED  702
#define WEIBEL  703
#define LOBE_GEOMETRIC 704
#define ALT_LOBE_GEOMETRIC 705

//options for linear solve
#define ITERATIVE    801
#define DIRECT       802

//options for initialisation
#define EMPTY 901
#define FULL 902      

//options for perturbations
#define PERT_AREA 1001
#define PERT_LENGTH 1002
#define PERT_ELASTICITY 1003

//options for volume/pressure input
#define VOLUME 1101
#define PRESSURE 1102

//options for dispersion
#define TAYLOR 1201
#define SCHERER 1202
#define NONE 1203

//Options for concentration file printout
#define VTK 1301
#define CSV 1302
#define VTK_CSV 1303

//Options for leak types
#define INSPIRATORY_LEAK 1401
#define EXPIRATORY_LEAK 1402
#define INSPIRATORY_AND_EXPIRATORY_LEAK 1403

/********************************************/
//Unit conversions
#define cmH20_Pa 100    //pressure conv from cm^2H20 to Pa
#define L_m3 0.001      //volume from litres to m^3
#define cm_m 0.01       //cm to metres
#define PEDLEY_C 1.85   //pedley constant

//Default values of parameters
#define NGENTOT 23              //Number of generations total
#define NGENDEF 15              //Number of generations before model terminates
#define VISCDEF 1.93E-07            //Air Viscosity (cmH20 s) (actual, 1.93E-07 at 37C) -- engineering toolbox
#define DENSITYDEF 1.138        //Air Density (kg m^-3) 1.138 at 37C -- engineering toolbox
#define EDEF 5.0              //Elastance (cmH20 L^-1)
#define RBDEF 0.2              //bag resistance
#define RMDEF 0.6              //mouth resistance
#define DIFFUSIONDEF 0.105       //Diffusion constant (cm2/s) (Helium is approx 0.71, Oxygen 0.19)
#define MAXPECDEF 10            //Max element peclet number
#define K0DEF 0                 //Uptake term
#define LDDEF 3.0                //Length to diameter ratio
#define LD2DEF 2.3                //Length to diameter ratio
#define VFRCDEF 3.0            //Functional residual capacity in L
#define VDDEF 0.12             //Conducting airway dead space
#define VDMDEF 0.05            //Volume of mouth cavity
#define VDUCTDEF 0.20            //fraction of acinus consisting of duct at FRC
#define TINDEF 2.5              //Duration of first breath in (secs)
#define P0DEF 1.0           //Max pressure applied at distal end (cmH20)
#define VTDEF 1.0
#define RUNDEF 40               //Total Simulation Time (number of breathing cycles)
#define LAMBDADEF 0.794          //Geometric parameter
#define LAMBDA2DEF 0.92         //Geom parameter acinus
#define STDEF 0               //Time source is left at opening (in number of breaths)
#define DEFPERTSIZE 0.1       //perturbation size
#define DEFPRINTERVAL 10       //time between prints (in terms of length of first breath)
#define DEFMINGENSIZE 1          //minimum number of points in a generation
#define MINGENSIZE_ACIN 4
#define DEFSTARTDELAY 10         //number of breaths to simulate before starting transport simulation
#define DFACTOR 0.2   //increase in eff area in acinus (default)

//Tolerances precisions etc.
#define PTOL 1E-12  //precision for stree calc
#define TOL 1E-12   //precision for ltree calc
#define OUTPUT_PRECISION_CONC 12   //precision for conc file output
#define OUTPUT_PRECISION_FLUX 12   //precision for flux file output
#define DEF_TS 0.01               //defined timestep
#define DEF_DX 0.025               //defined max spacestep

//define stree numbers for lobar models -- trachea is 0
#define BR 1    //right branch
#define BL 2    //left branch
#define BRML 3   //right middle/lower
#define BRU 4    //right upper
#define BLL 5    //left lower
#define BLU 6    //left upper
#define BRM 7    //right middle
#define BRL 8    //right lower
#define BRL1 9   //right lower major
#define BRL2 10  //right lower minor
#define BLL1 11  //left lower major
#define BLL2 12  //left lower minor

using namespace std;

class Defect    //class for storing ariway/acinar defect information 
{
public:
	unsigned jn;   //generation number
	unsigned long kstart, kend;   //branch number start and end
	unsigned type;       //defect type
	double mag;      //magnitude of defect
	Defect()
	{
		jn = 0;
		kstart = 0;
		kend = 0;
		type = 0;   //should return error if unchanged
		mag = 0;
	}
	inline bool operator==(Defect d)   //equality operator
	{
		if (type == BLOCKAGE) return(type == d.type);  //defects are equal if both blocked
		else return (type == d.type && ((int)(1000000 * mag)) == ((int)(1000000 * d.mag)));  //defects equal if same type and magnitude
	}
};

class Leak    //class for storing leak details - for simulating inspiratory and expiratory leaks at mouth
{
public:
	bool exists;
	double start, end, size;    //start time (breaths), end time (breaths), size (fraction)
	unsigned type;             //Insp only, exp only or both
	Leak()
	{
		exists = false;
		start = 0;
		end = 100000;
		size = 0;
		type = INSPIRATORY_AND_EXPIRATORY_LEAK;
	}
};

class Conversions    //class for storing unit conversions
{
public:
	double LL_to_cm;    //lung length (sim units) to cm
	double P_to_cmH20;   //pressure (sim units) to cmH20
	double t_to_s;      //time (sim units) to seconds
	double V_to_Litres;    //volume (sim units) to litres
};

class Options    //for storing input options
{
public:
	/***Determined by input file***/

	/**Overarching options**/
	unsigned TreeOp;         //Tree format
	unsigned BcOp;           //Boundary condition
	unsigned UpwindOp;       //Upwind interpolation scheme
	unsigned PressOp;        //Function for pressure profile
	unsigned RespFunc;       //Function for elastic resistance
	unsigned FlowType;       //Poiseuille or turbulent
	unsigned ShapeOp;        //Determines the shape and size of the tubes
	unsigned SolverOp;
	unsigned InitOp;        //Determines c initialisation (empty c=0, full c=1)
	unsigned InputOp;         //Determines whether volume or pressure is inputted
	unsigned TaylorDisp;    //True if Taylor Dispersion is considered
	unsigned OutputOp;   //output format option
	bool output_perts;  //whether to output pert files or not
	bool output_lung_volumes;
	Leak N2leak;

	/**Parameters**/
	unsigned Ngen;           //Number of generations before alveoli
	unsigned Ngen2;          //Total number of generations in tree
	unsigned MinGenSize, MinAcinGenSize;
	double dt, dxmax;                //timestep
	double Viscosity;        //Air viscosity
	double Density;          //Air density
	double Diffusion;      //Diffusion constant
	double k0;          //bare uptake rate
	double LDratio, LDratio2;          //L0 length of trachea
	double VFRC;        //total lung resting volume
	double Tin;         //period first inhalation
	double P0;          //Pressure on first inhalation (always negative)
	double VT, VD, VDM, Vductvol;          //Tidal volume
	double E, Rb, Rmouth, V0;           //Elastance of bag, (Total) Resistance of bag, volume of bag
	double RunTime;     //Total run time of simulation in breaths
	double PrerunTime;
	double lambda, lambda2;      //Geometric ratio of tree generations
	double StimTime;    //Time stimulation kept at inlet
	double printerval;  //Time between data outputs
	double MaxPeclet;   //Maximum element peclet number allowed
	double AcinAreaFactor;
	vector<vector<Defect>> def;     //Vector of large magnitude perturbations to symmetric tree
	//vector<Defect> pert;      //Vector of small magnitude perturbations to consider individually
	//vector<vector<char>> *tree_pert;
	unordered_map<unsigned long, vector<Defect>> def_map;   //map from position in tree to defects at that location
	unsigned long kmtot;   //total number of finite volumes

	/***Determined in code***/
	chrono::time_point<std::chrono::system_clock> tstart, tend;
	unsigned Ntrees;               //Number of subtrees - identifier for this run
	string SimID;
	double Peclet;        //Peclet number for branch 0
	/***Constructor - Default Values***/
	Options()
	{
		TreeOp = NOPERT;
		BcOp = NOFLUX;
		UpwindOp = FIRST_ORDER_UPWIND;
		PressOp = SINUSOIDAL;
		RespFunc = LINEAR_RESP;
		FlowType = POISEUILLE;
		ShapeOp = GEOMETRIC;
		SolverOp = ITERATIVE;
		InitOp = FULL;
		InputOp = VOLUME;
		TaylorDisp = SCHERER;
		output_perts = true;
		output_lung_volumes = true;
		OutputOp = CSV;
		Ngen = NGENDEF;
		Ngen2 = NGENTOT;
		MinGenSize = DEFMINGENSIZE;
		MinAcinGenSize = MINGENSIZE_ACIN;
		Viscosity = VISCDEF;
		Density = DENSITYDEF;
		E = EDEF;
		Rb = RBDEF;
		Rmouth = RMDEF;
		Diffusion = DIFFUSIONDEF;
		MaxPeclet = MAXPECDEF;
		k0 = K0DEF;
		LDratio = LDDEF;
		LDratio2 = LD2DEF;
		VFRC = VFRCDEF;
		VD = VDDEF;
		VDM = VDMDEF;
		Vductvol = VDUCTDEF;
		Tin = TINDEF;
		P0 = P0DEF;
		VT = VTDEF;
		RunTime = RUNDEF;
		PrerunTime = DEFSTARTDELAY;
		lambda = LAMBDADEF;
		lambda2 = LAMBDA2DEF;
		StimTime = STDEF;
		printerval = DEFPRINTERVAL;
		dt = DEF_TS;
		dxmax = DEF_DX;
		AcinAreaFactor = DFACTOR;
		def.resize(13);
	}
	//functions to return text for options
	string read_tree_option();
	string read_bc_option();
	string read_flow_option();
	string read_resp_option();
	string read_pressure_option();
	string read_shape_option();
	string read_upwind_option();
	string read_perttype();
	string read_taylor_option();
	string read_solver_option();
	string read_init_option();
	string read_input_option();
	string read_output_option();
	string read_leak_type();
};

struct node     //stores info at a single finite volume node
{
	double c, cold, Aold, Anew, DA, x, uc, y0;        //values of conc, previous conc, totarea, position, vel at element centre
	double ul, ur, al, ar, Dl, Dr;                             //velocity, diffusion and area on left and right of element
	unsigned long km;                                       //corresponding entrance in matrix
	vector<unsigned> iup[3], jup[3], kup[3];        //stores indices for neighbours n-1, n and n+1 respectively
	vector<double> dxup[3], dcfr[3], dcfl[3], ucfrpos[3], ucflpos[3], ucfrneg[3], ucflneg[3], sr, sul;   //vector of dx at n-1, n  and n+1 and cross-sections
	double sl, sdr;   //cross section left sl and right (for node to left) sdr 
	node()   //constructor
	{
		c = 0;
		cold = 0;
		Aold = 0;
		Anew = 0;
		DA = 0;
		x = 0;
		uc = 0;
		ul = 0;
		ur = 0;
		al = 0;
		ar = 0;
		Dl = 0;
		Dr = 0;
		km = 0;
		sl = 0;
		sdr = 0;
	}
	void set_quantities_zero()   //rest relevant quantities to zeros
	{
		c = 0;
		cold = 0;
		Aold = 0;
		Anew = 0;
		DA = 0;
		x = 0;
		uc = 0;
		ul = 0;
		ur = 0;
		al = 0;
		ar = 0;
		Dl = 0;
		Dr = 0;
		sl = 0;
		sdr = 0;
		km = 0;
		for (unsigned counter = 0; counter < 3; counter++)
		{
			for (unsigned counter2 = 0; counter2 < dxup[counter].size(); counter2++)
			{
				dxup[counter][counter2] = 0;
				dcfr[counter][counter2] = 0;
				dcfl[counter][counter2] = 0;
				ucfrpos[counter][counter2] = 0;
				ucflpos[counter][counter2] = 0;
				ucfrneg[counter][counter2] = 0;
				ucflneg[counter][counter2] = 0;
			}

		}
		for (unsigned counter2 = 0; counter2 < sr.size(); counter2++) sr[counter2] = 0;
		for (unsigned counter2 = 0; counter2 < sul.size(); counter2++) sul[counter2] = 0;
	}
};

struct gen    //stores a generation of finite volume nodes
{
	double dx;              //element length
	unsigned long Nb;        //No of branches in this gen
	double x0;         //start position
	vector<node> p;    //vector of nodes
	double Ap, Lp;     //defects/perturbations
	gen()    //Constructor
	{
		dx = 0;
		Nb = 0;
		x0 = 0;
	}
	void set_quantities_zero()
	{
		dx = 0;
		x0 = 0;
		Ap = 0;
		Lp = 0;
	}
};

struct subtree    //stores a subtree consisting of generations connected in series
{
	unsigned StartGen;           //generation where subtree starts
	unsigned long StartBranch;        //branch no. of tree base
	unsigned EndGen;             //generation where subtree ends
	unsigned Ncond, Ntot;         //Gen at which conducting airways terminate and alv airways terminate
	unsigned imeanpath;         //number of mean path
	double A0, L0, A1, L1;
	double Vold;            //equivalent total volume of bags at end at time t
	double Vnew;            //equivalent total volume of bags at end at time t+dt
	double Valv0, sValv0, V0, Vacinduct;           //length of all pipes in alveoli region
	double qend;
	vector<gen> gn;                 //generation
	double fluxin, fluxout; //instantantaneous flux at inlet and outlet
	double masstot, totIGvol;         //total mass in system
	double totIGvolold;         //total mass in system
	int treein;             //number of the tree feeding this one (-1 if root)
	int treeout[2];         //numbers of tree fed by this one (both -1 if base)
	double E, Ep, Rb, Rbp;              //tree elastance and bag resistance
	double dR;              //effective resistance change at start gen of tree
	double Rj;                  //effective resistance of subtree
	double Rtree;                 //resistance of this tree only
	double dVacinus;           //stores the volume change of the associated acinus (or acini)
	double z0;   //coordinate for plotting
	vector<unsigned> EndSubtrees;              //list of terminating subtrees connected to this subtree
	vector<unsigned> isub;                     //list of all subtrees
	bool blocked;     //true if blocked off, false if not
	unsigned i_stree;   //stores index of corresponding subtree in stree
	subtree()   //constructor
	{
		StartGen = 0;
		EndGen = NGENTOT;
		Vold = 0;
		Vnew = 0;
		Valv0 = 0;
		sValv0 = 0;
		Vacinduct = 0;
		L0 = 0;
		A0 = 0;
		L1 = 0;
		A1 = 0;
		V0 = 0;
		qend = 0;
		fluxin = 0;
		fluxout = 0;
		masstot = 0;
		totIGvol = 0;
		totIGvolold = 0;
		treein = -1;
		treeout[0] = -1;
		treeout[1] = -1;
		E = EDEF;
		Ep = 0;
		Rb = RBDEF;
		Rbp = 0;
		dR = 0;
		Rj = 0;
		Rtree = 0;
		z0 = 0;
		blocked = false;
		i_stree = 0;
		dVacinus = 0;
	}

	~subtree()
	{
		gn.clear();
	}
	void set_quantities_zero() //set relevant quantities to zero for perturbative model
	{
		Vold = 0;
		Vnew = 0;
		Valv0 = 0;
		sValv0 = 0;
		V0 = 0;
		qend = 0;
		fluxin = 0;
		fluxout = 0;
		masstot = 0;
		totIGvol = 0;
		totIGvolold = 0;
		EndSubtrees.clear();
		Ep = 0;
		Rbp = 0;
		dR = 0;
		Rj = 0;
	}
	void allocate_gn(unsigned size)
	{
		gn.resize(size);
	}
	void deallocate_gn()
	{
		gn.clear();
	}
	double length_function(unsigned j, Options &o);    //returns gen j length
	double area_function(double dx, unsigned j, Options &o);       //returns dx along gen j cross section  
	double alveolar_density(unsigned j, unsigned k, Options &o);    //returns alveolar density function
	void subtree_resistance(void);                               //computes subtree resistance
	void tot_area_function(double Vtot);                         //distributes sac volume to cross section A
	void tot_lp_area_function(subtree &sst, double DVtot);        //same as abouve for perturbed tree
	void velocity_calc(double qend, Options &o);                //calculates flow velocity from A
	void velocity_lp_calc(subtree &sst, double qend, Options &o);          //same as above for perturbed tree
};

class Tree                  //for storing whole airway tree data
{
public:
	vector<subtree> st;
	double masstot, masstotold, fluxin;            //total mass in system
	double Ppl, Pplold;                //pleural pressure
	unsigned pert_ij[2];          //stores st and generation of perturbation
	char pert_type;               //'K', 'A' or 'L' depending on which per
	//Eigen matrices and functions fro ventialtion and transport problems
	Eigen::HouseholderQR<Eigen::MatrixXd> *AmatQR;
	Eigen::PartialPivLU<Eigen::MatrixXd> *AmatLU;
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::IncompleteLUT<double>> *solver_iter;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> *solver_dir;
	Eigen::MatrixXd Amat, Bmat, Kmat;
	Eigen::VectorXd KVs, Vs;
	double VDfront, VDfrontold, VDfrontoldold, cmouth, cmouthold, Vlung, Vairways;  //volume of mouth dead space filled this breath
	double Vlungold;
	vector<double> ctop, Vtot;
	vector<unsigned> EndSubtrees, meanpath;         
	unordered_map<unsigned, unsigned> STno_to_ESTno;   //subtree number to end subtree number
	vector<vector<unsigned>> mpsubtrees;  //stores subtrees that are in a given mean path
	double R0;
	unsigned long kmtot;  //total number of points for transport equations
	Tree() //constructor
	{
		masstot = 0;
		fluxin = 0;
		Ppl = 0;
		Pplold = 0;
		pert_ij[0] = 0;
		pert_ij[1] = 0;
		pert_type = 'N';  //N for none
		cmouth = 1;
		cmouthold = 1;
		VDfront = 0;
		VDfrontold = 0;
		VDfrontoldold = 0;
		kmtot = 0;
		ctop.clear();
		Vtot.clear();
	}
	~Tree()
	{
		ctop.clear();
		Vtot.clear();
	}
	//output functions
	int printfunc_tree(string filename, Conversions cons);
	int printfunc_lp_tree(string filename, Tree &stree, Conversions cons);
	int printfunc_tree_csv(string filename, Conversions cons);
	int printfunc_lp_tree_csv(string filename, Tree &stree, Conversions cons);
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

//read and set simulation options
int parse_options(string infile, Options &o);
void check_option_consistency(Options &o);
int convert_units(Tree &stree, Conversions &cons, double LL, Options &o);

//initialisation
int initialise_system(Tree &stree, Conversions &cons, Options &o);
int initialise_defects(Tree &stree, Options &o);
int initialise_tree(Tree &stree, Options &o);
int pre_run_breaths(Tree &stree, vector<Tree> &ltree, double &time, Options &o);
int initialise_lung_volumes(Tree &stree, vector<Tree> &ltree, Options &o);

//tree building
double build_alternative_lobe_branches(Tree &stree, Options &o);
double build_symmetric_branches(Tree &stree, Options &o);
string read_lobe_no(unsigned i, Options &o);

//build perturbed trees
void setup_pert_trees(vector<Tree> &ltree, Tree &stree, Options &o);
void setup_lp_tree(vector<Tree> &ltree, Tree &stree, Options &o);
void calc_perts(vector<Tree> &ltree, Tree &stree);
void apply_pert(Tree &ltree, Tree &stree);

//build and solve ventilation problem
int update_flux(Tree &stree, vector<Tree> &ltree, double time, Options &o);
int tree_flux(Tree &stree, double time, Options &o);
int tree_lp_flux(Tree &ltree, Tree &stree, Options &o);

//build and solve transport problem
int update_conc(Tree &stree, vector<Tree> &ltree, double time, Options &o);
int update_c_stree(Tree &stree, double time, Options &o);
int update_c_ltree(Tree &ltree, Tree &stree, double time, Options &o);
void calc_gas_concs(Tree &stree, Eigen::VectorXd &XC, Options &o);
void fill_Ab(Tree &stree, Eigen::SparseMatrix<double, Eigen::RowMajor> &AC, Eigen::VectorXd &XC, Eigen::VectorXd &BC, unsigned long Ntot, double time, Options &o);
void calc_gas_concs_lp(Tree &ltree, Tree &stree, Eigen::VectorXd &XC, Options &o);
void fill_Ab_lp(Tree &ltree, Tree &stree, Eigen::SparseMatrix<double, Eigen::RowMajor> &AC, Eigen::VectorXd &XC, Eigen::VectorXd &BC, unsigned long Ntot, Options &o);

//output functions
int printfunc(Tree &stree, vector<Tree> &ltree, Options &o, Conversions cons, double t);
int append_masterout(string filename, double time, Tree &tree, Options &o, Conversions cons);
int print_simops(ofstream &summary_file, Options &o);
int print_params(ofstream &summary_file, Options &o);

//various input functions determining parameters
double pressure_func(double t, Options &o);
double c_stim(double t, Options &o);

double total_run_time(Options &o);
unsigned long calc_ts_number(Options &o);
unsigned lobe_gens(unsigned nlobe, Options &o);  //returns gen number of final branch
unsigned long ijk_index(unsigned i, unsigned j, unsigned long k, Options &o);
double weibel_length(unsigned j);
double weibel_area(unsigned j);

//calculating resistance matrices for ventilation problem
void subtree_resistance(Tree &tree, unsigned i, Options &o);
void subtree_lp_symm_resistance(subtree &st, subtree &symm, Options &o);
int full_tree_fluxcalc(Tree &tree, double time, Options &o);
int full_lp_tree_fluxcalc(Tree &ltree, Tree &stree, Options &o);
int build_resistance_matrix(Tree &tree, Options &o);
bool isparent(Tree &tree, unsigned ip, unsigned id);

//finite difference functions and solving routines
void upwind_coeffs(vector<double> ucfl[3], vector<double> ucfr[3], vector<double> dxh[3], double u0, Options &o);
void lp_upwind_coeffs(vector<double> ducfl[3], vector<double> ducfr[3], vector<double> dx0[3], vector<double> dxh[3], double u0, Options &o);
void fd_coeffs(vector<double> dcfl[3], vector<double> dcfr[3], vector<double> dxh[3]);
void lp_fd_coeffs(vector<double> ddcfl[3], vector<double> ddcfr[3], vector<double> dx0[3], vector<double> dxh[3]);
void tree_point_pos(Tree &tree, unsigned i, unsigned j, unsigned k, int n, vector<unsigned> &in, vector<unsigned> &jn, vector<unsigned> &kn);
int iterative_solver(Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::IncompleteLUT<double>> *solver, Eigen::SparseMatrix<double, Eigen::RowMajor> &AC, Eigen::VectorXd &x, Eigen::VectorXd &b);
int LU_solver(Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> *solver, Eigen::SparseMatrix<double> &A, Eigen::VectorXd &x, Eigen::VectorXd &b);

inline unsigned bitflip(unsigned m)
{
	if (m == 0) return 1;
	else return 0;
}

inline bool jcomp(Defect d1, Defect d2)  //compare defect by i, then j, then k, then type, then magntiude
{
	return ((d1.jn < d2.jn) || ((d1.jn == d2.jn) && (d1.kstart < d2.kstart)) || ((d1.jn == d2.jn) && (d1.kstart == d2.kstart) && (d1.type < d2.type)));
}
