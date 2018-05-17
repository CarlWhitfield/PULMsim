//
//  read_inputs.cpp
//  
//
//  Created by Carl Whitfield on 09/03/2017.
//
//

//--Checked and commented 05/07/17--//

#include "lung_model_discrete_branching.hpp"

int parse_options(string infile, Options &o)  //crunch through input file
{
	ifstream options;
	options.open(infile);
	char buff[1000];
	unsigned buffsize = 1000;
	string temp, nc;
	stringstream ss;
	bool done;
	unsigned i;
	
	if (!options.good())
	{
		cout << "Error, could not read input.\n";    //check file
		return 1;
	}
	while (options.good())                            //read
	{
		options.getline(buff, buffsize);             //read line from file
		ss.str(buff);
		temp.clear();
		ss >> skipws >> temp;                //parse line
		done = false;
		if (temp.size() == 0 || temp[0] == '%') done = true;       //ignore empty strings and comments
		else   //read in values
		{
			ss >> skipws >> nc;                     //read value into string
			//-----------Code Options (determine what simulation type to run-----------//
			
			if (temp == "TREE")    //tree type
			{
				switch (nc[0])  //set tree option
				{
					case 'n': o.TreeOp = NOPERT;
						break;
						
					case 'l': o.TreeOp = LINEARPERT;
						break;
						
					default:
					{
						cout << "Did not recognise TREE option: " << nc << '\n';
						cout << "Using TREE default option:" << o.read_tree_option() << '\n';
					}
				}
				done = true;
			}
			
			if (temp == "BC")  //boundary condition at end
			{
				switch (nc[0])
				{
					case 'b': o.BcOp = BAG;
						break;
						
					case 'n': o.BcOp = NOFLUX;
						break;
						
					case 's': o.BcOp = SINK;
						break;
						
					default:
					{
						cout << "Did not recognise BC option: " << nc << '\n';
						cout << "Using BC default option: " << o.read_bc_option() << '\n';
					}
				}
				done = true;
			}
			
			
			
			if (temp == "UPWIND")   //upwind interpolation
			{
				switch (nc[0])
				{
					case 'c': o.UpwindOp = CENTRAL_UPWIND;
						break;
						
					case 'f': o.UpwindOp = FIRST_ORDER_UPWIND;
						break;
						
					default:
					{
						cout << "Did not recognise UPWIND option: " << nc << '\n';
						cout << "Using UPWIND default option: " << o.read_upwind_option() << '\n';
					}
				}
				done = true;
			}
			
			
			
			if (temp == "PRESSURE")      //applied pressure function
			{
				switch (nc[0])
				{
					case 'd': o.PressOp = STEP_FUNCTION;
						break;
						
					case 's': o.PressOp = SIGMOIDAL;
						break;
						
					case 'w': o.PressOp = SINUSOIDAL;
						break;
						
					case 'l': o.PressOp = LINEAR;
						break;
						
					default:
					{
						cout << "Did not recognise PRESSURE option: " << nc << '\n';
						cout << "Using PRESSURE default option: " << o.read_pressure_option() << '\n';
					}
				}
				done = true;
			}
			
			
			
			if (temp == "RESP")             //elastic response function
			{
				switch (nc[0])
				{
					case 'l': o.RespFunc = LINEAR_RESP;
						break;
						
					case 'n': o.RespFunc = NONLINEAR_RESP;
						break;
						
					default:
					{
						cout << "Did not recognise RESP option: " << nc << '\n';
						cout << "Using RESP default option: " << o.read_resp_option() << '\n';
					}
				}
				done = true;
			}
			
			
			
			if (temp == "FLOW")            //flow pressure relationship
			{
				switch (nc[0])
				{
					case 'p': o.FlowType = POISEUILLE;
						break;
						
					case 't': o.FlowType = PEDLEY;
						break;
						
					default:
					{
						cout << "Did not recognise FLOW option: " << nc << '\n';
						cout << "Using FLOW default option: " << o.read_flow_option() << '\n';
					}
				}
				done = true;
			}
			
			
			
			if (temp == "SHAPE")         //area distribution function
			{
				switch (nc[0])
				{
					case 'g': o.ShapeOp = GEOMETRIC;
						break;
						
					case 'h': o.ShapeOp = HOMOGENISED;
						break;
						
					case 'w': o.ShapeOp = WEIBEL;
						break;
						
					case 'l': o.ShapeOp = LOBE_GEOMETRIC;
						break;

					case 'a': o.ShapeOp = ALT_LOBE_GEOMETRIC;
						break;
						
					default:
					{
						cout << "Did not recognise SHAPE option: " << nc << '\n';
						cout << "Using SHAPE default option: " << o.read_shape_option() << '\n';
					}
				}
				
				done = true;
			}
			
			if (temp == "SOLVER")         //linear solver option
			{
				switch (nc[0])
				{
					case 'i': o.SolverOp = ITERATIVE;
						break;
						
					case 'd': o.SolverOp = DIRECT;
						break;
						
					default:
					{
						cout << "Did not recognise SOLVER option: " << nc << '\n';
						cout << "Using SOLVER default option: " << o.read_solver_option() << '\n';
					}
				}
				
				done = true;
			}
			
			if (temp == "TAYLOR")         //Talyor dispersion?
			{
				switch (nc[0])
				{
					case 't': o.TaylorDisp = TAYLOR;
						break;
						
					case 'f': o.TaylorDisp = NONE;
						break;

					case 's':o.TaylorDisp = SCHERER;
						break;
						
					default:
					{
						cout << "Did not recognise TAYLOR option: " << nc << '\n';
						cout << "Using TAYLOR default option: " << o.read_taylor_option() << '\n';
					}
				}
				
				done = true;
			}
			
			if (temp == "INIT")         //linear solver option
			{
				switch (nc[0])
				{
					case 'e': o.InitOp = EMPTY;
						break;
						
					case 'f': o.InitOp = FULL;
						break;
						
					default:
					{
						cout << "Did not recognise INIT option: " << nc << '\n';
						cout << "Using INIT default option: " << o.read_init_option() << '\n';
					}
				}
				
				done = true;
			}
			
			if (temp == "INPUT")         //Volume or pressure input
			{
				switch (nc[0])
				{
					case 'v': o.InputOp = VOLUME;
						break;
						
					case 'p': o.InputOp = PRESSURE;
						break;
						
					default:
					{
						cout << "Did not recognise INPUT option: " << nc << '\n';
						cout << "Using INPUT default option: " << o.read_input_option() << '\n';
					}
				}
				done = true;
			}

			if (temp == "CPOUT")         //Volume or pressure input
			{
				switch (nc[0])
				{
				case 't': 
				case 'y':
					o.output_perts = true;
					break;

				case 'f':
				case 'n':
					o.output_perts = false;
					break;

				default:
				{
					cout << "Did not recognise CPOUT (output_perts) option: " << nc << '\n';
					cout << "Using CPOUT default option: " << o.output_perts << '\n';
				}
				}
				done = true;
			}
			
			if(temp == "LUOUT" || temp == "VOLUMESOUT" || temp == "UNITSOUT")         //Volume or pressure input
			{
				switch (nc[0])
				{
					case 't':
					case 'y':
					o.output_lung_volumes = true;
					break;
					
					case 'f':
					case 'n':
					o.output_lung_volumes = false;
					break;
					
					default:
					{
						cout << "Did not recognise LUOUT (output_lung_volumes) option: " << nc << '\n';
						cout << "Using LUOUT default option: " << o.output_perts << '\n';
					}
				}
				done = true;
			}
			
			if (temp == "OUTPUT")         //Volume or pressure input
			{
				switch (nc[0])
				{
				case 'v':
					o.OutputOp = VTK;
					break;

				case 'c':
					o.OutputOp = CSV;
					break;

				case 'b':
					o.OutputOp = VTK_CSV;
					break;

				default:
				{
					cout << "Did not recognise OUTPUT option: " << nc << '\n';
					cout << "Using OUTPUT default option: " << o.read_output_option() << '\n';
				}
				}
				done = true;
			}
			
			if (temp == "DEFECTS")
			{
				unsigned long Ndef = ((unsigned long)atoi(nc.c_str()));   //number of defects to read
				unsigned ih = 0;
				if(o.ShapeOp == LOBE_GEOMETRIC) o.def.resize(9);
				if (o.ShapeOp == ALT_LOBE_GEOMETRIC) o.def.resize(13);
				else o.def.resize(1);
				for (unsigned long n = 0; n < Ndef; n++)
				{
					Defect dh;   //defined inside loop (reconstructed each time)
					options.getline(buff, buffsize);             //read line from file
					ss.str("");
					ss.clear();				//throw away rest of line after
					ss.str(buff);
					i = 0;
					nc = "";              //reset nc
					ss >> skipws >> nc;                     //read value into string
					while (nc.size() != 0 && nc[0] != '%' && i < 6)                     //comments start with %
					{
						switch (i)
						{
							case 0: dh.jn = ((unsigned)atoi(nc.c_str()));     //first number is gen in mean path
								break;
								
							case 1: dh.kstart = ((unsigned)atoi(nc.c_str()));     //second number is start branch number in mean path
								break;
								
                            case 2: dh.kend = ((unsigned)atoi(nc.c_str()));     //third number is end branch number in mean path
                                break;
                                
							case 3:
							{
								switch (nc[0])
								{
									case 'b': dh.type = BLOCKAGE;
										break;
										
									case 'a': dh.type = AREA;
										break;
										
									case 'l': dh.type = LENGTH;
										break;
										
									case 'e': dh.type = ELASTICITY;
										break;
										
									case 'r': dh.type = BAG_RESISTANCE;
										break;
										
									default: dh.type = 0; //error
								}
							} break;
								
							case 4:
							{
								if (dh.type == BLOCKAGE) dh.mag = 1.;
								else dh.mag = atof(nc.c_str());
							}
								break;
								
							case 5:
							{
								if (o.ShapeOp == LOBE_GEOMETRIC || o.ShapeOp == ALT_LOBE_GEOMETRIC)
								{
									if(nc == "0" || nc == "T") ih = 0;
									else
									{
										if(nc == "R") ih = BR;
										else
										{
											if(nc == "L") ih = BL;
											else
											{
												if(nc == "RML") ih = BRML;
												else
												{
													if(nc == "RU") ih = BRU;
													else
													{
														if(nc == "RL") ih = BRL;
														else
														{
															if(nc == "RM") ih = BRM;
															else
															{
																if(nc == "LU") ih = BLU;
																else
																{
																	if(nc == "LL") ih = BLL;
																	else
																	{
																		if (nc == "LL1" && o.ShapeOp == ALT_LOBE_GEOMETRIC) ih = BLL1;
																		else
																		{
																			if (nc == "LL2" && o.ShapeOp == ALT_LOBE_GEOMETRIC) ih = BLL2;
																			else
																			{
																				if (nc == "RL1" && o.ShapeOp == ALT_LOBE_GEOMETRIC) ih = BRL1;
																				else
																				{
																					if (nc == "RL2" && o.ShapeOp == ALT_LOBE_GEOMETRIC) ih = BRL2;
																					else
																					{
																						cout << "Did not recognise lobe code. Assuming RU.\n";
																						ih = BRU;
																					}
																				}
																			}
																		}
																	}
																}
															}
														}
													}
												}
											}
										}
									}
								}
								else ih = 0;
								
							}
								break;
						}
						ss >> skipws >> nc;                     //read value into string
						i++;
					}
					o.def[ih].push_back(dh);   //add defect to list
				}
				done = true;
			}
			
			//-----------------------------------------------------------------------------//
			
			//------------------------Code Parameters------------------------------//
			
			if (temp == "Leak_size" || temp == "Leak")
			{
				o.N2leak.exists = true;
				o.N2leak.size = atof(nc.c_str());
				done = true;
			}

			if (temp == "Leak_start")
			{
				o.N2leak.exists = true;
				o.N2leak.start = atof(nc.c_str());
				done = true;
			}

			if (temp == "Leak_end")
			{
				o.N2leak.exists = true;
				o.N2leak.end = atof(nc.c_str());
				done = true;
			}

			if (temp == "Leak_type")
			{
				o.N2leak.exists = true;
				switch (nc[0])
				{
				case 'i':
				case 'I':
				{
					o.N2leak.type = INSPIRATORY_LEAK;
				}	break;

				case 'e':
				case 'E':
				{
					o.N2leak.type = EXPIRATORY_LEAK;
				}	break;

				case 'b':
				case 'B':
				{
					o.N2leak.type = INSPIRATORY_AND_EXPIRATORY_LEAK;
				}	break;

				default:
				{
					cout << "Did not recognise leak type: " << nc << ". Using default type: " << o.read_leak_type() << '\n';
				}
				}
				done = true;
			}
			
			if (temp == "Y" || temp == "AcinFactor" || temp == "AcinAreaFactor")
			{
				o.AcinAreaFactor = atof(nc.c_str());
				done = true;
			}
			
			if (temp == "dx" || temp == "SpaceStep" || temp == "Space_step")
			{
				o.dxmax = atof(nc.c_str());
				done = true;
			}
			
			if (temp == "dt" || temp == "TimeStep" || temp == "Time_step")
			{
				o.dt = atof(nc.c_str());
				done = true;
			}
			
			if (temp == "MinAcinGenSize")
			{
				o.MinAcinGenSize = atoi(nc.c_str());
				done = true;
			}

			if (temp == "MinGenSize")
			{
				o.MinGenSize = atoi(nc.c_str());
				done = true;
			}
			
			if (temp == "Ngen")
			{
				o.Ngen = atoi(nc.c_str());
				done = true;
			}
			
			if (temp == "Ngen2")
			{
				o.Ngen2 = atoi(nc.c_str());
				done = true;
			}
			
			if (temp == "Viscosity")
			{
				o.Viscosity = atof(nc.c_str());
				done = true;
			}
			
			if (temp == "Density")
			{
				o.Density = atof(nc.c_str());
				done = true;
			}
			
			if (temp == "Diffusion")
			{
				o.Diffusion = atof(nc.c_str());
				done = true;
			}
			
			if (temp == "MaxPeclet")
			{
				o.MaxPeclet = atof(nc.c_str());
				done = true;
			}
			
			if (temp == "K0")
			{
				o.k0 = atof(nc.c_str());
				done = true;
			}
			
			if (temp == "P0")
			{
				o.P0 = atof(nc.c_str());
				done = true;
			}
			
			if (temp == "Tin")
			{
				o.Tin = atof(nc.c_str());
				done = true;
			}
			
			if (temp == "LDratio")
			{
				o.LDratio = atof(nc.c_str());
				done = true;
			}

			if (temp == "LDratio2")
			{
				o.LDratio2 = atof(nc.c_str());
				done = true;
			}
			
			if (temp == "VFRC")
			{
				o.VFRC = atof(nc.c_str());
				done = true;
			}
			
			if (temp == "VD")
			{
				o.VD = atof(nc.c_str());
				done = true;
			}
			
			if (temp == "VDM")
			{
				o.VDM = atof(nc.c_str());
				done = true;
			}

			if (temp == "DuctFrac")
			{
				o.Vductvol = atof(nc.c_str());
				done = true;
			}
			
			if (temp == "E")
			{
				o.E = atof(nc.c_str());
				done = true;
			}
			
			if (temp == "Rb" || temp == "Rbag")
			{
				o.Rb = atof(nc.c_str());
				done = true;
			}

			if (temp == "Rm" || temp == "Rmouth")
			{
				o.Rmouth = atof(nc.c_str());
				done = true;
			}
			
			if (temp == "RunTime")
			{
				o.RunTime = atof(nc.c_str());
				done = true;
			}
			
			if (temp == "Lambda")
			{
				o.lambda = atof(nc.c_str());
				done = true;
			}
			
			if (temp == "Lambda2")
			{
				o.lambda2 = atof(nc.c_str());
				done = true;
			}
			
			if (temp == "StimTime")
			{
				o.StimTime = atof(nc.c_str());
				done = true;
			}
			
			if (temp == "PrerunTime")
			{
				o.PrerunTime = atof(nc.c_str());
				done = true;
			}
			
			if (temp == "Printerval")
			{
				o.printerval = atof(nc.c_str());
				done = true;
			}
		}
		if (!done)
		{
			cout << "Did not recognise parameter " << temp << ".\n";
		}
		
		ss.str("");
		ss.clear();				//throw away rest of line after
	}
	
	return 0;
}

void check_option_consistency(Options &o)   //function to check options are self consistent
{
	vector<unsigned long> errors;
	unsigned n=0,n2=0;
	unsigned NMP;
	if(o.ShapeOp == LOBE_GEOMETRIC) NMP = 9;
	if (o.ShapeOp == ALT_LOBE_GEOMETRIC) NMP = 13;
	else NMP = 1;
	
	for(unsigned i = 0; i < NMP; i++)
	{
		while(n < o.def[i].size()) //check defects are consistent
		{
			if (o.def[i][n].type == 0 || o.def[i][n].jn > lobe_gens(i,o) || o.def[i][n].kstart > pow(2, o.def[i][n].jn) - 1 || o.def[i][n].kend > pow(2, o.def[i][n].jn) - 1)
			{
				cout << "Error parsing defect " << n2 << ". Ignoring...\n";
				o.def.erase(o.def.begin() + n);
			}
			else n++;    //only counts remaining elements
			n2++;   //counts original number of elements in vector
		}
		
		
		sort(o.def[i].begin(), o.def[i].end(), jcomp);   	//sort defect vector by i, then j, then k, then type, then magntiude
		
		for (n = 0; n < o.def[i].size(); n++)    //register defects in unordered map (each branch can have more than one defect)
		{
			o.def_map[ijk_index(i,o.def[i][n].jn,o.def[i][n].kstart,o)].push_back(o.def[i][n]);
		}
	}
	
	o.VD = L_m3*o.VD;   //convert to m3
	o.VDM = L_m3*o.VDM;   //convert to m3
	o.VFRC = L_m3*o.VFRC;   //convert to m3
	
	//Parameter inconsistencies
	if (o.InputOp == VOLUME)
	{
		o.VT = L_m3*o.P0;     //if volume is inputted P0 parameter specifies VT
		o.P0 = 0;
	}
	else o.VT = 0;     //otherwise it specifies pressure magnitude
	
	if(o.E < 0)
	{
		cout << "Value of E inconsistent. Setting to default.\n";
		o.E = EDEF;
	}
	
	if(o.dt <= 0)
	{
		cout << "Value of dt inconsistent. Setting to default.\n";
		o.dt= DEF_TS;
	}
	
	if(o.dxmax <= 0)
	{
		cout << "Value of dx inconsistent. Setting to default.\n";
		o.dxmax = DEF_DX;
	}
	
	if(o.AcinAreaFactor < 0)
	{
		cout << "Value of Diffusice Acinar Area Factor (Y) inconsisten. Setting to default.\n";
		o.AcinAreaFactor = DFACTOR;
	}
	
	if(o.MinAcinGenSize <= 0)
	{
		cout << "Value of MinAcinGenSize inconsistent. Setting to default.\n";
		o.MinAcinGenSize = MINGENSIZE_ACIN;
	}
	
	if (o.Rb < 0)
	{
		cout << "Value of Rb inconsistent. Setting to default.\n";
		o.Rb = RBDEF;
	}

	if (o.Rmouth < 0)
	{
		cout << "Value of Rmouth inconsistent. Setting to default.\n";
		o.Rmouth = RMDEF;
	}
	
	if(o.Viscosity <= 0)
	{
		cout << "Value of Viscosity inconsistent. Setting to default.\n";
		o.Viscosity = VISCDEF;
	}
	
	if(o.Diffusion <= 0)
	{
		cout << "Value of Diffusion constant inconsistent. Setting to default.\n";
		o.Diffusion = DIFFUSIONDEF;
	}
	
	if(o.MaxPeclet <= 0)
	{
		cout << "Value of MaxPeclet inconsistent. Setting to default.\n";
		o.MaxPeclet = MAXPECDEF;
	}
	
	if ((o.Ngen <= 3 && o.ShapeOp == LOBE_GEOMETRIC) || (o.Ngen <= 4 && o.ShapeOp == ALT_LOBE_GEOMETRIC))
	{
		cout << "Value of Number of Generations inconsistent. Setting to default.\n";
		o.Ngen = NGENDEF;
	}
	
	if(o.Ngen2 <= o.Ngen)
	{
		cout << "Value of total Number of Generations inconsistent. Setting to default.\n";
		if(o.Ngen >= NGENTOT) o.Ngen2 = o.Ngen + 1;
		else o.Ngen2 = NGENTOT;
	}

	if (o.ShapeOp == WEIBEL)
	{
		o.Ngen2 = 23;
		o.Ngen = 15;
	}

	if (o.BcOp == BAG)
	{
		o.Ngen2 = o.Ngen;
	}
	
	if(o.k0 < 0)
	{
		cout << "Value of k0 inconsistent. Setting to default.\n";
		o.k0 = K0DEF;
	}
	
	if(o.LDratio <= 0)
	{
		cout << "Value of length-diameter ratio inconsistent. Setting to default.\n";
		o.LDratio = LDDEF;
	}

	if (o.LDratio2 <= 0)
	{
		cout << "Value of length-diameter ratio inconsistent. Setting to default.\n";
		o.LDratio2 = LD2DEF;
	}
	
	if (o.VFRC <= 0)
	{
		cout << "Value of VFRC inconsistent. Setting to default.\n";
		o.VFRC = VFRCDEF;
	}
	
	if (o.VD <= 0)
	{
		cout << "Value of VD inconsistent. Setting to default.\n";
		o.VD = VDDEF;
	}
	
	if (o.VDM < 0)
	{
		cout << "Value of VDM inconsistent. Setting to default.\n";
		o.VDM = VDMDEF;
	}

	if (o.Vductvol <= 0)
	{
		cout << "Value of DuctFrac inconsistent. Setting to default.\n";
		o.Vductvol = VDUCTDEF;
	}
	
	if(o.RunTime < 0)
	{
		cout << "Value of RunTime inconsistent. Setting to default.\n";
		o.RunTime = RUNDEF;
	}
	
	if(o.PrerunTime < 0)
	{
		cout << "Value of PrerunTime inconsistent. Setting to default.\n";
		o.PrerunTime = DEFSTARTDELAY;
	}
	if(o.lambda <= 0)
	{
		cout << "Value of lambda inconsistent. Setting to default.\n";
		o.lambda = LAMBDADEF;
	}
	if (o.lambda2 <= 0)
	{
		cout << "Value of lambda2 inconsistent. Setting to default.\n";
		o.lambda2 = LAMBDA2DEF;
	}
	if(o.StimTime < 0)
	{
		cout << "Value of stim time inconsistent. Setting to default.\n";
		o.StimTime = STDEF;
	}
	
	if(o.printerval < 0)
	{
		cout << "Value of printerval inconsistent. Setting to default.\n";
		o.printerval = DEFPRINTERVAL;
	}
	
	if(o.MinGenSize < 1)
	{
		cout << "Value of minimum gen size inconsistent. Setting to default.\n";
		o.MinGenSize = DEFMINGENSIZE;
	}
	
	switch (o.ShapeOp)
	{
		case GEOMETRIC:
		{
			
		}break;
			
	}
}

//----Functions for converting options to strings-----//
string Options::read_tree_option()
{
	switch (TreeOp)
	{
		case NOPERT:
			return "noperts";
			break;
			
		case LINEARPERT:
			return "linearpert";
			break;
			
		default:
			return "other";
	}
}

string Options::read_bc_option()
{
	switch (BcOp)
	{
		case BAG:
			return "bag";
			break;
			
		case NOFLUX:
			return "noflux";
			break;
			
		case SINK:
			return "sink";
			break;
			
		default:
			return "other";
	}
}

string Options::read_upwind_option()
{
	switch (UpwindOp)
	{
		case CENTRAL_UPWIND:
			return "central";
			break;
			
		case FIRST_ORDER_UPWIND:
			return "upwind";
			break;
			
		default:
			return "other";
	}
}

string Options::read_flow_option()
{
	switch (FlowType)
	{
		case POISEUILLE:
			return "poiseuille";
			break;
			
		case PEDLEY:
			return "pedley";
			break;
			
		default:
			return "other";
	}
}

string Options::read_resp_option()
{
	switch (RespFunc)
	{
		case LINEAR_RESP:
			return "linear";
			break;
			
		case NONLINEAR_RESP:
			return "non_linear";
			break;
			
		default:
			return "other";
	}
}

string Options::read_pressure_option()
{
	switch (PressOp)
	{
		case STEP_FUNCTION:
			return "step";
			break;
			
		case SIGMOIDAL:
			return "sigmoidal";
			break;
			
		case SINUSOIDAL:
			return "sinusoidal";
			break;
			
		case LINEAR:
			return "linear_increase";
			break;
			
		default:
			return "other";
	}
}

string Options::read_shape_option()
{
	switch (ShapeOp)
	{
		case GEOMETRIC:
			return "geometric";
			break;
			
		case HOMOGENISED:
			return "homogenised";
			break;
			
		case WEIBEL:
			return "weibel";
			break;
			
		case LOBE_GEOMETRIC:
			return "lobe";
			break;

		case ALT_LOBE_GEOMETRIC:
			return "alt_lobe";
			break;
			
		default:
			return "other";
	}
}

string Options::read_taylor_option()
{
	switch (TaylorDisp)
	{
	case SCHERER:
		return "Scherer";
		break;

	case TAYLOR:
		return "Taylor";
		break;

	case NONE:
		return "NoTaylor";
		break;

	default:
		return "other";
	}
}

string Options::read_solver_option()
{
	switch (SolverOp)
	{
		case ITERATIVE:
			return "Iterative";
			break;
			
		case DIRECT:
			return "Direct";
			break;
			
		default:
			return "other";
	}
}

string Options::read_init_option()
{
	switch (InitOp)
	{
		case EMPTY:
			return "Empty";
			break;
			
		case FULL:
			return "Full";
			break;
			
		default:
			return "other";
	}
}

string Options::read_input_option()
{
	switch (InputOp)
	{
		case VOLUME:
			return "Volume";
			break;
			
		case PRESSURE:
			return "Pressure";
			break;
			
		default:
			return "other";
	}
}

string Options::read_output_option()
{
	switch (OutputOp)
	{
	case VTK:
		return "VTK";
		break;

	case CSV:
		return "CSV";
		break;

	case VTK_CSV:
		return "VTK_and_CSV";
		break;

	default:
		return "other";
	}
}

string Options::read_leak_type()
{
	switch (N2leak.type)
	{
	case INSPIRATORY_LEAK:
		return "Inspiratory";
		break;

	case EXPIRATORY_LEAK:
		return "Expiratory";
		break;

	case INSPIRATORY_AND_EXPIRATORY_LEAK:
		return "Inspiratory_and_Expiratory";
		break;

	default:
		return "other";
	}
}


string read_lobe_no(unsigned i, Options &o)
{
	switch(o.ShapeOp)
	{
		case GEOMETRIC:
		case HOMOGENISED:
			return "symm";
			break;
			
		case ALT_LOBE_GEOMETRIC:
		case LOBE_GEOMETRIC:
		{
			switch(i)
			{
				case 0: return "TB";
					break;
					
				case BR: return "RB";
					break;
					
				case BL: return "LB";
					break;
					
				case BRML: return "RMLB";
					break;
					
				case BRU: return "RU";
					break;
					
				case BLL: return "LLB";
					break;
					
				case BLU: return "LU";
					break;
					
				case BRM: return "RM";
					break;
					
				case BRL: return "RLB";
					break;
					
				case BRL1: return "RLmaj";
					break;
				
				case BRL2: return "RLmin";
					break;
					
				case BLL1: return "LLmaj";
					break;
					
				case BLL2: return "LLmin";
					break;
			}
		}
	}
	stringstream retval;
	retval.str("");
	retval << i;
	return retval.str();
}

int convert_units(Tree &stree, Conversions &cons, double LL, Options &o)
{
	double M0;
	
    //volumes in m^-3
    
	cons.P_to_cmH20 = o.E*o.VFRC/L_m3;   //E*V0 = 1
	cons.t_to_s = o.Tin;       //store breath time in s
	cons.LL_to_cm = LL/cm_m;       //LL in m
	cons.V_to_Litres = LL*LL*LL/L_m3;   //Volume to litres
	
	M0 = (cmH20_Pa*cons.P_to_cmH20)*LL*o.Tin*o.Tin;   //Mass scale from P0 (kg)
	
	for(unsigned i = 0; i < stree.st.size(); i++)
	{
		stree.st[i].L0 = stree.st[i].L0/LL;          //convert lung length to unitless
		stree.st[i].A0 = stree.st[i].A0/(LL*LL);          //convert area to unitless
		stree.st[i].L1 = stree.st[i].L1 / LL;          //convert lung length to unitless
		stree.st[i].A1 = stree.st[i].A1 / (LL*LL);
		stree.st[i].V0 = stree.st[i].V0/(LL*LL*LL);
	}
	
	//convert to simulation units (length scale LL, time scale Tin, mass scale M0)
	o.Viscosity = (cmH20_Pa*o.Viscosity)*(LL*o.Tin / M0);
	o.Density = (o.Density*LL*LL*LL / M0);
	o.Diffusion = (cm_m*cm_m*o.Diffusion)*(o.Tin / (LL*LL));
	o.VFRC = o.VFRC / (LL*LL*LL);
	o.V0 = o.V0 / (LL*LL*LL);
	o.VD = o.VD / (LL*LL*LL);
	o.VDM = o.VDM / (LL*LL*LL);
	if (o.InputOp == PRESSURE) o.P0 = cmH20_Pa*o.P0*(LL*o.Tin*o.Tin/M0);     //convert pressure to kg m^(-1) s^(-2) then to unitless
	else o.VT = o.VT / (LL*LL*LL);
	o.Rb = (cmH20_Pa*o.Rb / L_m3)*(LL*LL*LL*LL*o.Tin / M0);                    //convert Rb to kg m^(-4) s^(-1)
	o.Rmouth = (cmH20_Pa*o.Rmouth / L_m3)*(LL*LL*LL*LL*o.Tin / M0);      //convert Rmouth to kg m^(-4) s^(-1)
	//--TO DO: UNits of k0--//
	
	//Normalised mass, length and timescales
	o.E = 1.0 / o.VFRC;
	o.Tin = 1.0;
	
	//Calculate run time etc.
	o.RunTime = total_run_time(o);
	
	return 0;
}

