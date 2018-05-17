#include "lung_model_discrete_branching.hpp"

//Checked and commented 16/05/18
double subtree::length_function(unsigned j, Options &o)  //branch length at subtree gen j (relative to start)
{
	double length = 0;

	//can add other options here
	switch (o.ShapeOp)
	{
	case LOBE_GEOMETRIC:
	case ALT_LOBE_GEOMETRIC:
	case GEOMETRIC:
	{
		if (j < Ncond) length = L0*pow(o.lambda, j);
		else length = L1*pow(o.lambda2, (j - Ncond));
	} break;

	case HOMOGENISED:
	{
		length = L0*pow(o.lambda, j);
	} break;

	case WEIBEL:
	{
		length = L0*weibel_length(j);
	} break;
	default:
	{
		cout << "Length option not recognised\n";
		length = 0;  //return 0 on error
	}
	}
	return length;
}

double subtree::area_function(double dx, unsigned j, Options &o) //returns tube cross sectional area at gen j, position dx
{
	double area = 0;
	//add an option here to interpolate areas linearly instead of by generation?

	switch (o.ShapeOp)
	{
	case GEOMETRIC:
	case LOBE_GEOMETRIC:
	case ALT_LOBE_GEOMETRIC:
	{
		if (j < Ncond) area = A0*pow(o.lambda, 2.0*j);
		else area = A1*pow(o.lambda2, 2.0*(j-Ncond));
	} break;

	case HOMOGENISED:
	{
		if (dx < 0) area = A0;
		else  area = A0*pow(o.lambda, 2.0*j)*pow((1 - dx*(1 - o.lambda) / (L0*pow(o.lambda, j))), 2 + log(2) / log(o.lambda));
	} break;

	case WEIBEL:
	{
		area = A0*weibel_area(j + StartGen);
	} break;

	default:
	{
		cout << "Area option not recognised\n";
		area = 0; //returns 0 on error
	}
	}

	return area;
}

double pressure_func(double t, Options &o)   //returns applied pressure as function of time or applied flow rate depending on input option
{
	double P;
	switch (o.PressOp)
	{
	case SIGMOIDAL:
	{
		int t0 = ((int) (t+0.5));
		double sigma = 50;
		if (o.InputOp == VOLUME) P = o.VT*((-1.0 + 2.0/(1.0 + exp(-sigma*(t-t0))))*(1-(t0%2)) + (1.0 - 2.0/(1.0 + exp(-sigma*(t-t0))))*(t0%2));
		else P = o.P0*(-(1 - (t0 % 2))*1.0 / (1.0 + exp(-sigma*(t - t0))) - (t0 % 2)*1.0 / (1.0 + exp(sigma*(t - t0))));
	} break;

	case STEP_FUNCTION:
	{
		if (((int)t) % 2)   //breathe out
		{
			if (o.InputOp == VOLUME) P = -o.VT;
			else P = 0.0;
		}
		else
		{
			if (o.InputOp == VOLUME) P = o.VT;
			else P = -o.P0;
		}
	} break;

	case SINUSOIDAL:
	{
		if (o.InputOp == VOLUME) P = 0.5*o.VT*M_PI*(sin(M_PI*t));
		else P = -0.5*o.P0*(1.0 + sin(M_PI*t));
	} break;

	case LINEAR:
	{
		if (o.InputOp == VOLUME) P = o.VT;
		else P = -o.P0*t;
	} break;

	default: P = 0;
	}

	return P;
}

double subtree::alveolar_density(unsigned j, unsigned k, Options &o)    //alveolar density in each gen for mean path subtree
{
	if (o.BcOp == BAG)   //bags only in final generation
	{
		if (j == Ncond)
		{
			return 1.0;
		}
		else return 0.0;
	}
	else   //weibel model
	{
		if (j < Ncond) return 0.0;
		if (j == Ncond) return 0.2;
		if (j == Ncond + 1) return 0.4;
		if (j == Ncond + 2) return 0.7;
		if (j > Ncond + 2) return 1.0;
	}
	cout << "Something went wrong in alveolar_density\n";
	return 0.0;
}

void subtree::tot_area_function(double Vtot) //reassigns tube + alveoli cross sectional area
{
	unsigned j, k;

	for (j = StartGen + Ncond; j <= EndGen; j++)
	{
		for (k = 0; k<gn[j].p.size(); k++)
		{
			gn[j].p[k].Aold = gn[j].p[k].Anew;
			gn[j].p[k].Anew = 0.5*(gn[j].p[k].al + gn[j].p[k].ar) + gn[j].p[k].DA*Vtot;
		}
	}
}

void subtree::tot_lp_area_function(subtree &sst, double DVtot) //reassigns tube + alveoli cross sectional area
{
	unsigned j, k;

	for (j = StartGen + Ncond; j <= EndGen; j++)
	{
		for (k = 0; k<gn[j].p.size(); k++)
		{
			gn[j].p[k].Aold = gn[j].p[k].Anew;
			gn[j].p[k].Anew = 0.5*(gn[j].p[k].al + gn[j].p[k].ar) 
				              + (gn[j].p[k].DA / sst.gn[j].p[k].DA)*(sst.gn[j].p[k].Anew - 0.5*(sst.gn[j].p[k].al + sst.gn[j].p[k].ar))
				              + sst.gn[j].p[k].DA*DVtot;
		}
	}
}

void subtree::velocity_calc(double qend, Options &o)   //calculates velocity on subtree from updated A
{
	double srtot, sul, uul;
	unsigned Nj;

	for (unsigned jc = 0; jc <= (EndGen - StartGen); jc++)   //moving up tree
	{
		unsigned j = EndGen - jc;     //counting up from end
		Nj = ((unsigned) gn[j].p.size());              //number of points in gen
		for (unsigned kc = 0; kc < Nj; kc++)   //moving up tree
		{
			unsigned k = Nj - 1 - kc;     //counting up from end
			srtot = 0;
			for (unsigned m = 0; m < gn[j].p[k].sr.size(); m++) srtot += gn[j].p[k].sr[m];
			if (k == Nj - 1 && j == EndGen)    //last cell, sul = sr (c-s of last point)
			{
				sul = srtot;
				uul = qend / srtot;  //velocity at left edge cell to the right
			}
			else
			{
				unsigned ju, ku;
				ju = gn[j].p[k].jup[2][0];
				ku = gn[j].p[k].kup[2][0];
				uul = gn[ju].p[ku].ul;                   //vel at left edge of k+1
				sul = gn[j].p[k].sul[0];    //c-s of left edge of k+1
			}
			if (srtot == 0.) gn[j].p[k].ur = 0;
			else gn[j].p[k].ur = sul*uul / srtot;         //ur diff to uul when there is a jump between cells
			if (gn[j].p[k].sl == 0.) gn[j].p[k].ul = 0;
			else gn[j].p[k].ul = (srtot*gn[j].p[k].ur / gn[j].p[k].sl) + gn[j].Nb*gn[j].dx*(gn[j].p[k].Anew - gn[j].p[k].Aold) / (o.dt*gn[j].p[k].sl);  //integrate over cell
			if (o.TaylorDisp==TAYLOR && j < StartGen + Ncond)
			{
				gn[j].p[k].Dl = o.Diffusion + gn[j].p[k].ul*gn[j].p[k].ul*gn[j].p[k].al / (192 * M_PI*o.Diffusion);
				gn[j].p[k].Dr = o.Diffusion + gn[j].p[k].ur*gn[j].p[k].ur*gn[j].p[k].ar / (192 * M_PI*o.Diffusion);
			}
			if (o.TaylorDisp == SCHERER && j < StartGen + Ncond)
			{
				double schcoeff;
				if (Vnew - Vold >= 0) schcoeff = 1.08;   //see Scherer 1975
				else schcoeff = 0.37;
				gn[j].p[k].Dl = o.Diffusion + schcoeff*fabs(gn[j].p[k].ul)*sqrt(gn[j].p[k].al*4.0/M_PI);
				gn[j].p[k].Dr = o.Diffusion + schcoeff*fabs(gn[j].p[k].ur)*sqrt(gn[j].p[k].ar*4.0/M_PI);
			}
		}
	}
}

void subtree::velocity_lp_calc(subtree &sst, double dqend, Options &o)  //calculates perturbed velocity on subtree from updated A
{
	double srtot, sul, uul;
	double dsrtot, dsul, duul;
	unsigned Nj;

	for (unsigned jc = 0; jc <= (EndGen - StartGen); jc++)   //moving up tree
	{
		unsigned j = EndGen - jc;     //counting up from end
		Nj = ((unsigned)gn[j].p.size());              //number of points in gen
		for (unsigned kc = 0; kc < Nj; kc++)   //moving up tree
		{
			unsigned k = Nj - 1 - kc;     //counting up from end
			srtot = sst.gn[j].p[k].ar*gn[j].Nb;    //cs area in symm tree
			dsrtot = 0;
			for (unsigned m = 0; m < gn[j].p[k].sr.size(); m++) dsrtot += gn[j].p[k].sr[m];   //cs area for pert tree
			if (k == Nj - 1 && j == EndGen)    //last cell, sul = sr (c-s of last point)
			{
				sul = srtot;
				dsul = dsrtot;
				uul = sst.gn[j].p[k].ur;
				duul = (dqend - dsrtot*sst.gn[j].p[k].ur) / srtot ;  //delta velocity at left edge cell to the right
			}
			else
			{
				unsigned ju, ku;
				ju = gn[j].p[k].jup[2][0];
				ku = gn[j].p[k].kup[2][0];
				uul = sst.gn[ju].p[ku].ul;                   //symm vel at left edge of k+1
				sul = sst.gn[ju].p[ku].al*gn[ju].Nb;                 //symm c-s of left edge of k+1
				duul = gn[ju].p[ku].ul;                   //delta vel at left edge of k+1
				dsul = gn[j].p[k].sul[0];                 //delta c-s of left edge of k+1
			}
			if (srtot == 0.)
			{
				if (dsrtot > 0) sst.gn[j].p[k].ur = (sul*duul + dsul*uul)/dsrtot;   //this is order epsilon -- always multiplied by dsr
				else sst.gn[j].p[k].ur = 0;
				gn[j].p[k].ur = 0;   //always multiplied by sr so value is unimportant
			}
			else gn[j].p[k].ur = (dsul*uul + sul*duul - (dsrtot/srtot)*sul*uul) / srtot;         //ur diff to uul when there is a jump between cells
			if (sst.gn[j].p[k].al == 0.)
			{
				if (gn[j].p[k].sl > 0) sst.gn[j].p[k].ul = (dsrtot*sst.gn[j].p[k].ur / gn[j].p[k].sl) + gn[j].Nb*sst.gn[j].dx*(gn[j].p[k].Anew - gn[j].p[k].Aold) / (o.dt*gn[j].p[k].sl);
				else sst.gn[j].p[k].ul = 0;
				gn[j].p[k].ul = 0;  //always multiplied by sl so value is unimportant
			}
			else gn[j].p[k].ul = (srtot*gn[j].p[k].ur + dsrtot*sst.gn[j].p[k].ur 
				                - srtot*sst.gn[j].p[k].ur*gn[j].p[k].sl / (gn[j].Nb*sst.gn[j].p[k].al)) / (gn[j].Nb*sst.gn[j].p[k].al)
								+ (sst.gn[j].dx*(gn[j].p[k].Anew - gn[j].p[k].Aold) 
								+ gn[j].dx*(sst.gn[j].p[k].Anew - sst.gn[j].p[k].Aold)
								- (gn[j].p[k].al / sst.gn[j].p[k].al)*sst.gn[j].dx*(sst.gn[j].p[k].Anew - sst.gn[j].p[k].Aold)) / (o.dt*sst.gn[j].p[k].al);  //integrate over cell
			if (o.TaylorDisp == TAYLOR && j < StartGen + Ncond)
			{
				gn[j].p[k].Dl = (2.0*gn[j].p[k].ul*sst.gn[j].p[k].ul*sst.gn[j].p[k].al
					             + sst.gn[j].p[k].ul*sst.gn[j].p[k].ul*gn[j].p[k].al) / (192 * M_PI*o.Diffusion);
				gn[j].p[k].Dr = (2.0*gn[j].p[k].ur*sst.gn[j].p[k].ur*sst.gn[j].p[k].ar 
					             + sst.gn[j].p[k].ur*sst.gn[j].p[k].ur*gn[j].p[k].ar) / (192 * M_PI*o.Diffusion);
			}
			if (o.TaylorDisp == SCHERER && j < StartGen + Ncond)
			{
				double schcoeff;
				if (sst.Vnew - sst.Vold >= 0) schcoeff = 1.08;   //see Scherer 1975
				else schcoeff = 0.37;
				//				gn[j].p[k].Dl = o.Diffusion + schcoeff*fabs(gn[j].p[k].ul)*sqrt(gn[j].p[k].al*4.0/M_PI);
				
				gn[j].p[k].Dl = schcoeff*sqrt(4.0 / M_PI)*(0.5*fabs(sst.gn[j].p[k].ul)*gn[j].p[k].al / sqrt(sst.gn[j].p[k].al)
					                                       + ((gn[j].p[k].ul*sst.gn[j].p[k].ul) / fabs(sst.gn[j].p[k].ul))*sqrt(sst.gn[j].p[k].al));
				gn[j].p[k].Dr = schcoeff*sqrt(4.0 / M_PI)*(0.5*fabs(sst.gn[j].p[k].ur)*gn[j].p[k].ar / sqrt(sst.gn[j].p[k].ar)
					                                       + ((gn[j].p[k].ur*sst.gn[j].p[k].ur) / fabs(sst.gn[j].p[k].ur))*sqrt(sst.gn[j].p[k].ar));
			}
		}
	}
}

double c_stim(double t, Options &o)   //inspiratory boundary condition for concentration
{
	if (t < o.StimTime) return 1.0;   //while connected to source
	if (o.N2leak.exists && o.N2leak.type != EXPIRATORY_LEAK)    //washout but leak exists on inspiration (N2)
	{
		if (t - o.StimTime >= o.N2leak.start && t - o.StimTime <= o.N2leak.end)
		{
			return o.N2leak.size;
		}
	}
	return 0.0;     //regular washout
}

double total_run_time(Options &o)
{
	return o.RunTime*o.Tin;   //want to add options for different breath profiles at some point
}

unsigned lobe_gens(unsigned nlobe, Options &o)  //returns gen number of final branch
{
	if(o.ShapeOp == LOBE_GEOMETRIC)
	{
		switch (nlobe)
		{
		case 0: return o.Ngen2 + 1;
			break;

		case BR: return o.Ngen2;
			break;

		case BL:
		case BRML:
		{
			return o.Ngen2 - 1;
		} break;

		case BRM: return o.Ngen2 - 3;
			break;

		default: return o.Ngen2 - 2;
			break;
		}
	}
	else
	{
		if (o.ShapeOp == ALT_LOBE_GEOMETRIC)
		{
			switch (nlobe)
			{
			case 0: return o.Ngen2 + 2;
				break;

			case BR: return o.Ngen2 + 1;
				break;

			case BL:
			case BRML:
			{
				return o.Ngen2;
			} break;

			case BRL:
			case BLL:
			{
				return o.Ngen2 - 1;
			} break;

			case BRM: return o.Ngen2 - 3;
				break;

			case BRL2: return o.Ngen2 - 4;
				break;

			case BLL2: return o.Ngen2 - 4;
				break;

			default: return o.Ngen2 - 2;
				break;
			}
		}
		else return o.Ngen2;
	}
}

unsigned long ijk_index(unsigned i, unsigned j, unsigned long k, Options &o)
{
	if (o.ShapeOp == LOBE_GEOMETRIC || o.ShapeOp == ALT_LOBE_GEOMETRIC)
	{
		unsigned im=0;
		unsigned index=0;
		while (im < i)
		{
			index += ((unsigned long) pow(2,lobe_gens(im,o)+1))-1;
			im++;
		}
		index += ((unsigned long) pow(2,j)) - 1 + k;
		return index;
	}
	else return (((unsigned long) pow(2,j)) - 1 + k);
	
}

double weibel_length(unsigned j) //returns branch length relative to L0
{
	double l = 0; //returns 0 on error

	switch (j)
	{
	case 0: l = 1;
		break;

	case 1: l = 0.3967;
		break;

	case 2: l = 0.1583;
		break;

	case 3: l = 0.06333;
		break;

	case 4: l = 0.1058;
		break;

	case 5: l = 0.08917;
		break;

	case 6: l = 0.07500;
		break;

	case 7: l = 0.06333;
		break;

	case 8: l = 0.05333;
		break;

	case 9: l = 0.04500;
		break;

	case 10: l = 0.03833;
		break;

	case 11: l = 0.03250;
		break;

	case 12: l = 0.02750;
		break;

	case 13: l = 0.02250;
		break;

	case 14: l = 0.01917;
		break;

	case 15: l = 0.01167;
		break;

	case 16: l = 0.01108;
		break;

	case 17: l = 0.01000;
		break;

	case 18: l = 0.007750;
		break;

	case 19: l = 0.006917;
		break;

	case 20:
	case 21:
	case 22:
	case 23: { l = 0.005833; }
			 break;
	}

	return l;
}

double weibel_area(unsigned j) //returns branch length relative to L0
{
	double a = 0;  //returns 0 on error

	switch (j)
	{
	case 0: a = 1;
		break;

	case 1: a = 0.4594;
		break;

	case 2: a = 0.2126;
		break;

	case 3: a = 0.09679;
		break;

	case 4: a = 0.06250;
		break;

	case 5: a = 0.03781;
		break;

	case 6: a = 0.02420;
		break;

	case 7: a = 0.01633;
		break;

	case 8: a = 0.01068;
		break;

	case 9: a = 0.007320;
		break;

	case 10: a = 0.005216;
		break;

	case 11: a = 0.003667;
		break;

	case 12: a = 0.002786;
		break;

	case 13: a = 0.002126;
		break;

	case 14: a = 0.001690;
		break;

	case 15: a = 0.0007716;
		break;

	case 16: a = 0.0007716;
		break;

	case 17: a = 0.0007410;
		break;

	case 18: a = 0.0004938;
		break;

	case 19: a = 0.0004457;
		break;

	case 20: a = 0.0004000;
		break;

	case 21: a = 0.0003568;
		break;

	case 22: a = 0.0002966;
		break;

	case 23: a = 0.0002596;
		break;
	}

	return a;
}


