#include "lung_model_discrete_branching.hpp"

void fill_Ab(Tree &stree, Eigen::SparseMatrix<double, Eigen::RowMajor> &AC, Eigen::VectorXd &XC, Eigen::VectorXd &BC, unsigned long Ntot, double time, Options &o)
{
	//functions to fill A and b in Ax=b for concentration update
	AC.reserve(4 + 3 * (Ntot - 2) + stree.st.size() - 1);

	//---------mouth bc-------------//
	stree.VDfrontoldold = stree.VDfrontold;
	stree.VDfrontold = stree.VDfront;  //at n*dt
	stree.VDfront += stree.st[0].Vnew - stree.st[0].Vold;   //at (n+1)*dt
	if (stree.VDfront < 0) stree.VDfront = 0;
	if (stree.VDfront > o.VDM) stree.VDfront = o.VDM;
	stree.Vtot.push_back(stree.st[0].Vnew);   //store old volume
	stree.cmouthold = stree.cmouth;          //store old conc
	double A0[2] = { 0, 0 }, b0 = 0, x0 = 0;
	unsigned t = ((unsigned)stree.Vtot.size() - 1);   //=n+1
	if (stree.st[0].Vnew - stree.st[0].Vold < 0)   //on exhalation
	{
		A0[0] = 1.0 / (stree.st[0].gn[0].p[0].Anew*((double)stree.st[0].gn[0].Nb));
		A0[1] = -1.0 / (stree.st[0].gn[0].p[1].Anew*((double)stree.st[0].gn[0].Nb));   //gradient condition
	}
	else   //on inhalation
	{
		A0[0] = 1.0;
		unsigned tc = t - 1;
		double ceff = 0;
		if (o.VDM > 1E-10)  //if mouth dead space is > 0
		{
			//work out characteristics - VDfront / o.VDM is proportion of VDM exhaled in this breath
			if ((1.0 - stree.VDfront / o.VDM) <= 1E-10)   //allow for rounding errors  -- assuming VT always > VDM, VDM filled
			{
				//count back to last time VD front was the same
				while (tc > 0 && (stree.Vtot[t] - (stree.Vtot[tc - 1] + o.VDM))* (stree.Vtot[t] - (stree.Vtot[tc] + o.VDM)) > 0) //when this > 0, both these terms have same sign, hence not the right point
				{
					tc--;
				}
				double ctc = (c_stim((tc - 1)*o.dt, o)* fabs(stree.Vtot[t] - (stree.Vtot[tc] + o.VDM)) + c_stim((tc - 1)*o.dt, o)* fabs(stree.Vtot[t] - (stree.Vtot[tc - 1] + o.VDM))) / (fabs(stree.Vtot[t] - (stree.Vtot[tc - 1] + o.VDM)) + fabs(stree.Vtot[t] - (stree.Vtot[tc] + o.VDM)));
				if ((1.0 - stree.VDfrontold / o.VDM) <= 1E-10)   //if the last time-point also had full VDM
				{
					ceff = ctc;
					if ((1.0 - stree.VDfrontoldold / o.VDM) > 1E-10)
					{
						tc = t - 1;
						while (tc > 0 && (stree.Vtot[t - 1] - (stree.Vtot[tc - 1] + o.VDM))* (stree.Vtot[t - 1] - (stree.Vtot[tc] + o.VDM)) > 0) //when this > 0, both these terms have same sign, hence not the right point
						{
							tc--;
						}
						double ctc1 = (c_stim((tc - 1)*o.dt, o)* fabs(stree.Vtot[t - 1] - (stree.Vtot[tc] + o.VDM)) + c_stim((tc - 1)*o.dt, o)* fabs(stree.Vtot[t - 1] - (stree.Vtot[tc - 1] + o.VDM))) / (fabs(stree.Vtot[t - 1] - (stree.Vtot[tc - 1] + o.VDM)) + fabs(stree.Vtot[t - 1] - (stree.Vtot[tc] + o.VDM)));
						stree.st[0].gn[0].p[0].cold = ctc1;  //change value of cold
					}
				}
				else  //need to compensate for missed advection 0.5*(Vnew - Vold)*(cnew + cold) = (Vnew - Vold)
				{
					ceff = (2.0*(o.VDM - stree.VDfrontold) / (stree.st[0].Vnew - stree.st[0].Vold) - 1.0)*stree.st[0].gn[0].p[0].cold
						+ 2.0*(1.0 - (o.VDM - stree.VDfrontold) / (stree.st[0].Vnew - stree.st[0].Vold))*ctc;
				}
			}
			else   //VDM already been exhaled
			{
				while (tc > 0 && (stree.Vtot[t] - stree.Vtot[tc - 1])*(stree.Vtot[t] - stree.Vtot[tc]) > 0) //when this > 0, both these terms have same sign, hence not the right point
				{
					tc--;
				}
				if (tc <= 1)
				{
					if (o.InitOp == FULL) ceff = 1.0;
					if (o.InitOp == EMPTY) ceff = 0.0;
				}
				else
				{
					ceff = (stree.ctop[tc - 1] * fabs(stree.Vtot[t] - stree.Vtot[tc]) + stree.ctop[tc] * fabs(stree.Vtot[t] - stree.Vtot[tc - 1])) / (fabs(stree.Vtot[t] - stree.Vtot[tc]) + fabs(stree.Vtot[t] - stree.Vtot[tc - 1]));
				}
			}
			b0 = ceff*stree.st[0].gn[0].p[0].Anew*((double)stree.st[0].gn[0].Nb);
			x0 = ceff*stree.st[0].gn[0].p[0].Anew*((double)stree.st[0].gn[0].Nb);
		}
		else //no mouth dead space included
		{
			b0 = c_stim(time + o.dt, o)*stree.st[0].gn[0].p[0].Anew*((double)stree.st[0].gn[0].Nb);
			x0 = c_stim(time + o.dt, o)*stree.st[0].gn[0].p[0].Anew*((double)stree.st[0].gn[0].Nb);
		}
		stree.cmouth = c_stim(time + o.dt, o); //cmouth at n+1
	}
	//--------End of mouth bc--------//

	for (unsigned i = 0; i < stree.st.size(); i++) //loop over subtrees
	{
		for (unsigned j = stree.st[i].StartGen; j <= stree.st[i].EndGen; j++)   //loop over generations
		{
			unsigned k0;
			if (j == 0) k0 = 1;  //start index (j=0 includes ghost point at 0)
			else k0 = 0;
			unsigned Nj = ((unsigned)stree.st[i].gn[j].p.size());  //no. of nodes in this gen
			unsigned ih, jh, kh;

			for (unsigned k = k0; k < Nj; k++)
			{
				vector<double> ucfr[3], ucfl[3], dcfr[3], dcfl[3];   //finite difference coeffs
				unsigned iright0 = stree.st[i].gn[j].p[k].iup[2][0];
				unsigned jright0 = stree.st[i].gn[j].p[k].jup[2][0];
				unsigned kright0 = stree.st[i].gn[j].p[k].kup[2][0];
				vector<double> cd0[3];
				unsigned nodes_right = ((unsigned)stree.st[i].gn[j].p[k].iup[2].size());
				for (unsigned n = 0; n < 3; n++)   //vector of positions of nearest 3 points
				{
					cd0[n].clear();
					ih = stree.st[i].gn[j].p[k].iup[n][0];  //indices of neighbouring points
					jh = stree.st[i].gn[j].p[k].jup[n][0];
					kh = stree.st[i].gn[j].p[k].kup[n][0];
					dcfr[n] = stree.st[i].gn[j].p[k].dcfr[n];
					dcfl[n] = stree.st[i].gn[j].p[k].dcfl[n];
					if (stree.st[i].gn[j].p[k].ul > 0) //ful determined by ul velocity
					{
						ucfl[n].push_back(stree.st[i].gn[j].p[k].ucflpos[n][0]);
					}
					else
					{
						ucfl[n].push_back(stree.st[i].gn[j].p[k].ucflneg[n][0]);
					}
					if (stree.st[iright0].gn[jright0].p[kright0].ul > 0)  //fur determined by ul at next point along (at bifuractions this can be different for m=0 and m=1
					{
						ucfr[n].push_back(stree.st[i].gn[j].p[k].ucfrpos[n][0]);    //upwind coefficients
					}
					else
					{
						ucfr[n].push_back(stree.st[i].gn[j].p[k].ucfrneg[n][0]);
					}
					cd0[n].push_back(stree.st[ih].gn[jh].p[kh].cold);     //conc at neighbour points
					if (n == 2 && nodes_right>1)
					{
						ih = stree.st[i].gn[j].p[k].iup[n][1];
						jh = stree.st[i].gn[j].p[k].jup[n][1];
						kh = stree.st[i].gn[j].p[k].kup[n][1];
						cd0[n].push_back(stree.st[ih].gn[jh].p[kh].cold);     //conc at neighbour points
					}
				}
				if(nodes_right>1)
				{
					for (unsigned n = 0; n < 3; n++)   //vector of positions of nearest 3 points
					{
						unsigned iright1 = stree.st[i].gn[j].p[k].iup[2][1];
						unsigned jright1 = stree.st[i].gn[j].p[k].jup[2][1];
						unsigned kright1 = stree.st[i].gn[j].p[k].kup[2][1];
						if (stree.st[i].gn[j].p[k].ul > 0)
						{
							ucfl[n].push_back(stree.st[i].gn[j].p[k].ucflpos[n][1]);
						}
						else
						{
							ucfl[n].push_back(stree.st[i].gn[j].p[k].ucflneg[n][1]);
						}
						if (stree.st[iright1].gn[jright1].p[kright1].ul > 0)  //fur determined by ul at next point along
						{
							ucfr[n].push_back(stree.st[i].gn[j].p[k].ucfrpos[n][1]);    //upwind coefficients
						}
						else
						{
							ucfr[n].push_back(stree.st[i].gn[j].p[k].ucfrneg[n][1]);
						}
					}
				}

				//---Work out left and right fluxes---//
				//flux at left boundary
				vector<double> Snew[3], Sold[3];     //total cs at right edge(k) and left edge(k+1)
				ih = stree.st[i].gn[j].p[k].iup[0][0];
				jh = stree.st[i].gn[j].p[k].jup[0][0];
				kh = stree.st[i].gn[j].p[k].kup[0][0];
				Snew[0].push_back(stree.st[ih].gn[jh].p[kh].Anew*((double)stree.st[ih].gn[jh].Nb));  //total cross-section of left element
				Snew[1].push_back(stree.st[i].gn[j].p[k].Anew*((double)stree.st[i].gn[j].Nb));  //total cross-section of centre element
				Sold[0].push_back(stree.st[ih].gn[jh].p[kh].Aold*((double)stree.st[ih].gn[jh].Nb));  //total cross-section of left element
				Sold[1].push_back(stree.st[i].gn[j].p[k].Aold*((double)stree.st[i].gn[j].Nb));  //total cross-section of centre element
				double ful, fdl;            //advective and diffusive flux at left edge(k)
				vector<double> fur, fdr;    //advective and diffusive flux at right edge(k)

				ful = stree.st[i].gn[j].p[k].ul*stree.st[i].gn[j].p[k].sl;

				double sdrtot = 0;
				double srtot = 0;
				for (unsigned m = 0; m < stree.st[ih].gn[jh].p[kh].sr.size(); m++) sdrtot += stree.st[ih].gn[jh].p[kh].sr[m];
				for (unsigned m = 0; m < nodes_right; m++) srtot += stree.st[i].gn[j].p[k].sr[m];
				double Seffdr;
				if (sdrtot> 0) Seffdr = stree.st[i].gn[j].p[k].sdr*(1.0 + o.AcinAreaFactor*(0.5*(Snew[0][0] + Sold[0][0]) - sdrtot) / sdrtot);   //work out effective radii increases
				else Seffdr = 0.5*o.AcinAreaFactor*(Snew[0][0] + Sold[0][0]) / (stree.st[ih].gn[jh].p[kh].sr.size());
				double Seffl = stree.st[i].gn[j].p[k].sl + o.AcinAreaFactor*(0.5*(Snew[1][0] + Sold[1][0]) - stree.st[i].gn[j].p[k].sl);
				if (Seffl <= Seffdr) fdl = -0.5*(stree.st[i].gn[j].p[k].Dl + stree.st[ih].gn[jh].p[kh].Dr)*Seffl;
				else fdl = -0.5*(stree.st[i].gn[j].p[k].Dl + stree.st[ih].gn[jh].p[kh].Dr)*Seffdr;
				//                fdl = -0.5*(stree.st[i].gn[j].p[k].Dl*stree.st[i].gn[j].p[k].sl + stree.st[ih].gn[jh].p[kh].Dr*stree.st[i].gn[j].p[k].sdr);

				vector<double> Seffr, Sefful;
				for (unsigned m = 0; m < nodes_right; m++) //flux(es) at right boundary
				{
					ih = stree.st[i].gn[j].p[k].iup[2][m];
					jh = stree.st[i].gn[j].p[k].jup[2][m];
					kh = stree.st[i].gn[j].p[k].kup[2][m];
					Snew[2].push_back(stree.st[ih].gn[jh].p[kh].Anew*((double)stree.st[ih].gn[jh].Nb)); //cs of right element
					Sold[2].push_back(stree.st[ih].gn[jh].p[kh].Aold*((double)stree.st[ih].gn[jh].Nb)); //cs of right element
					fur.push_back(stree.st[ih].gn[jh].p[kh].ul*stree.st[i].gn[j].p[k].sul[m]);
					if (srtot>0) Seffr.push_back(stree.st[i].gn[j].p[k].sr[m] * (1.0 + o.AcinAreaFactor*(0.5*(Snew[1][0] + Sold[1][0]) - srtot) / srtot));   //work out effective radii increases
					else Seffr.push_back(0.5*o.AcinAreaFactor*(Snew[1][0] + Sold[1][0]) / (stree.st[i].gn[j].p[k].sr.size())); //duct cs = 0, diffusion area just split evenly between downstream ducts
					Sefful.push_back(stree.st[i].gn[j].p[k].sul[m] + o.AcinAreaFactor*(0.5*(Snew[2][m] + Sold[2][m]) - stree.st[i].gn[j].p[k].sul[m]));
					if (Sefful[m] <= Seffr[m]) fdr.push_back(-0.5*(stree.st[ih].gn[jh].p[kh].Dl + stree.st[i].gn[j].p[k].Dr)*Sefful[m]);
					else fdr.push_back(-0.5*(stree.st[ih].gn[jh].p[kh].Dl + stree.st[i].gn[j].p[k].Dr)*Seffr[m]);
					//                    fdr.push_back(-0.5*(stree.st[i].gn[j].p[k].Dr*stree.st[i].gn[j].p[k].sr[m] + stree.st[ih].gn[jh].p[kh].Dl*stree.st[i].gn[j].p[k].sul[m]));
					if (k == Nj - 1 && j == stree.st[stree.st[i].imeanpath].StartGen + stree.st[stree.st[i].imeanpath].Ntot)   //last point in system
					{
						fur[m] = 0.0;
						fdr[m] = 0.0;
					}
				}

				unsigned long kmh = stree.st[i].gn[j].p[k].km;
				double Ah[4] = { 0, 0, 0, 0 }, bh = 0, xh = 0;
				if (i == 0 && j == 0 && k == 1) fdl = 0;   //no diffusion out of mouth
				for (unsigned n = 0; n < 2; n++)   //fill matrix for left and central term
				{
					if (n == 1)
					{
						Ah[n] += 1.0;
						bh += (((double)stree.st[i].gn[j].Nb)*stree.st[i].gn[j].p[k].Aold)*cd0[n][0];
						xh += (((double)stree.st[i].gn[j].Nb)*stree.st[i].gn[j].p[k].Aold)*cd0[n][0];
						for (unsigned m = 0; m < nodes_right; m++)    //contribution from right flux
						{
							if (Snew[n][0] > 0.) Ah[n] += -0.5*o.dt*((-ucfr[n][m] * fur[m]) + (-dcfr[n][m] * fdr[m])) / (Snew[n][0] * stree.st[i].gn[j].dx);

							bh += 0.5*o.dt*cd0[n][0] * ((-ucfr[n][m] * fur[m]) + (-dcfr[n][m] * fdr[m])) / (stree.st[i].gn[j].dx);

							xh += o.dt*cd0[n][0] * ((-ucfr[n][m] * fur[m]) + (-dcfr[n][m] * fdr[m])) / (stree.st[i].gn[j].dx);
						}
					}
					//contribution from ful
					if (Snew[n][0] > 0.) Ah[n] += -0.5*o.dt*((ucfl[n][0] * ful) + (dcfl[n][0] * fdl)) / (Snew[n][0] * stree.st[i].gn[j].dx);

					bh += 0.5*o.dt*cd0[n][0] * ((ucfl[n][0] * ful) + (dcfl[n][0] * fdl)) / (stree.st[i].gn[j].dx);

					xh += o.dt*cd0[n][0] * ((ucfl[n][0] * ful) + (dcfl[n][0] * fdl)) / (stree.st[i].gn[j].dx);
				}

				for (unsigned m = 0; m < nodes_right; m++)    //contribution from fur, fill in for right term
				{
					unsigned n = 2;
					if (Snew[n][m] > 0.) Ah[n + m] += -0.5*o.dt*((-ucfr[n][m] * fur[m]) + (-dcfr[n][m] * fdr[m])) / (Snew[n][m] * stree.st[i].gn[j].dx);

					bh += 0.5*o.dt*cd0[n][m] * ((-ucfr[n][m] * fur[m]) + (-dcfr[n][m] * fdr[m])) / (stree.st[i].gn[j].dx);

					xh += o.dt*cd0[n][m] * ((-ucfr[n][m] * fur[m]) + (-dcfr[n][m] * fdr[m])) / (stree.st[i].gn[j].dx);
				}

				if(o.BcOp == SINK) //set c = 0 BC (for testing only)
				{
					if (j == stree.st[stree.st[i].imeanpath].StartGen + stree.st[stree.st[i].imeanpath].Ntot)   //last point in system
					{
						Ah[0] = 0;
						Ah[1] = 1.0;
						Ah[2] = 0;
						Ah[3] = 0;
						xh = 0;
						bh = 0;
					}
				}

					//set matrix
				unsigned long kmdwn = stree.st[stree.st[i].gn[j].p[k].iup[0][0]].gn[stree.st[i].gn[j].p[k].jup[0][0]].p[stree.st[i].gn[j].p[k].kup[0][0]].km;
				unsigned long kmup0 = stree.st[stree.st[i].gn[j].p[k].iup[2][0]].gn[stree.st[i].gn[j].p[k].jup[2][0]].p[stree.st[i].gn[j].p[k].kup[2][0]].km;

				//Add this line of equation to matrix AC
				AC.insert(kmh, kmdwn) = Ah[0];
				AC.insert(kmh, kmh) = Ah[1];
				if (kmup0 != kmh)
				{
					AC.insert(kmh, kmup0) = Ah[2];
				}
				else
				{
					AC.coeffRef(kmh, kmh) += Ah[2];
				}
				if (nodes_right > 1)
				{
					unsigned long kmup1 = stree.st[stree.st[i].gn[j].p[k].iup[2][1]].gn[stree.st[i].gn[j].p[k].jup[2][1]].p[stree.st[i].gn[j].p[k].kup[2][1]].km;
					if (kmup1 != kmh)
					{
						AC.insert(kmh, kmup1) = Ah[3];
					}
					else
					{
						AC.coeffRef(kmh, kmh) += Ah[3];
					}
				}
				BC(kmh) = bh;
				XC(kmh) = xh;
			}
		}
	}
	//cout << AC << '\n';


	if (stree.st[0].Vnew - stree.st[0].Vold < 0)
	{
		x0 = XC(1);  //guess on exhalation
	}
	AC.insert(0, 0) = A0[0];
	AC.insert(0, 1) = A0[1];
	BC(0) = b0;
	XC(0) = x0;
}

void calc_gas_concs(Tree &stree, Eigen::VectorXd &XC, Options &o)
{
	//count up total volume of gas in tree
	stree.masstot = 0;
	unsigned t = ((unsigned)stree.Vtot.size() - 1);
	for (unsigned i = 0; i < stree.st.size(); i++)
	{
		stree.st[i].masstot = 0;
		for (unsigned j = stree.st[i].StartGen; j <= stree.st[i].EndGen; j++)
		{
			unsigned Nj = ((unsigned)stree.st[i].gn[j].p.size());
			for (unsigned k = 0; k < Nj; k++)
			{
				//if (stree.st[i].gn[j].p[k].Anew > 0) stree.st[i].gn[j].p[k].c = gsl_vector_get(x, stree.st[i].gn[j].p[k].km) / (stree.st[i].gn[j].Nb*stree.st[i].gn[j].p[k].Anew);
				if (stree.st[i].gn[j].p[k].Anew > 0) stree.st[i].gn[j].p[k].c = XC(stree.st[i].gn[j].p[k].km) / (((double)stree.st[i].gn[j].Nb)*stree.st[i].gn[j].p[k].Anew);
				else stree.st[i].gn[j].p[k].c = stree.st[i].gn[j].p[k].cold;
				if (j > 0 || k > 0) stree.st[i].masstot += ((double)stree.st[i].gn[j].Nb)*(stree.st[i].gn[j].p[k].Anew*stree.st[i].gn[j].p[k].c*stree.st[i].gn[j].dx);
			}
		}
		stree.masstot += stree.st[i].masstot;
	}

	//calculate vaule of cmouth (if exhaling)
	stree.ctop.push_back(stree.st[0].gn[0].p[0].c);   //trach value of c at top of tree
	if (stree.st[0].Vnew - stree.st[0].Vold < 0)
	{
		unsigned tc = t - 1;
		double ceff = 0;
		double Vnh = stree.Vtot[t];
		if (o.VDM > 1E-10)
		{
			if (stree.VDfront > 1E-10)   //front still in dead-space
			{
				while (tc > 0 && (Vnh - stree.Vtot[tc - 1])* (Vnh - stree.Vtot[tc]) > 0) //when this > 0, both these terms have same sign, hence not the right point
				{
					tc--;
				}
				ceff = (c_stim((tc - 1)*o.dt, o)* fabs(Vnh - stree.Vtot[tc]) + c_stim(tc*o.dt, o)* fabs(Vnh - stree.Vtot[tc - 1])) / (fabs(Vnh - stree.Vtot[tc - 1]) + fabs(Vnh - stree.Vtot[tc]));
			}
			else
			{
				while (tc > 0 && (Vnh - (stree.Vtot[tc - 1] - o.VDM))* (Vnh - (stree.Vtot[tc] - o.VDM)) > 0) //when this > 0, both these terms have same sign, hence not the right point
				{
					tc--;
				}
				if (tc <= 1)
				{
					if (o.InitOp == FULL) ceff = 1.0;
					if (o.InitOp == EMPTY) ceff = 0.0;
				}
				else
				{
					ceff = (stree.ctop[tc - 1] * fabs(Vnh - (stree.Vtot[tc] - o.VDM)) + stree.ctop[tc] * fabs(Vnh - (stree.Vtot[tc - 1] - o.VDM))) / (fabs(Vnh - (stree.Vtot[tc] - o.VDM)) + fabs(Vnh - (stree.Vtot[tc - 1] - o.VDM)));
				}
			}
			stree.cmouth = ceff; //cmouth at n+1
		}
		else stree.cmouth = stree.st[0].gn[0].p[0].c;
	}

	//calculate flux in at top of subtree and whole tree
	for (unsigned i = 0; i < stree.st.size(); i++)
	{
		unsigned j = stree.st[i].StartGen;
		unsigned ih, jh, kh;
		unsigned k0;
		if (j == 0) k0 = 1;
		else k0 = 0;
		vector<double> ucfr, ucfl, dcfr, dcfl;
		stree.st[i].fluxin = 0;
		for (unsigned n = 0; n < 3; n++)
		{
			ih = stree.st[i].gn[j].p[k0].iup[n][0];
			jh = stree.st[i].gn[j].p[k0].jup[n][0];
			kh = stree.st[i].gn[j].p[k0].kup[n][0];
			dcfr = stree.st[i].gn[j].p[k0].dcfr[n];
			dcfl = stree.st[i].gn[j].p[k0].dcfl[n];
			if (stree.st[i].gn[j].p[k0].ul > 0)
			{
				ucfr = stree.st[i].gn[j].p[k0].ucfrpos[n];
				ucfl = stree.st[i].gn[j].p[k0].ucflpos[n];
			}
			else
			{
				ucfr = stree.st[i].gn[j].p[k0].ucfrneg[n];
				ucfl = stree.st[i].gn[j].p[k0].ucflneg[n];
			}
			if (i == 0) dcfl[0] = 0; //turning off diffusion at mouth
			stree.st[i].fluxin += 0.5*o.dt*((double)stree.st[i].gn[j].Nb)*stree.st[i].gn[j].p[k0].al*(stree.st[ih].gn[jh].p[kh].cold + stree.st[ih].gn[jh].p[kh].c)*(ucfl[0] * stree.st[i].gn[j].p[k0].ul - dcfl[0] * stree.st[i].gn[j].p[k0].Dl);
		}
		//        double Check1 = stree.st[i].fluxin;
		//        double Check2 = (stree.st[i].masstot - mtotold[i]);
		//        if(fabs(Check1 - Check2) > 1E-11)
		//        {
		//            cout << "Here\n";
		//        }
	}
	stree.fluxin = stree.st[0].fluxin;

	//calculate inert gas vol in alveoli
	for (unsigned m = 0; m < stree.EndSubtrees.size(); m++)  //loop trees that reach acinus
	{
		unsigned i = stree.EndSubtrees[m];   //current subtree
		stree.st[i].totIGvol = 0;
		for (unsigned n = 0; n < stree.st[i].isub.size(); n++)  //loop over subtrees
		{
			if (i != stree.st[i].isub[n])stree.st[stree.st[i].isub[n]].totIGvol = 0;
			for (unsigned j = stree.st[stree.st[i].isub[n]].StartGen + stree.st[stree.st[i].isub[n]].Ncond; j <= stree.st[stree.st[i].isub[n]].EndGen; j++)
			{
				for (unsigned k = 0; k < stree.st[stree.st[i].isub[n]].gn[j].p.size(); k++)
				{
					stree.st[stree.st[i].isub[n]].totIGvol += stree.st[stree.st[i].isub[n]].gn[j].p[k].c*stree.st[stree.st[i].isub[n]].gn[j].p[k].Anew*stree.st[stree.st[i].isub[n]].gn[j].dx*((double)stree.st[stree.st[i].isub[n]].gn[j].Nb);
				}
			}
			if (i != stree.st[i].isub[n]) stree.st[i].totIGvol += stree.st[stree.st[i].isub[n]].totIGvol;
		}
	}
}


void fill_Ab_lp(Tree &ltree, Tree &stree, Eigen::SparseMatrix<double, Eigen::RowMajor> &AC, Eigen::VectorXd &XC, Eigen::VectorXd &BC, unsigned long Ntot, Options &o)
{
	//fill A matrix for concentration update on perturbed tree ltree
	AC.reserve(4 + 3 * (Ntot - 2) + ltree.st.size() - 1);
	double A0[2] = { 0, 0 }, b0 = 0, x0 = 0;
	unsigned t = ((unsigned)stree.Vtot.size() - 1);

	//---------mouth bc-------------//
	if (stree.st[0].Vnew - stree.st[0].Vold < 0)   //BC depends on stree flow
	{
		A0[0] = 1.0 / (stree.st[0].gn[0].p[0].Anew*((double)ltree.st[0].gn[0].Nb));
		A0[1] = -1.0 / (stree.st[0].gn[0].p[1].Anew*((double)ltree.st[0].gn[0].Nb));   //gradient condition
	}
	else
	{
		//position of volume front does not change -- flow rate the same and DS same as before
		A0[0] = 1.0;
		unsigned tc = t - 1;
		double dceff = 0;
		if (o.VDM > 1E-10)
		{
			if ((1.0 - stree.VDfront / o.VDM) <= 1E-10)   //allow for rounding errors  -- assuming VT always > VDM, VDM filled
			{
				//if mouth DS filled, pert same as stree value
				if ((1.0 - stree.VDfrontold / o.VDM) <= 1E-10) //if the last time-point also had full VDM
				{
					dceff = 0;
					if ((1.0 - stree.VDfrontoldold / o.VDM) > 1E-10) ltree.st[0].gn[0].p[0].cold = 0;  //change value of cold
				}
				else  //need to compensate for missed advection 0.5*(Vnew - Vold)*(cnew + cold) = (Vnew - Vold)
				{
					dceff = (2.0*(o.VDM - stree.VDfrontold) / (stree.st[0].Vnew - stree.st[0].Vold) - 1.0)*ltree.st[0].gn[0].p[0].cold;
				}
			}
			else
			{
				while (tc > 0 && (stree.Vtot[t] - stree.Vtot[tc - 1])*(stree.Vtot[t] - stree.Vtot[tc]) > 0) //when this > 0, both these terms have same sign, hence not the right point
				{
					tc--;
				}
				if (tc <= 1)
				{
					dceff = 0; //no change in c at time 0
				}
				else
				{
					dceff = (ltree.ctop[tc - 1] * fabs(stree.Vtot[t] - stree.Vtot[tc]) + ltree.ctop[tc] * fabs(stree.Vtot[t] - stree.Vtot[tc - 1])) / (fabs(stree.Vtot[t] - stree.Vtot[tc]) + fabs(stree.Vtot[t] - stree.Vtot[tc - 1]));
				}
			}
			b0 = dceff*stree.st[0].gn[0].p[0].Anew*((double)ltree.st[0].gn[0].Nb);
			x0 = dceff*stree.st[0].gn[0].p[0].Anew*((double)ltree.st[0].gn[0].Nb);
		}
		ltree.cmouth = 0;    //cmouth unchanged for inhalation
	}
	//--------End of mouth bc--------//
	for (unsigned i = 0; i < ltree.st.size(); i++)
	{
		unsigned istree = ltree.st[i].i_stree;   //corresponding st number in stree
		for (unsigned j = ltree.st[i].StartGen; j <= ltree.st[i].EndGen; j++)
		{
			unsigned k0;
			if (j == 0) k0 = 1;
			else k0 = 0;
			unsigned Nj = ((unsigned)ltree.st[i].gn[j].p.size());
			unsigned ih, ihs, jh, kh;
			for (unsigned k = k0; k < Nj; k++)
			{
				vector<double> dcd0[3], cd0[3], cd0old[3];
				vector<double> ucfr[3], ucfl[3], dcfr[3], dcfl[3];
				vector<double> ducfr[3], ducfl[3], ddcfr[3], ddcfl[3];
				unsigned istreeright0 = stree.st[istree].gn[j].p[k].iup[2][0];
				unsigned jright0 = ltree.st[i].gn[j].p[k].jup[2][0];
				unsigned kright0 = ltree.st[i].gn[j].p[k].kup[2][0];
				unsigned stree_nodes_right = ((unsigned) stree.st[istree].gn[j].p[k].iup[2].size());
				unsigned nodes_right = ((unsigned)ltree.st[i].gn[j].p[k].iup[2].size());
				for (unsigned n = 0; n < 3; n++)  //vector of positions of nearest 3 points
				{
					cd0[n].clear();
					cd0old[n].clear();
					dcd0[n].clear();
					for (unsigned m = 0; m < ltree.st[i].gn[j].p[k].iup[n].size(); m++)
					{
						ih = ltree.st[i].gn[j].p[k].iup[n][m];
						if (stree.st[istree].gn[j].p[k].iup[n].size() > 1) ihs = stree.st[istree].gn[j].p[k].iup[n][m];
						else ihs = stree.st[istree].gn[j].p[k].iup[n][0];
						jh = ltree.st[i].gn[j].p[k].jup[n][m];
						kh = ltree.st[i].gn[j].p[k].kup[n][m];
						cd0[n].push_back(stree.st[ihs].gn[jh].p[kh].c);     //stree conc at neighbour points
						cd0old[n].push_back(stree.st[ihs].gn[jh].p[kh].cold);     //old stree conc at neighbour points
						dcd0[n].push_back(ltree.st[ih].gn[jh].p[kh].cold);     //old ltree conc at neighbour points
					}
					dcfr[n] = stree.st[istree].gn[j].p[k].dcfr[n];
					dcfl[n] = stree.st[istree].gn[j].p[k].dcfl[n];
					ddcfr[n] = ltree.st[i].gn[j].p[k].dcfr[n];
					ddcfl[n] = ltree.st[i].gn[j].p[k].dcfl[n];
					if (stree.st[istree].gn[j].p[k].ul > 0) //upwind determined by ul velocity
					{
						ucfl[n].push_back(stree.st[istree].gn[j].p[k].ucflpos[n][0]);
						ducfl[n].push_back(ltree.st[i].gn[j].p[k].ucflpos[n][0]);
					}
					else
					{
						ucfl[n].push_back(stree.st[istree].gn[j].p[k].ucflneg[n][0]);
						ducfl[n].push_back(ltree.st[i].gn[j].p[k].ucflneg[n][0]);
					}
					if (stree.st[istreeright0].gn[jright0].p[kright0].ul > 0)  //fur determined by ul at next point along (at bifuractions this can be different for m=0 and m=1
					{
						ucfr[n].push_back(stree.st[istree].gn[j].p[k].ucfrpos[n][0]);    //upwind coefficients
						ducfr[n].push_back(ltree.st[i].gn[j].p[k].ucfrpos[n][0]);    //upwind coefficients
					}
					else
					{
						ucfr[n].push_back(stree.st[istree].gn[j].p[k].ucfrneg[n][0]);
						ducfr[n].push_back(ltree.st[i].gn[j].p[k].ucfrneg[n][0]);
					}
				}

				if (nodes_right > 1)   //re-size stree vectors to account for new tree structure
				{
					unsigned jright1 = ltree.st[i].gn[j].p[k].jup[2][1];
					unsigned kright1 = ltree.st[i].gn[j].p[k].kup[2][1];
					if (stree_nodes_right > 1)
					{
						unsigned istreeright1 = stree.st[istree].gn[j].p[k].iup[2][1];
						for (unsigned n = 0; n < 3; n++)
						{
							if (stree.st[istree].gn[j].p[k].ul > 0)
							{
								ucfl[n].push_back(stree.st[istree].gn[j].p[k].ucflpos[n][1]);
								ducfl[n].push_back(ltree.st[i].gn[j].p[k].ucflpos[n][1]);
							}
							else
							{
								ucfl[n].push_back(stree.st[istree].gn[j].p[k].ucflneg[n][1]);
								ducfl[n].push_back(ltree.st[i].gn[j].p[k].ucflneg[n][1]);
							}
							if (stree.st[istreeright1].gn[jright1].p[kright1].ul > 0)  //fur determined by ul at next point along
							{
								ucfr[n].push_back(stree.st[istree].gn[j].p[k].ucfrpos[n][1]);    //upwind coefficients
								ducfr[n].push_back(ltree.st[i].gn[j].p[k].ucfrpos[n][1]);    //upwind coefficients
							}
							else
							{
								ucfr[n].push_back(stree.st[istree].gn[j].p[k].ucfrneg[n][1]);
								ducfr[n].push_back(ltree.st[i].gn[j].p[k].ucfrneg[n][1]);
							}
						}
					}
					else
					{
						for (unsigned n = 0; n < 3; n++)
						{
							if (stree.st[istree].gn[j].p[k].ul > 0)
							{
								ducfl[n].push_back(ltree.st[i].gn[j].p[k].ucflpos[n][1]);
							}
							else
							{
								ducfl[n].push_back(ltree.st[i].gn[j].p[k].ucflneg[n][1]);
							}
							if (stree.st[istreeright0].gn[jright0].p[kright0].ul > 0)  //fur determined by ul at next point along
							{
								ducfr[n].push_back(ltree.st[i].gn[j].p[k].ucfrpos[n][1]);    //upwind coefficients
							}
							else
							{
								ducfr[n].push_back(ltree.st[i].gn[j].p[k].ucfrneg[n][1]);
							}
							ucfr[n].push_back(ucfr[n][0]);
							ucfl[n].push_back(ucfl[n][0]);
							dcfr[n].push_back(dcfr[n][0]);
							dcfl[n].push_back(dcfl[n][0]);
						}
					}
				}

				//flux at left boundary
				vector<double> Snew[3], dSnew[3];     //total cs at right edge(k) and left edge(k+1)
				ih = ltree.st[i].gn[j].p[k].iup[0][0];
				ihs = stree.st[istree].gn[j].p[k].iup[0][0];
				jh = ltree.st[i].gn[j].p[k].jup[0][0];
				kh = ltree.st[i].gn[j].p[k].kup[0][0];
				//total cross sec in stree
				Snew[0].push_back(stree.st[ihs].gn[jh].p[kh].Anew*((double)ltree.st[ih].gn[jh].Nb));  //total cross-section of left element
				Snew[1].push_back(stree.st[istree].gn[j].p[k].Anew*((double)ltree.st[i].gn[j].Nb));  //total cross-section of centre element
				//total cross sec in ltree
				dSnew[0].push_back(ltree.st[ih].gn[jh].p[kh].Anew*((double)ltree.st[ih].gn[jh].Nb));  //total cross-section of left element
				dSnew[1].push_back(ltree.st[i].gn[j].p[k].Anew*((double)ltree.st[i].gn[j].Nb));
				double ful, fdl, dful, dfdl;            //advective and diffusive flux at left edge(k)
				vector<double> fur, fdr, dfur, dfdr;    //advective and diffusive flux at right edge(k)

				//advective flux left
				ful = stree.st[istree].gn[j].p[k].ul*stree.st[istree].gn[j].p[k].al*((double)ltree.st[i].gn[j].Nb);
				dful = ltree.st[i].gn[j].p[k].ul*stree.st[istree].gn[j].p[k].al*((double)ltree.st[i].gn[j].Nb)
					+ stree.st[istree].gn[j].p[k].ul*ltree.st[i].gn[j].p[k].sl;

				//sr and sdr can be split in two, work out total
				double sdrtot = 0;
				double dsdrtot = 0;
				double srtot = 0;
				double dsrtot = 0;
				for (unsigned m = 0; m < stree.st[ihs].gn[jh].p[kh].sr.size(); m++)
				{
					sdrtot += (((double)ltree.st[ih].gn[jh].Nb) / ((double)stree.st[ihs].gn[jh].Nb))*stree.st[ihs].gn[jh].p[kh].sr[m];
				}
				for (unsigned m = 0; m < stree.st[istree].gn[j].p[k].sr.size(); m++)
				{
					srtot += (((double)ltree.st[i].gn[j].Nb) / ((double)stree.st[istree].gn[j].Nb))*stree.st[istree].gn[j].p[k].sr[m];
				}
				for (unsigned m = 0; m < ltree.st[ih].gn[jh].p[kh].sr.size(); m++)
				{
					dsdrtot += ltree.st[ih].gn[jh].p[kh].sr[m];
				}
				for (unsigned m = 0; m < nodes_right; m++)
				{
					dsrtot += ltree.st[i].gn[j].p[k].sr[m];
				}

				//Seff is effectiv diffusive cross section
				double Seffdr, dSeffdr;
				double sdrscale = ((double)ltree.st[ih].gn[jh].Nb*stree.st[ihs].gn[jh].p[kh].sr.size())
					/ ((double)stree.st[ihs].gn[jh].Nb*ltree.st[ih].gn[jh].p[kh].sr.size());  //converts stree sdr to ltree structure
				double slscale = ((double)ltree.st[i].gn[j].Nb) / ((double)stree.st[istree].gn[j].Nb);   //converts stree sl to ltree structure

				if (sdrtot> 0) Seffdr = sdrscale*stree.st[istree].gn[j].p[k].sdr*(1.0 + o.AcinAreaFactor*(0.5*(Snew[0][0] + stree.st[ihs].gn[jh].p[kh].Aold*((double)ltree.st[ih].gn[jh].Nb)) / sdrtot - 1.0));   //work out effective radii increases
				else Seffdr = 0.5*o.AcinAreaFactor*(Snew[0][0] + stree.st[ihs].gn[jh].p[kh].Aold*((double)ltree.st[ih].gn[jh].Nb))
					/ (ltree.st[ih].gn[jh].p[kh].sr.size());

				double Seffl = slscale*stree.st[istree].gn[j].p[k].sl + o.AcinAreaFactor*(0.5*(Snew[1][0] +
					stree.st[istree].gn[j].p[k].Aold*((double)ltree.st[i].gn[j].Nb))
					- slscale*stree.st[istree].gn[j].p[k].sl);

				if (sdrtot > 0)
				{
					dSeffdr = ltree.st[i].gn[j].p[k].sdr*(1.0 + o.AcinAreaFactor*(0.5*(Snew[0][0] + stree.st[ihs].gn[jh].p[kh].Aold*((double)ltree.st[ih].gn[jh].Nb)) / sdrtot - 1.0))
						+ 0.5*sdrscale*stree.st[istree].gn[j].p[k].sdr*o.AcinAreaFactor*(dSnew[0][0] + ltree.st[ih].gn[jh].p[kh].Aold*((double)ltree.st[ih].gn[jh].Nb)) / sdrtot
						- 0.5*sdrscale*stree.st[istree].gn[j].p[k].sdr*o.AcinAreaFactor*dsdrtot*(Snew[0][0] + stree.st[ihs].gn[jh].p[kh].Aold*((double)ltree.st[ih].gn[jh].Nb)) / (sdrtot*sdrtot);   //work out effective radii increases
				}
				else
				{
					dSeffdr = 0.5*o.AcinAreaFactor*(dSnew[0][0] + ltree.st[ih].gn[jh].p[kh].Aold*((double)ltree.st[ih].gn[jh].Nb)) / (ltree.st[ih].gn[jh].p[kh].sr.size());
				}

				double dSeffl = ltree.st[i].gn[j].p[k].sl + o.AcinAreaFactor*(0.5*(dSnew[1][0] + ltree.st[i].gn[j].p[k].Aold*((double)ltree.st[i].gn[j].Nb)) - ltree.st[i].gn[j].p[k].sl);

				//diffusive flux left
				if (Seffl <= Seffdr)
				{
					fdl = -0.5*(stree.st[istree].gn[j].p[k].Dl + stree.st[ihs].gn[jh].p[kh].Dr)*Seffl;
					dfdl = -0.5*(ltree.st[i].gn[j].p[k].Dl + ltree.st[ih].gn[jh].p[kh].Dr)*Seffl - 0.5*(stree.st[istree].gn[j].p[k].Dl + stree.st[ihs].gn[jh].p[kh].Dr)*dSeffl;
				}
				else
				{
					fdl = -0.5*(stree.st[ihs].gn[jh].p[kh].Dr + stree.st[istree].gn[j].p[k].Dl)*Seffdr;
					dfdl = -0.5*(ltree.st[ih].gn[jh].p[kh].Dr + ltree.st[i].gn[j].p[k].Dl)*Seffdr - 0.5*(stree.st[ihs].gn[jh].p[kh].Dr + stree.st[istree].gn[j].p[k].Dl)*dSeffdr;
				}

				vector<double> Seffr, dSeffr, Sefful, dSefful;
				//Fluxes on RHS of cell
				for (unsigned m = 0; m < nodes_right; m++) //flux(es) at right boundary
				{
					double sr, sul;
					ih = ltree.st[i].gn[j].p[k].iup[2][m];
					jh = ltree.st[i].gn[j].p[k].jup[2][m];
					kh = ltree.st[i].gn[j].p[k].kup[2][m];
					if (stree.st[istree].gn[j].p[k].iup[2].size() > 1)
					{
						ihs = stree.st[istree].gn[j].p[k].iup[2][m];
						fur.push_back(stree.st[ihs].gn[jh].p[kh].ul*stree.st[istree].gn[j].p[k].sul[m] * ((double)ltree.st[ih].gn[jh].Nb) / ((double)stree.st[ihs].gn[jh].Nb));
						dfur.push_back(ltree.st[ih].gn[jh].p[kh].ul*stree.st[istree].gn[j].p[k].sul[m] * ((double)ltree.st[ih].gn[jh].Nb) / ((double)stree.st[ihs].gn[jh].Nb)
							+ (stree.st[ihs].gn[jh].p[kh].ul*ltree.st[i].gn[j].p[k].sul[m]));
						sr = stree.st[istree].gn[j].p[k].sr[m] * ((double)ltree.st[ih].gn[jh].Nb) / ((double)stree.st[ihs].gn[jh].Nb);
						sul = stree.st[istree].gn[j].p[k].sul[m] * ((double)ltree.st[ih].gn[jh].Nb) / ((double)stree.st[ihs].gn[jh].Nb);
					}
					else
					{
						ihs = stree.st[istree].gn[j].p[k].iup[2][0];
						fur.push_back(stree.st[ihs].gn[jh].p[kh].ul*stree.st[istree].gn[j].p[k].sul[0] * ((double)ltree.st[ih].gn[jh].Nb) / ((double)stree.st[ihs].gn[jh].Nb));
						dfur.push_back(ltree.st[ih].gn[jh].p[kh].ul*stree.st[istree].gn[j].p[k].sul[0] * ((double)ltree.st[ih].gn[jh].Nb) / ((double)stree.st[ihs].gn[jh].Nb)
							+ stree.st[ihs].gn[jh].p[kh].ul*ltree.st[i].gn[j].p[k].sul[m]);
						sr = stree.st[istree].gn[j].p[k].sr[0] * ((double)ltree.st[ih].gn[jh].Nb) / ((double)stree.st[ihs].gn[jh].Nb);
						sul = stree.st[istree].gn[j].p[k].sul[0] * ((double)ltree.st[ih].gn[jh].Nb) / ((double)stree.st[ihs].gn[jh].Nb);
					}
					Snew[2].push_back(stree.st[ihs].gn[jh].p[kh].Anew*((double)ltree.st[ih].gn[jh].Nb)); //cs of right element
					dSnew[2].push_back(ltree.st[ih].gn[jh].p[kh].Anew*((double)ltree.st[ih].gn[jh].Nb)); //cs of right element
					if (srtot > 0)
					{
						Seffr.push_back(sr * (1.0 + o.AcinAreaFactor*(0.5*(Snew[1][0] + stree.st[istree].gn[j].p[k].Aold*((double)ltree.st[i].gn[j].Nb)) - srtot) / srtot));   //work out effective radii increases
						dSeffr.push_back(ltree.st[i].gn[j].p[k].sr[m] * (1.0 + o.AcinAreaFactor*(0.5*(Snew[1][0] + stree.st[istree].gn[j].p[k].Aold*((double)ltree.st[i].gn[j].Nb)) - srtot) / srtot)
							+ 0.5*sr*o.AcinAreaFactor*(dSnew[1][0] + ltree.st[i].gn[j].p[k].Aold*((double)ltree.st[i].gn[j].Nb)) / srtot
							- 0.5*sr*o.AcinAreaFactor*dsrtot*(Snew[1][0] + stree.st[istree].gn[j].p[k].Aold*((double)ltree.st[i].gn[j].Nb)) / (srtot*srtot));
					}
					else
					{
						Seffr.push_back(0.5*o.AcinAreaFactor*(Snew[1][0] + stree.st[istree].gn[j].p[k].Aold*((double)ltree.st[i].gn[j].Nb)) / ((double)ltree.st[i].gn[j].p[k].sr.size())); //duct cs = 0, diffusion area just split evenly between downstream ducts
						dSeffr.push_back(0.5*o.AcinAreaFactor*(dSnew[1][0] + ltree.st[i].gn[j].p[k].Aold*((double)ltree.st[i].gn[j].Nb)) / ((double)ltree.st[i].gn[j].p[k].sr.size())); //duct cs = 0, diffusion area just split evenly between downstream ducts
					}
					Sefful.push_back(sul + o.AcinAreaFactor*(0.5*(Snew[2][m] + stree.st[ihs].gn[jh].p[kh].Aold*((double)ltree.st[ih].gn[jh].Nb)) - sul));
					dSefful.push_back(ltree.st[i].gn[j].p[k].sul[m] + o.AcinAreaFactor*(0.5*(dSnew[2][m] + ltree.st[ih].gn[jh].p[kh].Aold*((double)ltree.st[ih].gn[jh].Nb))
						- ltree.st[i].gn[j].p[k].sul[m]));
					if (Sefful[m] <= Seffr[m])
					{
						fdr.push_back(-0.5*(stree.st[ihs].gn[jh].p[kh].Dl + stree.st[istree].gn[j].p[k].Dr)*Sefful[m]);
						dfdr.push_back(-0.5*(ltree.st[ih].gn[jh].p[kh].Dl + ltree.st[i].gn[j].p[k].Dr)*Sefful[m] - 0.5*(stree.st[ihs].gn[jh].p[kh].Dl + stree.st[istree].gn[j].p[k].Dr)*dSefful[m]);
					}
					else
					{
						fdr.push_back(-0.5*(stree.st[istree].gn[j].p[k].Dr + stree.st[ihs].gn[jh].p[kh].Dl)*Seffr[m]);
						dfdr.push_back(-0.5*(ltree.st[i].gn[j].p[k].Dr + ltree.st[ih].gn[jh].p[kh].Dl)*Seffr[m] - 0.5*(stree.st[istree].gn[j].p[k].Dr + stree.st[ihs].gn[jh].p[kh].Dl)*dSeffr[m]);
					}
					if (k == Nj - 1 && j == stree.st[stree.st[istree].imeanpath].StartGen + stree.st[stree.st[istree].imeanpath].Ntot)   //last point in system
					{
						fur[m] = 0.0;
						dfur[m] = 0.0;
						fdr[m] = 0.0;
						dfdr[m] = 0.0;
					}
				}

				//Fill in matrix
				unsigned long kmh = ltree.st[i].gn[j].p[k].km;
				if (i == 0 && j == 0 && k == 1)
				{
					fdl = 0;   //no diffusion out of mouth -- check this
					dfdl = 0;
				}
				double Ah[4] = { 0, 0, 0, 0 }, bh = 0, xh = 0;
				for (unsigned n = 0; n < 2; n++)   //fill matrix for left and central term
				{
					if (n == 1)
					{
						Ah[n] += 1.0;
						bh += ((double)ltree.st[i].gn[j].Nb)*(ltree.st[i].gn[j].p[k].Aold*cd0old[n][0] + stree.st[istree].gn[j].p[k].Aold*dcd0[n][0] - ltree.st[i].gn[j].p[k].Anew*cd0[n][0]);
						xh += ((double)ltree.st[i].gn[j].Nb)*(ltree.st[i].gn[j].p[k].Aold*cd0old[n][0] + stree.st[istree].gn[j].p[k].Aold*dcd0[n][0] - ltree.st[i].gn[j].p[k].Anew*cd0[n][0]);
						for (unsigned m = 0; m < nodes_right; m++)    //contribution from right flux
						{
							if (Snew[n][0] > 0.) Ah[n] += -0.5*o.dt*((-ucfr[n][m] * fur[m]) + (-dcfr[n][m] * fdr[m])) / (Snew[n][0] * stree.st[istree].gn[j].dx);       //stays the same

							bh += 0.5*o.dt*dcd0[n][0] * ((-ucfr[n][m] * fur[m]) + (-dcfr[n][m] * fdr[m])) / (stree.st[istree].gn[j].dx)
								+ 0.5*o.dt*(cd0[n][0] + cd0old[n][0])*((-ducfr[n][m] * fur[m] - ucfr[n][m] * dfur[m]) + (-ddcfr[n][m] * fdr[m] - dcfr[n][m] * dfdr[m])) / (stree.st[istree].gn[j].dx)
								- 0.5*o.dt*(cd0[n][0] + cd0old[n][0])*ltree.st[i].gn[j].dx*((-ucfr[n][m] * fur[m]) + (-dcfr[n][m] * fdr[m])) / (stree.st[istree].gn[j].dx*stree.st[istree].gn[j].dx);

							xh += o.dt*dcd0[n][0] * ((-ucfr[n][m] * fur[m]) + (-dcfr[n][m] * fdr[m])) / (stree.st[istree].gn[j].dx)
								+ 0.5*o.dt*(cd0[n][0] + cd0old[n][0])* ((-ducfr[n][m] * fur[m] - ucfr[n][m] * dfur[m]) + (-ddcfr[n][m] * fdr[m] - dcfr[n][m] * dfdr[m])) / (stree.st[istree].gn[j].dx)
								- 0.5*o.dt*(cd0[n][0] + cd0old[n][0])*ltree.st[i].gn[j].dx*((-ucfr[n][m] * fur[m]) + (-dcfr[n][m] * fdr[m])) / (stree.st[istree].gn[j].dx*stree.st[istree].gn[j].dx);
						}
					}
					if (Snew[n][0] > 0.) Ah[n] += -0.5*o.dt*((ucfl[n][0] * ful) + (dcfl[n][0] * fdl)) / (Snew[n][0] * stree.st[istree].gn[j].dx);

					bh += 0.5*o.dt*dcd0[n][0] * ((ucfl[n][0] * ful) + (dcfl[n][0] * fdl)) / (stree.st[istree].gn[j].dx)
						+ 0.5*o.dt*(cd0[n][0] + cd0old[n][0])*((ducfl[n][0] * ful + ucfl[n][0] * dful) + (ddcfl[n][0] * fdl + dcfl[n][0] * dfdl)) / (stree.st[istree].gn[j].dx)
						- 0.5*o.dt*(cd0[n][0] + cd0old[n][0])*ltree.st[i].gn[j].dx*((ucfl[n][0] * ful) + (dcfl[n][0] * fdl)) / (stree.st[istree].gn[j].dx*stree.st[istree].gn[j].dx);

					xh += o.dt*dcd0[n][0] * ((ucfl[n][0] * ful) + (dcfl[n][0] * fdl)) / (stree.st[istree].gn[j].dx)
						+ 0.5*o.dt*(cd0[n][0] + cd0old[n][0])*((ducfl[n][0] * ful + ucfl[n][0] * dful) + (ddcfl[n][0] * fdl + dcfl[n][0] * dfdl)) / (stree.st[istree].gn[j].dx)
						- 0.5*o.dt*(cd0[n][0] + cd0old[n][0])*ltree.st[i].gn[j].dx*((ucfl[n][0] * ful) + (dcfl[n][0] * fdl)) / (stree.st[istree].gn[j].dx*stree.st[istree].gn[j].dx);
				}

				{
					unsigned n = 2;
					for (unsigned m = 0; m < nodes_right; m++)    //contribution from fur, fill in for right term
					{
						if (Snew[n][m] > 0.) Ah[n + m] += -0.5*o.dt*((-ucfr[n][m] * fur[m]) + (-dcfr[n][m] * fdr[m])) / (Snew[n][m] * stree.st[istree].gn[j].dx);

						bh += 0.5*o.dt*dcd0[n][m] * ((-ucfr[n][m] * fur[m]) + (-dcfr[n][m] * fdr[m])) / (stree.st[istree].gn[j].dx)
							+ 0.5*o.dt*(cd0[n][m] + cd0old[n][m])*((-ducfr[n][m] * fur[m] - ucfr[n][m] * dfur[m]) + (-ddcfr[n][m] * fdr[m] - dcfr[n][m] * dfdr[m])) / (stree.st[istree].gn[j].dx)
							- 0.5*o.dt*(cd0[n][m] + cd0old[n][m])*ltree.st[i].gn[j].dx*((-ucfr[n][m] * fur[m]) + (-dcfr[n][m] * fdr[m])) / (stree.st[istree].gn[j].dx*stree.st[istree].gn[j].dx);

						xh += o.dt*dcd0[n][m] * ((-ucfr[n][m] * fur[m]) + (-dcfr[n][m] * fdr[m])) / (stree.st[istree].gn[j].dx)
							+ 0.5*o.dt*(cd0[n][m] + cd0old[n][m])*((-ducfr[n][m] * fur[m] - ucfr[n][m] * dfur[m]) + (-ddcfr[n][m] * fdr[m] - dcfr[n][m] * dfdr[m])) / (stree.st[istree].gn[j].dx)
							- 0.5*o.dt*(cd0[n][m] + cd0old[n][m])*ltree.st[i].gn[j].dx*((-ucfr[n][m] * fur[m]) + (-dcfr[n][m] * fdr[m])) / (stree.st[istree].gn[j].dx*stree.st[istree].gn[j].dx);
					}
				}

				//set matrix AC
				unsigned long kmdwn = ltree.st[ltree.st[i].gn[j].p[k].iup[0][0]].gn[ltree.st[i].gn[j].p[k].jup[0][0]].p[ltree.st[i].gn[j].p[k].kup[0][0]].km;
				unsigned long kmup0 = ltree.st[ltree.st[i].gn[j].p[k].iup[2][0]].gn[ltree.st[i].gn[j].p[k].jup[2][0]].p[ltree.st[i].gn[j].p[k].kup[2][0]].km;
				AC.insert(kmh, kmdwn) = Ah[0];
				AC.insert(kmh, kmh) = Ah[1];
				if (kmup0 != kmh)
				{
					AC.insert(kmh, kmup0) = Ah[2];
				}
				else
				{
					AC.coeffRef(kmh, kmh) += Ah[2];
				}
				if (nodes_right > 1)
				{
					unsigned long kmup1 = ltree.st[ltree.st[i].gn[j].p[k].iup[2][1]].gn[ltree.st[i].gn[j].p[k].jup[2][1]].p[ltree.st[i].gn[j].p[k].kup[2][1]].km;
					if (kmup1 != kmh)
					{
						AC.insert(kmh, kmup1) = Ah[3];
					}
					else
					{
						AC.coeffRef(kmh, kmh) += Ah[3];
					}
				}
				BC(kmh) = bh;
				XC(kmh) = xh;
			}
		}
	}

	if (stree.st[0].Vnew - stree.st[0].Vold < 0) x0 = XC(1);  //guess on exhalation
	AC.insert(0, 0) = A0[0];
	AC.insert(0, 1) = A0[1];
	BC(0) = b0;
	XC(0) = x0;
}


void calc_gas_concs_lp(Tree &ltree, Tree &stree, Eigen::VectorXd &XC, Options &o)
{
	//count up linear difference in gas vol on ltree versus stree
	unsigned t = ((unsigned)stree.Vtot.size() - 1);
	ltree.masstot = 0;
	vector<double> dmtotold;
	for (unsigned i = 0; i < ltree.st.size(); i++)
	{
		unsigned istree = ltree.st[i].i_stree;
		dmtotold.push_back(ltree.st[i].masstot);
		ltree.st[i].masstot = 0;
		for (unsigned j = ltree.st[i].StartGen; j <= ltree.st[i].EndGen; j++)
		{
			unsigned Nj = ((unsigned)ltree.st[i].gn[j].p.size());
			for (unsigned k = 0; k < Nj; k++)
			{
				if (stree.st[istree].gn[j].p[k].Anew > 0) ltree.st[i].gn[j].p[k].c = XC(ltree.st[i].gn[j].p[k].km) / (((double)ltree.st[i].gn[j].Nb)*stree.st[istree].gn[j].p[k].Anew);
				else ltree.st[i].gn[j].p[k].c = ltree.st[i].gn[j].p[k].cold;
				if (j>0 || k > 0) ltree.st[i].masstot += ((double)ltree.st[i].gn[j].Nb)*(stree.st[istree].gn[j].p[k].Anew*stree.st[istree].gn[j].p[k].c*ltree.st[i].gn[j].dx + stree.st[istree].gn[j].p[k].Anew*ltree.st[i].gn[j].p[k].c*stree.st[istree].gn[j].dx
					+ ltree.st[i].gn[j].p[k].Anew*stree.st[istree].gn[j].p[k].c*stree.st[istree].gn[j].dx);
			}
		}
		ltree.masstot += ltree.st[i].masstot;
	}

	//calc linear difference in mouth conc on ltree versus stree
	ltree.ctop.push_back(ltree.st[0].gn[0].p[0].c);   //trach value of c at top of tree
	if (stree.st[0].Vnew - stree.st[0].Vold < 0)   //exhalation
	{
		unsigned tc = t - 1;
		double dceff = 0;
		double Vnh = stree.Vtot[t];
		if (o.VDM > 1E-10)
		{
			if (stree.VDfront > 1E-10)   //front still in dead-space
			{
				dceff = 0;
			}
			else
			{
				while (tc > 0 && (Vnh - (stree.Vtot[tc - 1] - o.VDM))* (Vnh - (stree.Vtot[tc] - o.VDM)) > 0) //when this > 0, both these terms have same sign, hence not the right point
				{
					tc--;
				}
				if (tc <= 1) dceff = 0;
				else
				{
					dceff = (ltree.ctop[tc - 1] * fabs(Vnh - (stree.Vtot[tc] - o.VDM)) + ltree.ctop[tc] * fabs(Vnh - (stree.Vtot[tc - 1] - o.VDM))) / (fabs(Vnh - (stree.Vtot[tc] - o.VDM)) + fabs(Vnh - (stree.Vtot[tc - 1] - o.VDM)));
				}
			}
			ltree.cmouth = dceff; //cmouth at n+1
		}
		else ltree.cmouth = ltree.st[0].gn[0].p[0].c;
	}

	//calc linear difference in flux into subtrees on ltree versus stree
	for (unsigned i = 0; i < ltree.st.size(); i++)
	{
		unsigned istree = ltree.st[i].i_stree;
		unsigned j = ltree.st[i].StartGen;
		unsigned ih, jh, kh, ihs;
		unsigned k0;
		if (j == 0) k0 = 1;
		else k0 = 0;
		vector<double> ucfr[3], ucfl[3], dcfr[3], dcfl[3], ducfr[3], ducfl[3], ddcfr[3], ddcfl[3];
		ltree.st[i].fluxin = 0;
		for (unsigned n = 0; n < 3; n++)
		{
			ih = ltree.st[i].gn[j].p[k0].iup[n][0];
			ihs = stree.st[istree].gn[j].p[k0].iup[n][0];
			jh = ltree.st[i].gn[j].p[k0].jup[n][0];
			kh = ltree.st[i].gn[j].p[k0].kup[n][0];
			dcfr[n] = stree.st[istree].gn[j].p[k0].dcfr[n];
			dcfl[n] = stree.st[istree].gn[j].p[k0].dcfl[n];
			ddcfr[n] = ltree.st[i].gn[j].p[k0].dcfr[n];
			ddcfl[n] = ltree.st[i].gn[j].p[k0].dcfl[n];
			if (stree.st[istree].gn[j].p[k0].ul > 0)
			{
				ucfr[n] = stree.st[istree].gn[j].p[k0].ucfrpos[n];
				ucfl[n] = stree.st[istree].gn[j].p[k0].ucflpos[n];
				ducfr[n] = ltree.st[i].gn[j].p[k0].ucfrpos[n];
				ducfl[n] = ltree.st[i].gn[j].p[k0].ucflpos[n];
			}
			else
			{
				ucfr[n] = stree.st[istree].gn[j].p[k0].ucfrneg[n];
				ucfl[n] = stree.st[istree].gn[j].p[k0].ucflneg[n];
				ducfr[n] = ltree.st[i].gn[j].p[k0].ucfrneg[n];
				ducfl[n] = ltree.st[i].gn[j].p[k0].ucflneg[n];
			}
			ltree.st[i].fluxin += 0.5*o.dt*((double)ltree.st[i].gn[j].Nb)*
				(ltree.st[i].gn[j].p[k0].al*ucfl[n][0] * stree.st[istree].gn[j].p[k0].ul*(stree.st[ihs].gn[jh].p[kh].cold + stree.st[ihs].gn[jh].p[kh].c) +
				stree.st[istree].gn[j].p[k0].al*ducfl[n][0] * stree.st[istree].gn[j].p[k0].ul*(stree.st[ihs].gn[jh].p[kh].cold + stree.st[ihs].gn[jh].p[kh].c) +
				stree.st[istree].gn[j].p[k0].al*ucfl[n][0] * ltree.st[i].gn[j].p[k0].ul*(stree.st[ihs].gn[jh].p[kh].cold + stree.st[ihs].gn[jh].p[kh].c) +
				stree.st[istree].gn[j].p[k0].al*ucfl[n][0] * stree.st[istree].gn[j].p[k0].ul*(ltree.st[ih].gn[jh].p[kh].cold + ltree.st[ih].gn[jh].p[kh].c));
		}
		//        double Check1 = stree.st[i].fluxin;
		//        double Check2 = (stree.st[i].masstot - mtotold[i]);
		//        if(fabs(Check1 - Check2) > 1E-11)
		//        {
		//            cout << "Here\n";
		//        }
	}
	ltree.fluxin = ltree.st[0].fluxin;

	//count up difference in gas vol in each acinus of ltree versus stree
	for (unsigned m = 0; m < ltree.EndSubtrees.size(); m++)  //loop trees that reach acinus
	{
		unsigned i = ltree.EndSubtrees[m];   //current subtree
		ltree.st[i].totIGvol = 0;
		for (unsigned n = 0; n < ltree.st[i].isub.size(); n++)  //loop over subtrees
		{
			if (i != ltree.st[i].isub[n]) ltree.st[ltree.st[i].isub[n]].totIGvol = 0;
			unsigned isst = ltree.st[ltree.st[i].isub[n]].i_stree;
			for (unsigned j = ltree.st[ltree.st[i].isub[n]].StartGen + ltree.st[ltree.st[i].isub[n]].Ncond; j <= ltree.st[ltree.st[i].isub[n]].EndGen; j++)
			{
				for (unsigned k = 0; k < ltree.st[ltree.st[i].isub[n]].gn[j].p.size(); k++)
				{
					ltree.st[ltree.st[i].isub[n]].totIGvol += ((double)ltree.st[ltree.st[i].isub[n]].gn[j].Nb)
						*(ltree.st[ltree.st[i].isub[n]].gn[j].p[k].c*stree.st[isst].gn[j].p[k].Anew*stree.st[isst].gn[j].dx
						+ stree.st[isst].gn[j].p[k].c*ltree.st[ltree.st[i].isub[n]].gn[j].p[k].Anew*stree.st[isst].gn[j].dx
						+ stree.st[isst].gn[j].p[k].c*stree.st[isst].gn[j].p[k].Anew*ltree.st[ltree.st[i].isub[n]].gn[j].dx);
				}
			}
			if (i != ltree.st[i].isub[n]) ltree.st[i].totIGvol += ltree.st[ltree.st[i].isub[n]].totIGvol;
		}
	}
}
