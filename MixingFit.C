#include "MixingFit.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cassert> 
#include "TMinuit.h"
#include "TRandom.h" 
#include <algorithm> 
#include <iomanip> 
//The commenting that was missed
TMatrixTSym<matrixfp> MixingResult::sigma; //this is the covariance matrix for inversion for the chi2
TMatrixTSym<matrixfp> MixingResult::invsigma; //the inverse covariance matrix
std::vector<double> MixingResult::epsilon; //the error for a single measurement for epsilon_i W_ij epsilon_j
std::vector<MixingResult*> MixingResult::allResults; //All the results we load
std::map<std::string, MixingResult*> MixingResult::mapResults;//map the results from a string to a mixing result
bool MixingResult::initialised = false; //Class helper. See if the result is initialized
const int MixingResult::nParams = MixingResult::NUMSENSE; //number of parameters
TMinuit* MixingResult::minuit = new TMinuit(MixingResult::nParams);//initialize minuit
MixingResult* MixingResult::fitResult = new MixingResult("fitresult", "plainx", "none", 0, 0, 0, 0, 1);  //plain initializer for the mixing fit
bool* MixingResult::isSensitive = new bool[MixingResult::nParams]; //important. return sensitivity for an individual mixing result
bool debugCorrelations = false; //do you want to debug correlations?
bool printChisq = false; //do you want to print the chi2 at every step (dump this into a file, otherwise you'll be sorry)
bool use_hfag_convention = true; //hfag convention for the sign of the primes?
bool special_invert = false; //do rolf's inversion?
MixingResult::CpvAllowed allowcpv = MixingResult::NOCPV; //the type of cpv allowed i.e. direct, indirect, etc
MixingResult::ConstraintType fit_for_which = MixingResult::BOTH_FREE; // Fit for phi, |q/p|, or both. 
bool special_alex_fit = false; // If true, we fit for x_12, y_12, phi_12. 
bool skipGraphics = false; 
char strbuffer[1000]; 

double fiveSigma  = 28.74;
double fourSigma  = 19.33;
double threeSigma = 11.83;
double twoSigma   =  6.18;
double oneSigma   =  2.30;


MixingResult::MixingResult (const char* n, 
			    const char* rtype, 
			    const char* ptype, 
			    double res,
			    double stat,
			    double syst,
			    double modl,
			    double coef) 
  : name(n)//fill the shiz for the class
  , measurement(res)
  , error(0)
  , correlations()
  , index(-1)
  , myType(getResultType(rtype))
  , myPrime(getPrimeType(ptype))
  , active(false) 
  , coefficient(coef)
{

  if (name == "fitresult") active = false;
  allResults.push_back(this); 
  mapResults[name] = this; 
  error = sqrt(pow(stat, 2) + pow(syst, 2) + pow(modl, 2)); 

  
  std::cout << name << " : "
	    << measurement << " +/- "
	    << error << " ("
	    << stat << ", "
	    << syst << ", "
	    << modl << ") "
            << coef << " "
	    << rtype << " "
	    << myType << " "
	    << ptype << " " 
	    << myPrime << " " 
	    << std::endl; 

}

MixingResult::~MixingResult () {//destructor
  mapResults[name] = 0; 
  std::vector<MixingResult*>::iterator me = std::find(allResults.begin(), allResults.end(), this); 
  if (me != allResults.end()) allResults.erase(me); 
}

void MixingResult::randomise () {//randomizer
  static TRandom donram(42); 
  std::cout << "Randomising " << getName() << "; old value "
	    << measurement;
  double delta = donram.Gaus()*error;
  measurement += delta;
  std::cout << " delta " << delta << " new value " << measurement << std::endl; 
}

double getMixX (double x12, double y12, double phi12) {//
  x12 *= x12;
  y12 *= y12; 

  double ret = sqrt(x12 - y12 + sqrt(pow(x12+y12, 2) - x12*y12*pow(2*sin(phi12), 2)));
  if (cos(phi12) < 0) ret *= -1; 
  ret *= sqrt(0.5); 

  return ret;
}

double getMixY (double x12, double y12, double phi12) {
  x12 *= x12;
  y12 *= y12; 
  double ret = sqrt(y12 - x12 + sqrt(pow(x12+y12, 2) - x12*y12*pow(2*sin(phi12), 2)));
  assert(ret == ret);
  ret *= sqrt(0.5); 
  
  return ret;
}

double getQoverP (double x12, double y12, double phi12) {
  double ret =  x12*x12 + y12*y12 + 2*x12*y12*sin(phi12);
  ret       /= (x12*x12 + y12*y12 - 2*x12*y12*sin(phi12));
  ret = pow(ret, 0.25);
  assert(ret == ret);
  return ret; 
}


double getPhi (double x12, double y12, double phi12) {
  if (fabs(x12) <= fabs(y12)*1e-10) return 0; 

  double ret = -sin(2*phi12);
  ret /= (cos(2*phi12) + (y12*y12/(x12*x12)));
  ret = atan(ret); 
  assert(ret == ret);
  ret *= 0.5; 
  return ret;  
}

int sign (double x) {
  if (x < 0) return -1;
  return 1; 
}

double MixingResult::calcEpsilon (double x, 
				  double y, 
				  double deltaKpi, 
				  double deltaKpipi, 
				  double phi, 
				  double rsubdm, 
				  double rsubdp,
				  double qoverp) {

  if (special_alex_fit) {
    // In this case, interpret x, y, phi as underlying x_12, y_12, phi_12,
    // and calculate x, y, phi, q/p from them.
    double mixx = getMixX(x, y, phi);
    double mixy = getMixY(x, y, phi); 
    double qovp = getQoverP(x, y, phi);
    double rphi = getPhi(x, y, phi);
    x = mixx;
    y = mixy;
    qoverp = qovp;
    phi = rphi; 
  }
  else {
    switch (fit_for_which) {
    case MixingResult::PHI_FREE:
      //qoverp = 1 - (y/x)*tan(phi);
      qoverp = sqrt((x-y*tan(phi))/(x+y*tan(phi))); 
      // This is eqn 20 of Grossman/Nir/Perez. It assumes small CPV, ie |sin(phi_12)|<<1.
      // NB! This function's 'phi' is not the same as 'phi_12'! 
      break;
    case QP_FREE:
      if (fabs(x) > fabs(y)*1e12) phi = sign(x)*sign(y)*M_PI_2; 
      // ie, y approaches zero, tan phi approaches infinity, phi is pi/2.
      else phi = atan(((1-qoverp*qoverp)/(1+qoverp*qoverp))*(x/y)); 
      break;
    default:
      break;
    }
  }



  static int primeSign = use_hfag_convention ? 1 : -1; 
  double ret = 0; 
  double prime = 0; 
  double poverq = 1.0 / qoverp;
  double a_m = (qoverp*qoverp - poverq*poverq);
  a_m /= (qoverp*qoverp + poverq*poverq);
  switch (myType) {
  case PLAINX: ret = (x - measurement); break;
  case XSQUARE: ret = (x*x - measurement); break;
  case PLAINY: 
  case YCP:    
    ret  = (qoverp + poverq)*y*cos(phi);
    ret -= (qoverp - poverq)*x*sin(phi);
    ret *= 0.5;
    ret -= measurement; 
    break;
  case AGAMMA:
    ret  = (qoverp - poverq)*y*cos(phi);
    ret -= (qoverp + poverq)*x*sin(phi);
    ret *= 0.5;
    ret -= measurement; 
    break;
  case XPRIME: 
    switch (myPrime) {
    case KPI:    prime = pow(x*cos(deltaKpi)   + primeSign*y*sin(deltaKpi), 2); break; // NB squaring - Kpi result is assumed x'^2. 
    case KPIPI0: prime =     x*cos(deltaKpipi) + primeSign*y*sin(deltaKpipi); break;
    default:
    case NOPRIME:
      std::cout << "Error: Measurement " << name << " is not a primed quantity." << std::endl; 
      assert(false); 
      break;
    }
    ret = prime - measurement; 
    break;
  case XPRIME_P:
    switch (myPrime) {
    case KPI:    
      prime  = (x*cos(deltaKpi)   + primeSign*y*sin(deltaKpi)) * cos(phi); 
      prime += (y*cos(deltaKpi)   - primeSign*x*sin(deltaKpi)) * sin(phi); 
      prime *= pow((1+a_m)/(1-a_m), 0.25);
      prime *= prime;
      break;
    case KPIPI0: 
    default:
    case NOPRIME:
      std::cout << "Error: Measurement " << name << " is not a primed quantity." << std::endl; 
      assert(false); 
      break;
    }
    ret = prime - measurement; 
    break;
  case XPRIME_M:
    switch (myPrime) {
    case KPI:    
      prime  = (x*cos(deltaKpi)   + primeSign*y*sin(deltaKpi)) * cos(phi); 
      prime -= (y*cos(deltaKpi)   - primeSign*x*sin(deltaKpi)) * sin(phi); 
      prime *= pow((1-a_m)/(1+a_m), 0.25);
      prime *= prime;
      break;
    case KPIPI0: 
    default:
    case NOPRIME:
      std::cout << "Error: Measurement " << name << " is not a primed quantity." << std::endl; 
      assert(false); 
      break;
    }
    ret = prime - measurement; 
    break;


  case YPRIME:
    switch (myPrime) { 
    case KPI:    prime = y*cos(deltaKpi)   - primeSign*x*sin(deltaKpi); break;
    case KPIPI0: prime = y*cos(deltaKpipi) - primeSign*x*sin(deltaKpipi); break;
    default:
    case NOPRIME:
      std::cout << "Error: Measurement " << name << " is not a primed quantity." << std::endl; 
      assert(false); 
      break;
    }
    ret = prime - measurement; 
    break;
  case YPRIME_P:
    switch (myPrime) {
    case KPI:    
      prime  = (y*cos(deltaKpi)   - primeSign*x*sin(deltaKpi)) * cos(phi); 
      prime -= (x*cos(deltaKpi)   + primeSign*y*sin(deltaKpi)) * sin(phi); 
      prime *= pow((1+a_m)/(1-a_m), 0.25);
      break;
    case KPIPI0: 
    default:
    case NOPRIME:
      std::cout << "Error: Measurement " << name << " is not a primed quantity." << std::endl; 
      assert(false); 
      break;
    }
    ret = prime - measurement; 
    break;
  case YPRIME_M:
    switch (myPrime) {
    case KPI:    
      prime  = (y*cos(deltaKpi)   - primeSign*x*sin(deltaKpi)) * cos(phi); 
      prime += (x*cos(deltaKpi)   + primeSign*y*sin(deltaKpi)) * sin(phi); 
      prime *= pow((1-a_m)/(1+a_m), 0.25);
      break;
    case KPIPI0: 
    default:
    case NOPRIME:
      std::cout << "Error: Measurement " << name << " is not a primed quantity." << std::endl; 
      assert(false); 
      break;
    }
    ret = prime - measurement; 
    break;
  case COSANGLE:
    switch (myPrime) {
    case KPI:    ret = cos(deltaKpi) - measurement; break;
    case KPIPI0: ret = cos(deltaKpipi) - measurement; break;
    default:
      std::cout << "Error: Measurement " << name << " is not a primed quantity." << std::endl;
      assert(false);
      break;
    }
    break;
  case SINANGLE:
    switch (myPrime) {
    case KPI:    ret = sin(deltaKpi) - measurement; break;
    case KPIPI0: ret = sin(deltaKpipi) - measurement; break;
    default:
      std::cout << "Error: Measurement " << name << " is not a primed quantity." << std::endl;
      assert(false);
      break;
    }
    break;
    
  case ANGLE:
    switch (myPrime) {
    case KPI:    ret = deltaKpi - measurement; break;
    case KPIPI0: ret = deltaKpipi - measurement; break; 
    default:
    case NOPRIME:
      std::cout << "Error: Measurement " << name << " is not a primed quantity." << std::endl;
      assert(false);
      break;
    }
    break;

  case NODCPV:
    switch (myPrime) {
      // Eqns 14, 15, 48, 52 in Kagan/Sokoloff Phys Rev D80, 0.6008 (2009)
      // Mixing x in terms of underlying x12, y12, phi12
      // NB 14 and 15 corrected by factor 1/sqrt(2) - email from Alex
      // Note that *formal parameters* x, y are now *variables* x12, y12 - confusing! 
    case KPI: // Mixing x
      ret = getMixX(x, y, phi) - measurement;
      break; 
    case KPIPI0: // y
      ret = getMixY(x, y, phi) - measurement;
      break;
    case QOVERP: // q over p
      ret = getQoverP(x, y, phi) - measurement; 
      break;
    case NOPRIME: // phi
      ret = getPhi(x, y, phi) - measurement; 
      break;
    default:
      std::cout << "Result " << getName() << " has bad prime type; check the input.\n";
      assert(false);
      break;
    }
    break; 

  case RD:  ret = (ALL_CPV == allowcpv ? 0.5*(rsubdm+rsubdp): rsubdm) - measurement; break;
  case RDM: ret = rsubdm - measurement; break;
  case RDP: ret = rsubdp - measurement; break;
  case QP:  ret = qoverp - measurement; break;
  case PHI: ret = phi  - measurement; break;

  case NUMTYPES: 
  default: 
    std::cout << "Error: Unknown result type." << std::endl; 
    assert(false); 
    break;
  }

  if (active) {
    assert(index < (int) epsilon.size()); 
    epsilon[index] = ret;
  }
  return ret; 
}

double MixingResult::correlation (MixingResult* dat) {
  if (dat == this) return 1; 
  return correlations[dat]; 
}

void MixingResult::setCorrelation (MixingResult* dat, double corr) {
  assert(dat); 
  assert(fabs(corr) <= 1); 
  initialised = false; 
  correlations[dat] = corr; 
}


void MixingResult::initialise () {
  if (initialised) return;
  initialised = true; 
  //std::cout << "Entering initialisation method" << std::endl; 

  int numActive = 0;
  std::vector<MixingResult*> actives; 
  for (ResultIterator i = begin(); i != end(); ++i) {
    if (!(*i)->active) continue;
    (*i)->index = numActive++; 
    actives.push_back(*i); 
    //std::cout << "Active " << numActive << " is " << (*i)->name << std::endl;
  }

  sigma.ResizeTo(numActive, numActive);
  invsigma.ResizeTo(numActive, numActive);
  epsilon.clear();
  epsilon.resize(numActive); 
 
  static const matrixfp mintol = 2.35e-16;

  for (int ii = 0; ii < numActive; ++ii) {
    //std::cout << actives[ii]->name << " " << actives[ii]->error << std::endl; 
    for (int jj = 0; jj < numActive; ++jj) {
      sigma(ii, jj) = actives[ii]->error * actives[jj]->error * actives[ii]->correlation(actives[jj]); 
      invsigma(ii,jj) = sigma(ii, jj);
    } 

    // Avoid inversion tolerance - kludgy.
    if (sigma(ii, ii) < mintol) {
      std::cout << "Low value " << ii << " " << (double) sigma(ii, ii) << std::endl;
      sigma(ii, ii) = mintol;
      invsigma(ii, ii) =  mintol;
    }
  } 

  invsigma.Invert(0);

  if (special_invert) {
    for (int j = 0; j < numActive; ++j) {
      for (int i = 0; i < numActive; ++i) {
	invsigma(i, j) = 0;
      }
    }
    invsigma(3, 3) = 1.0 / sigma(3, 3);
    invsigma(6, 6) = 1.0 / sigma(6, 6);
    invsigma(7, 7) = 1.0 / sigma(7, 7);
    double det = sigma(4, 4)*sigma(5, 5) - sigma(4, 5)*sigma(5, 4);
    det = 1.0 / det;
    invsigma(4, 4) = sigma(5, 5)*det;
    invsigma(5, 5) = sigma(4, 4)*det;
    invsigma(4, 5) = -sigma(4, 5)*det;
    invsigma(5, 4) = -sigma(5, 4)*det;
    
    TMatrixTSym<double> threeblock;
    threeblock.ResizeTo(3, 3); 
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
	threeblock(i, j) = sigma(i, j);
      }
    }
    threeblock.Invert(0);
    
    for (int j = 0; j < 3; ++j) {
      for (int i = 0; i < 3; ++i) {
	invsigma(i, j) = threeblock(i, j); 
      }
    }
  }
  /*
  std::cout << "Inversions:\n"; 
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      std::cout << std::setprecision(20) << threeblock(i, j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  for (int j = 0; j < numActive; ++j) {
    for (int i = 0; i < numActive; ++i) {
      std::cout << std::setprecision(20) << invsigma(i, j) << " ";
    }
    std::cout << std::endl;
  }
  */

  /*
  std::cout << " sigma matrix = " << std::endl;
  sigma.Print();
  std::cout << " invsigma matrix = " << std::endl;
  invsigma.Print();
  TMatrix prod = invsigma*sigma;
  std::cout << " prod matrix = " << std::endl;
  prod.Print();
  */

}

bool MixingResult::isSensitiveTo (int param) const {
  assert(param >= 0);
  assert(param < nParams);

  switch (param) {
  case EKS: 
    if (YCP      == myType) return (ALL_CPV == allowcpv); 
    if (AGAMMA   == myType) return ((ALL_CPV == allowcpv) || ((INDIRECT_CPV == allowcpv) && (MixingResult::PHI_FREE == fit_for_which))); 
    if (PLAINY   == myType) return false; 
    if (RDM      == myType) return false;
    if (RDP      == myType) return false;
    if (COSANGLE == myType) return false;
    if (SINANGLE == myType) return false;
    if (ANGLE    == myType) return false;
    return true; 
  case WYE: 
    if (AGAMMA   == myType) return ((ALL_CPV == allowcpv) || ((INDIRECT_CPV == allowcpv) && (QP_FREE == fit_for_which))); 
    if (PLAINX   == myType) return false;
    if (XSQUARE  == myType) return false;
    if (RDM      == myType) return false; 
    if (RDP      == myType) return false; 
    if (COSANGLE == myType) return false;
    if (SINANGLE == myType) return false;
    if (ANGLE    == myType) return false;
    return true;
  case DELTAKPI: 
    if (NODCPV != myType) return (KPI == myPrime);
    return false;
  case DELTAKPIPI: 
    if (NODCPV != myType) return (KPIPI0 == myPrime); 
    return false;
  case PHI_12:
    if ((INDIRECT_CPV == allowcpv) && (QP_FREE == fit_for_which)) return false; 
    if (NOCPV == allowcpv) return false; 
    if (AGAMMA   == myType) return true;
    if (YCP      == myType) return true;
    if (NODCPV   == myType) return true; 
    if (YPRIME_P == myType) return true; 
    if (YPRIME_M == myType) return true; 
    if (XPRIME_P == myType) return true; 
    if (XPRIME_M == myType) return true; 
    return false;
  case RSUBDP:
    if ((RD == myType) && (ALL_CPV == allowcpv)) return true; 
    return (RDP == myType);
  case RSUBDM:
    if (RD == myType) return true; 
    return (RDM == myType);    
  case MAGQOVERP:
    if (NOCPV == allowcpv) return false; 
    if ((INDIRECT_CPV == allowcpv) && (MixingResult::PHI_FREE == fit_for_which)) return false; 
    if (special_alex_fit)   return false; 
    if (AGAMMA   == myType) return true;
    if (YCP      == myType) return true;
    if (YPRIME_P == myType) return true;
    if (YPRIME_M == myType) return true;
    if (XPRIME_P == myType) return true;
    if (XPRIME_M == myType) return true;
    if (NODCPV == myType)   return true; 
    return (QP == myType); 
  default: assert(false); 
  }
  return false; 
}

MixingResult::ResultType MixingResult::getResultType (std::string n) {
  if (n == "plainx")   return PLAINX;
  if (n == "plainy")   return PLAINY;
  if (n == "xsquare")   return XSQUARE; 
  if (n == "xprime")   return XPRIME;
  if (n == "xprimem")  return XPRIME_M;
  if (n == "xprimep")  return XPRIME_P;
  if (n == "yprime")   return YPRIME;
  if (n == "yprimem")  return YPRIME_M;
  if (n == "yprimep")  return YPRIME_P;
  if (n == "angle")    return ANGLE;
  if (n == "cosangle") return COSANGLE;
  if (n == "sinangle") return SINANGLE;
  if (n == "ycp")      return YCP;
  if (n == "nodcpv")   return NODCPV; 
  if (n == "rsubd")    return RD;
  if (n == "rsubdm")   return RDM;
  if (n == "rsubdp")   return RDP;
  if (n == "qoverp")   return QP;
  if (n == "phi")      return PHI;
  if (n == "agamma")   return AGAMMA; 
  std::cout << "Error: Unknown result type \"" << n << "\"\n"; 
  assert(false); 
  return PLAINX; 
}

MixingResult::PrimeType MixingResult::getPrimeType (std::string n) {
  if ((n == "kpi") || (n == "xmix"))   return KPI;
  if ((n == "kpipi") || (n == "ymix")) return KPIPI0;
  if  (n == "qoverp")                  return QOVERP; 
  if ((n == "none") || (n == "phi"))   return NOPRIME; 
  std::cout << "Warning: Unknown prime type" << n << std::endl;
  return NOPRIME; 
}

void MixChisqFcn (int& npar, double* deriv, double& func, double param[], int flag) {
  if (!MixingResult::initialised) {
    MixingResult::initialise(); 
  }

  for (std::vector<MixingResult*>::iterator i = MixingResult::allResults.begin(); i != MixingResult::allResults.end(); ++i) {
    if (!(*i)->active) continue;
    (*i)->calcEpsilon(param[0], param[1], param[2], param[3], param[4], param[5], param[6], param[7]); 
  }

  if (debugCorrelations) {
    std::cout << "Call with x="  << param[0] 
	      << ", y="          << param[1] 
	      << ", deltaKpi="   << param[2] 
	      << ", deltaKpipi=" << param[3] 
	      << ", phi="        << param[4] 
	      << ", R_D_m="      << param[5]
	      << ", R_D_p="      << param[6] 
	      << ", |q/p|="      << param[7] 
	      << std::endl;
  }

  double chisq = 0;
  double extra = 0; 
  for (std::vector<MixingResult*>::iterator i = MixingResult::allResults.begin(); i != MixingResult::allResults.end(); ++i) {
    if (!(*i)->active) continue;
    int idx1 = (*i)->index; 
    double currChisq = 0; 
    for (std::vector<MixingResult*>::iterator j = MixingResult::allResults.begin(); j != MixingResult::allResults.end(); ++j) {
      if (!(*j)->active) continue;
      int idx2 = (*j)->index; 
      double term = MixingResult::epsilon[idx1] * MixingResult::invsigma(idx1, idx2) * MixingResult::epsilon[idx2];
      extra += term; 
      currChisq += term; 
      chisq += MixingResult::epsilon[idx1] * MixingResult::invsigma(idx1, idx2) * MixingResult::epsilon[idx2];

      if ((debugCorrelations) && (fabs(term) > 0.0001)) {
	std::cout << (*i)->name << ", " << (*j)->name << " contribution: " 
		  << term << " from "
		  << MixingResult::epsilon[idx1] << " " 
		  << (double) MixingResult::invsigma(idx1, idx2) << " "
		  << MixingResult::epsilon[idx2] << " "
	  //<< extra << " " << chisq << " " << (extra-chisq) 
		  << std::endl; 
      }
    }
    if (printChisq) {
      std::cout << (*i)->name << ":\t " << currChisq << " " << extra << std::endl; 
    }
  }
  //func = chisq;
  func = extra; 

  if (debugCorrelations) {
    std::cout << "Total chisquare: " << func << std::endl << std::endl; 
  }

  //sprintf(strbuffer, "%.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f", param[0], param[1], param[2], param[3], param[4], param[5], param[6], func);
    //sprintf(strbuffer, "%.12f %.12f %.12f %.12f %.12f", param[0], param[1], param[2], param[5], func);



    //std::cout << strbuffer << std::endl;
}


int main (int argc, char** argv) {
  if (argc < 3) {
    std::cout << "Usage: mixingfit fitfile graphicsfile" << std::endl;
    return 1; 
  }
  
  std::ifstream fitOpts(argv[1]); 
  if (fitOpts.fail()) {
    std::cout << "Failed to open fitting file " << argv[1] << std::endl;
    return 1; 
  }

  std::ifstream grpOpts(argv[2]); 
  if (grpOpts.fail()) {
    std::cout << "Failed to open graphics file " << argv[2] << std::endl;
    return 1; 
  }

  int numActive = 0; 
  char token[10000]; 
  MixingResult* res = 0;
  std::ofstream my_out_file;
  std::string my_string = std::string(argv[1])+"_out.tex";
  //  std::cout<<"my_string = "<<my_string<<std::endl;
  my_out_file.open(my_string.c_str());
  my_out_file<<"\\begin{table}[htdp]"<<std::endl
	     <<"\\caption{All CPV allowed inputs}"<<std::endl
	     <<"\\begin{center}"<<std::endl
	     <<"\\begin{tabular}{|c|c|c|}"<<std::endl
	     <<"\\hline"<<std::endl
	     <<"Result & Value & Correlation Coefficients \\\\"<<std::endl
	     <<"\\hline \\hline"<<std::endl;


  bool which_override = false; 
  while (!fitOpts.eof()) {
    fitOpts >> token; 
    if (token[0] == '#') {
      fitOpts.getline(token, 10000); 
      continue; 
    }
    
    std::string lineType(token); 
    if (lineType == "result") {
      char name[1000];
      char rtype[1000];
      char ptype[1000];
      double result = 0;
      double stat = 0;
      double syst = 0;
      double modl = 0;
      double coef = 0; 
      fitOpts >> name >> rtype >> ptype >> result >> stat >> syst >> modl >> coef; 
      res = new MixingResult(name, rtype, ptype, result, stat, syst, modl, coef); 
      my_out_file<<"$" <<name <<"("<<ptype<<")"<<"$"<<" & $"<<result<<"\\pm"<<stat<<"\\pm"<<syst<<"$ & \\\\"<<std::endl;
    }
    else if (lineType == "correlation") {
      char numOne[1000];
      char numTwo[1000];
      double corr = 0;
      fitOpts >> numOne >> numTwo >> corr;
      MixingResult* one = MixingResult::getByName(numOne);
      MixingResult* two = MixingResult::getByName(numTwo);
      assert(one);
      assert(two); 
      one->setCorrelation(two, corr); 
      two->setCorrelation(one, corr); 
    }   
    else if (lineType == "fitlist") {
      std::string name;
      while (true) {
	fitOpts >> name;
	if (name == "endline") break;
	MixingResult* one = MixingResult::getByName(name);
	one->setActive(true); 
	std::cout << "Activating result: " << name << std::endl; 
	numActive++; 
      }
    }
    else if (lineType == "silent") {
      MixingResult::minuit->SetPrintLevel(-1); 
    }
    else if (lineType == "allow_indirect_cpv") {
      allowcpv = MixingResult::INDIRECT_CPV; 
      if (!which_override) fit_for_which = MixingResult::QP_FREE; 
    }
    else if (lineType == "allow_all_cpv") allowcpv = MixingResult::ALL_CPV; 
    else if (lineType == "use_hfag_convention") use_hfag_convention = true;
    else if (lineType == "fit_for_qp") {fit_for_which = MixingResult::QP_FREE; which_override = true;}
    else if (lineType == "fit_for_phi") {fit_for_which = MixingResult::PHI_FREE; which_override = true;} 
    else if (lineType == "special_alex_fit") special_alex_fit = true; 
    else if (lineType == "skipGraphics") skipGraphics = true; 
    else if (lineType == "use_block_diag") {
      std::string name;
      fitOpts >> name;
      special_invert = (name == "yes"); 
    }
    else if (lineType == "randomise") {
      std::string name;
      fitOpts >> name;
      MixingResult* one = MixingResult::getByName(name);
      assert(one); 
      one->randomise(); 
    }
    else if (lineType == "endfile") break; 
  }
  my_out_file<<"\\hline"<<std::endl
	     <<"\\end{tabular}"<<std::endl
	     <<"\\end{center}"<<std::endl
	     <<"\\label{table:"<<std::string(argv[1]).c_str()<<"}"<<std::endl
	     <<"\\end{table}"<<std::endl;

  my_out_file.close();
  if (0 == numActive) {
    std::cout << "Warning: No measurements specified in fitlist! No fit will be done." << std::endl; 
  }
  else {
    std::cout << "Following results will be used in fit: \n";
    for (MixingResult::ResultIterator res = MixingResult::begin(); res != MixingResult::end(); ++res) {
      if (!(*res)->isActive()) continue;
      std::cout << "  " << (*res)->getName() << std::endl; 
    }
  }

  if (3 < argc) {
    if (std::string(argv[3]) == "qp")  fit_for_which = MixingResult::QP_FREE;
    if (std::string(argv[3]) == "phi") fit_for_which = MixingResult::PHI_FREE;
  }

  double par[MixingResult::nParams];               
  double stepSize[MixingResult::nParams];          
  double minVal[MixingResult::nParams];            
  double maxVal[MixingResult::nParams];            
  double outpar[MixingResult::nParams];            
  double minerr[MixingResult::nParams];            
  std::string parName[MixingResult::nParams];
  par[0] = 0.0012;            
  stepSize[0] = 0.0001;       
  minVal[0] = -0.1;         
  maxVal[0] = 0.15;
  parName[0] = special_alex_fit ? "x_12" : "x";
  par[1] = 0.0070;          
  stepSize[1] = 0.0001;     
  minVal[1] = -0.1;         
  maxVal[1] = 0.15;
  parName[1] = special_alex_fit ? "y_12" : "y";
  par[2]      = 0.3;       
  stepSize[2] = 0.01;      
  minVal[2]   = -1.6;      
  maxVal[2]   = 2.6;
  parName[2]  = "deltaKpi";
  par[3]      = 0.39; 
  stepSize[3] = 0.01;
  minVal[3]   = -3;
  maxVal[3]   = 2;
  parName[3]  = "deltaKpipi0";
  par[4]      = 0; 
  stepSize[4] = 0.001;
  minVal[4]   = -3.15;
  maxVal[4]   = 3.15;
  parName[4]  = special_alex_fit ? "phi_12" : "phi"; 
  par[5]      = 0.0034; 
  stepSize[5] = 0.0001;
  minVal[5]   = -0.01;
  maxVal[5]   = 0.01;
  parName[5]  = "rsubdm";
  par[6]      = 0.0034; 
  stepSize[6] = 0.0001;
  minVal[6]   = -0.01;
  maxVal[6]   = 0.01; 
  parName[6]  = "rsubdp"; 

  par[7]      = 1.00; 
  stepSize[7] = 0.01;
  minVal[7]   = 0.25;
  maxVal[7]   = 1.9; 
  parName[7]  = "qoverp"; 
  
  /*
  par[0] = 4.81799e-03;            
  stepSize[0] = 0.0001;       
  minVal[0] = -0.1;         
  maxVal[0] = 0.1;
  parName[0] = "x";

  par[1] = 6.86772e-03;          
  stepSize[1] = 0.0001;     
  minVal[1] = -0.1;         
  maxVal[1] = 0.1;
  parName[1] = "y";

  par[2]      = 3.24594e-01;       
  stepSize[2] = 0.01;      
  minVal[2]   = -1.6;      
  maxVal[2]   = 2.6;
  parName[2]  = "deltaKpi";

  par[3]      = 3.90000e-01; 
  stepSize[3] = 0.01;
  minVal[3]   = -3;
  maxVal[3]   = 2;
  parName[3]  = "deltaKpipi0";

  par[4]      = -6.23232e-02; 
  stepSize[4] = 0.001;
  minVal[4]   = -3.15;
  maxVal[4]   = 3.15;
  parName[4]  = "phi"; 

  par[5]      = 3.56763e-03; 
  stepSize[5] = 0.0001;
  minVal[5]   = -0.01;
  maxVal[5]   = 0.01;
  parName[5]  = "rsubdm";

  par[6]      = 3.54714e-03; 
  stepSize[6] = 0.0001;
  minVal[6]   = -0.01;
  maxVal[6]   = 0.01; 
  parName[6]  = "rsubdp"; 

  par[7]      = 9.51249e-01; 
  stepSize[7] = 0.01;
  minVal[7]   = 0.75;
  maxVal[7]   = 1.25; 
  parName[7]  = "qoverp"; 
  */
  for (int i = 0; i < MixingResult::nParams; i++) {
    MixingResult::isSensitive[i] = false; 

    for (MixingResult::ResultIterator r = MixingResult::begin(); r != MixingResult::end(); ++r) {
      if (!(*r)->isActive()) continue; 
      if (!(*r)->isSensitiveTo(i)) continue;
      std::cout << (*r)->getName() << " is sensitive to " << parName[i] << std::endl; 
      MixingResult::isSensitive[i] = true;
      break; 
    }
  }

  //MixingResult::isSensitive[2] = false; 
  
  for (int i = 0; i < MixingResult::nParams; i++) {
    MixingResult::minuit->DefineParameter(i, parName[i].c_str(), par[i], stepSize[i], minVal[i], maxVal[i]);
    //MixingResult::minuit->FixParameter(i);
    if (!MixingResult::isSensitive[i]) MixingResult::minuit->FixParameter(i); 
  }

  MixingResult::minuit->SetFCN(MixChisqFcn);
  int dummy = 0; 
  MixingResult::minuit->mncomd("SET STR 2", dummy); 
  assert(0 == dummy); 

  if (0 < numActive) {
    MixingResult::minuit->Migrad();
  }
  for (int i=0; i<MixingResult::nParams; i++){
    MixingResult::minuit->GetParameter(i,outpar[i],minerr[i]);
    std::cout << " i, outpar[i], err[i] =  "
              << i << ":  "
              << outpar[i] << " +- "
              << minerr[i]
              << std::endl;
  }

  double fitpar[MixingResult::nParams];            
  double fiterr[MixingResult::nParams];            
  double chisq;
  int npars = MixingResult::nParams;

  for (int i = 0; i < MixingResult::nParams; ++i) {
    MixingResult::minuit->GetParameter(i, fitpar[i], fiterr[i]);
  }  

  printChisq = false;
  MixChisqFcn(npars, 0, chisq, fitpar, 4);

  //#define DEBUGCORR 1
#ifdef DEBUGCORR

  debugCorrelations = true;
  MixChisqFcn(npars, 0, chisq, fitpar, 4);

  double step = 0.1;

  for (double xdiff = 0*step*fiterr[0]; xdiff <= 5*step*fiterr[0]; xdiff += step*fiterr[0]) {
    for (double ydiff = 0*step*fiterr[1]; ydiff < 5*step*fiterr[1]; ydiff += step*fiterr[1]) {
      MixingResult::minuit->GetParameter(0, fitpar[0], fiterr[0]);
      MixingResult::minuit->GetParameter(1, fitpar[1], fiterr[1]);
      fitpar[0] += xdiff;
      fitpar[1] += ydiff; 

      MixChisqFcn(npars, 0, chisq, fitpar, 4);
    }
  }

#endif 

  if (skipGraphics) return 0; 

#define DRAWSTUFF 1
#ifdef DRAWSTUFF

  debugCorrelations = false;
  MixDrawer* drawer = new MixDrawer(); 

  DrawOptions* prevOpt = 0; 
  while (!grpOpts.eof()) {
    grpOpts >> token; 
    if (token[0] == '#') {
      grpOpts.getline(token, 10000); 
      continue; 
    }
    std::string lineType(token); 
    if (lineType == "size") {
      grpOpts >> drawer->xmin >> drawer->ymin >> drawer->xmax >> drawer->ymax;
    }    
    else if (lineType == "directory") {
      grpOpts >> drawer->directory;
    }
    else if (lineType == "include") {
      char name[1000];
      int col; 
      grpOpts >> name;
      std::cout<<"name = "<<name<<std::endl;//ad 8/18/13
      std::string sname(name); 
      MixingResult* toDraw = MixingResult::getByName(sname); 
      std::cout<<"toDraw->getName() = "<<toDraw->getName()<<std::endl;//AD 8/18/13
      if (!toDraw) {
	std::cout << "Could not find result " << name << std::endl;
	assert(toDraw); 
      }

      prevOpt = new DrawOptions();
      grpOpts >> name;
      std::cout<<"name = "<<name<<std::endl;//ad 8/18/13
      MixingResult* second = MixingResult::getByName(name);
      //std::cout<<"second->getName() = "<<second->getName()<<std::endl;//AD 8/18/13
      if (second) prevOpt->drawWith.push_back(second);
      grpOpts >> name;
      prevOpt->drawType = DrawOptions::getDrawType(name); 
      grpOpts >> name >> col >> prevOpt->numContours;
      std::cout<<"prevOpt->numContours = "<<prevOpt->numContours<<std::endl;//ad 8/18/13
      prevOpt->setColor(name, col); 
      drawer->addResult(toDraw, prevOpt); 
      
    }
    else if (lineType == "xaxis") {
      grpOpts >> MixDrawer::graphicsXIndex >> drawer->xaxisTitle;
      std::cout<<"MixDrawer::graphicsXIndex = "<<MixDrawer::graphicsXIndex<<std::endl;//ad 8/19/13
    }
    else if (lineType == "yaxis") {
      grpOpts >> MixDrawer::graphicsYIndex >> drawer->yaxisTitle; 
    }
    else if (lineType == "ymult") {
      grpOpts >> MixDrawer::ymult; 
    }    
    else if (lineType == "addtoset") {
      assert(prevOpt); 
      char name[1000];
      grpOpts >> name; 
      MixingResult* second = MixingResult::getByName(name);
      if (second) prevOpt->drawWith.push_back(second);
      //std::cout<<"second->getName() = "<<second->getName()<<std::endl;//ad 8/18/13
    }
    else if (lineType == "pointsPerContour") {
      grpOpts >> MixDrawer::pointsPerContour; 
    }
    else if (lineType == "contours") {
      grpOpts >> fiveSigma >> fourSigma >> threeSigma >> twoSigma >> oneSigma; 
    }
    else if (lineType == "endfile") break; 
  }
  
  drawer->draw(); 
#endif 
  return 0; 
}
