#ifndef MIXINGFIT_HH
#define MIXINGFIT_HH

#include <string>
#include <vector>
#include <cassert> 
#include <map>
#include "TMatrix.h"
#include "TMatrixTSym.h"
#include "TMinuit.h" 
#include <iostream> 
#include "TH2F.h"
class TCanvas; 
class TGraph; 

//typedef __float128 matrixfp;  // I wish.
typedef double matrixfp;  

extern double fiveSigma;
extern double fourSigma;
extern double threeSigma;
extern double twoSigma;
extern double oneSigma;

void MixChisqFcn (int& npar, double* deriv, double& func, double param[], int flag);

class MixingResult {
  friend void MixChisqFcn (int& npar, double* deriv, double& func, double param[], int flag);
  friend class MixDrawer;
public:
  MixingResult (const char* n, 
		const char* rtype, 
		const char* ptype, 
		double res,
		double stat,
		double syst,
		double modl,
		double coef); 
  ~MixingResult ();

  typedef std::vector<MixingResult*>::iterator ResultIterator; 
  enum ResultType {PLAINX, PLAINY, YCP, AGAMMA, XPRIME, XPRIME_P, XPRIME_M, YPRIME, YPRIME_P, YPRIME_M, COSANGLE, SINANGLE, ANGLE, NODCPV, RDM, RDP, RD, AD, QP, PHI, XSQUARE, NUMTYPES}; 
  enum PrimeType {KPI, KPIPI0, QOVERP, NOPRIME, NUMPRIMES};
  enum Sensitivity {EKS = 0, WYE, DELTAKPI, DELTAKPIPI, PHI_12, RSUBDM, RSUBDP, MAGQOVERP, NUMSENSE}; 
  enum CpvAllowed {NOCPV = 0, INDIRECT_CPV, ALL_CPV, NUMCPV}; 
  enum ConstraintType {BOTH_FREE = 0, PHI_FREE, QP_FREE, DELTAGAMMA_FREE, NUMCONSTRAINTS}; 

  std::string getName () const {return name;} 
  double calcEpsilon (double x, 
		      double y, 
		      double deltaKpi, 
		      double deltaKpipi, 
		      double phi12, 
		      double rsubdm, 
		      double rsubdp,
		      double qoverp); 
  double correlation (MixingResult* dat); 
  void setCorrelation (MixingResult* dat, double corr); 
  void setActive (bool t) {if (t != active) initialised = false; active = t;} 
  bool isActive () const {return active;} 
  bool isSensitiveTo (int param) const; 
  void randomise (); 

  static ResultIterator begin () {return allResults.begin();}
  static ResultIterator end   () {return allResults.end();} 

  static ResultType getResultType (std::string n); 
  static PrimeType  getPrimeType  (std::string n); 
  static MixingResult* getByName  (std::string n) {return mapResults[n];} 
  static bool initialised; 
  static void initialise (); 
  static const int nParams; 
  static TMinuit* minuit; 
  static bool* isSensitive; 

private:
  std::string name; 
  double measurement;
  double error;
  std::map<MixingResult*, double> correlations; 
  int index; 
  ResultType myType;
  PrimeType myPrime; 
  bool active; 
  double coefficient; // For use with bands, so that measurement = x + coefficient*y. 

  static TMatrixTSym<matrixfp> sigma; 
  static TMatrixTSym<matrixfp> invsigma; 
  static std::vector<double> epsilon; 
  static std::vector<MixingResult*> allResults; 
  static std::map<std::string, MixingResult*> mapResults;
  static MixingResult* fitResult; 
};

struct DrawOptions {
  enum DrawingTypes {Ellipse = 0, Ybar, Xbar, Band, Annulus, AnnuKidney, Special, Unknown};

  std::vector<MixingResult*> drawWith; 
  DrawingTypes drawType; 
  int colour; 
  int numContours; 
  
  void setColor (std::string base, int mod); 
  static DrawingTypes getDrawType (std::string name); 
  
}; 

class MixDrawer {
public:
  MixDrawer (); 

  std::string directory; 
  double xmin;
  double xmax;

  double ymin;
  double ymax; 
  std::string xaxisTitle;
  std::string yaxisTitle; 

  double extraXMultiplier;
  double extraYMultiplier; 

  void addResult (MixingResult* dat, DrawOptions* dis);
  void draw (); 
  static int pointsPerContour; 
  static int graphicsXIndex; // Stores the Minuit parameter index to use as x axis in plots. 
  static int graphicsYIndex; // Same for y axis. 
  static int ymult;          // Zoom level of y axis. 

private:
  void drawResult (MixingResult* dat, DrawOptions* dis, TCanvas* foo);
  TGraph* getEllipse (double errorDef); 
  void findPoint (TGraph* ret, int idx, double angle, double errorDef, int par1, int par2);//ad 8/19/13 
  std::pair<TGraph*, TGraph*> drawEllipse (DrawOptions* dis, TCanvas* foo);
  std::vector<TGraph*> drawEllipse3 (DrawOptions* dis, TCanvas* foo);//ad 8/18/13
  void drawEllipseForce (DrawOptions* dis, TCanvas* foo); 
  void drawEllipseForce_Adam (DrawOptions* dis, TCanvas* foo); 
  void drawYbar (MixingResult* dat, DrawOptions* dis, TCanvas* foo);
  void drawBand (MixingResult* dat, DrawOptions* dis, TCanvas* foo);
  void drawAnnulus (MixingResult* dat, DrawOptions* dis, TCanvas* foo);
  std::map<MixingResult*, DrawOptions*> drawmap; 

  double runFit (); 
};

#endif
