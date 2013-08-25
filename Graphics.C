#include "MixingFit.h" 
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TF1.h"
#include "TPolyLine.h" 
#include "TEllipse.h" 
#include "TBox.h" 
#include "TAxis.h"
#include "TColor.h"
#include "TGaxis.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TPaveText.h"
#include <cmath> 
#include <vector> 
#include <map>

int MixDrawer::pointsPerContour = 3600; 
int MixDrawer::graphicsXIndex = 0; 
int MixDrawer::graphicsYIndex = 1; 
int MixDrawer::ymult = 1; 

MixDrawer::MixDrawer ()
  : directory("testing")
  , xmin(-0.02)
  , xmax(0.02)
  , ymin(0)
  , ymax(0.045)
  , xaxisTitle("x")
  , yaxisTitle("y")
{}


void MixDrawer::draw () {
  //gROOT->ProcessLine(".L ./lhcb_style.cc");

  // Use times new roman, precision 2 
  Int_t lhcbFont        = 132;  // Old LHCb style: 62;
  // Line thickness
  Double_t lhcbWidth    = 2.00; // Old LHCb style: 3.00;
  // Text size
  Double_t lhcbTSize    = 0.06; 
  
  // use plain black on white colors
  gROOT->SetStyle("Plain"); 
  TStyle *lhcbStyle= new TStyle("lhcbStyle","LHCb plots style");
  
  //lhcbStyle->SetErrorX(0); //  don't suppress the error bar along X

  lhcbStyle->SetFillColor(1);
  lhcbStyle->SetFillStyle(1001);   // solid
  lhcbStyle->SetFrameFillColor(0);
  lhcbStyle->SetFrameBorderMode(0);
  lhcbStyle->SetPadBorderMode(0);
  lhcbStyle->SetPadColor(0);
  lhcbStyle->SetCanvasBorderMode(0);
  lhcbStyle->SetCanvasColor(0);
  lhcbStyle->SetStatColor(0);
  lhcbStyle->SetLegendBorderSize(0);

  // If you want the usual gradient palette (blue -> red)
  lhcbStyle->SetPalette(1);
  // If you want colors that correspond to gray scale in black and white:
  int colors[8] = {0,5,7,3,6,2,4,1};
  lhcbStyle->SetPalette(8,colors);

  // set the paper & margin sizes
  lhcbStyle->SetPaperSize(20,26);
  lhcbStyle->SetPadTopMargin(0.05);
  lhcbStyle->SetPadRightMargin(0.09); // increase for colz plots
  lhcbStyle->SetPadBottomMargin(0.16);
  lhcbStyle->SetPadLeftMargin(0.20);
  
  // use large fonts
  lhcbStyle->SetTextFont(lhcbFont);
  lhcbStyle->SetTextSize(lhcbTSize);
  lhcbStyle->SetLabelFont(lhcbFont,"x");
  lhcbStyle->SetLabelFont(lhcbFont,"y");
  lhcbStyle->SetLabelFont(lhcbFont,"z");
  lhcbStyle->SetLabelSize(lhcbTSize,"x");
  lhcbStyle->SetLabelSize(lhcbTSize,"y");
  lhcbStyle->SetLabelSize(lhcbTSize,"z");
  lhcbStyle->SetTitleFont(lhcbFont);
  lhcbStyle->SetTitleFont(lhcbFont,"x");
  lhcbStyle->SetTitleFont(lhcbFont,"y");
  lhcbStyle->SetTitleFont(lhcbFont,"z");
  lhcbStyle->SetTitleSize(1.2*lhcbTSize,"x");
  lhcbStyle->SetTitleSize(1.2*lhcbTSize,"y");
  lhcbStyle->SetTitleSize(1.2*lhcbTSize,"z");

  // use medium bold lines and thick markers
  lhcbStyle->SetLineWidth(lhcbWidth);
  lhcbStyle->SetFrameLineWidth(lhcbWidth);
  lhcbStyle->SetHistLineWidth(lhcbWidth);
  lhcbStyle->SetFuncWidth(lhcbWidth);
  lhcbStyle->SetGridWidth(lhcbWidth);
  lhcbStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  lhcbStyle->SetMarkerStyle(20);
  lhcbStyle->SetMarkerSize(1.0);

  // label offsets
  lhcbStyle->SetLabelOffset(0.010,"X");
  lhcbStyle->SetLabelOffset(0.010,"Y");

  // by default, do not display histogram decorations:
  lhcbStyle->SetOptStat(0);  
  //lhcbStyle->SetOptStat("emr");  // show only nent -e , mean - m , rms -r
  // full opts at http://root.cern.ch/root/html/TStyle.html#TStyle:SetOptStat
  lhcbStyle->SetStatFormat("6.3g"); // specified as c printf options
  lhcbStyle->SetOptTitle(0);
  lhcbStyle->SetOptFit(0);
  //lhcbStyle->SetOptFit(1011); // order is probability, Chi2, errors, parameters
  //titles
  lhcbStyle->SetTitleOffset(1.05,"X");
  lhcbStyle->SetTitleOffset(1.55,"Y");
  lhcbStyle->SetTitleOffset(1.2,"Z");
  lhcbStyle->SetTitleFillColor(0);
  lhcbStyle->SetTitleStyle(0);
  lhcbStyle->SetTitleBorderSize(0);
  lhcbStyle->SetTitleFont(lhcbFont,"title");
  lhcbStyle->SetTitleX(0.0);
  lhcbStyle->SetTitleY(1.0); 
  lhcbStyle->SetTitleW(1.0);
  lhcbStyle->SetTitleH(0.05);
  
  // look of the statistics box:
  lhcbStyle->SetStatBorderSize(0);
  lhcbStyle->SetStatFont(lhcbFont);
  lhcbStyle->SetStatFontSize(0.05);
  lhcbStyle->SetStatX(0.9);
  lhcbStyle->SetStatY(0.9);
  lhcbStyle->SetStatW(0.25);
  lhcbStyle->SetStatH(0.15);

  // put tick marks on top and RHS of plots
  lhcbStyle->SetPadTickX(1);
  lhcbStyle->SetPadTickY(1);

  // histogram divisions: only 5 in x to avoid label overlaps
  lhcbStyle->SetNdivisions(505,"x");
  lhcbStyle->SetNdivisions(510,"y");
  
  gROOT->SetStyle("lhcbStyle");
  gROOT->ForceStyle();

  // add LHCb label
  TPaveText *lhcbName = new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
                           0.87 - gStyle->GetPadTopMargin(),
                           gStyle->GetPadLeftMargin() + 0.20,
                           0.95 - gStyle->GetPadTopMargin(),
                           "BRNDC");
  lhcbName->AddText("LHCb");
  lhcbName->SetFillColor(0);
  lhcbName->SetTextAlign(12);
  lhcbName->SetBorderSize(0);

  TText *lhcbLabel = new TText();
  lhcbLabel->SetTextFont(lhcbFont);
  lhcbLabel->SetTextColor(1);
  lhcbLabel->SetTextSize(lhcbTSize);
  lhcbLabel->SetTextAlign(12);

  TLatex *lhcbLatex = new TLatex();
  lhcbLatex->SetTextFont(lhcbFont);
  lhcbLatex->SetTextColor(1);
  lhcbLatex->SetTextSize(lhcbTSize);
  lhcbLatex->SetTextAlign(12);

  std::cout << "-------------------------" << std::endl;  
  std::cout << "Set LHCb Style - Feb 2012" << std::endl;
  std::cout << "-------------------------" << std::endl;  

  /*  
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleColor(1);
  gStyle->SetStatColor(0);
  gStyle->SetFillColor(0);
  gStyle->SetFuncWidth(1);
  gStyle->SetLineWidth(1);
  gStyle->SetLineColor(1);
  gStyle->SetPalette(1, 0);
  gStyle->SetOptStat(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleSize(0.04);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadRightMargin(0.13);
  gStyle->SetPadTopMargin(0.13);
  gStyle->SetLineStyleString(11, "0 5000");
  gStyle->SetLineStyleString(12, "100 50");
  TGaxis::SetMaxDigits(2);*/
  TCanvas* foo = new TCanvas();
  foo->Size(40,40);
  foo->SetFillColor(0);
  foo->SetBorderMode(0);     
  foo->SetFrameBorderMode(0);

  // Kludge to get nice axis endpoints. 
  double XMIN = (xmax+11*xmin)/12;
  double XMAX = (xmin+11*xmax)/12;
  double YMIN = (ymax+11*ymin)/12;
  double YMAX = (ymin+11*ymax)/12;
  double xpo[4] = {XMIN, XMIN, XMAX, XMAX};
  double ypo[4] = {YMIN, YMAX, YMIN, YMAX};
  TGraph* gr23 = new TGraph(4, xpo, ypo);
  gr23->SetTitle("");
  gr23->SetFillColor(38+150);
  XMIN = gr23->GetXaxis()->GetXmin();
  XMAX = gr23->GetXaxis()->GetXmax();
  YMIN = gr23->GetYaxis()->GetXmin();
  YMAX = gr23->GetYaxis()->GetXmax();
  gr23->GetXaxis()->SetTitle(xaxisTitle.c_str());
  gr23->GetYaxis()->SetTitle(yaxisTitle.c_str());
  gr23->Draw("alf");
  
  TF1* multXBy1k = new TF1("multXBy1k", "x", xmin, xmax);
  TF1* multYBy1k = new TF1("multYBy1k", "1000*x", ymin*ymult, ymax*ymult);
  /*  TGaxis* xaxis = new TGaxis(xmin,ymin,xmax,ymin,"multXBy1k",510,"+");
  TGaxis* yaxis = new TGaxis(xmin,ymin,xmin,ymax,"multYBy1k",510,"-");
  
  xaxis->SetTitle((xaxisTitle+ " #times 10^{3}").c_str()); 
  xaxis->CenterTitle();
  xaxis->SetTitleSize(0.06);
  xaxis->SetLabelSize(0.06);
  if (1000 == ymult) yaxis->SetTitle((yaxisTitle + " #times 10^{3}").c_str()); 
  else yaxis->SetTitle(yaxisTitle.c_str()); 
  yaxis->CenterTitle();
  yaxis->SetTitleSize(0.06);
  yaxis->SetLabelSize(0.06);
  yaxis->SetTitleOffset(1.1); 
  xaxis->SetTitleOffset(0.9); 
  gr23->GetXaxis()->SetLabelColor(kWhite);
  gr23->GetXaxis()->SetTitleColor(kWhite);
  gr23->GetXaxis()->SetAxisColor(kWhite);
  gr23->GetYaxis()->SetLabelColor(kWhite);
  gr23->GetYaxis()->SetTitleColor(kWhite);
  gr23->GetYaxis()->SetAxisColor(kWhite);

  xaxis->Draw();
  yaxis->Draw();*/
  foo->SaveAs((directory + "/axisplot.png").c_str()); 
  delete multXBy1k;
  delete multYBy1k;

  for (std::map<MixingResult*, DrawOptions*>::iterator i = drawmap.begin(); i != drawmap.end(); ++i) {
    if ((*i).first == MixingResult::fitResult) continue; 
    drawResult((*i).first, (*i).second, foo); 
  }

  if (drawmap[MixingResult::fitResult]) {
    if (MixingResult::isSensitive[MixingResult::WYE]) {
      if (MixingResult::isSensitive[MixingResult::EKS]) {
	//drawEllipse3(drawmap[MixingResult::fitResult], foo); 
	drawEllipseForce(drawmap[MixingResult::fitResult], foo); 
	//drawEllipse(drawmap[MixingResult::fitResult], foo); 
      }
      else {
	MixingResult::minuit->GetParameter(MixingResult::WYE, MixingResult::fitResult->measurement, MixingResult::fitResult->error);
	drawYbar(MixingResult::fitResult, drawmap[MixingResult::fitResult], foo);
      }
    }
  }

  foo->SaveAs((directory + "/finalplot.png").c_str()); 
  //foo->SaveAs((directory + "/finalplot.pdf").c_str()); 
  delete foo; 
}


void MixDrawer::drawAnnulus (MixingResult* dat, DrawOptions* dis, TCanvas* foo) {
  assert(1 == dis->drawWith.size()); 
  MixingResult* xres = dat;
  MixingResult* yres = dis->drawWith[0];
  switch (xres->myType) {
  case MixingResult::XPRIME:
    break;
  case MixingResult::YPRIME:
    xres = yres;
    yres = dat;
    break;
  default:
    std::cout << "Error: Cannot draw annulus with non-prime result " << dat->getName() << std::endl; 
    assert(false);
  }
  if (yres->myType != MixingResult::YPRIME) {
    std::cout << "Error: Cannot draw annulus with non-prime result " << dat->getName() << std::endl;
    assert(false);
  }


  // Use separate fitter to avoid side effects on main result. 
  static TMinuit* annuMin = 0;
  if (0 == annuMin) {
    annuMin = new TMinuit(MixingResult::nParams); 
    annuMin->DefineParameter(0, "special_x",      0.005, 0.0001, -0.10, 0.1);
    annuMin->DefineParameter(1, "special_y",      0.005, 0.0001, -0.10, 0.1);
    annuMin->DefineParameter(2, "special_dkpi",   0.800, 0.0100, -0.05, 6.5); 
    annuMin->DefineParameter(3, "special_dkpipi", 0.800, 0.0100, -0.05, 6.5); 
    annuMin->SetFCN(MixChisqFcn); 
  }
  
  //annuMin->SetPrintLevel(-1);
  char cmdbuf[1000];
  int cmdres = 0; 
  int didx = -1; 
  switch (xres->myPrime) {
  case MixingResult::KPI:    didx = 3; break;
  case MixingResult::KPIPI0: didx = 4; break; // NB, these are for use with MINUIT which is Fortran indexed. 
  default: break;
  }

  if (-1 == didx) {
    std::cout << "Error: Cannot draw annulus with non-prime result " << xres->getName() << " " << xres->myPrime << std::endl; 
    assert(false); 
  }
  
  annuMin->FixParameter(2); 
  annuMin->FixParameter(3); 

  TMinuit* temp = MixingResult::minuit;
  MixingResult::minuit = annuMin; 

  double anglestep = 0.001; 
  for (double delta = 0; delta < 6.281; delta += anglestep) {
    // For this value of delta, what (x, y) are compatible with
    // the measured x', y'? 

    sprintf(cmdbuf, "SET PAR %i %f", didx, delta); 
    annuMin->mncomd(cmdbuf, cmdres);
    if (0 != cmdres) std::cout << "Problem with command, returned " << cmdres << std::endl; 
    TGraph* curr = getEllipse(2.3);
    curr->SetFillColor(dis->colour);
    curr->Draw("if");
  }

  if (DrawOptions::AnnuKidney == dis->drawType) {
    double delta = 0;
    double delerr = 0;
    MixingResult::minuit->GetParameter(didx-1, delta, delerr); // Back to C indexing
    sprintf(cmdbuf, "SET PAR %i %f", didx, delta); 
    annuMin->mncomd(cmdbuf, cmdres);
    if (0 != cmdres) std::cout << "Problem with command, returned " << cmdres << std::endl; 
    annuMin->mnfree(didx); 
    TGraph* curr = getEllipse(2.3);
    curr->SetFillColor(dis->colour+4);
    curr->Draw("if");
  }

  MixingResult::minuit = temp; 

}

/*
void MixDrawer::drawAnnulus (MixingResult* dat, DrawOptions* dis, TCanvas* foo) {
  assert(1 == dis->drawWith.size()); 
  MixingResult* xres = dat;
  MixingResult* yres = dis->drawWith[0];
  switch (xres->myType) {
  case MixingResult::XPRIME:
    break;
  case MixingResult::YPRIME:
    xres = yres;
    yres = dat;
    break;
  default:
    std::cout << "Error: Cannot draw annulus with non-prime result " << dat->getName() << std::endl; 
    assert(false);
  }
  if (yres->myType != MixingResult::YPRIME) {
    std::cout << "Error: Cannot draw annulus with non-prime result " << dat->getName() << std::endl;
    assert(false);
  }

  double bestx = xres->measurement;
  double besty = yres->measurement; 
  if ((MixingResult::KPI == xres->myPrime) || (bestx < 0)) {
    // Move result into physically allowed region
    double sigmas = -1*(bestx / xres->error); 
    bestx = 0; 
    besty += sigmas * yres->error * yres->correlation(xres); 
  }

  if (MixingResult::KPI == xres->myPrime) bestx = sqrt(bestx); 

  double fitpar[MixingResult::nParams];
  double fiterr[MixingResult::nParams];
  double chisq;
  int npars = MixingResult::nParams;

  fitpar[0] = bestx;
  fitpar[1] = besty; 
  for (int i = 0; i < MixingResult::nParams; ++i) {
    MixingResult::minuit->GetParameter(i, fitpar[i], fiterr[i]);
  }
  MixChisqFcn(npars, 0, chisq, fitpar, 4);
  double minimumChisq = chisq; 
  double normDist  = sqrt(pow(bestx, 2) + pow(besty, 2));
  double maxLoDist = 1;
  double minLoDist = 0;
  double minHiDist = 1; 
  double maxHiDist = 1000;

  fitpar[0] = 0;
  fitpar[1] = 0;
  MixChisqFcn(npars, 0, chisq, fitpar, 4);
  const double tolerance = 0.01; 
  if (chisq - minimumChisq < 2.3) {
    maxLoDist = 0; 
  }
  else {
    for (int i = 0; i < 1000; ++i) {
      double currDist = minLoDist + 0.5*(maxLoDist - minLoDist);
      fitpar[0] = bestx * currDist;
      fitpar[1] = besty * currDist;
      MixChisqFcn(npars, 0, chisq, fitpar, 4);
      if (chisq - minimumChisq < 2.3) maxLoDist = currDist;
      else minLoDist = currDist; 

      std::cout << "Annulus minimum [" << minLoDist << ", " << maxLoDist << "] " << (chisq - minimumChisq) << std::endl;
      if (maxLoDist - minLoDist < tolerance) break;
    }
  }

  for (int i = 0; i < 1000; ++i) {
    double currDist = minHiDist + 0.5*(maxHiDist - minHiDist);
    fitpar[0] = bestx * currDist;
    fitpar[1] = besty * currDist;
    MixChisqFcn(npars, 0, chisq, fitpar, 4);
    if (chisq - minimumChisq < 2.3) minHiDist = currDist;
    else maxHiDist = currDist;
    std::cout << "Annulus maximum [" << minHiDist << ", " << maxHiDist << "] " << (chisq - minimumChisq) << std::endl;
    if (maxHiDist - minHiDist < tolerance) break;
  }

  double minRadius = 0.5*(minLoDist + maxLoDist)*normDist;
  double maxRadius = 0.5*(minHiDist + maxHiDist)*normDist;

  std::cout << normDist << " " << bestx << " " << besty << " " << minRadius << " " << maxRadius << std::endl; 
  double step = (maxRadius - minRadius)*0.005;
  for (double rad = maxRadius; rad >= minRadius; rad -= step) {
    std::cout << "Ellipse radius " << rad << std::endl; 
    TEllipse* currEllipse = new TEllipse(0, 0, rad); 
    currEllipse->SetLineColor(dis->colour); 
    currEllipse->Draw(); 
  }  
}
*/

void MixDrawer::drawResult (MixingResult* dat, DrawOptions* dis, TCanvas* foo) {
  std::vector<MixingResult*> actives; 
  for (MixingResult::ResultIterator i = MixingResult::begin(); i != MixingResult::end(); ++i) {
    if((*i)->isActive()){ actives.push_back(*i);
      std::cout<<"i->getName() = "<<(*i)->getName()<<std::endl;//ad 8/19/13
    }
    (*i)->setActive(false); 
  }
  dat->setActive(true); 
  for (std::vector<MixingResult*>::iterator i = dis->drawWith.begin(); i != dis->drawWith.end(); ++i) {
    (*i)->setActive(true); 
  }

  switch (dis->drawType) {
  case DrawOptions::Ellipse:
    {
      assert(dis->drawWith.size() > 0); 
      //std::pair<TGraph*, TGraph*> ells = drawEllipse(dis, foo); 
      std::vector<TGraph*> ells = drawEllipse3(dis, foo); 
    }
    break; 
  case DrawOptions::Ybar:
    drawYbar(dat, dis, foo); 
    break;
  case DrawOptions::Xbar:
  case DrawOptions::Band:
    drawBand(dat, dis, foo); 
    break; 
  case DrawOptions::Annulus:
  case DrawOptions::AnnuKidney:
    drawAnnulus(dat, dis, foo); 
    break; 
  default:
  case DrawOptions::Unknown: 
    std::cout << "Warning: Could not draw result " << dat->name << "; unknown drawing type." << std::endl; 
    break;
  }


  for (std::vector<MixingResult*>::iterator i = actives.begin(); i != actives.end(); ++i) {
    (*i)->setActive(true); 
  }
}

void MixDrawer::findPoint (TGraph* ret, int idx, double angle, double errorDef, int par1, int par2) {
  //we're finding a point at a certain angle here, but we don't know where R is, so we have to iterate to find it.
  double fitpar[MixingResult::nParams];//first get the parameters
  double fiterr[MixingResult::nParams];//next get the error
  double chisq;//initialize the chi2
  int npars = MixingResult::nParams;//get the number of parameters
  
  for (int i = 0; i < MixingResult::nParams; ++i) {
    MixingResult::minuit->GetParameter(i, fitpar[i], fiterr[i]);//get all the results
    //std::cout<<"i = "<<i<<", fitpar["<<i<<"] = "<<fitpar[i]<<", fiterr["<<i<<"] = "<<fiterr[i]<<std::endl;
  }
  
  MixChisqFcn(npars, 0, chisq, fitpar, 4);//get the chi2. Derivative (0) isn't used in the FCN, so it doesn't matter. What is flag=4?

  double initX = fitpar[par1];//this is the initial x position
  double initY = fitpar[par2]; //this is the initial y position
  double xincrement = cos(angle);//this is the cos of the angle we want
  double yincrement = sin(angle); //sin of the angle we want
  double minimum = 0;//minimum chi2
  double chiAtMinimum = chisq; //chi2 at the minimum
  double chiAtSolution = chisq; //chi2 at the current solution
  //  std::cout<<"initX = "<<initX<<", initY = "<<initY<<", xincrement = "<<xincrement<<", yincrement = "<<yincrement<<std::endl;//ad 8/22/13
  double maximum = 0.001; //this is the radius that we want to test
  double chiAtMaximum = chisq; 
  while (chiAtMaximum < chiAtSolution + errorDef) {//break out of the loop as soon as we exceed the chi2. The error def is then the delta chi2
    maximum *= 2; 
    fitpar[0] = initX + maximum*xincrement;//x+(r cos(angle))
    fitpar[1] = initY + maximum*yincrement;//y+(r sin(angle)
    MixChisqFcn(npars, 0, chisq, fitpar, 4);//check the fit
    chiAtMaximum = chisq; 
  }
  //we have the maximum
  //ret->SetPoint(idx, fitpar[0],fitpar[1]); 
  // zero in for higher tolerance
  static const double tolerance = 0.0001; //what is this?
  for (int i = 0; i < 1000; ++i) {//iterate 1000 times
    double currDist = minimum + (maximum - minimum)*0.5; //take mean of max and min
    fitpar[0] = initX + currDist*xincrement;//x= x_0 +currDist*cos(angle)
    fitpar[1] = initY + currDist*yincrement;//y= y_0+ currDist*sin(angle)
    MixChisqFcn(npars, 0, chisq, fitpar, 4);
    if (chisq > chiAtSolution + errorDef) {//if we found a better delta chi2
      maximum = currDist;
      chiAtMaximum = chisq;
    }
    else {//otherwise update the minimum
      minimum = currDist;
      chiAtMinimum = chisq;
    }
  
    if (1 == idx) std::cout << "Scan: currDist =  " << currDist << ", delta Chi2 =  " << (chisq - chiAtSolution) << ",chi2 at solution =  " << chiAtSolution << std::endl; 
    assert(chiAtMaximum > chiAtMinimum);
    if (chiAtMaximum - chiAtMinimum < tolerance) break; 
  }

  double dist = minimum + 0.5*(maximum-minimum); 
  std::cout << "Found point " << initX + dist*xincrement << ", " << initY + dist*yincrement << std::endl;
  ret->SetPoint(idx, initX + dist*xincrement, initY + dist*yincrement); 
  
}

bool intersect (TLine& one, TLine& two) {
  // Here I have special-case knowledge that one is always horizontal. 
  if ((two.IsHorizontal()) && (one.GetY1() != two.GetY1())) return false;
  // So now they cannot be parallel, they must meet somewhere. 

  double x1 = one.GetX1();
  double x2 = one.GetX2();
  double x3 = two.GetX1();
  double x4 = two.GetX2();

  double y1 = one.GetY1();
  double y2 = one.GetY2();
  double y3 = two.GetY1();
  double y4 = two.GetY2();

  double det = (x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4);
  if (fabs(det) < 1e-12) return false; // Almost parallel

  double xInter = (x1*y2 - y1*x2)*(x3 - x4) - (x1 - x2)*(x3*y4 - y3*x4);
  double yInter = (x1*y2 - y1*x2)*(y3 - y4) - (y1 - y2)*(x3*y4 - y3*x4);

  xInter /= det;
  yInter /= det;

  // Check that intersection occurs on line segment one. 
  if (xInter < x1) return false;
  if (xInter > x2) return false;
  if (yInter < y1) return false;
  if (yInter > y2) return false;  

  // And segment two. 
  if (xInter < x3) return false;
  if (xInter > x4) return false;
  if (yInter < y3) return false;
  if (yInter > y4) return false;  

  /*
  std::cout << "Intersect of (" 
	    << x1 << ", " << y1 << ") to (" << x2 << ", " << y2 << ") with ("
	    << x3 << ", " << y3 << ") to (" << x4 << ", " << y4 << ") at point ("
	    << xInter << ", " << yInter << ")\n"; 
  */ 
  return true; 
}

bool isWithin (double xp, double yp, TGraph* poly, int npoints) {
  // Returns true if the point (xp, yp) is within the polygon
  // described by the first npoints points of poly. 

  if (npoints < 3) return false; 
  int intersects = 0;
  TLine one(xp, yp, xp+100000, yp);
  for (int i = 1; i < npoints; ++i) {
    double p1, p2, p3, p4;
    poly->GetPoint(i-1, p1, p2);
    poly->GetPoint(i,   p3, p4);
    TLine two(p1, p2, p3, p4);
    if (intersect(one, two)) intersects++; 
  }
  return (1 == intersects % 2);
}

TGraph* MixDrawer::getEllipse (double errorDef) {
  MixingResult::minuit->SetErrorDef(errorDef);
  TGraph* ret = (TGraph*) MixingResult::minuit->Contour(pointsPerContour, graphicsXIndex, graphicsYIndex);
  if (!ret) return 0; 

  std::cout<<"calling rolf's find point with errordef " << errorDef << std::endl;
  if(graphicsXIndex==0){
    if ((!ret) || (ret->GetN() < pointsPerContour)) {
      double xp = 0;
      double yp = 0;
      // Discard doubly-wrapped points. 
      bool endOfLine = false; 
      for (int i = 1; i < ret->GetN(); ++i) {
	if (endOfLine) {
	  ret->RemovePoint(i);
	  i--;
	  continue; 
	}
	ret->GetPoint(i, xp, yp);
	//std::cout << "Testing point " << i << " (" << xp << ", " << yp << ") ... ";
	if (isWithin(xp, yp, ret, i-1)) {
	  //std::cout << " discarding\n"; 
	  ret->RemovePoint(i);
	  i--;
	  endOfLine = true; 
	}
	//else std::cout << " keeping\n"; 
      }
    }
  }
  
  return ret; 
}

std::pair<TGraph*, TGraph*> MixDrawer::drawEllipse (DrawOptions* dis, TCanvas* foo) {
  std::pair<TGraph*, TGraph*> ret(0, 0); 

  /*
  for (int i = 2; i < MixingResult::NUMSENSE; ++i) {
    bool sensitive = false; 
    for (MixingResult::ResultIterator m = MixingResult::begin(); m != MixingResult::end(); ++m) {
      if (!(*m)->isActive()) continue;
      if (!(*m)->isSensitiveTo(i)) continue;
      sensitive = true;
      break;
    }
    if (!sensitive) MixingResult::minuit->FixParameter(i);
  }
  MixingResult::minuit->FixParameter(2);
  MixingResult::minuit->FixParameter(3);
  */
  MixingResult::initialised = false; 

  std::cout << "Drawing " << dis->numContours << " ellipses\n";

  //  In 2 dimensions, an error ellipse nominally
  //  providing 68.27% coverage [equivalent to
  //  1 sigma in 1 dimension] has chisq=2.30
  MixingResult::minuit->Migrad(); 
  ret.second = getEllipse(6.18);
  if (ret.second) ret.second->SetFillColor(dis->colour + 2);
  ret.first = getEllipse(2.30);
  if (ret.first) ret.first->SetFillColor(dis->colour);

  if ((dis->numContours > 1) && (ret.second)) {
    ret.second->Draw("if9");
    //std::cout << "Outer ellipse drawn\n"; 
  }
  if (ret.first) {
    ret.first->Draw("if9");
    //std::cout << "Inner ellipse drawn\n"; 
  }

  MixingResult::minuit->SetErrorDef(1.00);
  return ret; 
}

double MixDrawer::runFit () {
  static double fitpar[10];
  static double fiterr[10];
  static int npar = MixingResult::nParams;
  MixingResult::initialised = false; 
  MixingResult::minuit->Migrad(); 
  for (int p = 0; p < MixingResult::nParams; ++p) {
    MixingResult::minuit->GetParameter(p, fitpar[p], fiterr[p]);
  }
  double chisq = 0; 
  MixChisqFcn(npar, 0, chisq, fitpar, 4);
  return chisq;
}

void MixDrawer::drawEllipseForce (DrawOptions* dis, TCanvas* foo) {
  double fitpar[MixingResult::nParams];
  double orgpar[MixingResult::nParams];
  double wrkpar[MixingResult::nParams];
  double fiterr[MixingResult::nParams];
  for (int p = 0; p < MixingResult::nParams; ++p) {
    MixingResult::minuit->GetParameter(p, fitpar[p], fiterr[p]);
    orgpar[p] = fitpar[p]; 
    wrkpar[p] = fitpar[p]; 
  }

  double XMIN = (xmax+11*xmin)/12;
  double XMAX = (xmin+11*xmax)/12;
  double YMIN = (ymax+11*ymin)/12;
  double YMAX = (ymin+11*ymax)/12;
  TH2F* histogram = new TH2F("hist", "", 1000, XMIN, XMAX, 1000, YMIN, YMAX); 
  histogram->SetStats(false); 
  TH2F* edmhist = new TH2F("edmhist", "", 1000, XMIN, XMAX, 1000, YMIN, YMAX); 
  edmhist->SetStats(false); 
  TH1F edm1d("edm1d", "", 1000, 0, 0.00001); 
  

  char strbuf[100];
  int dummy = MixingResult::nParams;
  double chisq = 0;
  double maxEdm = 0; 
  MixChisqFcn(dummy, 0, chisq, fitpar, 4);
  double centralChisq = chisq; 
  MixingResult::minuit->FixParameter(0);
  MixingResult::minuit->FixParameter(1);
  MixingResult::minuit->SetErrorDef(1); 
  MixingResult::minuit->SetPrintLevel(-1); 

  for (int i = 1; i <= 1000; ++i) {
    //for (int i = 100; i <= 400; ++i) {
    double prevChisq = 1000; 
    double xval = xmin + (i + 0.5)*(xmax - xmin)*0.001;
    sprintf(strbuf, "SET PAR 1 %f", xval);
    MixingResult::minuit->mncomd(strbuf, dummy);
    assert(0 == dummy);
    if (0 == i%25) std::cout << "Fitting line " << i << std::endl; 
    //for (int j = 200; j <= 1000; ++j) {
      for (int j = 2; j <= 1000; ++j) {
      double yval = ymin + (j + 0.5)*(ymax - ymin)*0.001;
      sprintf(strbuf, "SET PAR 2 %f", yval);
      MixingResult::minuit->mncomd(strbuf, dummy);
      assert(0 == dummy);
      chisq = runFit() - centralChisq;

      if (chisq > fiveSigma) {
	std::cout << "Bad fit for (" << i << ", " << j << ") " << chisq << ", retrying with last working set.\n";
	for (int p = 3; p <= 8; ++p) {
	  sprintf(strbuf, "SET PAR %i %f", p, wrkpar[p-1]);
	  MixingResult::minuit->mncomd(strbuf, dummy);
	  assert(0 == dummy);
	}
	chisq = runFit() - centralChisq;
	if (chisq > 2*fiveSigma) {
	  std::cout << "  Still bad, " << chisq << ", retrying with central values.\n";
	  for (int p = 3; p <= 8; ++p) {
	    sprintf(strbuf, "SET PAR %i %f", p, orgpar[p-1]);
	    MixingResult::minuit->mncomd(strbuf, dummy);
	    assert(0 == dummy);
	  }
	  chisq = runFit() - centralChisq;
	  if (chisq > 3*fiveSigma) {
	    std::cout << "  Final fit is " << chisq << ", giving up.\n"; 
	    chisq = 0; // Give up
	  }
	  else for (int p = 0; p < MixingResult::nParams; ++p) MixingResult::minuit->GetParameter(p, wrkpar[p], fiterr[p]);
	}
	else for (int p = 0; p < MixingResult::nParams; ++p) MixingResult::minuit->GetParameter(p, wrkpar[p], fiterr[p]);
      }
      else for (int p = 0; p < MixingResult::nParams; ++p) MixingResult::minuit->GetParameter(p, wrkpar[p], fiterr[p]);
      
      double fmin, fedm, errdef;
      int nparx, istat; 
      MixingResult::minuit->mnstat(fmin, fedm, errdef, dummy, nparx, istat);
      if (3 > istat) chisq = 0; 
      //if (fedm > 0.001) chisq = 0; 
      //if (chisq >= 1.5*prevChisq) chisq = 0;

      if (chisq > 0) prevChisq = chisq; 
      histogram->SetBinContent(i, j, chisq); 
      edmhist->SetBinContent(i, j, fedm);
      edm1d.Fill(fedm); 
      maxEdm = (fedm > maxEdm ? fedm : maxEdm); 
      //break; 
    }
    //break; 
  }

  /*
  for (int i = 2; i < 1000; ++i) {
    for (int j = 2; j < 1000; ++j) {
      if (0 < histogram->GetBinContent(i, j)) continue;
      double avg = 0; 
      //avg       += histogram->GetBinContent(i+0, j+1);
      //avg       += histogram->GetBinContent(i+0, j-1);
      avg       += histogram->GetBinContent(i+1, j+0);
      avg       += histogram->GetBinContent(i-1, j+0);
      avg *= 0.5;
      histogram->SetBinContent(i, j, avg); 
    }
  }
  */

  for (int i = 2; i < 1000; ++i) {
    for (int j = 2; j < 1000; ++j) {
      double curr = histogram->GetBinContent(i, j);
      if (0 == curr) continue;
      // Values chosen for five-colour plot. 
      if      (curr < oneSigma)   curr = 1;
      else if (curr < twoSigma)   curr = 7;
      else if (curr < threeSigma) curr = 15;
      else if (curr < fourSigma)  curr = 20; 
      else if (curr < fiveSigma)  curr = fiveSigma-1; 
      
      histogram->SetBinContent(i, j, curr); 
    }
  }

  // Search for bad fits no found by earlier methods. For each bin, is it an outlier in its region?
  /*
  for (int i = 2; i < 1000; ++i) {
    for (int j = 2; j < 1000; ++j) {
      double curr = histogram->GetBinContent(i, j);
      if (0 == curr) continue;
      std::map<double, int> surrounds; 
      for (int xp = -1; xp <= 1; ++xp) {
	for (int yp = -1; yp <= 1; ++yp) {
	  if ((0 == xp) && (0 == yp)) continue; 
	  double area = histogram->GetBinContent(i+xp, j+yp);
	  if (0 == area) continue; 
	  surrounds[area]++;
	}
      }
      if (0 < surrounds[curr]) continue; 
      
      histogram->SetBinContent(i, j, 0); 
    }
    }*/

  // Fix areas of bad fits. Take modal value of surrounding bins. 
  for (int i = 2; i < 1000; ++i) {
    for (int j = 2; j < 1000; ++j) {
      if (0 < histogram->GetBinContent(i, j)) continue;
      // Right edge detection
      int numNonzero = 0;
      for (int xp = 1; xp <= 10; ++xp) {
	if (0 == histogram->GetBinContent(i+xp, j)) continue;
	numNonzero++;
      }
      if (0 == numNonzero) continue; // We're at a right edge. 
	
      std::map<double, int> surrounds; 
      for (int xp = -1; xp <= 1; ++xp) {
	for (int yp = -1; yp <= 1; ++yp) {
	  if ((0 == xp) && (0 == yp)) continue; 
	  double area = histogram->GetBinContent(i+xp, j+yp);
	  if (0 == area) continue; 
	  surrounds[area]++;
	}
      }
      std::map<double, int>::iterator best = surrounds.begin(); 
      for (std::map<double, int>::iterator c = surrounds.begin(); c != surrounds.end(); ++c) {
	if ((*c).second < (*best).second) continue;
	best = c; 
      }

      histogram->SetBinContent(i, j, (*best).first); 
    }
  }

  edmhist->GetZaxis()->SetRangeUser(0, 0.00001); 
  edmhist->Draw("colz"); 
  foo->SaveAs("edms.png"); 

  edm1d.Draw();
  foo->SaveAs("edm1d.png"); 

  int colors[5];
  colors[0] = kBlue;
  colors[1] = kCyan-7;
  colors[2] = 8;
  colors[3] = 42;
  colors[4] = 46;
  TColor::SetPalette(5, colors); 

  histogram->Draw("col"); 
}


std::vector<TGraph*> MixDrawer::drawEllipse3 (DrawOptions* dis, TCanvas* foo) {
  double fitpar[MixingResult::nParams];
  double fiterr[MixingResult::nParams];
  for (int i = 0; i < MixingResult::nParams; ++i) {
    MixingResult::minuit->GetParameter(i, fitpar[i], fiterr[i]);
  }

  std::vector<TGraph*> ret(dis->numContours); 
  MixingResult::initialised = false; 

  std::cout << "Drawing " << dis->numContours << " ellipses\n";

  //  In 2 dimensions, an error ellipse nominally
  //  providing 68.27% coverage [equivalent to
  //  1 sigma in 1 dimension] has chisq=2.30
  //MixingResult::minuit->Migrad(); 
  if(dis->numContours>4){
    ret[4] = getEllipse(fiveSigma);
    if(ret[4]) ret[4]->SetFillColor(46);
  }

  //char strbuf[100];
  //int dummy = 0; 

  if(dis->numContours>3){
    ret[3] = getEllipse(fourSigma); //4 sigma contour
    if(ret[3]) ret[3]->SetFillColor(42);
  }

  /*
  for (int i = 0; i < MixingResult::nParams; ++i) {
    sprintf(strbuf, "SET PAR %i %f", i, fitpar[i]); 
    MixingResult::minuit->mncomd(strbuf, dummy);
    assert(0 == dummy); 
  }
  */ 

  if(dis->numContours>2){
    ret[2] = getEllipse(threeSigma);//3 sigma in 2d
    if (ret[2]) ret[2]->SetFillColor(8);
  }
  ret[1] = getEllipse(twoSigma);
  if (ret[1]) ret[1]->SetFillColor(kCyan-7);
  ret[0] = getEllipse(oneSigma);
  if (ret[0]) ret[0]->SetFillColor(kBlue);

  //ret[1]->RemovePoint(11);

  for(unsigned int ell_num=(ret.size())-1; ell_num >= 1; --ell_num) {
    std::cout<<"drawing contour"<<ell_num<<std::endl;
    if (ret[ell_num]) ret[ell_num]->Draw("if9");
  }
  ret[0]->Draw("if9");
  MixingResult::minuit->SetErrorDef(1.00);
  return ret; 
}

bool lineIntersectionHelper (double xcoef, double ycoef, 
			     double con, double val, 
			     bool isXline, 
			     double min, double max, 
			     std::pair<double, double>& ret) {
  // Returns intersection of line xcoef*x + ycoef*y = con with line (x or y) = val. 

  if (isXline) {
    ret.second = val;
    ret.first = (con - ycoef*val) / xcoef; 
    if (ret.first < min) return false;
    if (ret.first > max) return false;
  }
  else {
    ret.first = val;
    ret.second = (con - xcoef*val) / ycoef; 
    if (ret.second < min) return false;
    if (ret.second > max) return false;
  }

  return true;
}

void MixDrawer::drawBand (MixingResult* dat, DrawOptions* dis, TCanvas* foo) {
  double xCoefficient = 0;
  double yCoefficient = 0;
  switch (dis->drawType) {
  case DrawOptions::Xbar: 
    xCoefficient = 1;
    yCoefficient = 0;
    break;
  case DrawOptions::Ybar:
    xCoefficient = 0;
    yCoefficient = 1;
    break;
  case DrawOptions::Band:
    xCoefficient = 1/dat->coefficient;
    yCoefficient = 1; 
    break;
  default:
    assert(false);
    break;
  }

  // Line consists of points satisfying x*xcoef + y*ycoef = result.
  // Adding errors, we get x*xcoef + y*ycoef = result +- error.
  std::pair<double, double> upLeft(0, 0);
  std::pair<double, double> upRigt(0, 0);
  std::pair<double, double> dnLeft(0, 0);
  std::pair<double, double> dnRigt(0, 0);

  // First try left y-axis. 
  bool isGood = lineIntersectionHelper(xCoefficient, yCoefficient, 
					dat->measurement + dat->error, 
					xmin, 
					false, 
					ymin, ymax,
					upLeft);
  if (isGood) {
    // Other intersection can be any of the other three lines. 
    // First try right y-axis. 
    isGood = lineIntersectionHelper(xCoefficient, yCoefficient, 
				    dat->measurement + dat->error, 
				    xmax, 
				    false, 
				    ymin, ymax,
				    upRigt);
    // Next try upper x-axis.
    if (!isGood) isGood = lineIntersectionHelper(xCoefficient, yCoefficient, 
						 dat->measurement + dat->error, 
						 ymax, 
						 true, 
						 xmin, xmax,
						 upRigt);
    // Finally lower y-axis.
    if (!isGood) isGood = lineIntersectionHelper(xCoefficient, yCoefficient, 
						 dat->measurement + dat->error, 
						 ymin, 
						 true, 
						 xmin, xmax,
						 upRigt);
    // That had better work! 
    assert(isGood); 
  }
  else {
    // Ok, try upper y-axis. 
    isGood = lineIntersectionHelper(xCoefficient, yCoefficient, 
				    dat->measurement + dat->error, 
				    ymax, 
				    true, 
				    xmin, xmax,
				    upLeft);
    if (isGood) {
      // Other point may be lower y-axis, or right x-axis. 
      isGood = lineIntersectionHelper(xCoefficient, yCoefficient, 
				      dat->measurement + dat->error, 
				      ymin, 
				      true, 
				      xmin, xmax,
				      upRigt);
      if (!isGood) isGood = lineIntersectionHelper(xCoefficient, yCoefficient, 
						   dat->measurement + dat->error, 
						   xmax, 
						   false, 
						   ymin, ymax,
						   upRigt);
      assert(isGood); 
    }
    else {
      // No? Then the remaining possibility is lower y-axis and right x-axis. 
      isGood = lineIntersectionHelper(xCoefficient, yCoefficient, 
				      dat->measurement + dat->error, 
				      ymin, 
				      true, 
				      xmin, xmax,
				      upLeft);
      assert(isGood); 
      if (!isGood) isGood = lineIntersectionHelper(xCoefficient, yCoefficient, 
						   dat->measurement + dat->error, 
						   xmax, 
						   false, 
						   ymin, ymax,
						   upRigt);
      assert(isGood); 
    }
  }
  
  // Repeat for other line. 
  isGood = lineIntersectionHelper(xCoefficient, yCoefficient, 
				  dat->measurement - dat->error, 
				  xmin, 
				  false, 
				  ymin, ymax,
				  dnLeft);
  if (isGood) {
    // Other intersection can be any of the other three lines. 
    // First try right y-axis. 
    isGood = lineIntersectionHelper(xCoefficient, yCoefficient, 
				    dat->measurement - dat->error, 
				    xmax, 
				    false, 
				    ymin, ymax,
				    dnRigt);
    // Next try upper x-axis.
    if (!isGood) isGood = lineIntersectionHelper(xCoefficient, yCoefficient, 
						 dat->measurement - dat->error, 
						 ymax, 
						 true, 
						 xmin, xmax,
						 dnRigt);
    // Finally lower y-axis.
    if (!isGood) isGood = lineIntersectionHelper(xCoefficient, yCoefficient, 
						 dat->measurement - dat->error, 
						 ymin, 
						 true, 
						 xmin, xmax,
						 dnRigt);
    // That had better work! 
    assert(isGood); 
  }
  else {
    // Ok, try upper y-axis. 
    isGood = lineIntersectionHelper(xCoefficient, yCoefficient, 
				    dat->measurement - dat->error, 
				    ymax, 
				    true, 
				    xmin, xmax,
				    dnLeft);
    if (isGood) {
      // Other point may be lower y-axis, or right x-axis. 
      isGood = lineIntersectionHelper(xCoefficient, yCoefficient, 
				      dat->measurement - dat->error, 
				      ymin, 
				      true, 
				      xmin, xmax,
				      dnRigt);
      if (!isGood) isGood = lineIntersectionHelper(xCoefficient, yCoefficient, 
						   dat->measurement - dat->error, 
						   xmax, 
						   false, 
						   ymin, ymax,
						   dnRigt);
      assert(isGood); 
    }
    else {
      // No? Then the remaining possibility is lower y-axis and right x-axis. 
      isGood = lineIntersectionHelper(xCoefficient, yCoefficient, 
				      dat->measurement - dat->error, 
				      ymin, 
				      true, 
				      xmin, xmax,
				      dnLeft);
      assert(isGood); 
      if (!isGood) isGood = lineIntersectionHelper(xCoefficient, yCoefficient, 
						   dat->measurement - dat->error, 
						   xmax, 
						   false, 
						   ymin, ymax,
						   dnRigt);
      assert(isGood); 
    }
  }

  TPolyLine* poly = new TPolyLine(5); 
  poly->SetPoint(0, upLeft.first, upLeft.second);
  poly->SetPoint(1, upRigt.first, upRigt.second);
  poly->SetPoint(2, dnRigt.first, dnRigt.second);
  poly->SetPoint(3, dnLeft.first, dnLeft.second);
  poly->SetPoint(4, upLeft.first, upLeft.second);

  std::cout << "Drawing polygon " << dat->name << " (" 
	    << upLeft.first << ", " << upLeft.second << ") ("
	    << upRigt.first << ", " << upRigt.second << ") ("
	    << dnRigt.first << ", " << dnRigt.second << ") ("
	    << dnLeft.first << ", " << dnLeft.second << ") \n"; 
   
  poly->SetFillColor(dis->colour); 
  poly->SetLineStyle(11);
  poly->Draw("f"); 

  if ((fabs(upLeft.first - dnLeft.first) > 0.00001) &&
      (fabs(upLeft.second - dnLeft.second) > 0.00001)) {
    TPolyLine* upTriangle = new TPolyLine(4);
    upTriangle->SetPoint(0, upLeft.first, upLeft.second);
    upTriangle->SetPoint(1, dnLeft.first, dnLeft.second);
    if (upLeft.second > dnLeft.second) 
      upTriangle->SetPoint(2, dnLeft.first, upLeft.second);
    else 
      upTriangle->SetPoint(2, upLeft.first, dnLeft.second);
    upTriangle->SetPoint(3, upLeft.first, upLeft.second);
    upTriangle->SetFillColor(dis->colour); 
    upTriangle->SetLineStyle(11);
    upTriangle->Draw("f"); 
  }

  if ((fabs(upRigt.first - dnRigt.first) > 0.00001) &&
      (fabs(upRigt.second - dnRigt.second) > 0.00001)) {
    TPolyLine* dnTriangle = new TPolyLine(4);
    dnTriangle->SetPoint(0, upRigt.first, upRigt.second);
    dnTriangle->SetPoint(1, dnRigt.first, dnRigt.second);
    if (upRigt.second < dnRigt.second) 
      dnTriangle->SetPoint(2, dnRigt.first, upRigt.second);
    else 
      dnTriangle->SetPoint(2, upRigt.first, dnRigt.second);
    dnTriangle->SetPoint(3, upRigt.first, upRigt.second);
    dnTriangle->SetFillColor(dis->colour); 
    dnTriangle->SetLineStyle(11);
    dnTriangle->Draw("f"); 
  }
}

void MixDrawer::drawYbar (MixingResult* dat, DrawOptions* dis, TCanvas* foo) {
  for (int i = dis->numContours; i > 0; --i) {
    std::cout << "Drawing contour " << i << " of " << dat->name << std::endl; 
    double maxy = std::min(dat->measurement + i*dat->error, ymax);
    double miny = std::max(dat->measurement - i*dat->error, ymin);
    
    TBox* ycpbox = new TBox(xmin, miny, xmax, maxy); 
    ycpbox->SetFillColor(dis->colour - 1 + 2*i);
    ycpbox->SetLineStyle(11);
    ycpbox->Draw();
  }
}

void MixDrawer::addResult (MixingResult* dat, DrawOptions* dis) {
  drawmap[dat] = dis;  
}

DrawOptions::DrawingTypes DrawOptions::getDrawType (std::string name) {
  if (name == "ellipse")    return Ellipse; 
  if (name == "ybar")       return Ybar;
  if (name == "xbar")       return Xbar;
  if (name == "band")       return Band; 
  if (name == "annulus")    return Annulus; 
  if (name == "annukidney") return AnnuKidney; 
  if (name == "special")    return Special; 
  return Unknown; 
}

void DrawOptions::setColor (std::string base, int mod) {
  if (base == "red") colour = kRed;
  else if (base == "blue") colour = kBlue;
  else if (base == "violet") colour = kViolet;
  else if (base == "cyan") colour = kCyan;
  else if (base == "green") colour = kGreen;
  else if (base == "yellow") colour = kYellow;
  else if (base == "black") colour = kBlack;
  else if (base == "orange") colour = kOrange;
  else {
    std::cout << "Unknown base colour " 
	      << base 
	      << ", using black."
	      << std::endl;
    colour = kBlack;
  }
  colour += mod; 
}



/*
void MixDrawer::findPoint_adam (TGraph* ret, int idx, double angle, double errorDef, int par1, int par2) {
  //we're finding a point at a certain angle here, but we don't know where R is, so we have to iterate to find it.
  double fitpar[MixingResult::nParams];//first get the parameters
  double fiterr[MixingResult::nParams];//next get the error
  double chisq;//initialize the chi2
  int npars = MixingResult::nParams;//get the number of parameters
  
  for (int i = 0; i < MixingResult::nParams; ++i) {
    MixingResult::minuit->GetParameter(i, fitpar[i], fiterr[i]);//get all the results
    //std::cout<<"i = "<<i<<", fitpar["<<i<<"] = "<<fitpar[i]<<", fiterr["<<i<<"] = "<<fiterr[i]<<std::endl;
  }
  
  MixChisqFcn(npars, 0, chisq, fitpar, 4);//get the chi2. Derivative (0) isn't used in the FCN, so it doesn't matter. What is flag=4?

  std::cout<<"FCN called. chi2 set to "<<chisq<<std::endl<<std::endl<<std::endl;
  std::cout<<"initializing points"<<std::endl;
  double initX = fitpar[par1];//this is the initial x position
  double initY = fitpar[par2]; //this is the initial y position

  double xincrement = cos(angle);//this is the cos of the angle we want
  double yincrement = sin(angle); //sin of the angle we want
  double minimum = 0;//minimum chi2
  double chiAtMinimum = chisq; //chi2 at the minimum
  double chiAtSolution = chisq; //chi2 at the current solution
  //  std::cout<<"initX = "<<initX<<", initY = "<<initY<<", xincrement = "<<xincrement<<", yincrement = "<<yincrement<<std::endl;//ad 8/22/13
  double maximum = TMath::Sqrt(initX*initX+initY*initY)+0.001; //initialize the radius plus epsilon
  double chiAtMaximum = chisq; 
  std::cout<<"initX = "<<initX<<", initY = "<<initY<<endl;
  std::cout<<"cos(angle) = "<<xincrement<<std::endl;
  std::cout<<"sin(angle) = "<<xincrement<<std::endl;
  std::cout<<""
    
  while (chiAtMaximum-chiAtSolution< errorDef) {//break out of the loop as soon as we exceed the delta chi2. The error def is then the delta chi2
    maximum *= 2; 
    fitpar[0] = initX + maximum*xincrement;//(r+dr) cos(angle)
    fitpar[1] = initY + maximum*yincrement;//(r+dr) sin(angle)
    MixChisqFcn(npars, 0, chisq, fitpar, 4);//check the fit
    
    chiAtMaximum = chisq; 
  }
  //we have the maximum
  

  static const double tolerance = 0.001; //what is this?
  for (int i = 0; i < 1000; ++i) {//iterate 1000 times
    double currDist = minimum + (maximum - minimum)*0.5; //take mean of max and min
    fitpar[0] = initX + currDist*xincrement;//x= x_0 +currDist*cos(angle)
    fitpar[1] = initY + currDist*yincrement;//y= y_0+ currDist*sin(angle)
    MixChisqFcn(npars, 0, chisq, fitpar, 4);
    if (chisq > chiAtSolution + errorDef) {
      maximum = currDist;
      chiAtMaximum = chisq;
    }
    else {
      minimum = currDist;
      chiAtMinimum = chisq;
    }
  
    //if (0 == idx) std::cout << "Scan " << currDist << " " << (chisq - chiAtSolution) << " " << chiAtSolution << std::endl; 
    assert(chiAtMaximum > chiAtMinimum);
    if (chiAtMaximum - chiAtMinimum < tolerance) break; 
  }

  double dist = minimum + 0.5*(maximum-minimum); 
  //std::cout << "Found point " << initX + dist*xincrement << ", " << initY + dist*yincrement << std::endl;
  //whole point here is to find from the point on the contour where the FCN is a minimum. I'm going to steal the code from mncont
  
  ret->SetPoint(idx, initX + dist*xincrement, initY + dist*yincrement); 
}
*/
