#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TASImage.h"
#include <string>


// Global variables

TString cmsText     = "CMS";
float cmsTextFont   = 61;  // default is helvetic-bold
float cmsX = 0.2;
float cmsY = 0.91;
float lumiY = 0.91;
bool writeExtraText = false;
TString extraText   = "Simulation Preliminary";
float extraTextFont = 52;  // default is helvetica-italics

float cmsFontSize = 0.055;   // Static text size
float lumiFontSize = 0.045;
float extraFontSize = 0.045;
float extraXOffset = 0.125;  
// Absolute font sizes (not relative to top margin)
float lumiTextSize     = 0.045;  // ~4.5% of canvas height
float lumiTextOffset   = 0.2;
float cmsTextSize      = 0.055;  // ~5.5% of canvas height
float cmsTextOffset    = 0.1;

float relPosX    = 0.045;
float relPosY    = 0.035;
float relExtraDY = 1.2;

// ratio of "CMS" and extra text size
float extraOverCmsTextSize  = 0.76;

TString lumi_13p6TeV = "34.74 fb^{-1}";
TString lumi_13TeV = "20.1 fb^{-1}";
TString lumi_8TeV  = "19.7 fb^{-1}";
TString lumi_7TeV  = "5.1 fb^{-1}";
TString lumi_sqrtS = "";

bool drawLogo      = false;

void CMS_lumi( TPad* pad, int iPeriod=3, int iPosX=10, TString eet="" );
