#include "CMS_lumi.h"
#include <iostream>

void CMS_lumi(TPad *pad, int iPeriod, int iPosX, TString extraExtraText)
{

  bool outOfFrame    = false;
  if( iPosX/10==0 ) 
  {
    outOfFrame = true;
  }
  int alignY_=3;
  int alignX_=2;
  if( iPosX/10==0 ) alignX_=1;
  if( iPosX==0    ) alignX_=1;
  if( iPosX==0    ) alignY_=1;
  if( iPosX/10==1 ) alignX_=1;
  if( iPosX/10==2 ) alignX_=2;
  if( iPosX/10==3 ) alignX_=3;
  //if( iPosX == 0  ) relPosX = 0.12;
  int align_ = 10*alignX_ + alignY_;

  float H = pad->GetWh();
  float W = pad->GetWw();
  float l = pad->GetLeftMargin();
  float t = pad->GetTopMargin();
  float r = pad->GetRightMargin();
  float b = pad->GetBottomMargin();
  //  float e = 0.025;

  pad->cd();

  // Construct the luminosity text
  TString lumiText;
  if( iPeriod==1 )
  {
    lumiText += lumi_7TeV;
    lumiText += " (7 TeV)";
  }
  else if ( iPeriod==2 )
  {
    lumiText += lumi_8TeV;
    lumiText += " (8 TeV)";
  }
  else if( iPeriod==3 ) 
  {
    lumiText = lumi_8TeV; 
    lumiText += " (8 TeV)";
    lumiText += " + ";
    lumiText += lumi_7TeV;
    lumiText += " (7 TeV)";
  }
  else if ( iPeriod==4 )
  {
    lumiText += lumi_13TeV;
    lumiText += " (13 TeV)";
  }
  else if ( iPeriod==7 )
  { 
    if( outOfFrame ) lumiText += "#scale[0.85]{";
    lumiText += lumi_13TeV; 
    lumiText += " (13 TeV)";
    lumiText += " + ";
    lumiText += lumi_8TeV; 
    lumiText += " (8 TeV)";
    lumiText += " + ";
    lumiText += lumi_7TeV;
    lumiText += " (7 TeV)";
    if( outOfFrame) lumiText += "}";
  }
  else if ( iPeriod==12 )
  {
    lumiText += "8 TeV";
  }
  else if (iPeriod == 2022)
  {
    lumiText += lumi_13p6TeV;
    lumiText += " (13.6 TeV)";
  }
  else if (iPeriod == 0)
  {
    lumiText += lumi_sqrtS;
  }

  std::cout << "[CMS_lumi] lumiText = " << lumiText << std::endl;

  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);


  // Set fixed positions and sizes (independent of margins)
  float cmsX = 0.162;   // Distance from left edge
  float cmsY = 0.905;   // Distance from bottom edge
  float lumiX = 0.97;  // Right-aligned near top-right
  float lumiY = 0.905;

  float cmsFontSize = 0.05;   // Static text size
  float lumiFontSize = 0.035;
  float extraFontSize = 0.035;
  float extraXOffset = 0.125;  // Downshift for "Preliminary"

  // Draw Lumi text (right-aligned)
  latex.SetTextFont(42);
  latex.SetTextAlign(31); // Right aligned
  latex.SetTextSize(lumiFontSize);
  latex.DrawLatex(lumiX, lumiY, lumiText);

  // Draw "CMS" on the left
  latex.SetTextFont(cmsTextFont);
  latex.SetTextAlign(11); // Left aligned
  latex.SetTextSize(cmsFontSize);
  latex.DrawLatex(cmsX, cmsY, cmsText);

  // Draw "Preliminary" or other extra text
  if (writeExtraText) {
    latex.SetTextFont(extraTextFont);
    //latex.SetTextAlign(11); // Left aligned
    latex.SetTextSize(extraFontSize);
    latex.DrawLatex(cmsX + extraXOffset, cmsY, extraText + " " + extraExtraText);
  }
  return;
}
