{
  TStyle *sysStyle  = new TStyle("sysStyle","Sys Style");

  // from ROOT plain style
  sysStyle->SetCanvasBorderMode(0);
  sysStyle->SetPadBorderMode(0);
  sysStyle->SetPadColor(0);
  sysStyle->SetCanvasColor(0);
  sysStyle->SetStatColor(0);
  sysStyle->SetLegendFillColor(0);
  sysStyle->SetPalette(1);
  sysStyle->SetStatBorderSize(1);
  sysStyle->SetStatFont(132);
  sysStyle->SetLegendBorderSize(1);
  sysStyle->SetLegendFont(132);

  sysStyle->SetCanvasDefH(600);
  sysStyle->SetCanvasDefW(750);
  sysStyle->SetPadBottomMargin(0.11);
  sysStyle->SetPadTopMargin(0.09);
  sysStyle->SetPadLeftMargin(0.14);
  sysStyle->SetPadRightMargin(0.22);
  sysStyle->SetOptFit(1111);
  sysStyle->SetPadTickX(1);
  sysStyle->SetPadTickY(1);

  sysStyle->SetFrameBorderMode(0);

  sysStyle->SetMarkerStyle(20);
  sysStyle->SetMarkerSize(0.6);
  // Fonts
  int font = 132;
  sysStyle->SetLabelFont(font, "xyz");
  sysStyle->SetTitleFont(font, "xyz");
  sysStyle->SetTitleFont(font);
  sysStyle->SetTitleOffset(1.4, "y");
  sysStyle->SetLegendFont(font);

  gROOT->SetStyle("sysStyle");

  // Color pallette
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

  gROOT->ForceStyle();
}

