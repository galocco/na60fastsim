#include "KMCPixelPlane.h"
#include <iostream>

#ifdef _NOVTPLANEGAPS_
const float KMCPixelPlane::DeadXShort = 0;// 0.15;  // short end cap
const float KMCPixelPlane::DeadXLong = 0;//0.45;   // long end cap   (readout)
const float KMCPixelPlane::DeadYBottom = 0;//0.0525; // bottom dead zone
const float KMCPixelPlane::DeadYTop = 0;//0.0525; // top dead zone
#else
const float KMCPixelPlane::DeadXShort = 0.15;  // short end cap
const float KMCPixelPlane::DeadXLong = 0.45;   // long end cap   (readout)
const float KMCPixelPlane::DeadYBottom = 0.0525; // bottom dead zone
const float KMCPixelPlane::DeadYTop = 0.0525; // top dead zone
#endif



const float KMCPixelPlane::DeadTopBotHalves = 0.012; // dead space between top and bottom halves
const float KMCPixelPlane::DeadXTile = 0.002; // dead space between tiles in X (apart from DeadXDataBackbone after each 3 tiles)
const float KMCPixelPlane::DeadXDataBackbone = 0.006; // dead space between segments
const float KMCPixelPlane::XSizeTot = 12.9996 + DeadXShort +  DeadXLong;// 13.5996;
const float KMCPixelPlane::YSizeTot = 13.5898 + DeadYBottom + DeadYTop; // 13.6948; // readout side
const int   KMCPixelPlane::NXTiles = 36;
const int   KMCPixelPlane::NXSegments = 12; // group of 3 tiles
const int   KMCPixelPlane::NYSensors = 7;
const float KMCPixelPlane::DYSens = YSizeTot / NYSensors;
const float KMCPixelPlane::DXSegment = (XSizeTot - DeadXShort - DeadXLong) / NXSegments;
const float KMCPixelPlane::DXTile = (DXSegment - DeadXDataBackbone) / 3;


KMCPixelPlane::KMCPixelPlane(const char *name, float zpos, float thickness,
			     float radL, float density, NaMaterial* mat, // substrate
			     float radLS, float densityS, NaMaterial* matS, // sensor
			     float offsX_, float offsY_, float interChipGap_, // offset of the inner corner of the 1st quadrant chip
			     float xRes, Float_t yRes, Float_t eff
			     ) : KMCPolyLayer(name)
{
  // subtrate settings

  std::cout<<"KMCPixelPlane "<<name<<"\n";
  std::cout<<"thickness*density/radL: "<<thickness*density/radL<<std::endl;
  std::cout<<"radL: "<<radL<<std::endl;
  std::cout<<"thickness: "<<thickness<<std::endl;
  std::cout<<"density: "<<density<<std::endl;
  SetZ(zpos);
  SetThickness(thickness);
  SetX2X0( radL>0 ? thickness*density/radL : 0);
  SetXTimesRho(thickness*density);
  SetMaterial(mat);
  //
  // sensors resolution and efficiency
  SetXRes(xRes);
  SetYRes(yRes);
  SetLayerEff(eff);
  SetDead(false);
  //
  sensMat = matS;
  
  // define 2-chip contour
  offsX = offsX_;
  offsY = offsY_;
  interChipGap = interChipGap_;
  
  float contourX[8] = {}, contourY[8] = {};
  contourX[0] = XSizeTot + interChipGap/2 + offsX;
  contourY[0] = YSizeTot + offsY + interChipGap/2;
  contourX[1] = contourX[0];
  contourY[1] = contourY[0] - YSizeTot;
  contourX[2] = contourX[1] - XSizeTot;
  contourY[2] = contourY[1];
  contourX[3] = contourX[2];
  contourY[3] = contourY[2] - offsY - offsY;

  contourX[4] = contourX[3] - (XSizeTot + interChipGap);
  contourY[4] = contourY[3];
  contourX[5] = contourX[4];
  contourY[5] = contourY[4] + YSizeTot;
  contourX[6] = contourX[5] + XSizeTot;
  contourY[6] = contourY[5];
  contourX[7] = contourX[6];
  contourY[7] = contourY[6] + offsY + offsY;
  //  contourX[8] = contourX[0];   contourY[8] = contourY[0];
   
  addPolygon(8, contourX, contourY, (radLS>0 && densityS>0) ? thickness/(radLS/densityS) : 0, thickness*densityS);
  setNSectorsPhiStart(2, 0.0); // define full set of chips as rotated copies of the 1st one
}

int KMCPixelPlane::isInAcc(float x, float y, float) const
{
  int res = 0;
  int pol = getPolygonID(x,y);
  if (pol<0) return res;
  res = -1; // in physicall acceptance
  bool rot = false;
  if (pol!=0) { // rotate to base poligon
    x = -x;
    y = -y;
    rot = true;
  }
  int module = 0;
  x -= offsX;
  if (x < 0.) {
    module = 1;
    y += offsY - interChipGap/2;
    x += XSizeTot + interChipGap/2; // relate to chip local left/bottom corner
  } else {    
    y -= offsY + interChipGap/2;
    y = YSizeTot - y;
    x = (XSizeTot + interChipGap/2) - x; // relate to chip local left/bottom corner
  }
  if (x < DeadXLong || x > XSizeTot - DeadXShort || y < DeadYBottom || y > YSizeTot - DeadYTop) return res;
  x -= DeadXLong;
  //  y -= DeadYBottom;
  int sens = int(y/DYSens);
  y -= sens*DYSens;
  if (y < DeadYBottom || y > DYSens - DeadYTop) return res;
  int segm = int(x/DXSegment);
  x -= segm*DXSegment;
  if (x > DXSegment - DeadXDataBackbone) return res;
  int tile = int(x/DXTile);
  x -= tile*DXTile;
  if (x<DeadXTile) return res;

  res = 1;
  
  return res;
}
  
