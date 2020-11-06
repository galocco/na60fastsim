void loadLibs(){

  gROOT->ProcessLine(".L ../KMCUtils.cxx+");
  gROOT->ProcessLine(".L ../NaMaterial.cxx+");
  gROOT->ProcessLine(".L ../KMCProbeFwd.cxx+");
  gROOT->ProcessLine(".L ../KMCClusterFwd.cxx+");
  gROOT->ProcessLine(".L ../KMCLayerFwd.cxx+");
  gROOT->ProcessLine(".L ../KMCPolyLayer.cxx+"); 
  gROOT->ProcessLine(".L ../KMCMagnetBuilder.cxx+");
  gROOT->ProcessLine(".L ../KMCDetectorFwd.cxx+");
 
}
