#ifndef _LIB_H		
#define _LIB_H


#define _CRT_SECURE_NO_WARNINGS


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <time.h>



#include <fstream>
#include <iostream>
using namespace std;

#include <opencv2/opencv.hpp>
#include <opencv2/nonfree/nonfree.hpp>
#include <opencv2/nonfree/features2d.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <cv.h>
#include <cxcore.h>
#include <highgui.h>

using namespace cv;



#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TProfile.h>



#pragma comment(lib, "opencv_highgui249.lib")  
#pragma comment( lib, "opencv_core249.lib" )     
#pragma comment( lib, "opencv_imgproc249.lib" )    
#pragma comment( lib, "opencv_nonfree249.lib" )    
#pragma comment( lib, "opencv_features2d249.lib" )    
#pragma comment( lib, "opencv_calib3d249.lib" )    
#pragma comment( lib, "opencv_video249.lib" )    
#pragma comment( lib, "opencv_legacy249.lib" )    

#pragma comment(lib, "complexDict.lib" ) 
#pragma comment(lib, "dequeDict.lib" ) 
#pragma comment(lib, "freetype.lib" ) 
#pragma comment(lib, "libAfterImage.lib" ) 
#pragma comment(lib, "libASImage.lib" ) 
#pragma comment(lib, "libASImageGui.lib" ) 
#pragma comment(lib, "libCint.lib" ) 
#pragma comment(lib, "libCintex.lib" ) 
#pragma comment(lib, "libCore.lib" ) 
#pragma comment(lib, "libEG.lib" ) 
#pragma comment(lib, "libEGPythia8.lib" ) 
#pragma comment(lib, "libEve.lib" ) 
#pragma comment(lib, "libFitPanel.lib" ) 
#pragma comment(lib, "libFoam.lib" ) 
#pragma comment(lib, "libFTGL.lib" ) 
#pragma comment(lib, "libFumili.lib" ) 
#pragma comment(lib, "libGdml.lib" ) 
#pragma comment(lib, "libGed.lib" ) 
#pragma comment(lib, "libGenetic.lib" ) 
#pragma comment(lib, "libGenVector.lib" ) 
#pragma comment(lib, "libGeom.lib" ) 
#pragma comment(lib, "libGeomBuilder.lib" ) 
#pragma comment(lib, "libGeomPainter.lib" ) 
#pragma comment(lib, "libGLEW.lib" ) 
#pragma comment(lib, "libGpad.lib" ) 
#pragma comment(lib, "libGraf.lib" ) 
#pragma comment(lib, "libGraf3d.lib" ) 
#pragma comment(lib, "libGui.lib" ) 
#pragma comment(lib, "libGuiBld.lib" ) 
#pragma comment(lib, "libGuiHtml.lib" ) 
#pragma comment(lib, "libGviz3d.lib" ) 
#pragma comment(lib, "libHist.lib" ) 
#pragma comment(lib, "libHistPainter.lib" ) 
#pragma comment(lib, "libHtml.lib" ) 
#pragma comment(lib, "libMathCore.lib" ) 
#pragma comment(lib, "libMathMore.lib" ) 
#pragma comment(lib, "libMatrix.lib" ) 
#pragma comment(lib, "libMinuit.lib" ) 
#pragma comment(lib, "libMinuit2.lib" ) 
#pragma comment(lib, "libMLP.lib" ) 
#pragma comment(lib, "libNet.lib" ) 
#pragma comment(lib, "libPhysics.lib" ) 
#pragma comment(lib, "libPostscript.lib" ) 
#pragma comment(lib, "libProof.lib" ) 
#pragma comment(lib, "libProofDraw.lib" ) 
#pragma comment(lib, "libProofPlayer.lib" ) 
#pragma comment(lib, "libPyROOT.lib" ) 
#pragma comment(lib, "libQuadp.lib" ) 
#pragma comment(lib, "libRecorder.lib" ) 
#pragma comment(lib, "libReflex.lib" ) 
#pragma comment(lib, "libReflexDict.lib" ) 
#pragma comment(lib, "libRGL.lib" ) 
#pragma comment(lib, "libRHTTP.lib" ) 
#pragma comment(lib, "libRint.lib" ) 
#pragma comment(lib, "libRIO.lib" ) 
#pragma comment(lib, "libRODBC.lib" ) 
#pragma comment(lib, "libRooFit.lib" ) 
#pragma comment(lib, "libRooFitCore.lib" ) 
#pragma comment(lib, "libRooStats.lib" ) 
#pragma comment(lib, "libRootAuth.lib" ) 
#pragma comment(lib, "libSessionViewer.lib" ) 
#pragma comment(lib, "libSmatrix.lib" ) 
#pragma comment(lib, "libSpectrum.lib" ) 
#pragma comment(lib, "libSpectrumPainter.lib" ) 
#pragma comment(lib, "libSPlot.lib" ) 
#pragma comment(lib, "libSQLIO.lib" ) 
#pragma comment(lib, "libTable.lib" ) 
#pragma comment(lib, "libThread.lib" ) 
#pragma comment(lib, "libTMVA.lib" ) 
#pragma comment(lib, "libTree.lib" ) 
#pragma comment(lib, "libTreePlayer.lib" ) 
#pragma comment(lib, "libTreeViewer.lib" ) 
#pragma comment(lib, "libUnuran.lib" ) 
#pragma comment(lib, "libVMC.lib" ) 
#pragma comment(lib, "libWin32gdk.lib" ) 
#pragma comment(lib, "libXMLIO.lib" ) 
#pragma comment(lib, "listDict.lib" ) 
#pragma comment(lib, "map2Dict.lib" ) 
#pragma comment(lib, "mapDict.lib" ) 
#pragma comment(lib, "mathtext.lib" ) 
#pragma comment(lib, "multimap2Dict.lib" ) 
#pragma comment(lib, "multimapDict.lib" ) 
#pragma comment(lib, "multisetDict.lib" ) 
#pragma comment(lib, "setDict.lib" ) 
#pragma comment(lib, "vectorDict.lib" ) 






#define PI 3.1415926536

#include "CTime.h"

#include "spline.h"
#include "CFrame.h"
#include "CPointToPeak.h"
#include "CEdgeHit.h"
#include "CLine.h"
#include "CEdgeHitChain.h"
#include "CLinePara.h"
#include "paras.h"
#include "CProfile.h"
#include "CColorCut.h"
#include "CPicToBin.h"
#include "CBinToProfile.h"
#include "CCubicBSpline.h"
#include "CIntProfileToFloat.h"
#include "CFloatMatch.h"
#include "CPicToEdge.h"
#include "CPointToNorm.h"
#include "CPointToCircle.h"
#include "CNormToPeak.h"
#include "CCircleToPeak.h"
#include "CEdgeToFloat.h"
#include "CEdgeToCircleChain.h"
#include "CLBP.h"
#include "CSixSword.h"
#include "CBmpDiff.h"
#include "CCurveAngle.h"
#include "CPointProjection.h"
#include "CRegionMatch.h"
#include "CChainToPoly.h"
#include "CPolyTag.h"
//#include "CChainMatch.h"
#include "CTagMatch.h"
#include "CSeparateChain.h"

#include"detectline.h"

#endif	