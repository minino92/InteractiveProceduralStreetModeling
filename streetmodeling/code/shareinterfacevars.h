

//### shareinterfacevars.h
#include "afxcmn.h"
#include "afxwin.h"

typedef struct SharedInterfaceVars
{
	int tenElemType;
	bool TensorDesignOn;  

	bool MoveElemOn;
	bool EditElementOn;
	bool RemoveElemOn;
	bool BrushInterfaceOn;
	bool CombineWithOtherPartsOn;

	bool ShowSingularitiesOn;
	bool ShowRegElemOn;
	bool ShowTensorLinesOn;
	bool ShowStreetGraphOn;
	bool ShowRoadMapOn;
	bool ShowTheMapOn;
	bool ShowExtractedBoundsOn;
	bool MeanFilterOn;
	bool RegionSmoothOn;
	bool GenTenFromExtractedBoundsOn;
	bool UseAllBoundsOn;
	bool DesignGridOn;
	bool ConnectDeadEndsOn;
	bool ShowRegionBlocksOn;
	bool ShowStreetGraphEndPointsOn;
	bool EnableSketchBasedDesign;
	bool CombinePopDensityOn;
	bool GenTenFromSketchesOn;
	bool ShowSketchesOn;
	bool CloseLoopOn;
	bool ShowPopDensityMapOn;
	bool ShowIBFVOn;
	bool EditStreetNetOn;
	bool ShowIntersectsOn;
	bool ShowStreetUseNetworkOn;
	bool SelStreetRegToEditOn;
	bool ShowLineStyleStreetsOn;

	bool UseBoundsAsSketchesOn;
	bool UseBoundsAsRoadsOn;
	bool UseMajRoadsAsSketchesOn;
	bool ShowScalarFieldOn;
	bool EnableHeighfieldDesignOn;
	bool ShowMajRoadsOn;
	bool ShowInitSeedsOn;
	bool ShowSegmentGraphOn;
	bool ShowVegMapOn;

	bool ApplyAsymFldOn;
	bool RemoveMajDeadEndsOn;

	bool ShowMajRoadGoogleStyleOn;

	bool AntiAliasingOn;

	bool BrushLikeRegSelOn;

	/*  in the major road setting dialog  */
	bool AllowMajRoadCrossRiverOn;
	bool AllowMajRoadFollowBoundaryOn;
	bool ShowMajRoadNetworkOn;
	bool AllowCrossSingularitiesOn;
	bool JobardMethodOn;

	/*  in the minor road setting dialog  */
	bool AllowMinCloseToMajOn;
	bool RemDeadEndsTraceOn;
	bool UseNewMedRemDeadEndOn;

	bool rdScalarSingularElem;
	bool rdSketchMajorHighways;

	/*  in the save project setting dialog  */
	unsigned char rdSaveProjDesignElem;
	unsigned char rdSaveProjSketches;
	unsigned char rdSaveProjBrushes;
	unsigned char rdSaveProjMajRoadNetwork;
	unsigned char rdSaveProjOtherSetting;
    unsigned char rdSaveProjStreetNetwork;

}ShareInterfaceVars;