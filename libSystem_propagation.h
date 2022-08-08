#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include <libEDM_library.h>
#include <libEDM_propagation.h>
#include <libEDM_random.h>

// forward declarations
class Interferer;
class NodeB;
class Position;
class Station;
class WrappedStation;

class PathgainType {
public:
	      double    pathgain;
    const Position *wrappedPosition;
	      metres    distance;
	      double    wrappedAntennaGain;
	      double    stationAntennaGain;
	      double    propagationGain;
          double    shadowFading;
          bool      isMaxPathgainLimited;         // Indicates that the maxPathgain limit has been applied to the pathgain
          bool      isLOSLimited;                 // Indicates that the LOS limit has been applied to the propagation gain
          bool      isLOSLimitedAfterShadowing;   // Indicates that the LOS limit was applied after shadowing

    PathgainType () : isMaxPathgainLimited(false), isLOSLimited(false), isLOSLimitedAfterShadowing(false), distance(metres(0.0)) {}
};

class FadeParams
{
public:
    const dB shadowFadingStdDev;
    const dB buildingFadingMedian;
    const dB buildingFadingStdDev;

    FadeParams( const dB &shadowFadingStdDev, const dB &buildingFadingMedian, const dB &buildingFadingStdDev ) : shadowFadingStdDev(shadowFadingStdDev), buildingFadingMedian(buildingFadingMedian), buildingFadingStdDev(buildingFadingStdDev)
    { assert( buildingFadingMedian <= 0.0 );}
};

class NetworkPropagation : public Propagation {
public:
    enum Model {LOS, PowerLaw, Hata, ExtendedHata, NoModel};

    class NetworkPropagationInit {
    public:
        Model       nodeBToStationModel;
		Model       interfererToStationModel;
        MHz         frequency_MHz;
		dB          defaultkdB;
		double      defaultGamma;
		MHz         systemBandwidth_MHz;
		dB          minCouplingLoss;
		ClutterType defaultClutterType;
		double      orthfactor;
		metres      minDistance;
        double      minBuildingPenetrationGain;

        NetworkPropagationInit() : nodeBToStationModel        (NoModel),
                                   interfererToStationModel   (NoModel),
                                   frequency_MHz              (MHz(0.0)),
                                   defaultkdB                 (dB(0.0)),
                                   defaultGamma               (0.0),
                                   systemBandwidth_MHz        (MHz(0.0)),
                                   minCouplingLoss            (dB(0.0)),
                                   defaultClutterType         (Urban),
                                   orthfactor                 (1.0),
                                   minDistance                (metres(0.0)) {}
    };

public:
    const double     maxCouplingGain;
	const double     orthfactor;
	const MHz        systemBandwidth;
    const FadeParams fadeParams;

	const Model  interfererToStationModel, nodeBToStationModel;

private:
    ShadowFading _shadowFading;
    ShadowFading _buildingFading;

    const double minBuildingFadeGain, maxBuildingFadeGain;

public:
    NetworkPropagation (const NetworkPropagationInit &propagationParams, const FadeParams &fadeParams)
		: Propagation                (propagationParams.frequency_MHz, propagationParams.defaultkdB, propagationParams.defaultGamma, propagationParams.defaultClutterType, propagationParams.minDistance),
		  nodeBToStationModel        (propagationParams.nodeBToStationModel),
          interfererToStationModel   (propagationParams.interfererToStationModel),
		  systemBandwidth            (propagationParams.systemBandwidth_MHz),
		  maxCouplingGain            (dB2linear(-propagationParams.minCouplingLoss)),
		  orthfactor                 (propagationParams.orthfactor),
          fadeParams                 (fadeParams),
          minBuildingFadeGain        (propagationParams.minBuildingPenetrationGain),
          maxBuildingFadeGain        (1.0),
          _shadowFading              (fadeParams.shadowFadingStdDev),
          _buildingFading            (fadeParams.buildingFadingStdDev, fadeParams.buildingFadingMedian)
    {assert( minBuildingFadeGain <= maxBuildingFadeGain );};

    ~NetworkPropagation () {}

    PathgainType pathgain (const Interferer &interferer, const Station &station);
    PathgainType pathgain (const NodeB      &nodeB,      const Station &station, const bool newSite);

	PathgainType pathgain (const WrappedStation &wrappedStation, const Station &station, const Model model, const bool newSite);
};