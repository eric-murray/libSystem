#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <limits>

#include <libEDM_library.h>

#include <libSystem_propagation.h>
#include <libSystem_stations.h>

using std::cerr;
using std::cos;
using std::max;
using std::min;
using std::norm;
using std::polar;

//
// class NetworkPropagation
//

PathgainType NetworkPropagation::pathgain (const Interferer &interferer, const Station &station)
{
    return pathgain(interferer, station, interfererToStationModel, true);
}

PathgainType NetworkPropagation::pathgain (const NodeB &nodeB, const Station &station, const bool newSite)
{
    return pathgain(nodeB, station, nodeBToStationModel, newSite);
}

PathgainType NetworkPropagation::pathgain (const WrappedStation &wrappedStation, const Station &station, const Model model, const bool newSite)
{
	// find position of interferer relative to station with highest path + antenna gain
    // Note: UE antenna gain excluded, as pointing direction unknown at this stage

	PathgainType data;

	if ( model != LOS )
	{
        double wrappedShadowFading = wrappedStation.position.shadowFading();
        double stationShadowFading = station.position.shadowFading();

        if ( newSite )
        {
            _shadowFading.reset();
            _buildingFading.reset();
        }

        double channelBuildingFading = _buildingFading * wrappedStation.position.buildingFading() * station.position.buildingFading();
        channelBuildingFading = std::max(minBuildingFadeGain, channelBuildingFading);
        channelBuildingFading = std::min(maxBuildingFadeGain, channelBuildingFading);

        data.shadowFading = wrappedShadowFading * _shadowFading * channelBuildingFading * stationShadowFading;
	}
    else
        data.shadowFading = 1.0;
    
    double highestPathgain = 0.0;
    bool   isLOSLimited    = false;
    for (size_t i=0; i<wrappedStation.positions.size(); i++)
    {
		const Position &position = wrappedStation.positions.position(i);

		// compute antenna gain from wrappedStation to station
		const radians angle              = station.angleFrom(position.point);
		const double  wrappedAntennaGain = wrappedStation.gain(angle);

        // compute LOS pathgain as this is always required
		const metres distance = station.distance(position.point);

		// compute pathgain between two positions
        double pathgain;
		switch ( model )
		{
		case LOS:
			pathgain = pathgain_LOS(distance);
			break;

		case PowerLaw:
			pathgain = pathgain_power_law(distance);
			break;

		case Hata:
		case ExtendedHata:
            pathgain = pathgain_ExtendedHata(distance, wrappedStation.height, station.height, defaultClutterType, isLOSLimited);
			break;

		default:
			error("Propagation model not defined");
			break;
		}

        const double combinedPathgain = wrappedAntennaGain * pathgain;
		if ( combinedPathgain > highestPathgain )
		{
			highestPathgain         = combinedPathgain;
            data.isLOSLimited       = isLOSLimited;
			data.wrappedPosition    = &position;
			data.distance           = distance;
			data.wrappedAntennaGain = wrappedAntennaGain;
			data.propagationGain    = pathgain;
		}
	}

    if ( data.isLOSLimited )
        // Restrict shadow fading to building penetration fading if propagation is LOS limited
        data.shadowFading = _buildingFading;

    data.pathgain = data.propagationGain * data.shadowFading;

    const double losPathgain = pathgain_LOS(data.distance);
    if ( data.pathgain > losPathgain )
    {
        // shadow fading has pushed pathgain higher than LOS
        data.pathgain = losPathgain;

        if ( !data.isLOSLimited)
            data.isLOSLimitedAfterShadowing = true;
    }

    // epdate pathgain with BS antenna gain
    data.pathgain = data.wrappedAntennaGain * data.pathgain;

    // ensure pathgain is no greater than specified maxPathgain
    if ( data.pathgain > maxCouplingGain )
    {
        data.isMaxPathgainLimited = true;
        data.pathgain             = maxCouplingGain;
    }

	return data;
}