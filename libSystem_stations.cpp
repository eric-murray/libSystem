#define _USE_MATH_DEFINES

#include <cmath>
#include <iomanip>
#include <iostream>

#include <itpp\base\math\log_exp.h>

#include <libEDM_library.h>

#include <libSystem_interference.h>
#include <libSystem_propagation.h>
#include <libSystem_networks.h>
#include <libSystem_sites.h>
#include <libSystem_stations.h>

using std::cout;
using std::endl;
using std::fixed;
using std::min;

//
// class Antenna
//

void Antenna::print(const degrees stepAngle, const degrees startAngle, const degrees stopAngle)
{
    cout << fixed << setprecision(2);
    cout << "Angle(rad) Angle(deg)   Gain(dB)" << endl;
    for (radians angle = radians(startAngle); angle <= radians(stopAngle); angle += radians(stepAngle))
    {
        cout << setw(10) << angle << " ";
        cout << setw(10) << degrees(angle) << " ";
        cout << setw(10) << linear2dB(gain(angle)) << " ";
        cout << endl;
    }
    cout << endl;
}

//
// class GaussianAntenna
//

double GaussianAntenna::gain (const radians azimuth) const
{
    // M_PI / 180.0 should be combined in multiplier
    return max(boresightGain * exp(sqr(azimuth()) * _multiplier), _minimumGain);
}


//
// class Interferer
//

const OmniAntenna *Interferer::_antenna = new OmniAntenna();

//
// class InterfererList
//

void InterfererList::print(ostream &file, const bool banner, const size_t leadingSpaces) const
{
    if ( !empty() && banner )
        file << " Id       X       Y  TxPower" << endl;

    file << fixed << setprecision(1);
    for (size_t i=0; i<this->size(); i++)
    {
        Interferer *interferer = this->at(i);

        file << setw(3) << interferer->id << " ";
        file << setw(7) << interferer->x() << " " << setw(7) << interferer->y() << " ";
        file << setw(8) << watts2dBm(interferer->txPower()) << " ";
        file << endl;
    }

    if ( !empty() && banner )
        file << endl;
}

//
// class NodeB
//

bool NodeB::check_powers() const
{
	double computedTxPower = commonChannelTxPower + sharedChannelTxPower;
    for (size_t i=0; i<_connectedUEs.size(); i++)
        computedTxPower += _connectedUEs.at(i)->transmitSignalPower();

    cout << setw( 7) << id << " ";
    cout << setw( 8) << watts2dBm(txPower(),       -999.0) << " ";
    cout << setw(16) << watts2dBm(computedTxPower, -999.0) << endl;

    if (fabs((computedTxPower/txPower()) - 1.0) > 0.0001)
        return false;
    else
        return true;
}

void NodeB::disconnect_last()
{
	UE *ue = _connectedUEs.back();
	_connectedUEs.pop_back();

	if ( (_connectedRCUEs.size() > 0) && (ue == _connectedRCUEs.back()) )
		_connectedRCUEs.pop_back();
	else
		_connectedPCUEs.pop_back();
}

Position NodeB::generate_random_position (const FadeParams &fadeParams) const
{
    return position.site->generate_random_position(sectorId, fadeParams);
}

void NodeB::reset()
{
	_connectedPCUEs.clear();
	_connectedRCUEs.clear();
	_connectedUEs.clearAndDelete();

	setTrafficChannelTxPower(0.0);
	setTxPower(commonChannelTxPower + sharedChannelTxPower);
}

//
// class NodeBList
//

const string NodeBList::nodeBBanner = "NodeBId       X       Y  Azimuth  Fading  CommonChannelTxPower  SharedChannelTxPower  TrafficChannelTxPower  TotalTxPower  MaxTrafficChannelTxPower  MaxTxPower  NumConnectedUEs";

bool NodeBList::check_powers(const bool banner) const
{
    if ( !empty() && banner )
        cout << "NodeBId  TxPower  ComputedTxPower" << endl;

    cout << setprecision(3);
    for (size_t i=0; i<this->size(); i++)
        if (!this->at(i)->check_powers())
            return false;

    if ( banner )
        cout << endl;

    return true;
}

void NodeBList::print(ostream &file, const bool banner, const size_t leadingSpaces) const
{
    print(file, banner, leadingSpaces, false);
}

void NodeBList::print(ostream &file, const bool banner, const size_t leadingSpaces, const bool printTessellated) const
{
    if ( banner )
        file << nodeBBanner.c_str() << endl;

    file << fixed << setprecision(1);
    for (size_t i=0; i<this->size(); i++)
    {
        NodeB *nodeB = this->at(i);

        size_t numPositionsToPrint = 1;
        if ( printTessellated )
            numPositionsToPrint = nodeB->positions.size();

        for (size_t positionId = 0; positionId < numPositionsToPrint; positionId++)
        {
            if ( !banner && ((i > 0) || (positionId > 0)) )
                // insert leading spaces
                file << string(leadingSpaces, ' ').c_str();

            file << setw( 7) << nodeB->id << " ";
            file << setw( 7) << nodeB->positions.at(positionId)->point.x << " ";
            file << setw( 7) << nodeB->positions.at(positionId)->point.y << " ";
            file << setw( 8) << degrees(nodeB->azimuth()) << " ";
            file << setw( 7) << linear2dB(nodeB->shadowFading(positionId))           << " ";
		    file << setw(21) << watts2dBm(nodeB->commonChannelTxPower,       -999.0) << " ";
		    file << setw(21) << watts2dBm(nodeB->sharedChannelTxPower,       -999.0) << " ";
            file << setw(22) << watts2dBm(nodeB->trafficChannelTxPower(),    -999.0) << " ";
            file << setw(13) << watts2dBm(nodeB->txPower(),                  -999.0) << " ";
            file << setw(25) << watts2dBm(nodeB->maxTrafficChannelTxPower,   -999.0) << " ";
            file << setw(11) << watts2dBm(nodeB->maxTxPower,                 -999.0) << " ";
            file << setw(16) << nodeB->numConnectedUEs();
            file << endl;
        }
    }

    if ( banner )
        file << endl;
}

void NodeBList::reset()
{
    for (size_t i=0; i<this->size(); i++)
        this->at(i)->reset();
}

//
// class UE
//

ulong UE::_numCreatedUEs = 0;

UE::~UE()
{
    delete _connectedNodeBData;
    delete _closestNodeBData;

    _connectedNodeBData = NULL;
    _closestNodeBData   = NULL;
}

bool UE::affiliate (Network &network, const AffiliationType affiliationType, const size_t debug_level, const size_t stats_level)
{
    // initialise _pathgains and _fadegains
    _pathgains.resize(network.numNodeBs());
    _fadegains.resize(network.numNodeBs(), 1.0);

    // shuffle order UE sorts through network sites
    network.shuffle();

    // compute and store pathgains
    double highestPathgain = 0.0;
    double closestDistance = numeric_limits<double>::infinity();

    vector<PathgainType> pathgainData(network.numNodeBs());
    for (size_t i=0; i<network.numSites(); i++)
    {
        Site *site   = network.site(i);
        bool newSite = true;

        for (size_t j=0; j<site->size(); j++)
        {
            NodeB *nodeB = site->nodeB(j);

            PathgainType data = network.propagation.pathgain(*nodeB,*this,newSite);
			pathgainData[nodeB->id] = data;
            _pathgains  [nodeB->id] = data.pathgain;

			_adjacentNodeBs.insert(pair<double,NodeB*>(_pathgains[nodeB->id],nodeB));

            if ( _pathgains[nodeB->id] > highestPathgain * handoverMargin )
            {
                highestPathgain = _pathgains[nodeB->id];
                _connectedNodeB = nodeB;
            }

            if ( data.distance < closestDistance )
            {
                closestDistance = data.distance;
                _closestNodeB   = nodeB;
            }

			// set newSite false so that shadow fading values are correlated at a site
            newSite = false;
        }
    }

    switch ( affiliationType )
    {
    case BestServer:
        // _connectedNodeB is correct
        break;

    case ClosestServer:
        _connectedNodeB = _closestNodeB;
        break;

    case BestServerIfClosest:
        if ( _connectedNodeB != _closestNodeB )
            return false;
        break;

    default:
        error("Invalid affiliation type specified");
    }

    if ( _connectedNodeB == NULL )
        // problem finding a nodeB to connect to
        error("Cannot find NodeB to connect to");

	if ( _connectedNodeB->numConnectedUEs() == _connectedNodeB->maxConnectedCalls )
		// NodeB is full, so reject call
		return false;

    _connectedNodeBData = new PathgainType(pathgainData[_connectedNodeB->id]);
    _closestNodeBData   = new PathgainType(pathgainData[_closestNodeB->id]);

    // set UE azimuth to point to best server
    _azimuth = this->angleTo(_connectedNodeBData->wrappedPosition->point);

	// add UE antenna gain to computed pathgains
    for (size_t nodeBId = 0; nodeBId < _pathgains.size(); nodeBId++)
    {
        // update interfering link limiting status
        /*if ( _connectedNodeB->id != nodeBId )
        {
            if ( pathgainData[nodeBId].isLOSLimited )
                _interferingLinkData.numInterferingLinksLOSLimited++;

            if ( pathgainData[nodeBId].isLOSLimitedAfterShadowing )
                _interferingLinkData.numInterferingLinksLOSLimitedAfterShadowing++;

            if ( pathgainData[nodeBId].isMaxPathgainLimited )
                _interferingLinkData.numInterferingLinksMaxPathgainLimited++;
        }*/

	    const radians angle         = this->angleTo(pathgainData[nodeBId].wrappedPosition->point);
	    const double  ueAntennaGain = this->gain(angle);

        _pathgains[nodeBId] *= ueAntennaGain;

        if ( _connectedNodeB->id == nodeBId )
            _connectedNodeBData->stationAntennaGain = ueAntennaGain;
        if ( _closestNodeB->id == nodeBId )
            _closestNodeBData->stationAntennaGain = ueAntennaGain;
    }

	if ( debug_level > 0 )
	{
		cout << "Adjancent NodeB list for UE " << id << endl;
		multimap<double,NodeB*,greater<double>>::iterator i;
		for (i = _adjacentNodeBs.begin(); i != _adjacentNodeBs.end(); i++)
			cout << "NodeB " << i->second->id << " has pathgain of " << linear2dB(i->first) << endl;
		cout << "Affiliated NodeB is NodeB " << _connectedNodeB->id << endl << endl;
	}

    // randomise active nodeBs
    network.randomise_active_nodeBs(*_connectedNodeB);

    // Compute effective noise floor of UE
	for (size_t i=0; i<network.numSites(); i++)
        for (size_t j=0; j<network.site(i)->numNodeBs(); j++)
        {
            NodeB *nodeB = network.site(i)->nodeB(j);

            if ( nodeB->isActive && nodeB->carrierId == _connectedNodeB->carrierId )
            {
			    const double orthfact = (nodeB->id == _connectedNodeB->id) ? Interference::orthfactor() : 1.0;

                if ( nodeB->commonChannelTxPower != 0.0 )
                    _totalReceivedCommonChannelPower += nodeB->commonChannelTxPower * pathgain(nodeB->id) * orthfact;
                if ( nodeB->sharedChannelTxPower != 0.0 )
    			    _totalReceivedSharedChannelPower += nodeB->sharedChannelTxPower * pathgain(nodeB->id) * orthfact;
            }
        }

    // Compute interference power received by UE
    setInterfererNoise(network.compute_interference(*this));

	return true;
}

double UE::commonChannelSINR () const
{
	return commonChannelRcvdPwr() / (_effectiveReceivedPower - commonChannelRcvdPwr() * Interference::orthfactor());
}

double UE::SINR () const
{
    const double signalPower = receivedSignalPower();
	const double SINR =  signalPower / (_effectiveReceivedPower - signalPower * Interference::orthfactor());
	return SINR;
}

double UE::commonChannelRcvdPwr() const
{
	return _connectedNodeB->commonChannelTxPower * pathgain();
}

void UE::computeEffectiveReceivedPower()
{
	_effectiveReceivedPower = noiseFloor + interfererNoise();

	for (size_t m=0; m<Interference::network()->numSites(); m++)
		for (size_t n=0; n<Interference::network()->site(m)->numNodeBs(); n++)
		{
			const NodeB *nodeB = Interference::network()->site(m)->nodeB(n);

            if ( nodeB->isActive && _connectedNodeB->carrierId == nodeB->carrierId )
            {
			    if ( _connectedNodeB->id == nodeB->id )
				    // NodeB is that to which UE is connected, so orthfactor applies
				    _effectiveReceivedPower += nodeB->txPower() * Interference::orthfactor() * pathgain();
			    else
				    _effectiveReceivedPower += nodeB->txPower() * pathgain(nodeB->id);
            }
		}
}

double UE::intercellInterference() const
{
	double intercellInterference = _effectiveReceivedPower - noiseFloor - interfererNoise() - receivedNodeBPower();

    return std::max(intercellInterference, 0.0);
}

double UE::intracellInterference() const
{
	return receivedNodeBPower() - receivedSignalPower();
}

double UE::interToIntraRatio() const
{
    double inter = 0.0;
    for (size_t i=0; i<_pathgains.size(); i++)
        if (i != _connectedNodeB->id)
            inter += pathgain(i);

    return inter / pathgain();
}

double UE::pathgain() const
{
	return pathgain(_connectedNodeB->id);
}

void UE::print(ostream &file) const
{
    file << setw( 4) << id << " ";
    file << setw( 7) << x() << " " << setw(7) << y() << " ";
    file << setw( 6) << connectedNodeB()->id << " ";
    file << setw( 9) << -linear2dB(pathgain()) << " ";
    file << setw(11) << watts2dBm(noiseFloor) << " ";
    file << setw(16) << watts2dBm(interfererNoise(),        -999.0) << " ";
    file << setw( 8) << watts2dBm(commonChannelNoise(),     -999.0) << " ";
    file << setw(12) << watts2dBm(effectiveNoiseFloor(),    -999.0) << " ";
    file << setw( 8) << watts2dBm(effectiveReceivedPower(), -999.0) << " ";
    file << setw(12) << watts2dBm(transmitSignalPower(),    -999.0) << " ";
    file << setw(12) << watts2dBm(receivedSignalPower(),    -999.0) << " ";
    file << setw( 6) << linear2dB(SINR(),                   -999.0) << " ";
	file << setw(10) << watts2dBm(commonChannelRcvdPwr(),   -999.0) << " ";
	file << setw( 6) << linear2dB(commonChannelSINR(),      -999.0) << " ";
}

double UE::receivedNodeBPower() const
{
	return _connectedNodeB->txPower() * Interference::orthfactor() * pathgain();
}

double UE::throughput() const
{
    const double sinr = SINR();

    double bitsPerSymbol, throughput;
    switch ( _throughputMapping )
    {
    case ThroughputTable:
        bitsPerSymbol = _throughputTable->y(linear2dB(sinr));
        throughput    = _numStreams * _throughputBandwidthScaleFactor * Interference::systemBandwidth() * 1.0E6 * bitsPerSymbol;
        break;

    case ThroughputShannon:
        bitsPerSymbol = min(log2(1.0 + _throughputSINRScaleFactor * sinr), _maxBitsPerSymbol);
        throughput    = _numStreams * _throughputBandwidthScaleFactor * Interference::systemBandwidth() * 1.0E6 * bitsPerSymbol;
        break;

    case ThroughputCurveFit:
        throughput = 1.0E6 * _throughputMax / (1.0 + exp(_throughputShape * (_throughputOffset - linear2dB(sinr))));
        break;
    }

	return throughput;
}

void UE::updateFastFading()
{
    switch ( fastFadingModel )
    {
    case Rayleigh:
        for (size_t i=0; i<_fadegains.size(); i++)
            _fadegains[i] = globalRandom.rayleigh();

   		computeEffectiveReceivedPower();
        break;

    case this->NoFastFading:
        break;
    }
}

//
// class PowerControlledUE
//

bool PowerControlledUE::affiliate (Network &network, const AffiliationType affiliationType, const size_t debug_level, const size_t stats_level)
{
	if ( ! UE::affiliate(network, affiliationType, debug_level, stats_level) )
		// Cannot find a NodeB to connect to
		return false;

	_connectedNodeB->connect(this);
	_crossCouplingMultiplier = 1.0 / (pathgain() * (Interference::orthfactor() + gamma()));

    // Attempt to compute and update NodeB transmit powers
    // Return false if cannot compute powers and disconnect UE from nodeB
    if ( !Interference::update(*this, debug_level > 0) )
    {
		_connectedNodeB->disconnect_last();
        return false;
    }
    else
        return true;
}

void PowerControlledUE::print(ostream &file) const
{
	UE::print(file);

    file << setw( 6) << linear2dB(ebNo(),                   -999.0) << " ";
    file << endl;
}

double PowerControlledUE::receivedSignalPower() const
{
	return _effectiveReceivedPower / (Interference::orthfactor() + gamma());
}

double PowerControlledUE::receivedTrafficPower() const
{
	return _connectedNodeB->trafficChannelTxPower() * pathgain();
}

//
// class RateControlledUE
//

bool RateControlledUE::affiliate (Network &network, const AffiliationType affiliationType, const size_t debug_level, const size_t stats_level)
{
	if (! UE::affiliate(network, affiliationType, debug_level, stats_level) )
		// Cannot connect to NodeB
		return false;
	else
	{
		_connectedNodeB->connect(this);

        if ( _connectedNodeB->trafficChannelTxPower() == 0.0 )
            // All nodeB transmit power is shared or common
            _effectiveReceivedPower = effectiveNoiseFloor();
        else
		    computeEffectiveReceivedPower();

		return true;
	}
}
	
double RateControlledUE::sharedChannelSINR () const
{
	double signalPower = sharedChannelRcvdPwr();
	double interferencePower = _effectiveReceivedPower - signalPower * Interference::orthfactor();

	assert( signalPower >= 0.0 );
	assert( interferencePower >= 0.0 );

	return signalPower / interferencePower;
}

double RateControlledUE::sharedChannelRcvdPwr() const
{
	// Temporary fix
	// _connectedNodeB may not be top of list due to handover margin
	// need to think more about how to handle hard and soft handoever
	if ( macroDiversityBranches == 1 )
		return _connectedNodeB->sharedChannelTxPower * pathgain();
	else
	{
		double sharedChannelRcvdPwr = 0.0;
		multimap<double,NodeB*,greater<double>>::const_iterator iterator = _adjacentNodeBs.begin();
		for (size_t nodeBId = 0; nodeBId < macroDiversityBranches; nodeBId++)
		{
			sharedChannelRcvdPwr += iterator->second->sharedChannelTxPower * iterator->first;
			iterator++;
		}
		return sharedChannelRcvdPwr;
	}
}