#pragma once

#define _USE_MATH_DEFINES 
#include <cmath>
#include <cstdlib>
#include <map>
#include <vector>

#include <boost/lexical_cast.hpp>

#include <libEDM_library.h>
#include <libEDM_random.h>

#include <libSystem_position.h>

using std::greater;

// forward declarations
class Network;
class NodeB;
class PowerControlledUE;
class RateControlledUE;
class UE;
class Site;
class WrappedStation;
class PathgainType;

//
// class Antenna
//

class Antenna {
public:
    const bool   omni;
    const double boresightGain;

    virtual double gain (const radians azimuth) const = 0;

    void print(const degrees stepAngle = oneDegree, const degrees startAngle = oneDegree * -180.0, const degrees stopAngle = oneDegree * 180.0);

    Antenna (const dBi boresightGain_dBi = dBi(0.0), const bool omni = false) : boresightGain(dB2linear(boresightGain_dBi)), omni(omni) {}
};

//
// class OmniAntenna
//

class OmniAntenna : public Antenna {
public:
    double gain (const radians azimuth) const {return boresightGain;}

    OmniAntenna (const dBi boresightGain_dBi = dBi(0.0)) : Antenna(boresightGain_dBi, true) {}
};

//
// class GaussianAntenna
//

class GaussianAntenna : public Antenna {
public:
    double gain (const radians azimuth) const;

    GaussianAntenna (const radians beamwidth, const dBi boresightGain_dBi, const dB frontToBackRatio_dB = dB(dINFINITY))
		: Antenna(boresightGain_dBi, false),
		  _beamwidth(beamwidth),
		  _minimumGain(dB2linear(boresightGain_dBi - frontToBackRatio_dB)),
		  _multiplier(-2.772588722 / sqr(beamwidth())) {}

private:
    const double _beamwidth, _multiplier, _minimumGain;
};

//
// class Service
//

class Service {
public:
    const double targetEbNo;
    const double userBitRate;
    const MHz    systemBandwidth;
    const double processingGain;
    const double gamma;

    Service (const dB targetEbNo_dB, const double userBitRate, const MHz systemBandwidth) :
        targetEbNo     (dB2linear(targetEbNo_dB)),
        userBitRate    (userBitRate),
        systemBandwidth(systemBandwidth),
        processingGain (systemBandwidth * 1.0E6 / userBitRate),
        gamma          (processingGain / targetEbNo)
        {}
};

//
// class Station
//

class Station {
public:
    // attributes
	const metres   height;
    const ulong    id;
    const double   noiseFloor;
    const Position position;

	// function attributes
    radians azimuth        ()                            const {return _azimuth;}
    metres  x              ()                            const {return position.point.x;}
    metres  y              ()                            const {return position.point.y;}
    string  name           ()                            const {return _type + "_" + boost::lexical_cast<string,ulong>(id);}
    double  txPower        ()                            const {return _txPower;}
    double  shadowFading   (const size_t positionId = 0) const {return position.shadowFading();}
    double  interfererNoise()                            const {return _interfererNoise;}

    // modifiers
    void setTxPower(double txPower) {_txPower = txPower;}
    void setInterfererNoise (double interfererPower) {_interfererNoise = interfererPower;}

	// angle and gain computations
	double  gain      (const radians &azimuth) const {return _antenna.omni ? _antenna.boresightGain : _antenna.gain(azimuth - this->azimuth());}
    radians angleFrom (const Point   &point)   const {return absolute_angle(point, this->position.point);}
    radians angleTo   (const Point   &point)   const {return absolute_angle(this->position.point, point);}

	// computes distance from station to a given position or station
    metres distance (const Point &point) const {return sqrt(sqr(x() - point.x) + sqr(y() - point.y));}

    // constructors
    Station (const ulong id, const string &type, const Position &position, const metres &height, const dBm &noiseFloor, const dBm &txPower, const Antenna &antenna, const radians &azimuth = globalRandom.angle())
        : id              (id),
          position        (position),
          height          (height),
          noiseFloor      (dBm2Watts(noiseFloor)),
          _azimuth        (azimuth),
          _antenna        (antenna),
          _interfererNoise(0.0),
          _type           (type),
          _txPower        (dBm2Watts(txPower))
          {}

    ~Station () {}

protected:
    radians _azimuth;

private:
    const Antenna &_antenna;
    const string   _type;
          double   _txPower;  // transmit power in Watts
          double   _interfererNoise;
};

//
// class StationList
//

template <class Type>
class StationList : public PointerVector<Type> {
public:
    StationList(const bool deleteObjectsWhenFinished = true) : PointerVector<Type>(deleteObjectsWhenFinished) {}
    virtual void print(ostream &file = cout, const bool banner = true, const size_t leadingSpaces = 0) const = 0;
};

//
// class NodeBList
//

class NodeBList : protected StationList<NodeB> {
public:
    using StationList<NodeB>::size;
	using StationList<NodeB>::operator[];

    void print(ostream &file = cout, const bool banner = true, const size_t leadingSpaces = 0) const;
    void print(ostream &file, const bool banner = true, const size_t leadingSpaces = 0, const bool printTessellated = false) const;
    void reset();
    bool check_powers(const bool banner = true) const;

    static const string nodeBBanner;
};

//
// class UEListTemplate
//

template <class Type>
class UEListTemplate : public StationList<Type> {
public:
    UEListTemplate(const bool deleteObjectsWhenFinished = true) : StationList(deleteObjectsWhenFinished) {}

    void clearAndDelete();
    void print(ostream &file, const bool banner = true, const size_t leadingSpaces = 0) const;
    void print_positions(ostream &file) const;
};

template <class Type>
void UEListTemplate<Type>::clearAndDelete()
{
    StationList<Type>::clearAndDelete();
    UE::reset();
}

template <class Type>
void UEListTemplate<Type>::print(ostream &file, const bool banner, const size_t leadingSpaces) const
{
    if ( !empty() && banner )
        file << " Id       X       Y  NodeB  Pathloss  NoiseFloor  InterfererNoise  CCNoise  EffectiveNF  RcvdPwr  TxSignalPwr  RxSignalPwr    C/I  Eb/No  CCRcvdPwr  CCC/I" << endl;

    file << fixed << setprecision(1);
    for (size_t i=0; i<this->size(); i++)
    {
        Type *ue = this->at(i);
        ue->print(file);
    }

    if ( banner )
        file << endl;
}

template <class Type>
void UEListTemplate<Type>::print_positions(ostream &file) const
{
    for (size_t i=0; i<this->size(); i++)
    {
        Type *ue = this->at(i);
        file << ue->x() << ", " << ue->y() << endl;
    }
}

//
// class UEList
//

class UEList : public UEListTemplate<UE> {
public:
    UEList(const bool deleteObjectsWhenFinished = true) : UEListTemplate(deleteObjectsWhenFinished) {}
};

//
// class PowerControlledUEList
//

class PowerControlledUEList : public UEListTemplate<PowerControlledUE> {
};

//
// class RateControlledUEList
//

class RateControlledUEList : public UEListTemplate<RateControlledUE> {
};

//
// class NodeBParams
//

class NodeBParams {
public:
    const Antenna    &antenna;
	const metres      height;
    const size_t      maxConnectedCalls;
    const dBm         noiseFloor;
    const size_t      numSectors;
	const bool        orthogonalSectorCarriers;
    const bool        rotateSectorCarriers;
    const FadeParams  fadeParams;

    const double   maxTxPower;
    const double   maxTrafficChannelTxPower;
    const double   commonChannelTxPower;
	const double   sharedChannelTxPower;

    NodeBParams(const size_t      maxConnectedCalls,
                const dBm        &maxTxPower_dBm,
                const double      commonChannelTxPowerFraction,
                const double      sharedChannelTxPowerFraction,
				const metres     &height,
                const dBm        &noiseFloor,
                const Antenna    &antenna,
                const FadeParams &fadeParams,
                const size_t      numSectors = 1,
				const bool        orthogonalSectorCarriers = false,
                const bool        rotateSectorCarriers = false)
        : maxConnectedCalls        (maxConnectedCalls),
          maxTxPower               (dBm2Watts(maxTxPower_dBm)),
          maxTrafficChannelTxPower ((1.0 - commonChannelTxPowerFraction - sharedChannelTxPowerFraction) * maxTxPower),
          commonChannelTxPower     (commonChannelTxPowerFraction * maxTxPower),
          sharedChannelTxPower     (sharedChannelTxPowerFraction * maxTxPower),
		  height                   (height),
          noiseFloor               (noiseFloor),
          antenna                  (antenna),
          fadeParams               (fadeParams),
          numSectors               (numSectors),
		  orthogonalSectorCarriers (orthogonalSectorCarriers),
          rotateSectorCarriers     (rotateSectorCarriers)
	{assert(commonChannelTxPowerFraction + sharedChannelTxPowerFraction <= 1.0);}
};

//
// class WrappedStation
//

class WrappedStation : public Station {
public:
    const PositionList positions;

    double shadowFading (const size_t positionId = 0) const {return positions.at(positionId)->shadowFading();}

    WrappedStation (const ulong id, const string &type, const PositionList &positions, const metres &height, const dBm &noiseFloor, const dBm &txPower, const Antenna &antenna, const radians &azimuth = globalRandom.angle())
        : Station(id, type, *positions.at(0), height, noiseFloor, txPower, antenna, azimuth), positions(positions) {}
};

//
// class UE
//

class UE : public Station {
private:
    struct InterferingLinkData {
    public:
        size_t numInterferingLinksLOSLimited;
        size_t numInterferingLinksLOSLimitedAfterShadowing;
        size_t numInterferingLinksMaxPathgainLimited;

        InterferingLinkData () : numInterferingLinksLOSLimited(0), numInterferingLinksLOSLimitedAfterShadowing(0), numInterferingLinksMaxPathgainLimited(0) {}
    };
public:
    enum FastFadingModel  {NoFastFading, Rayleigh};
    enum AffiliationType  {ClosestServer, BestServer, BestServerIfClosest};
    enum ThroughputMapping{ThroughputTable, ThroughputCurveFit, ThroughputShannon};

    const double            handoverMargin;
	const size_t            macroDiversityBranches;
    const FastFadingModel   fastFadingModel;

    static void reset() {_numCreatedUEs = 0;}

    size_t numInterferingLinksLOSLimited               () const {return _interferingLinkData.numInterferingLinksLOSLimited;}
    size_t numInterferingLinksLOSLimitedAfterShadowing () const {return _interferingLinkData.numInterferingLinksLOSLimitedAfterShadowing;}
    size_t numInterferingLinksMaxPathgainLimited       () const {return _interferingLinkData.numInterferingLinksMaxPathgainLimited;}

    NodeB  *closestNodeB()            const {return _closestNodeB;}
    NodeB  *connectedNodeB()          const {return _connectedNodeB;}
    double  commonChannelNoise()      const {return _totalReceivedCommonChannelPower;}
    double  effectiveNoiseFloor()     const {return noiseFloor + interfererNoise() + _totalReceivedCommonChannelPower + _totalReceivedSharedChannelPower;}
	double  geometry()                const {return receivedNodeBPower() / (intercellInterference() + noiseFloor + interfererNoise());}
    double  pathgain()                const;
    double  pathgain(size_t nodeBId)  const {return _pathgains[nodeBId] * _fadegains[nodeBId];}
    double  effectiveReceivedPower()  const {return _effectiveReceivedPower;}
	double  throughput()              const;

	PathgainType *connectedNodeBData() const {return _connectedNodeBData;}
	PathgainType *closestNodeBData()   const {return _closestNodeBData;}

	double  commonChannelSINR()     const;
	double  commonChannelRcvdPwr()  const;
	double  SINR()                  const;
	double  receivedNodeBPower()    const;
	double  intercellInterference() const;
	double  intracellInterference() const;
    double  transmitSignalPower()   const {return receivedSignalPower() / pathgain();}

    double interToIntraRatio() const;

    bool affiliate                     (Network &network, const AffiliationType affiliationType, const size_t debug_level, const size_t stats_level);
	void computeEffectiveReceivedPower ();
    void updateFastFading              ();

	// Constructors
	UE (const Position       &position,
        const metres         &height,
        const Antenna        &antenna,
        const dBm            &noiseFloor,
        const size_t          macroDiversityBranches,
        const dB             &handoverMargin,
        const size_t          numStreams,
        const double          throughputBandwidthScaleFactor,
        const double          throughputSINRScaleFactor,
        const double          maxBitsPerSymbol,
        const FastFadingModel fastFadingModel)
        :
        Station                         (_numCreatedUEs++, "UE", position, height, noiseFloor, -dBm::infinity(), antenna),
        handoverMargin                  (dB2linear(handoverMargin)),
		macroDiversityBranches          (macroDiversityBranches),
        _totalReceivedCommonChannelPower(0.0),
		_totalReceivedSharedChannelPower(0.0),
        _connectedNodeB                 (NULL),
        _numStreams                     (numStreams),
        _throughputMapping              (ThroughputShannon),
        _throughputBandwidthScaleFactor (throughputBandwidthScaleFactor),
		_throughputSINRScaleFactor      (throughputSINRScaleFactor),
        _maxBitsPerSymbol               (maxBitsPerSymbol),
        _throughputTable                (NULL),
        _throughputMax                  (1.0),
        _throughputShape                (1.0),
        _throughputOffset               (1.0),
        fastFadingModel                 (fastFadingModel)
        {}

	UE (const Position       &position,
        const metres         &height,
        const Antenna        &antenna,
        const dBm            &noiseFloor,
        const size_t          macroDiversityBranches,
        const dB             &handoverMargin,
        const double          throughputMax,
        const double          throughputShape,
        const double          throughputOffset,
        const FastFadingModel fastFadingModel)
        :
        Station                         (_numCreatedUEs++, "UE", position, height, noiseFloor, -dBm::infinity(), antenna),
        handoverMargin                  (dB2linear(handoverMargin)),
		macroDiversityBranches          (macroDiversityBranches),
        _totalReceivedCommonChannelPower(0.0),
		_totalReceivedSharedChannelPower(0.0),
        _connectedNodeB                 (NULL),
        _numStreams                     (1),
        _throughputMapping              (ThroughputCurveFit),
        _throughputBandwidthScaleFactor (1.0),
		_throughputSINRScaleFactor      (1.0),
        _maxBitsPerSymbol               (numeric_limits<double>::infinity()),
        _throughputTable                (NULL),
        _throughputMax                  (throughputMax),
        _throughputShape                (throughputShape),
        _throughputOffset               (throughputOffset),
        fastFadingModel                 (fastFadingModel)
        {}

   	UE (const Position       &position,
        const metres         &height,
        const Antenna        &antenna,
        const dBm            &noiseFloor,
        const size_t          macroDiversityBranches,
        const dB             &handoverMargin,
        const LookUpTable    *throughputTable,
        const double          throughputBandwidthScaleFactor,
        const FastFadingModel fastFadingModel)
        :
        Station                         (_numCreatedUEs++, "UE", position, height, noiseFloor, -dBm::infinity(), antenna),
        handoverMargin                  (dB2linear(handoverMargin)),
		macroDiversityBranches          (macroDiversityBranches),
        _totalReceivedCommonChannelPower(0.0),
		_totalReceivedSharedChannelPower(0.0),
        _connectedNodeB                 (NULL),
        _numStreams                     (1),
        _throughputMapping              (ThroughputTable),
        _throughputBandwidthScaleFactor (throughputBandwidthScaleFactor),
		_throughputSINRScaleFactor      (1.0),
        _maxBitsPerSymbol               (0.0),
        _throughputTable                (throughputTable),
        _throughputMax                  (1.0),
        _throughputShape                (1.0),
        _throughputOffset               (1.0),
        fastFadingModel                 (fastFadingModel)
        {}

    // Destructor
    ~UE ();

	// virtual functions
	virtual double receivedSignalPower()                     const = 0;
    virtual void   print              (ostream &file = cout) const;

protected:
    double                                      _totalReceivedCommonChannelPower, _totalReceivedSharedChannelPower;
    double                                      _effectiveReceivedPower;
    NodeB                                      *_connectedNodeB, *_closestNodeB;
    vector<double>                              _pathgains, _fadegains;
    multimap<double, NodeB*, greater<double>>  _adjacentNodeBs;

private:
    static ulong _numCreatedUEs;

    const ThroughputMapping _throughputMapping;

    const double _throughputBandwidthScaleFactor;
	const double _throughputSINRScaleFactor;
    const double _throughputMax;
	const double _throughputShape;
	const double _throughputOffset;
    const double _maxBitsPerSymbol;
    const size_t _numStreams;

    const LookUpTable *_throughputTable;

    PathgainType *_connectedNodeBData, *_closestNodeBData;

    InterferingLinkData _interferingLinkData;
};


//
// class PowerControlledUE
//

class PowerControlledUE : public UE {
public:
    double  crossCouplingMultiplier() const {return _crossCouplingMultiplier;}
    double  ebNo()                    const {return SINR() * _service->processingGain;}
    double  gamma()                   const {return _service->gamma;}
	double  receivedSignalPower()     const;
    double  receivedTrafficPower()    const;

    bool affiliate (Network &network, AffiliationType affiliationType, const size_t debug_level, const size_t stats_level);

	void print(ostream &file = cout) const;

    PowerControlledUE (const Position       &position,
                       const metres          height,
                       const Antenna        &antenna,
                       const dBm             noiseFloor,
                       const Service        *service,
                       const size_t          macroDiversityBranches = 1,
                       const dB              handoverMargin = dB(0.0),
                       const size_t          numStreams = 1,
                       const double          throughputBandwidthScaleFactor = 1.0,
                       const double          throughputSINRScaleFactor = 1.0,
                       const double          maxBitsPerSymbol = numeric_limits<double>::infinity(),
                       const FastFadingModel fastFadingModel = NoFastFading)
		: UE(position, height, antenna, noiseFloor, macroDiversityBranches, handoverMargin, numStreams, throughputBandwidthScaleFactor, throughputSINRScaleFactor, maxBitsPerSymbol, fastFadingModel), _service(service) {}

private:
    const Service *const _service;
	      double         _crossCouplingMultiplier;
};

//
// class RateControlledUE
//

class RateControlledUE : public UE {
public:
	double receivedSignalPower()  const {return sharedChannelRcvdPwr();}
	double sharedChannelSINR()    const;
	double sharedChannelRcvdPwr() const;
	
    bool affiliate (Network &network, AffiliationType affiliationType, const size_t debug_level, const size_t stats_level);

	RateControlledUE (const Position       &position,
                      const metres         &height,
                      const Antenna        &antenna,
                      const dBm            &noiseFloor,
                      const size_t          numStreams                     = 1,
                      const double          throughputBandwidthScaleFactor = 1.0,
                      const double          throughputSINRScaleFactor      = 1.0,
                      const double          maxBitsPerSymbol               = numeric_limits<double>::infinity(),
                      const size_t          macroDiversityBranches         = 1,
                      const dB             &handoverMargin                 = dB(0.0),
                      const FastFadingModel fastFadingModel                = NoFastFading)
		: UE(position, height, antenna, noiseFloor, macroDiversityBranches, handoverMargin, numStreams, throughputBandwidthScaleFactor, throughputSINRScaleFactor, maxBitsPerSymbol, fastFadingModel) {}

	RateControlledUE (const Position       &position,
                      const metres         &height,
                      const Antenna        &antenna,
                      const dBm            &noiseFloor,
                      const double          throughputMax,
                      const double          throughputShape,
                      const double          throughputOffset,
                      const size_t          macroDiversityBranches = 1,
                      const dB             &handoverMargin         = dB(0.0),
                      const FastFadingModel fastFadingModel        = NoFastFading)
		: UE(position, height, antenna, noiseFloor, macroDiversityBranches, handoverMargin, throughputMax, throughputShape, throughputOffset, fastFadingModel) {}

	RateControlledUE (const Position       &position,
                      const metres         &height,
                      const Antenna        &antenna,
                      const dBm            &noiseFloor,
                      const LookUpTable    *throughputTable,
                      const double          throughputBandwidthScaleFactor = 1.0,
                      const size_t          macroDiversityBranches         = 1,
                      const dB             &handoverMargin                 = dB(0.0),
                      const FastFadingModel fastFadingModel                = NoFastFading)
		: UE(position, height, antenna, noiseFloor, macroDiversityBranches, handoverMargin, throughputTable, throughputBandwidthScaleFactor, fastFadingModel) {if ( throughputTable == NULL ) error("NULL throughput table passed to UE constructor");}

private:
};

//
// class NodeB
//

class NodeB : public WrappedStation { 
public:
    bool isActive;

	const size_t carrierId;
    const double commonChannelTxPower, sharedChannelTxPower;
    const size_t maxConnectedCalls;
    const double maxTxPower;
    const double maxTrafficChannelTxPower;
    const size_t sectorId;

    size_t numConnectedUEs       () const {return _connectedUEs.size();}
    size_t numConnectedPCUEs     () const {return _connectedPCUEs.size();}
    double trafficChannelTxPower () const {return _trafficChannelTxPower;}

    double receivedUEPower (const UE &ue) {return ue.txPower() * ue.pathgain(id);}
    double receivedUESINR  (const UE &ue) {return receivedUEPower(ue) / (noiseFloor + interfererNoise());}

    const UEList                *connectedUEs      ()                  const {return &_connectedUEs;}
    const PowerControlledUEList *connectedPCUEs    ()                  const {return &_connectedPCUEs;}
    const RateControlledUEList  *connectedRCUEs    ()                  const {return &_connectedRCUEs;}
          UE                    *ue                (const size_t ueId) const {return  _connectedUEs.at(ueId);}
          PowerControlledUE     *powerControlledUE (const size_t ueId) const {return  _connectedPCUEs.at(ueId);}

    void setTrafficChannelTxPower (const double trafficChannelTxPower) {_trafficChannelTxPower = trafficChannelTxPower;}

    void connect (PowerControlledUE *ue) {_connectedPCUEs.push_back(ue); _connectedUEs.push_back(ue);}
    void connect (RateControlledUE  *ue) {_connectedRCUEs.push_back(ue); _connectedUEs.push_back(ue);}

	void disconnect_last();
    void reset          ();

    Position generate_random_position (const FadeParams &fadeParams) const;

    bool check_powers() const;

    NodeB (const ulong id, const size_t sectorId, Position &position, const NodeBParams &nodeBParams, const size_t carrierId, const radians azimuth = radians(), const string namePrefix = "NodeB") :
        WrappedStation           (id, namePrefix, PositionList(position), nodeBParams.height, nodeBParams.noiseFloor, watts2dBm(nodeBParams.commonChannelTxPower + nodeBParams.sharedChannelTxPower), nodeBParams.antenna, azimuth),
        sectorId                 (sectorId),
		carrierId                (carrierId),
        maxConnectedCalls        (nodeBParams.maxConnectedCalls),
        maxTxPower               (nodeBParams.maxTxPower),
        maxTrafficChannelTxPower (nodeBParams.maxTrafficChannelTxPower),
        commonChannelTxPower     (nodeBParams.commonChannelTxPower),
        sharedChannelTxPower     (nodeBParams.sharedChannelTxPower),
        _trafficChannelTxPower   (0.0),
        isActive                 (true)
        {}

    NodeB (const ulong id, const size_t sectorId, const PositionList &positions, const NodeBParams &nodeBParams, const size_t carrierId, const radians azimuth = radians(), const string namePrefix = "NodeB") :
        WrappedStation           (id, namePrefix, positions, nodeBParams.height, nodeBParams.noiseFloor, watts2dBm(nodeBParams.commonChannelTxPower + nodeBParams.sharedChannelTxPower), nodeBParams.antenna, azimuth),
        sectorId                 (sectorId),
		carrierId                (carrierId),
        maxConnectedCalls        (nodeBParams.maxConnectedCalls),
        maxTxPower               (nodeBParams.maxTxPower),
        maxTrafficChannelTxPower (nodeBParams.maxTrafficChannelTxPower),
        commonChannelTxPower     (nodeBParams.commonChannelTxPower),
        sharedChannelTxPower     (nodeBParams.sharedChannelTxPower),
        _trafficChannelTxPower   (0.0),
        isActive                 (true)
        {}

private:
    double                _trafficChannelTxPower;
    UEList                _connectedUEs;
    PowerControlledUEList _connectedPCUEs;
    RateControlledUEList  _connectedRCUEs;
};

//
// class Interferer
//

class Interferer : public WrappedStation {
public:
	Interferer (const ulong id, Position &position, const metres height, const dBm txPower) :
        WrappedStation(id, "Interferer", PositionList(position), height, dBm(0.0), txPower, *_antenna, radians())
        {}

    Interferer (const ulong id, PositionList &positions, const metres height, const dBm txPower) :
        WrappedStation(id, "Interferer", positions, height, dBm(0.0), txPower, *_antenna, radians())
		{}

private:
    static const OmniAntenna *_antenna;
};

//
// class InterfererList
//

class InterfererList : public StationList<Interferer> {
public:
    void print(ostream &file = cout, const bool banner = true, const size_t leadingSpaces = 0) const;
};