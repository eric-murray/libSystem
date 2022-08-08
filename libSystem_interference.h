#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include <libEDM_matrix.h>

class Network;
class PowerControlledUE;
class UE;
class UEList;

class Interference {
public:
	static Network *network        () {return _network;}
	static double   orthfactor     ();
	static MHz      systemBandwidth();

	static void initialise(Network *network);
	static void recompute ();
    static void reset     ();
    static bool update    (const PowerControlledUE &ue, const bool check = false);

private:
	static Network *_network;

    static dMatrix _crossCoupling;
    static dVector _noise;
	static dMatrix _previousCrossCoupling;
	static dVector _previousNoise;

	static dMatrix _inverse;
	static dVector _z, _v;

    static dMatrix compute_crosscoupling_matrix  ();
    static bool    compute_nodeb_transmit_powers (const bool check);
    static dVector compute_noise_vector          ();
    static void    compute_ue_receive_powers     ();
    static void    update_crosscoupling_matrix   (const PowerControlledUE &ue, const bool check);
    static void    update_noise_vector           (const PowerControlledUE &ue, const bool check);

};