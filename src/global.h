/* 
 * File:   global.h
 * Author: ale
 *
 * Created on September 27, 2013, 6:21 PM
 */
#include "qdp.h"

#ifndef GLOBAL_H
#define	GLOBAL_H
using namespace QDP;

#if BASE_PRECISION == 32
const Real fuzz = 1.0e-5;
#elif BASE_PRECISION == 64
const Real fuzz = 1.0e-10;
#endif

#define TIMEDIR 3

struct Cfg_t {
    int cfg_type; // storage order for stored gauge configuration
    string cfg_file;
    multi1d<int> nrow;

};

struct MCControl {
    QDP::Seed rng_seed;
    unsigned long start_update_num;
    unsigned long n_warm_up_updates;
    unsigned long n_production_updates;
    unsigned int n_updates_this_run;
    unsigned int save_interval;
    std::string save_prefix;
};

struct HBParams {
    Real beta;
    int nOver;
    /**************************************************
     * number of maximum HB tries for Creutz or KP a_0, 
     * negative or zero value - update every single link
     * (try infinitely long)
     **************************************************/
    int NmaxHB;
};


//! Main struct from input params and output restarts

struct HBControl {
    HBParams hb_params;
    MCControl mc_control;
    Cfg_t cfg;
    std::string inline_measurement_xml;
};


//! Main struct from input params and output restarts

extern Set rbOpen;
extern Set rbStapleCoeffOpen;
extern XMLReader xml_in;
extern Real s_normplaq, t_normplaq;
#endif	/* GLOBAL_H */

