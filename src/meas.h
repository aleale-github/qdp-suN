/* 
 * File:   meas.h
 * Author: ale
 *
 * Created on September 29, 2013, 7:22 PM
 */

#ifndef MEAS_H
#define	MEAS_H
#include "qdp.h"
#include "global.h"

using namespace QDP;

void plaquette(multi1d<LatticeColorMatrix>& u, Real& w_plaq, Real& s_plaq, Real& t_plaq);

void doMeas(multi1d<LatticeColorMatrix>& u, unsigned long cur_update);


#endif	/* MEAS_H */

