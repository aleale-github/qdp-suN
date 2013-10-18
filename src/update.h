/* 
 * File:   newfile.h
 * Author: ale
 *
 * Created on September 29, 2013, 7:20 PM
 */

#ifndef NEWFILE_H
#define	NEWFILE_H

#include "qdp.h"
using namespace QDP;

void staple(LatticeColorMatrix& u_mu_staple, multi1d<LatticeColorMatrix>& u,
        int mu, Set& gaugeset, int cb, const Real& beta);

void update(multi1d<LatticeColorMatrix>& u, const HBParams& hbp);


#endif	/* NEWFILE_H */

