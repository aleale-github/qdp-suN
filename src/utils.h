/* 
 * File:   utils.h
 * Author: ale
 *
 * Created on September 29, 2013, 6:56 PM
 */

#ifndef UTILS_H
#define	UTILS_H
#include "qdp.h"
using namespace QDP;

void expm12(LatticeColorMatrix& a);

void taproj(LatticeColorMatrix& a);


void su3over(LatticeColorMatrix& u,
        const LatticeColorMatrix& w,
        const Subset& sub);

void su2_a_0_kp(const LatticeReal& weight, LatticeReal& a_0,
        const Subset& sub, const int NmaxHB,
        LatticeBoolean& lAccept);

void print_field(const LatticeReal& a0);

void su2_a_0(const LatticeReal& weight, LatticeReal& a_0,
        const Subset& sub, const int NmaxHB,
        LatticeBoolean& lAccept);

inline void su2Extract(multi1d<LatticeReal>& r, const LatticeColorMatrix& source,
        int su2_index, const Subset& s);

void su2_hb_update(LatticeColorMatrix& u_mu, const LatticeColorMatrix&
        u_mu_staple, Real BetaMC,
        const Subset& sub, const int NmaxHB);

inline void sunFill(LatticeColorMatrix& dest,
        const multi1d<LatticeReal>& r,
        int su2_index,
        const Subset & s);

void reunit(LatticeColorMatrix & xa);


#endif	/* UTILS_H */

