/* 
 * File:   reunit.h
 * Author: ale
 *
 * Created on October 16, 2013, 6:32 PM
 */

#ifndef REUNIT_H
#define	REUNIT_H


#include "qdp.h"
using namespace QDP;

void reunit(LatticeColorMatrix& xa, const Subset& mstag);
void reunit_check_n(LatticeColorMatrix& xa);

#endif	/* REUNIT_H */

