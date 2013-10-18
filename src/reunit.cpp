// -*- C++ -*-
// $Id: reunit.cc,v 3.3 2007-02-22 21:11:49 bjoo Exp $

/*! \file
 *  \ingroup gauge
 *  \author Subsetting added by A. Hart
 *  \param[in,out] xa  The descriptor of matrices to be reunitarized.
 *            Must be of type LatticeColorMatrix
 *  \param[in] bad Descriptor of flags indicating sites violating unitarity.
 *            Only used if ruflag = REUNITARIZE_LABEL or
 *            REUNITARIZE_ERROR.
 *  \param[in] ruflag Can also be REUNITARIZE in which case the
 *            matrices are reunitarized but no complaints are made.
 *  \param[out] numbad Total number of matrices violating unitarity.
 *            ONLY USED IF ruflag is testing for ERROR or LABEL. 
 *  \param[in] mstag  An (un)ordered subset of sites
 *  \brief Reunitarize in place a color matrix to SU(N)
 *
 *  Reunitarize (to a SU(N)) inplace the matrix XA under some option
 */

#include "qdp.h"
#include "global.h"
#include "utils.h"

void reunit(LatticeColorMatrix& xa,
        const Subset& mstag) {


    QDP::StopWatch swatch;
    swatch.reset();
    swatch.start();

    multi2d<LatticeComplex> a(Nc, Nc);
    multi2d<LatticeComplex> b(Nc, Nc);
    LatticeReal t1;
    LatticeComplex t2;
    LatticeReal t3;

    // Extract initial components 
    for (int i = 0; i < Nc; ++i)
        for (int j = 0; j < Nc; ++j)
            (a[i][j])[mstag] = peekColor(xa, i, j);



    /* normalise the first row */
    /* t1 = sqrt(u^t . u) */
    t1[mstag] = localNorm2(a[0][0]);
    for (int c = 1; c < Nc; ++c)
        t1[mstag] += localNorm2(a[c][0]);
    t1[mstag] = sqrt(t1);


    /* overwrite the first row with the rescaled value */
    /* u <- u/t1 */
    t3[mstag] = 1 / t1;
    for (int c = 0; c < Nc; ++c)
        (a[c][0])[mstag] *= t3;

    /* Do Gramm-Schmidt on the remaining rows */
    for (int j = 1; j < Nc; j++) {
        for (int i = 0; i < j; i++) {
            /* t2 <- u^t.v */
            t2[mstag] = adj(a[0][i]) * a[0][j];
            for (int c = 1; c < Nc; ++c) {
                t2[mstag] += adj(a[c][i]) * a[c][j];
            }

            /* orthogonalize the j-th row relative to the i-th row */
            /* v <- v - t2*u */
            for (int c = 0; c < Nc; ++c) {
                (a[c][j])[mstag] -= t2 * a[c][i];
            }
        }

        /* normalise the j-th row */
        /* t1 = sqrt(v^t . v) */
        t1[mstag] = localNorm2(a[0][j]);
        for (int c = 1; c < Nc; ++c)
            t1[mstag] += localNorm2(a[c][j]);
        t1[mstag] = sqrt(t1);

        /* overwrite the j-th row with the rescaled value */
        /* v <- v/t1 */
        t3[mstag] = 1 / t1;
        for (int c = 0; c < Nc; ++c)
            (a[c][j])[mstag] *= t3;
    }

    /* Now we have a unitary matrix. We need to multiply the last
       row with a phase to make the determinant 1. */
    /* We compute the determinant by LU decomposition */
    for (int j = 0; j < Nc; j++)
        for (int i = 0; i < Nc; i++)
            (b[j][i])[mstag] = a[j][i];

    for (int j = 0; j < Nc; j++) {
        for (int i = 0; i <= j; i++) {
            t2[mstag] = b[j][i];
            for (int c = 0; c < i; c++)
                t2[mstag] -= b[c][i] * b[j][c];

            (b[j][i])[mstag] = t2;
        }

        for (int i = (j + 1); i < Nc; i++) {
            t2[mstag] = b[j][i];
            for (int c = 0; c < j; c++)
                t2[mstag] -= b[c][i] * b[j][c];

            (b[j][i])[mstag] = adj(b[j][j]) * t2 / localNorm2(b[j][j]);
        }
    }

    /* The determinant */
    t2[mstag] = b[0][0] * b[1][1];
    for (int c = 2; c < Nc; c++)
        t2[mstag] *= b[c][c];

    /* The phase of the determinant */
    t2[mstag] = conj(t2);
    for (int c = 0; c < Nc; ++c)
        (a[c][Nc - 1])[mstag] *= t2;



    // Insert final reunitarized components 
    for (int i = 0; i < Nc; ++i)
        for (int j = 0; j < Nc; ++j)
            pokeColor(xa[mstag], a[i][j], i, j);

    swatch.stop();
    QDPIO::cerr << "Time - Reunit: "
            << swatch.getTimeInSeconds() << " sec" << endl;
}



///*Reunit check**************************************************************/
//
//void printcolor(ColorMatrix& a) {
//    Real rex, imx;
//    QDPIO::cout << "{";
//    for (int ic = 0; ic < Nc; ic++) {
//        QDPIO::cout << "{";
//        for (int jc = 0; jc < Nc; jc++) {
//            rex = real(peekColor(a, ic, jc));
//            imx = imag(peekColor(a, ic, jc));
//            QDPIO::cout << rex << " + " << imx << " I";
//            if (jc != Nc - 1) QDPIO::cout << ",";
//        }
//        QDPIO::cout << "}";
//        if (ic != Nc - 1) QDPIO::cout << ",";
//    }
//    QDPIO::cout << "}";
//    QDPIO::cout << endl << endl;
//    ;
//
//
//}
//
//
#if 0
void reunit_check_n(LatticeColorMatrix& xa) {

    multi2d<LatticeComplex> a(Nc, Nc);
    multi2d<LatticeComplex> b(Nc, Nc);
    LatticeComplex t2;
    Real fuzz = 1.0e-14;

    int numbad;

    // Extract initial components 
    for (int i = 0; i < Nc; ++i)
        for (int j = 0; j < Nc; ++j)
            (a[i][j]) = peekColor(xa, i, j);

    for (int j = 0; j < Nc; j++)
        for (int i = 0; i < Nc; i++)
            (b[j][i]) = a[j][i];

    for (int j = 0; j < Nc; j++) {
        for (int i = 0; i <= j; i++) {
            t2 = b[j][i];
            for (int c = 0; c < i; c++)
                t2 -= b[c][i] * b[j][c];

            (b[j][i]) = t2;
        }

        for (int i = (j + 1); i < Nc; i++) {
            t2 = b[j][i];
            for (int c = 0; c < j; c++)
                t2 -= b[c][i] * b[j][c];

            (b[j][i]) = adj(b[j][j]) * t2 / localNorm2(b[j][j]);
        }
    }

    /* The determinant */
    t2 = b[0][0] * b[1][1];
    for (int c = 2; c < Nc; c++)
        t2 *= b[c][c];

    numbad = toInt(sum(where(
            ((real(t2) - Real(1.0000000000000000)) > fuzz || imag(t2) > fuzz),
            LatticeInteger(1), LatticeInteger(0))));

    QDPIO::cerr << "Matrices violating Sunitarity = " << numbad
            << "Sum=" << sum(t2) - cmplx(Real(Layout::vol()), Real(0.0000000000000000))<< endl;


}
#endif
/*Reunit check**************************************************************/


