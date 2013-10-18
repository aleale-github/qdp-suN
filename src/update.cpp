#include "qdp.h"
#include "global.h"
#include "utils.h"
#include "reunit.h"
#include <string.h> 

using namespace QDP;

#ifdef OPENBC

const Set& getGaugeSet(const int& mu) {
    if (mu != TIMEDIR) {
        return rb;
    } else {
        return rbOpen;
    }
}
#else

const Set& getGaugeSet(const int&) {
    return rb;
}
#endif

void staple(LatticeColorMatrix& u_mu_staple,
        multi1d<LatticeColorMatrix>& u,
        int mu,
        const Set& gaugeset,
        int cb,
        const Real& beta) {


    u_mu_staple = zero;
    LatticeColorMatrix tmp;
    StopWatch swatch;

    swatch.start();

    for (int nu = 0; nu < Nd; ++nu) {

        if (mu == nu) continue;

        // +forward staple
        tmp[gaugeset[cb]] =
                shift(u[nu], FORWARD, mu) *
                adj(shift(u[mu], FORWARD, nu)) *
                adj(u[nu]) *
                beta;

        tmp[gaugeset[cb]] +=
                shift(shift(adj(u[nu]), FORWARD, mu), BACKWARD, nu) *
                shift(adj(u[mu]), BACKWARD, nu) *
                shift(u[nu], BACKWARD, nu) *
                beta;
#ifdef OPENBC
        if ((nu != TIMEDIR) && (mu != TIMEDIR))
            tmp[rbStapleCoeffOpen[cb]] *= Real(0.5);
#endif
        u_mu_staple[gaugeset[cb]] += tmp;
    }
    swatch.stop();
    QDPIO::cerr << "Time - Staples: "
            << swatch.getTimeInSeconds() << " sec" << endl;

}

void update(multi1d<LatticeColorMatrix>& u, const HBParams& hbp) {

    LatticeColorMatrix u_mu_staple;

    const int num_subsets = 2;
    
    for ( int cb = 0; cb < num_subsets; ++cb) {
        for ( int mu = 0; mu < Nd; ++mu) {

            /* Overrelaxation Step*/
            for (int i = 0; i < hbp.nOver; i++) {
                // Calculate the staple
                staple(u_mu_staple, u, mu, getGaugeSet(mu), cb, hbp.beta);
                // Overrelaxation
                su3over(u[mu], u_mu_staple, getGaugeSet(mu)[cb]);
            }

            /* Heatbath step */
            // Ri-Calculate the staple
            staple(u_mu_staple, u, mu, getGaugeSet(mu), cb, hbp.beta);
            // heatbath update of su2 subgroups
            su2_hb_update(u[mu], u_mu_staple,
                    Real(2.0 / Nc), getGaugeSet(mu)[cb], hbp.NmaxHB);

            /* Reunitarize */
            reunit(u[mu],getGaugeSet(mu)[cb]);
//            reunit_check_n(u[mu]);

        } // closes mu loop            

    } // closes cb loop

}
