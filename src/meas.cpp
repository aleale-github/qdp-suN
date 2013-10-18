// $Id: mesplq.cc,v 3.1 2006-08-24 21:04:31 edwards Exp $
/*! \file
 *  \brief Plaquette measurement
 */

#include "qdp.h"
#include "global.h"
using namespace QDP;


void plaquette(multi1d<LatticeColorMatrix>& u, Real& w_plaq, Real& s_plaq, Real& t_plaq) {

    multi2d<Real> plane_plaq(Nd, Nd);

    // Compute the average plaquettes
    for (int mu = 1; mu < Nd; ++mu) {
        for (int nu = 0; nu < mu; ++nu) {

            Real tmp =
                    sum(real(trace(u[mu] * shift(u[nu], FORWARD, mu) * adj(shift(u[mu], FORWARD, nu)) * adj(u[nu]))));
            plane_plaq[mu][nu] = tmp;
        }
    }

    // Normalize the planes
    for (int mu = 1; mu < Nd; ++mu)
        for (int nu = 0; nu < mu; ++nu) {
            plane_plaq[nu][mu] = plane_plaq[mu][nu];
        }
    //Layout::lattSize()[j_decay]
    // Compute basic plaquettes
    w_plaq = s_plaq = t_plaq = zero;

    for (int mu = 1; mu < Nd; ++mu) {
        for (int nu = 0; nu < mu; ++nu) {
            Real tmp = plane_plaq[mu][nu];

            if (mu == TIMEDIR || nu == TIMEDIR)
                t_plaq += tmp;
            else
                s_plaq += tmp;
        }
    }

    // Normalize
//    QDPIO::cout << s_normplaq << " " << t_normplaq << endl;
    s_plaq *= s_normplaq;
    t_plaq *= t_normplaq;
    w_plaq = (s_plaq + t_plaq) / Real(2.0);



}

void doMeas(multi1d<LatticeColorMatrix>& u,
        unsigned long cur_update) {

    std::string xmlfilename;
    std::ostringstream os;
    os << "plaquette/plaq" << "." << cur_update << ".xml";
    xmlfilename = os.str();

    XMLFileWriter xmlPlaqOut(xmlfilename);

    push(xmlPlaqOut, "Plaquette");
    write(xmlPlaqOut, "update_no", cur_update);

    Real w_plaq, s_plaq, t_plaq;

    plaquette(u, w_plaq, s_plaq, t_plaq);

    write(xmlPlaqOut, "w_plaq", w_plaq);
    write(xmlPlaqOut, "s_plaq", s_plaq);
    write(xmlPlaqOut, "t_plaq", t_plaq);

    pop(xmlPlaqOut);

    xmlPlaqOut.close();


    // Reset the default gauge field


}

