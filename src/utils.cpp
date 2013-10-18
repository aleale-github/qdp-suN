#include "qdp.h"
#include "global.h"
using namespace QDP;

inline void sunFill(LatticeColorMatrix& dest,
        const multi1d<LatticeReal>& r,
        int su2_index,
        const Subset & s) {

    /* Determine the SU(N) indices corresponding to the SU(2) indices */
    /* of the SU(2) subgroup $3 */
    int i1, i2;
    int found = 0;
    int del_i = 0;
    int index = -1;

    while (del_i < (Nc - 1) && found == 0) {
        del_i++;
        for (i1 = 0; i1 < (Nc - del_i); i1++) {
            index++;
            if (index == su2_index) {
                found = 1;
                break;
            }
        }
    }
    i2 = i1 + del_i;

    if (found == 0) {
        QDPIO::cerr << __func__ << ": trouble with SU2 subgroup index" << endl;
        QDP_abort(1);
    }

    /* 
     * Insert the b(k) of A_SU(2) = b0 + i sum_k bk sigma_k 
     * back into the SU(N) matrix
     */
    dest[s] = 1.0;

    pokeColor(dest[s], cmplx(r[0], r[3]), i1, i1);
    pokeColor(dest[s], cmplx(r[2], r[1]), i1, i2);
    pokeColor(dest[s], cmplx(-r[2], r[1]), i2, i1);
    pokeColor(dest[s], cmplx(r[0], -r[3]), i2, i2);


}

inline void su2Extract(multi1d<LatticeReal>& r, const LatticeColorMatrix& source,
        int su2_index, const Subset& s) {


    if (r.size() != 4) {
        QDPIO::cerr << "su2Extract: return result invalid size" << endl;
        QDP_abort(1);
    }

    /* Determine the SU(N) indices corresponding to the SU(2) indices */
    /* of the SU(2) subgroup $3 */
    int i1, i2;
    int found = 0;
    int del_i = 0;
    int index = -1;

    while (del_i < (Nc - 1) && found == 0) {
        del_i++;
        for (i1 = 0; i1 < (Nc - del_i); i1++) {
            index++;
            if (index == su2_index) {
                found = 1;
                break;
            }
        }
    }
    i2 = i1 + del_i;

    if (found == 0) {
        QDPIO::cerr << __func__ << ": trouble with SU2 subgroup index" << endl;
        QDP_abort(1);
    }

    /* Compute the b(k) of A_SU(2) = b0 + i sum_k bk sigma_k */
    r[0][s] = real(peekColor(source, i1, i1)) + real(peekColor(source, i2, i2));
    r[1][s] = imag(peekColor(source, i1, i2)) + imag(peekColor(source, i2, i1));
    r[2][s] = real(peekColor(source, i1, i2)) - real(peekColor(source, i2, i1));
    r[3][s] = imag(peekColor(source, i1, i1)) - imag(peekColor(source, i2, i2));


}

void expm12(LatticeColorMatrix& a) {


    // aux1 = aux2 = a;  a = ONE + a 
    LatticeColorMatrix aux1 = a;
    LatticeColorMatrix aux2 = a;
    LatticeColorMatrix aux3;
    a += 1;

    // Do a 12th order exponentiation
    for (int i = 2; i <= 12; ++i) {
        Real dummy = Real(1) / Real(i);

        aux3 = aux2 * aux1;
        aux2 = aux3 * dummy;
        a += aux2;
    }


}



void taproj(LatticeColorMatrix& a) {

    // a = a - a_dagger  --- a -> antihermitian matrix
    LatticeColorMatrix aux_1 = a;
    a -= adj(aux_1);

    if (Nc > 1) {
        // tmp = Im Tr[ a ]
        LatticeReal tmp = imag(trace(a));

        // a = a - (1/Nc) * Im Tr[ a] = a - (1/Nc)*tmp
        tmp *= (Real(1) / Real(Nc));

        // Is this a fill or a UnitMatrix*I?
        LatticeColorMatrix aux = cmplx(0, tmp);
        a -= aux;
    }

    // Normalisation to make taproj idempotent
    a *= (Real(1) / Real(2));


}

void su3over(LatticeColorMatrix& u,
        const LatticeColorMatrix& w,
        const Subset& sub) {


    /* V = U*W */
    LatticeColorMatrix v;
    LatticeColorMatrix tmp;
    /* Extract components r_k proportional to SU(2) submatrix su2_index */
    /* from the "SU(Nc)" matrix V. The SU(2) matrix is parametrized in the */
    /* sigma matrix basis. */
    multi1d<LatticeReal> r(4);
    LatticeBoolean lbtmp;
    LatticeReal lftmp;
    multi1d<LatticeReal> a(4);
    LatticeReal r_l;
    StopWatch swatch;

    swatch.start();

    /*# Loop over SU(2) subgroup index */
    for (int su2_index = 0; su2_index < Nc * (Nc - 1) / 2; ++su2_index) {

        v[sub] = u * w;
        su2Extract(r, v, su2_index, sub);

        /*
         * Now project onto SU(2)
         */
        r_l[sub] = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2] + r[3] * r[3]);

        // Normalize
        lbtmp[sub] = r_l > fuzz;
        lftmp[sub] = 1.0 / where(lbtmp, r_l, LatticeReal(1));

        // Fill   (r[0]/r_l, -r[1]/r_l, -r[2]/r_l, -r[3]/r_l) for r_l > fuzz
        //  and   (1,0,0,0)  for sites with r_l < fuzz
        a[0][sub] = where(lbtmp, r[0] * lftmp, LatticeReal(1));
        a[1][sub] = where(lbtmp, -(r[1] * lftmp), LatticeReal(0));
        a[2][sub] = where(lbtmp, -(r[2] * lftmp), LatticeReal(0));
        a[3][sub] = where(lbtmp, -(r[3] * lftmp), LatticeReal(0));

        /* Microcanonical updating matrix is the square of this */
        r[0][sub] = a[0] * a[0] - a[1] * a[1] - a[2] * a[2] - a[3] * a[3];

        a[0][sub] *= 2;
        r[1][sub] = a[0] * a[1];
        r[2][sub] = a[0] * a[2];
        r[3][sub] = a[0] * a[3];

        /*
         * Now fill an SU(3) matrix V with the SU(2) submatrix su2_index
         * paramtrized by a_k in the sigma matrix basis.
         */
        sunFill(v, r, su2_index, sub);

        // U = V*U
        tmp[sub] = v * u;
        u[sub] = tmp;
    }

    swatch.stop();
    QDPIO::cerr << "Time - OverRelaxation: "
            << swatch.getTimeInSeconds() << " sec" << endl;
}

void su2_a_0_kp(const LatticeReal& weight, LatticeReal& a_0,
        const Subset& sub, const int NmaxHB,
        LatticeBoolean& lAccept) {
    // received weight=SqDet*BetaMC
    int vol_cb = Layout::vol() >> 1; //volume of cb sublattice
    lAccept = (1 < 0);
    LatticeInt ilbtmp = 0;
    int vol_accept;
    LatticeReal xr1, xr2, xr3, xr4;
    int n_runs = 0;
    do {
        n_runs++;
        random(xr1);
        random(xr2);
        random(xr3);
        random(xr4);
        xr1 = -(log(xr1) / weight);
        xr3 = cos(2.0l * M_PI * xr3);
        xr3 = xr3*xr3;
        xr2 = -(log(xr2) / weight);
        //a_0=xr2+xr1*xr3;
        a_0[sub] = where(lAccept, a_0, xr2 + xr1 * xr3);
        lAccept[sub] = where(lAccept, (1 > 0), (xr4 * xr4) < (1.0l - a_0 / 2.0l));
        ilbtmp[sub] = where(lAccept, 1, 0); //convert to 1/0
        vol_accept = toInt(sum(ilbtmp));
        //QDPIO::cout<<vol_accept<<endl;

    } while ((vol_accept < vol_cb) && ((NmaxHB <= 0) || (n_runs < NmaxHB)));
    a_0 = 1.0l - a_0;
}

void print_field(const LatticeReal& a0) {
    // single node only !!!
    double temp, fuzz = 1e-44;
    const int nodeSites = Layout::sitesOnNode();
    const int nodeNumber = Layout::nodeNumber();
    for (int i = 0; i < nodeSites; i++) {
        temp = toDouble(peekSite(a0, Layout::siteCoords(nodeNumber, i)));
        if (temp > fuzz || temp<-fuzz) cout << temp << " ";
        //if(temp > -0.1 && temp <0.01) cout<<temp<<" ";
        //cout<<temp<<" ";
    }
    cout << endl;

}

void su2_a_0(const LatticeReal& weight, LatticeReal& a_0,
        const Subset& sub, const int NmaxHB,
        LatticeBoolean& lAccept) {
    // received weight=SqDet*BetaMC
    int vol_cb = Layout::vol() >> 1; //volume of cb sublattice
    int vol_accept;
    LatticeInt ilbtmp = 0;
    lAccept = (1 < 0);
    LatticeReal w_exp;
    w_exp[sub] = exp(-2.0 * weight); //too small, need to avoid
    LatticeReal x; //container for random numbers
    int n_runs = 0;
    Real RDummy;
    do {
        n_runs++;
        //random(x[sub]);
        random(x, sub);
        random(RDummy);
        //a_0[sub]=where(lAccept,a_0,1.0+log(x*(1.0-w_exp)+w_exp)/weight);
        a_0[sub] = where(lAccept, a_0, 1.0 + log(w_exp * (1 - x) + x) / weight);
        //print_field(a_0);exit(1);
        //random(x[sub]);
        random(x, sub);
        random(RDummy);
        //lAccept[sub] = where(lAccept,(1 > 0),((x*x) < (1.0-a_0*a_0)));
        //x=1.0l-x;
        x = x*x;
        LatticeReal a_0trial;
        a_0trial = a_0*a_0;
        a_0trial = 1 - a_0trial;
        lAccept[sub] = where(lAccept, (1 > 0), (x < a_0trial));
        //lAccept[sub] = where(lAccept,(1 > 0),(x > (1.0-sqrt(1.0-a_0*a_0))));
        ilbtmp[sub] = where(lAccept, 1, 0); //convert to 1/0
        vol_accept = toInt(sum(ilbtmp));
    } while ((vol_accept < vol_cb) && ((NmaxHB <= 0) || (n_runs < NmaxHB)));
}

void su2_hb_update(LatticeColorMatrix& u_mu, const LatticeColorMatrix&
        u_mu_staple, Real BetaMC,
        const Subset& sub, const int NmaxHB) {
    // ****************************************** 
    //              Parameters
    Double fuzz = 1e-16; //the smallest number for devision (det)
    // ******************************************
    //              Temp Storages
    LatticeBoolean lbtmp; // storage for lattice booleans
    //LatticeReal lftmp; // storage for lattice float tmps
    // ******************************************
    //V=U*U_staple

    multi1d<LatticeReal> r(4);
    multi1d<LatticeReal> b(4);
    multi1d<LatticeReal> a(4);

    LatticeColorMatrix v;
    LatticeReal a_0, CosTheta, Phi, a_abs, a_r, SqDet;
    LatticeBoolean lAccept;
    //    LatticeInt lWarning;

    //    int iWarning;
    Real half(0.5), RDummy;
    StopWatch swatch;

    swatch.start();

    /*# Loop over SU(2) subgroup index */
    for (int su2_index = 0; su2_index < Nc * (Nc - 1) / 2; ++su2_index) {


        v[sub] = u_mu*u_mu_staple;

        // convert to Pauli matrices parametrization
        su2Extract(r, v, su2_index, sub);
        //compensate for extra 2 !!! su(2)
        r[0][sub] *= half;
        r[1][sub] *= half;
        r[2][sub] *= half;
        r[3][sub] *= half;
        //SqDet=-1;
        SqDet[sub] = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2] + r[3] * r[3]);
        //print_field(SqDet);
        // Normalize
        //lftmp[sub]=1.0/where(lbtmp,SqDet,LatticeReal(1)); //not needed???
        // Inverse matrix u^-1
        r[0][sub] = r[0] / SqDet;
        r[1][sub] = -r[1] / SqDet;
        r[2][sub] = -r[2] / SqDet;
        r[3][sub] = -r[3] / SqDet;

        su2_a_0(BetaMC*SqDet, a_0, sub, NmaxHB, lAccept);
        //su2_a_0_kp(BetaMC*SqDet, a_0, sub, NmaxHB,lAccept);
        //print_field(a_0);
        a[0][sub] = a_0;
        //other a components

        a_abs[sub] = 1.0 - a[0] * a[0];
        lbtmp = (1 > 0);
        lbtmp[sub] = (a_abs >= 0);
        //        lWarning = where(lbtmp, 0, 1);
        //        iWarning = toInt(sum(lWarning));
        //        if (iWarning > 0) QDPIO::cerr << "wrong a_0!!!" << endl;
        lbtmp = (1 > 0);
        lbtmp[sub] = (a_abs > fuzz);
        //        lWarning = where(lbtmp, 0, 1);
        //        iWarning = toInt(sum(lWarning));
        //        if (iWarning > 0) QDPIO::cerr << "large a_0!!!" << endl;
        a_r[sub] = sqrt(a_abs);
        random(CosTheta, sub);
        random(RDummy);
        CosTheta[sub] = 1.0 - 2.0 * CosTheta;
        a[3][sub] = a_r*CosTheta;
        //LatticeReal pr_a;
        //pr_a[sub]=a[3];
        //print_field(pr_a);
        CosTheta[sub] = (1 - CosTheta * CosTheta);
        CosTheta[sub] = sqrt(CosTheta); //SinTheta
        random(Phi, sub);
        random(RDummy);
        Phi[sub] *= 8.0 * atan(1.0);
        a_r[sub] *= CosTheta; //a_r*SinTheta
        a[1][sub] = a_r * cos(Phi);
        a[2][sub] = a_r * sin(Phi);
        //u'=uu^-1 -> b = a*r
        b[0][sub] = a[0] * r[0] - a[1] * r[1] - a[2] * r[2] - a[3] * r[3];
        b[1][sub] = a[0] * r[1] + a[1] * r[0] - a[2] * r[3] + a[3] * r[2];
        b[2][sub] = a[0] * r[2] + a[2] * r[0] - a[3] * r[1] + a[1] * r[3];
        b[3][sub] = a[0] * r[3] + a[3] * r[0] - a[1] * r[2] + a[2] * r[1];
        sunFill(v, b, su2_index, sub);
        u_mu[sub] = where(lAccept, v*u_mu, u_mu);

    }

    swatch.stop();
    QDPIO::cerr << "Time - Heatbath: "
            << swatch.getTimeInSeconds() << " sec" << endl;

}

