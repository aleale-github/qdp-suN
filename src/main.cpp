/* 
 * File:   main.cpp
 * Author: ale
 *
 * Created on September 27, 2013, 4:36 PM
 */


#include "qdp.h"
#include "global.h"
#include "utils.h"
#include "update.h"
#include "meas.h"
#include "reunit.h"
#include <time.h>
#include <stdio.h> 

using namespace QDP;


HBControl hb_control;
XMLReader xml_in;
Real s_normplaq, t_normplaq;
Set rbOpen, rbStapleCoeffOpen;

class SetRBFuncOpen : public SetFunc {
public:

    int operator() (const multi1d<int>& coordinate) const {
        int sum = 0;
        for (int m = 0; m < coordinate.size(); ++m)
            sum += coordinate[m];
        if (coordinate[TIMEDIR] < Layout::lattSize()[TIMEDIR] - 1) {
            return sum & 1;
        } else {
            return 2;
        };

    }

    int numSubsets() const {
        return 3;
    }
};

class SetRBStapleCoeffOpen : public SetFunc {
public:

    int operator() (const multi1d<int>& coordinate) const {
        int sum = 0;
        for (int m = 0; m < coordinate.size(); ++m)
            sum += coordinate[m];
        if ((coordinate[TIMEDIR] == Layout::lattSize()[TIMEDIR] - 1) ||
                (coordinate[TIMEDIR] == 0)) {
            return sum & 1;
        } else {
            return 2;
        };

    }

    int numSubsets() const {
        return 3;
    }
};

void proginfo(XMLWriter& xml) {

    push(xml, "ProgramInfo");

    push(xml, "code_version");
    write(xml, "qdp", QDP_PACKAGE_VERSION);
    write(xml, "basePrecision", BASE_PRECISION);
    pop(xml);
    time_t now;

    if (time(&now) == -1) {
        QDPIO::cerr << "proginfo: Cannot get the time.\n";
        return;
    }
    tm *tp = localtime(&now);

    char date[64];
    strftime(date, 63, "%d %b %y %X %Z", tp);
    write(xml, "run_date", date);

    push(xml, "Setgeom");
    write(xml, "latt_size", Layout::lattSize());
    write(xml, "logical_size", Layout::logicalSize());
    write(xml, "subgrid_size", Layout::subgridLattSize());
    write(xml, "total_volume", Layout::vol());
    write(xml, "subgrid_volume", Layout::sitesOnNode());
    pop(xml);
    pop(xml);


}

//struct Cfg_t {
//    CfgType cfg_type; // storage order for stored gauge configuration
//    string cfg_file;
//    multi1d<int> nrow;
//
//};


// Configuration input

void read(XMLReader& xml, const string& path, Cfg_t& input) {
    XMLReader inputtop(xml, path);

    read(inputtop, "cfg_type", input.cfg_type);
    read(inputtop, "cfg_file", input.cfg_file);
    read(inputtop, "nrow", input.nrow);
}

// Write a config struct

void write(XMLWriter& xml, const string& path, const Cfg_t& cfg) {
    push(xml, "Cfg");
    write(xml, "cfg_type", cfg.cfg_type);
    write(xml, "cfg_file", cfg.cfg_file);
    write(xml, "nrow", cfg.nrow);
    pop(xml);
}

//! Holds gauge action


//! Reader
//struct HBParams {
//    Real beta;
//    int nOver;
//    /**************************************************
//     * number of maximum HB tries for Creutz or KP a_0, 
//     * negative or zero value - update every single link
//     * (try infinitely long)
//     **************************************************/
//    int NmaxHB;
//};

void read(XMLReader& xml, const std::string& path, HBParams& p) {
    try {
        XMLReader paramtop(xml, path);
        read(paramtop, "NmaxHB", p.NmaxHB);
        read(paramtop, "nOver", p.nOver);
        read(paramtop, "beta", p.beta);
    } catch (const std::string& e) {
        QDPIO::cerr << "Caught Exception reading HBParams: " << e << endl;
        QDP_abort(1);
    }
}

//! Writer

void write(XMLWriter& xml, const std::string& path, const HBParams& p) {
    push(xml, path);

    write(xml, "NmaxHB", p.NmaxHB);
    write(xml, "nOver", p.nOver);
    write(xml, "beta", p.beta);
    pop(xml);
}

//struct MCControl {
//    QDP::Seed rng_seed;
//    unsigned long start_update_num;
//    unsigned long n_warm_up_updates;
//    unsigned long n_production_updates;
//    unsigned int n_updates_this_run;
//    unsigned int save_interval;
//    std::string save_prefix;
//};


//! Params controlling running of monte carlo

void read(XMLReader& xml, const std::string& path, MCControl& p) {
    try {
        XMLReader paramtop(xml, path);
        read(paramtop, "./RNG", p.rng_seed);
        read(paramtop, "./StartUpdateNum", p.start_update_num);
        read(paramtop, "./NWarmUpUpdates", p.n_warm_up_updates);
        read(paramtop, "./NProductionUpdates", p.n_production_updates);
        read(paramtop, "./NUpdatesThisRun", p.n_updates_this_run);
        read(paramtop, "./SaveInterval", p.save_interval);
        read(paramtop, "./SavePrefix", p.save_prefix);
        if (p.n_updates_this_run % p.save_interval != 0)
            throw string("UpdateThisRun not a multiple of SaveInterval");
    } catch (const std::string& e) {
        QDPIO::cerr << "Caught Exception reading MCControl: " << e << endl;
        QDP_abort(1);
    }
}

void write(XMLWriter& xml, const std::string& path, const MCControl& p) {
    push(xml, path);

    write(xml, "RNG", p.rng_seed);
    write(xml, "StartUpdateNum", p.start_update_num);
    write(xml, "NWarmUpUpdates", p.n_warm_up_updates);
    write(xml, "NProductionUpdates", p.n_production_updates);
    write(xml, "NUpdatesThisRun", p.n_updates_this_run);
    write(xml, "SaveInterval", p.save_interval);
    write(xml, "SavePrefix", p.save_prefix);

    pop(xml);
}

//struct HBControl {
//    HBParams hb_params;
//    MCControl mc_control;
//    Cfg_t cfg;
//    std::string inline_measurement_xml;
//};


//! Reader

void read(XMLReader& xml_in, const std::string& path, HBControl& p) {
    try {
        XMLReader paramtop(xml_in, path);

        read(paramtop, "HB", p.hb_params);
        read(paramtop, "MCControl", p.mc_control);
        read(paramtop, "Cfg", p.cfg);

        if (paramtop.count("./InlineMeasurements") == 0) {
            XMLBufferWriter dummy;
            push(dummy, "InlineMeasurements");
            pop(dummy); // InlineMeasurements
            p.inline_measurement_xml = dummy.printCurrentContext();
        } else {
            XMLReader measurements_xml(paramtop, "./InlineMeasurements");
            std::ostringstream inline_os;
            measurements_xml.print(inline_os);
            p.inline_measurement_xml = inline_os.str();
            QDPIO::cout << "InlineMeasurements are: " << endl;
            QDPIO::cout << p.inline_measurement_xml << endl;
        }
    } catch (const std::string& e) {
        QDPIO::cerr << "Caught Exception reading HBControl: " << e << endl;
        QDP_abort(1);
    }
}


//! Writer

void write(XMLWriter& xml, const std::string& path, const HBControl& p) {
    push(xml, path);

    write(xml, "Cfg", p.cfg);
    write(xml, "MCControl", p.mc_control);
    xml << p.inline_measurement_xml;
    write(xml, "HB", p.hb_params);

    pop(xml);
}

//--------------------------------------------------------------------------
// Specialise

MCControl newMCHeader(const MCControl& mc_control,
        unsigned long update_no) {


    // Copy old params
    MCControl p_new = mc_control;

    // Get Current RNG Seed
    QDP::RNG::savern(p_new.rng_seed);

    // Set the current traj number
    p_new.start_update_num = update_no;

    // Reset the warmups
    p_new.n_warm_up_updates = 0;

    // Set the num_updates_this_run
    unsigned long total = mc_control.n_production_updates;

    if (total < mc_control.n_updates_this_run + update_no) {
        p_new.n_updates_this_run = total - update_no;
    }



    return p_new;
}

void gaugeStartup(multi1d<LatticeColorMatrix>& u,
        Cfg_t& cfg) {

    XMLReader gauge_file_xml;
    XMLReader gauge_xml;

    switch (cfg.cfg_type) {

        case 0:
        {

            QDPFileReader to(gauge_file_xml, cfg.cfg_file, QDPIO_SERIAL);

            read(to, gauge_xml, u);
            if (to.bad()) {
                QDPIO::cerr << __func__ << ": error reading file " << cfg.cfg_file << endl;
                QDP_abort(1);
            }
            close(to);
        }
            break;

        case 1:
        {
            QDPIO::cout << "Starting up disordered (random/hot) config" << endl;
            for (int mu = 0; mu < Nd; ++mu) {
                gaussian(u[mu]); // Gaussian fill
                taproj(u[mu]); // Traceless anti-hermitian projection
                expm12(u[mu]); // Exponentiate         
                reunit(u[mu], all);
            }
            XMLBufferWriter file_xml, record_xml;
            push(file_xml, "gauge");
            write(file_xml, "id", int(0));
            pop(file_xml);
            push(record_xml, "disordered");
            pop(record_xml);

            gauge_file_xml.open(file_xml);
            gauge_xml.open(record_xml);
#ifdef OPENBC
            u[TIMEDIR][rbOpen[2]] = zero;
#endif  
        }
            break;

        case 2:
        {
            QDPIO::cout << "Starting up unit gauge (free) config" << endl;
            u = 1;

            XMLBufferWriter file_xml, record_xml;
            push(file_xml, "gauge");
            write(file_xml, "id", int(0));
            pop(file_xml);
            push(record_xml, "unit");
            pop(record_xml);

            gauge_file_xml.open(file_xml);
            gauge_xml.open(record_xml);
#ifdef OPENBC
            u[TIMEDIR][rbOpen[2]] = zero;
#endif  
        }
            break;

        default:
            QDPIO::cerr << __func__ << ": Configuration type is unsupported." << endl;
            QDP_abort(1);
            break;
    }

    // Reunitarize gauge field (eg for Double prec?)
    // This should get rid of those pesky big Delta H's 
    // going from single to double and they can't possibly hurt
    // in either single or double

    for (int mu = 0; mu < Nd; mu++) {
        reunit(u[mu], all);
    }

#ifdef OPENBC
    checkOpenBC(u[TIMEDIR]);
#endif  


}


// Write a Gauge field in QIO format

/*
 * \param file_xml    xml reader holding config info ( Modify )
 * \param record_xml  xml reader holding config info ( Modify )
 * \param u           gauge configuration ( Modify )
 * \param file        path ( Read )
 * \param serpar      either QDPIO_SERIAL, QDPIO_PARALLEL ( Read )
 */
void writeGauge(XMLBufferWriter& file_xml,
        XMLBufferWriter& record_xml,
        const multi1d<LatticeColorMatrix>& u,
        const string& file,
        QDP_serialparallel_t serpar) {

    QDPFileWriter to(file_xml, file, QDPIO_SINGLEFILE, serpar, QDPIO_OPEN);
    if (to.bad()) {
        QDPIO::cerr << __func__ << ": error writing file " << file << endl;
        QDP_abort(1);
    }

    write(to, record_xml, u); // Write in native precision
    close(to);
}

void saveState(const HBParams& update_params,
        MCControl& mc_control,
        unsigned long update_no,
        const string& inline_measurement_xml,
        const multi1d<LatticeColorMatrix>& u) {

    MCControl mc_new = newMCHeader(mc_control, update_no);

    // Files
    std::ostringstream restart_data_filename;
    std::ostringstream restart_config_filename;

    unsigned long save_num = update_no / mc_control.save_interval;
    restart_data_filename << mc_control.save_prefix << ".ini.xml" << save_num;
    restart_config_filename << mc_control.save_prefix << ".lime" << save_num;

    {
        HBControl hb;
        hb.hb_params = update_params;
        hb.mc_control = mc_new;
        hb.inline_measurement_xml = inline_measurement_xml;

        // Set the name of the restart file
        hb.cfg.cfg_file = restart_config_filename.str();

        // Hijack this for now and assumes it means what I want it to mean
        hb.cfg.cfg_type = 0;

        // Write a restart DATA file from the buffer XML
        XMLFileWriter restart_xml(restart_data_filename.str().c_str());
        write(restart_xml, "purgaug", hb);
        restart_xml.close();
    }

    {
        // Save the config

        // some dummy header for the file
        XMLBufferWriter file_xml;
        push(file_xml, "HB");
        proginfo(file_xml);
        pop(file_xml);

        XMLBufferWriter config_xml;
        push(config_xml, "ChromaHB");
        write(config_xml, "MCControl", mc_new);
        write(config_xml, "HBItr", update_params);
        pop(config_xml);

        // Save the config
        writeGauge(file_xml,
                config_xml,
                u,
                restart_config_filename.str(),
                QDPIO_SERIAL);
    }


}

//--------------------------------------------------------------------------

void doWarmUp(multi1d<LatticeColorMatrix>& u, HBControl& hb_control) {

    // Set the update number
    unsigned long cur_update = 0;

    // Compute how many updates to do
    unsigned long to_do = hb_control.mc_control.n_warm_up_updates;

    QDPIO::cout << "WarmUp Control: About to do " << to_do << " updates" << endl;

    for (unsigned int i = 0; i < to_do; i++) {

        // Increase current update counter
        cur_update++;

        // Do the update, but with no measurements
        update(u, hb_control.hb_params); //one hb sweep

        // Do measurements
        doMeas(u, cur_update);
    }
}


//--------------------------------------------------------------------------

void doProd(multi1d<LatticeColorMatrix>& u, HBControl& hb_control) {


    // Set the update number
    unsigned long cur_update = hb_control.mc_control.start_update_num;

    // Compute how many updates to do
    unsigned long total_updates = hb_control.mc_control.n_production_updates;

    unsigned long to_do = 0;
    if (total_updates > hb_control.mc_control.n_updates_this_run + cur_update + 1) {
        to_do = hb_control.mc_control.n_updates_this_run;
    } else {
        to_do = total_updates - cur_update;
    }

    QDPIO::cout << "MC Control: About to do " << to_do << " updates" << endl;
    QDPIO::cout << "Current update num = " << cur_update << endl;
    StopWatch swatch;

    // XML Output

    for (int i = 0; i < to_do; i++) {

        // Increase current update counter
        cur_update++;

        // Decide if the next update is a warm up or not
        QDPIO::cout << "Doing Update: " << cur_update << " warm_up_p = " << false << endl;


        // Do the update
        swatch.start();
        update(u, hb_control.hb_params); //one hb sweep
        swatch.stop();

        QDPIO::cerr << "*** Total Time - update no " << i << " : "
                << swatch.getTimeInSeconds() << " sec ***" << endl;
        swatch.reset();

        // Do measurements
        doMeas(u, cur_update);

        // Save if needed
        if (cur_update % hb_control.mc_control.save_interval == 0) {
            saveState(hb_control.hb_params,
                    hb_control.mc_control,
                    cur_update,
                    hb_control.inline_measurement_xml, u);
        }

    }

}

void init(int* argc, char ***argv, HBControl& hb_control) {
    if (!QDP_isInitialized())
        QDP_initialize(argc, argv);

    string input_filename;

    for (int i = 0; i < *argc; i++) {
        // Get argv[i] into a string
        string argv_i = string((*argv)[i]);

        // Search for -i or --chroma-i
        if (argv_i == string("-h") || argv_i == string("--help")) {
            QDPIO::cerr << "Usage: " << (*argv)[0] << "  <options>" << endl
                    << "   -h           help\n"
                    << "   --help       help\n"
                    << "   -i           [xml input file name]\n"
                    << endl;
            QDP_abort(0);
        }

        // Search for -i or --chroma-i
        if (argv_i == string("-i")) {
            if (i + 1 < *argc) {
                input_filename = string((*argv)[i + 1]);
                // Skip over next
                i++;
            } else {
                // i + 1 is too big
                QDPIO::cerr << "Error: dangling -i specified. " << endl;
                QDP_abort(1);
            }
        }
    }

    try {
        xml_in.open(input_filename);

        read(xml_in, "/purgaug", hb_control);

    } catch (const std::string& e) {
        QDPIO::cerr << "Caught Exception reading input XML: " << e << endl;
        QDP_abort(1);
    } catch (std::exception& e) {
        QDPIO::cerr << "Caught standard library exception: " << e.what() << endl;
        QDP_abort(1);
    } catch (...) {
        QDPIO::cerr << "Caught unknown exception " << endl;
        QDP_abort(1);
    }

    Layout::setLattSize(hb_control.cfg.nrow);
    Layout::create();
#ifdef OPENBC    
    // Normalization due to open boundary condition (temporal links outside [0,T]
    // are set to zero
    s_normplaq = 2.0 / Real((Nd - 1)*(Nd - 2) * Nc * Layout::vol());
    t_normplaq = 1.0 / Real((Nd - 1) * Nc * (Layout::vol() - Layout::vol() / Layout::lattSize()[TIMEDIR]));
    // Initialize the red/black checkerboard for open BC
    rbOpen.make(SetRBFuncOpen());
    rbStapleCoeffOpen.make(SetRBStapleCoeffOpen());
#else
    s_normplaq = 2.0 / Real((Nd - 1)*(Nd - 2) * Nc * Layout::vol());
    t_normplaq = 1.0 / Real((Nd - 1) * Nc * Layout::vol());
#endif

    //Initialize the Random Generator
    QDP::RNG::setrn(hb_control.mc_control.rng_seed);


}

#ifdef OPENBC

void checkOpenBC(LatticeColorMatrix& utime) {
    LatticeInt lWarning;
    LatticeBoolean lbtmp = true;
    int iError;
    ColorMatrix zerobc = zero;
    lbtmp = (utime[rbOpen[2]] == zerobc);
    lWarning = where(lbtmp, 0, 1);
    iError = toInt(sum(lWarning));
    if (iError > 0) {
        QDPIO::cerr << "Error: OpenBC Selected but temporal "
                "links on the boundaries non null!!!" << endl;
        QDP_abort(1);
    }

}
#endif

int main(int argc, char *argv[]) {

    init(&argc, &argv, hb_control);

    // Start up the config
    multi1d<LatticeColorMatrix> u(Nd);

    gaugeStartup(u, hb_control.cfg);
    
    doMeas(u, 0);

    // Run
    try {
        // If warmups are required, do them first
        if (hb_control.mc_control.n_warm_up_updates > 0) {
            doWarmUp(u, hb_control);
            hb_control.mc_control.n_warm_up_updates = 0; // reset
        }
        // Do the production updates
        doProd(u, hb_control);
    } catch (const std::string& e) {
        QDPIO::cerr << "Caught Exception: " << e << endl;
        QDP_abort(1);
    }

    QDP_finalize();

    exit(0);
}


//    {
//    XMLFileWriter outxml("gaugebad.xml");
//    write(outxml, "gaugebad0", u[0]);
//    outxml.close();
//    }
//    {
//    XMLFileWriter outxml("gaugebad.xml");
//    write(outxml, "gaugebad1", u[1]);
//    outxml.close();
//    }
//    {
//    XMLFileWriter outxml("gaugebad.xml");
//    write(outxml, "gaugebad2", u[2]);
//    outxml.close();
//    }
//    {
//    XMLFileWriter outxml("gaugebad.xml");
//    write(outxml, "gaugebad3", u[3]);
//    outxml.close();
//    }

    /*Test***********************************************************/
    //    LatticeColorMatrix u2;
    //    int mui=0;
    //    u2=u[mui];
    //    
    //    reunita(u[mui]);
    //    reunitx(u2);
    //    LatticeColorMatrix u_diff=u[mui]-u2;
    //        
    //    //Setting the gauge field to 1
    //    LatticeInt lWarning;
    //    LatticeBoolean lbtmp = false;
    //    int iError = 0;
    //    multi2d<LatticeReal> a(Nc, Nc);
    //
    //    // Extract initial components 
    //    for (int i = 0; i < Nc; ++i) {
    //        for (int j = 0; j < Nc; ++j) {
    //            (a[i][j]) = real(peekColor(u_diff, i, j)*conj(peekColor(u_diff, i, j)));
    //        }
    //    }
    //    for (int i = 0; i < Nc; ++i) {
    //        for (int j = 0; j < Nc; ++j) {
    //            lbtmp = (a[i][j] < 1.0e-15 );
    //            lWarning = where(lbtmp, 0, 1);
    //            iError += toInt(sum(lWarning));
    //        }
    //    }
    //
    //    QDPIO::cout << "Error is " << iError << endl;
    /*Test***********************************************************/


/*Test*******************************************************************/
//    int mu = 0;
//    LatticeColorMatrix u_mu_staple1, u_mu_staple2,staple_diff;
//    staplex(u_mu_staple1, u, mu, rb, 0, hb_control.hb_params.beta);
//    staplechromax(u_mu_staple2, u, mu, 0, hb_control.hb_params.beta);
//    staple_diff=u_mu_staple1 - u_mu_staple2;
//    //Setting the gauge field to 1
//    LatticeInt lWarning;
//    LatticeBoolean lbtmp = true;
//    int iError = 0;
//    multi2d<LatticeReal> a(Nc, Nc);
//
//
//    // Extract initial components 
//    for (int i = 0; i < Nc; ++i) {
//        for (int j = 0; j < Nc; ++j) {
//            (a[i][j]) = real(peekColor(staple_diff, i, j)*conj(peekColor(staple_diff, i, j)));
//        }
//    }
//    for (int i = 0; i < Nc; ++i) {
//        for (int j = 0; j < Nc; ++j) {
//            lbtmp = (a[i][j] < 1.0e-15 );
//            lWarning = where(lbtmp, 0, 1);
//            iError += toInt(sum(lWarning));
//        }
//    }
//
//
//    QDPIO::cout << "Error is " << iError << endl;
//
//    //    XMLFileWriter stapleout("staple_diff.xml");
//    //    write(stapleout, "staple", u_mu_staple1-u_mu_staple2);
//    //    stapleout.close();
//
//    QDP_finalize();
//
//    exit(0);
//    /*Test*******************************************************************/
