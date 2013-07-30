#include "erpel.hpp"
#include "mpisupp.hpp"
#include "TinyCmd.hpp"
#include "writevtk.hpp"
#include "tinyconf.h"

using std::string;
using boost::multi_array;

static bool there_were_failures = false;

static void run_mode_and_store (const matrix &mode, int number) {
    fprintf (stderr, "[driver::run_mode_and_store] setting strain mode %i\n", number);
    erpel::sysmat_initialise (mode);
    fprintf (stderr, "[driver::run_mode_and_store] entering conjugate gradient solver\n");
    bool converged = erpel::cg_loop ();
    if (!converged)
        there_were_failures = true;
    fprintf (stderr, "[driver::run_mode_and_store] failed = %i, storing hull %i\n",
             int (!converged), number);
    erpel::store_hull (number);
}


//
// pre-processing of the structure
//

namespace {

struct PhaseInfo {
    PhaseInfo (std::string name_) {
        name = name_;
        l = k = g = 0.;
        type_of_info = INFO_VOID;
    }

    std::string name;
    enum {
        INFO_VOID,
        INFO_KG,
        INFO_LG
    }
    type_of_info;
    double l, k, g;
};

std::vector <PhaseInfo> phases_info;

void check_nonneg_modulus (double x) {
    if (! (x >= 0.)) {
        die ("[check_nonneg_modulus] Modulus %f must be >= 0.\n",
             x);
    }
}

void read_phases_file (std::string phases_filename) {

    fprintf (stderr, "    reading phases file...\n");

    phases_info = std::vector <PhaseInfo> (1, PhaseInfo ("void phase"));
    Configuration cfg (phases_filename);

    // check that phase zero is not defined in the phases file
    // (renumbering logic just isn't there, and 0 is illegal in .bin
    // anyway)
    if (cfg.has_section ("phase0"))
        die ("[read_phases_file] The phases file contains information about "
             "phase zero, which is not allowed.\n");

    for (;;) {
        int id = phases_info.size ();
        PhaseInfo ph = PhaseInfo (cpp_sprintf ("phase %i", id));

        // is there another phase?
        std::string section = cpp_sprintf ("phase%i", id);
        if (!cfg.has_section (section))
            break;

        // get the name, defaulting to phase%i
        ph.name = cfg.string (section, "name", ph.name);

        // get the shear modulus
        if (cfg.has_key (section, "shearmod"))
            ph.g = cfg.floating (section, "shearmod");
        else
            die ("[read_phases_file] No shear modulus for phase %i.\n", id);
        check_nonneg_modulus (ph.g);

        // second modulus, whichever that is
        if (cfg.has_key (section, "bulkmod")) {
            ph.k = cfg.floating (section, "bulkmod");
            ph.type_of_info = PhaseInfo::INFO_KG;
            check_nonneg_modulus (ph.k);
        } else if (cfg.has_key (section, "lamemod")) {
            ph.l = cfg.floating (section, "lamemod");
            ph.type_of_info = PhaseInfo::INFO_LG;
            check_nonneg_modulus (ph.l);
        } else {
            die ("[read_phases_file] Need a second modulus (either Lamé "
                 "or bulk) for phase %i.\n", id);
        }

        phases_info.push_back (ph);
    }
}

namespace {
    std::vector <erpel::phase_index> phase_map;

    // a small wrapper class around erpel's bin reader which allows us
    // to a) map phases  b) intercept nvox[0]
    struct translating_reader : public erpel::structure_reader {
        translating_reader (structure_reader *mother) {
            mom = mother;
        }

        virtual ~translating_reader () {
        }

        virtual void read_header (int *nvox) {
            mom->read_header (nvox);
            this->nvox[0] = nvox[0];
            this->nvox[1] = nvox[1];
            this->nvox[2] = nvox[2];
        }

        virtual void readplane (int i, erpel::structure_t *dst, int dsti) {
            mom->readplane (i, dst, dsti);
            for (int j = 0; j != nvox[1]; ++j)
            for (int k = 0; k != nvox[2]; ++k) {
                if (unsigned ((*dst)[dsti][j][k]) > phase_map.size ())
                    die ("I don't have information on phase %u, which is present in the input file.",
                         unsigned ((*dst)[dsti][j][k]));
                if (unsigned ((*dst)[dsti][j][k]) == 0u)
                    die ("Phase 0 is present in the input file");
                (*dst)[dsti][j][k] = phase_map[(*dst)[dsti][j][k]];
            }
        }

        structure_reader *mom;
        int nvox[3];
    };

    struct solid_reader : public erpel::structure_reader {
        solid_reader (int nvox0, int nvox1, int nvox2) {
            nvox[0] = nvox0;
            nvox[1] = nvox1;
            nvox[2] = nvox2;
        }

        virtual ~solid_reader () {
        }

        virtual void read_header (int *nvox_) {
            std::copy (nvox, nvox+3, nvox_);
        }

        virtual void readplane (int, erpel::structure_t *dst, int dsti) {
            for (int j = 0; j != nvox[1]; ++j)
            for (int k = 0; k != nvox[2]; ++k)
                (*dst)[dsti][j][k] = 1;
        }

        int nvox[3];
    };
}

void translate_phases () {
    // map void phases to phase zero
    phase_map.clear ();
    for (int q = 0; q != (int)phases_info.size (); ++q) {
        bool voidph = q == 0;
        switch (phases_info[q].type_of_info) {
        case PhaseInfo::INFO_LG:
            if (phases_info[q].l == 0. && phases_info[q].g == 0.)
                voidph = true;
            break;
        case PhaseInfo::INFO_KG:
            if (phases_info[q].k == 0. && phases_info[q].g == 0.)
        case PhaseInfo::INFO_VOID:
                voidph = true;
            break;
        }

        phase_map.push_back (voidph ? 0 : q);
    }

    for (int i = 0; i != (int)phases_info.size (); ++i) {
        erpel::set_phase_name (i, phases_info[i].name);
        switch (phases_info[i].type_of_info) {
        case PhaseInfo::INFO_LG:
            erpel::set_phase_lame_shear_moduli (i, phases_info[i].l, phases_info[i].g);
            break;
        case PhaseInfo::INFO_KG:
            erpel::set_phase_bulk_shear_moduli (i, phases_info[i].k, phases_info[i].g);
            break;
        case PhaseInfo::INFO_VOID:
            erpel::set_phase_bulk_shear_moduli (i, 0., 0.);
            break;
        default:
            die ("WTF?");
        }
    }
}

// obsolete command-line options
const char *obsolete_options[] = {
    "solid-phase-is",
    "kg",
    0
};

void check_for_obsolete_options (CommandLine &cmdline) {
    for (int i = 0; obsolete_options[i]; ++i) {
        if (cmdline.was_given (obsolete_options[i]))
            die ("[check_for_obsolete_options] %s was given\n",
                 obsolete_options[i]);
    }
}


}

int main (int argc, char **argv) {
    mpiutil::initialise (argc, argv);

#ifdef HAVE_LINUX_SETAFFINITY_NP
    mpiutil::affinity_auto_pin (stderr);
    mpiutil::affinity_report (stderr);
#endif

    CommandLine cmdline (argc, argv);
    outfile_prefix = cmdline.string ("prefix", "");
    mpiutil::redirect_stderr (make_filename ("logfile%i.txt", world->rank ()));
#ifdef HAVE_LINUX_SETAFFINITY_NP
    mpiutil::affinity_report (stderr);
#endif

    check_for_obsolete_options (cmdline);
    time_t time_begin = time (0);

    // set nice level -- we're running nice 10 by default
    bool be_rude = cmdline.boolean ("rude", false);
    if (!be_rude) set_nice ();

    // prescribed precision target
    double mprec = cmdline.floating ("precision", 1e-6);
    erpel::set_moduli_precision (mprec);

    erpel::set_problem (cmdline.string ("problem", "periodic_linear_elasticity"));

    // open logfile
    if_master_node {
        string cvg_filename = make_filename ("convergence.log");
        erpel::open_convergence_log (cvg_filename);
        cvg_filename = make_filename ("distortion.log");
        erpel::open_distortion_log (cvg_filename);
    }

    fprintf (stderr, "[main] reading input\n");

    // get moduli
    if (!cmdline.was_given ("phases"))
        die ("[main] Option --phases is required.\n");
    read_phases_file (cmdline.string ("phases"));

    erpel::itermax = cmdline.integer ("itermax", erpel::itermax);

    // set up moduli
    translate_phases ();
        
    // read the structure
    int nvox0;
    {
        if (cmdline.was_given ("structure")) {
            // load  a structure from a bin.bz2 file
            mpiutil::block_sigs guard;
            erpel::structure_reader *str = erpel::open_bin_file (cmdline.string ("structure"));
            translating_reader trns (str);
            erpel::structure_initialise (&trns);
            nvox0 = trns.nvox[0];
            delete str;

        } else if (cmdline.was_given ("solid-structure")) {
            // set a completely solid structure
            solid_reader *str = new solid_reader (85, 41, 20);
            erpel::structure_initialise (str);
            nvox0 = str->nvox[0];
            delete str;

        } else {
            die ("missing a structure");
        }
    }

    // open the output file for --data
    FILE *dataf = 0;
    if_master_node
    if (cmdline.was_given ("data")) {
        string datafilename = cmdline.string ("data");
        dataf = fopen (datafilename.c_str (), "wt");
        if (!dataf)
            die ("error opening %s", datafilename.c_str ());
    }

    fprintf (stderr, "[main] initialisation complete, waiting for peers...\n");
    world->barrier ();
    fprintf (stderr, "[main] past barrier, beginning computations.\n");

    string mode = cmdline.string ("mode", "mehrabadi");
    if (mode == "mehrabadi") {
        fprintf (stderr, "FULL COWIN-MEHRABADI MATRIX MODE\n");
        if (cmdline.was_given ("save-energy") ||
            cmdline.was_given ("save-dispx")  ||
            cmdline.was_given ("save-dispy")  ||
            cmdline.was_given ("save-dispz")) {

            die ("fields cannot be saved in mehrabadi mode");
        }
        run_mode_and_store (erpel::XX_MODE, 0);
        run_mode_and_store (erpel::YY_MODE, 1);
        run_mode_and_store (erpel::ZZ_MODE, 2);
        run_mode_and_store (erpel::YZ_MODE, 3);
        run_mode_and_store (erpel::ZX_MODE, 4);
        run_mode_and_store (erpel::XY_MODE, 5);

        if_master_node {
            fprintf (stderr, "Cowin-Mehrabadi matrix:\n");
            for (int i = 0; i != 6; ++i) {
                for (int j = 0; j != 6; ++j) {
                    fprintf (stderr, "%6.9f ", erpel::hull_product (i, j));
                    if (dataf)
                    fprintf (dataf,  "%6.9f ", erpel::hull_product (i, j));
                }
                fprintf (stderr, "\n");
                if (dataf)
                fprintf (dataf, "\n");
            }
        }

    } else if (mode == "cubic") {
        fprintf (stderr, "CUBIC MODULI MODE\n");
        if (cmdline.was_given ("save-energy") ||
            cmdline.was_given ("save-dispx")  ||
            cmdline.was_given ("save-dispy")  ||
            cmdline.was_given ("save-dispz")) {

            die ("fields cannot be saved in cubic mode");
        }
        run_mode_and_store (erpel::BULKMOD_MODE, 0);
        run_mode_and_store (erpel::SHEARMOD_MODE, 1);
        run_mode_and_store (erpel::CUB2ND_MODE, 2);
        if_master_node {
            fprintf (stderr, "# bulkmod shearmod shear2mod\n%.12f %.12f %.12f %i\n",
                erpel::hull_product (0, 0), erpel::hull_product (1, 1),
                erpel::hull_product (2, 2), nvox0);
            if (dataf)
            fprintf (dataf, "# bulkmod shearmod shear2mod\n%.12f %.12f %.12f\n",
                erpel::hull_product (0, 0), erpel::hull_product (1, 1),
                erpel::hull_product (2, 2));
        }
    } else if (mode == "isotropic") {
        fprintf (stderr, "ISOTROPIC MODULI MODE\n");
        if (cmdline.was_given ("save-energy") ||
            cmdline.was_given ("save-dispx")  ||
            cmdline.was_given ("save-dispy")  ||
            cmdline.was_given ("save-dispz")) {

            die ("fields cannot be saved in isotropic mode");
        }
        run_mode_and_store (erpel::BULKMOD_MODE, 0);
        run_mode_and_store (erpel::SHEARMOD_MODE, 1);
        if_master_node {
            fprintf (stderr, "# bulkmod shearmod\n%.12f %.12f\n",
                erpel::hull_product (0, 0), erpel::hull_product (1, 1));
            if (dataf)
            fprintf (dataf, "# bulkmod shearmod\n%.12f %.12f\n",
                erpel::hull_product (0, 0), erpel::hull_product (1, 1));
        }
    } else if (mode == "bulkmod") {
        fprintf (stderr, "BULK MODULUS MODE\n");
        run_mode_and_store (erpel::BULKMOD_MODE, 0);
        if_master_node {
            fprintf (stderr, "# bulkmod\n%.12f\n", erpel::hull_product (0, 0));
            if (dataf)
            fprintf (dataf, "# bulkmod\n%.12f\n", erpel::hull_product (0, 0));
        }
    } else if (mode == "shearmod") {
        fprintf (stderr, "ISOTROPIC SHEAR MODULUS MODE\n");
        run_mode_and_store (erpel::SHEARMOD_MODE, 0);
        if_master_node {
            fprintf (stderr, "# shearmod\n%.12f\n", erpel::hull_product (0, 0));
            if (dataf)
            fprintf (dataf, "# shearmod\n%.12f\n", erpel::hull_product (0, 0));
        }
    } else if (mode == "shear2mod") {
        fprintf (stderr, "CUBIC SECOND SHEAR MODULUS MODE\n");
        run_mode_and_store (erpel::CUB2ND_MODE, 0);
        if_master_node {
            fprintf (stderr, "# shear2mod\n%.12f\n", erpel::hull_product (0, 0));
            if (dataf)
            fprintf (dataf, "# shear2mod\n%.12f\n", erpel::hull_product (0, 0));
        }
    } else {
        die ("invalid mode %s", mode.c_str ());
    }

    multi_array <double, 3> energies;
    if (cmdline.was_given ("save-energy")) {
        erpel::produce_energy_density_field (&energies);
        if_master_node
            write_to_vtk (make_filename (cmdline.string ("save-energy")),
                "energy density field", energies);
    }
    if (cmdline.was_given ("save-dispx")) {
        erpel::produce_displacement_field_0 (&energies);
        if_master_node
            write_to_vtk (make_filename (cmdline.string ("save-dispx")),
                "displacement x component", energies);
    }
    if (cmdline.was_given ("save-dispy")) {
        erpel::produce_displacement_field_1 (&energies);
        if_master_node
            write_to_vtk (make_filename (cmdline.string ("save-dispy")),
                "displacement y component", energies);
    }
    if (cmdline.was_given ("save-dispz")) {
        erpel::produce_displacement_field_2 (&energies);
        if_master_node
            write_to_vtk (make_filename (cmdline.string ("save-dispz")),
                "displacement z component", energies);
    }

    fprintf (stderr, "[main] computation took %f minutes.\n",
             double (time (0) - time_begin) / 60.);
    if (there_were_failures)
        fprintf (stderr, "SOME OF THE MODES DID NOT CONVERGE TO PRESCRIBED PRECISION [see failed = 1]\n");
    erpel::close_convergence_log ();
    erpel::close_distortion_log ();
    if (dataf)
        fclose (dataf);
    fprintf (stderr, "[main] computations complete -- waiting for peers to arrive.\n");
    world->barrier ();
    fprintf (stderr, "[main] terminating MPI.\n");
    mpiutil::terminate ();
    fprintf (stderr, "[main] MPI down, dumping timers\n");
    erpel::dump_perf_counters (stderr);
    fprintf (stderr, "[main] exiting.\n");
    return 0;
}
