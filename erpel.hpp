#ifndef ERPEL_HPP_INCLUDED
#define ERPEL_HPP_INCLUDED

#include "shared.hpp"




namespace erpel {
    typedef std::string string;
    typedef boost::multi_array <uint8_t,3> structure_t;

    extern int itermax;
    extern matrix 
        XX_MODE, YY_MODE, ZZ_MODE,
        XY_MODE, ZX_MODE, YZ_MODE,
        BULKMOD_MODE, SHEARMOD_MODE, CUB2ND_MODE;

    // structure reader class -- following assumptions can be made:
    // 1. read_header is called exactly once by the loader code, before any
    //    planes are read; in should store three integers specifying the
    //    dataset size
    // 2. planes are requested in sequence, though not all planes will be read
    // 3. after the last plane is read, the object will be destroyed
    struct structure_reader {
        virtual ~structure_reader ();
        virtual void read_header (int *) = 0;
        virtual void readplane (int pl_no, structure_t *dst, int dst_plane) = 0;
    };

    struct phase_index {
        typedef uint8_t base_t;

        phase_index (base_t val_) : val (val_) {}
        phase_index (int val_) {
            assert (val_ >= 0);
            assert (val_ < 256);
            val = base_t (val_);
        }

        operator base_t () const { return val; }

    private:
        base_t val;
    };

    structure_reader *open_bin_file (std::string filename);
    structure_reader *reduplicating_reader (int repx, int repy, int repz);
    structure_reader *full_solid_reader (int nvoxx, int nvoxy, int nvoxz, phase_index);

    void structure_initialise (structure_reader *, int dof_per_vertex = 3);
    void sysmat_initialise (const matrix &mode);
    bool cg_loop ();

    void produce_energy_density_field (boost::multi_array <double, 3> *);
    void produce_displacement_field_0 (boost::multi_array <double, 3> *);
    void produce_displacement_field_1 (boost::multi_array <double, 3> *);
    void produce_displacement_field_2 (boost::multi_array <double, 3> *);

    void open_convergence_log (string filename);
    void close_convergence_log ();
    void open_distortion_log (string filename);
    void close_distortion_log ();

    void dump_perf_counters (FILE *);

    void set_phase_lame_shear_moduli (int phase, double lame, double shear);
    void set_phase_bulk_shear_moduli (int phase, double bulk, double shear);
    void set_phase_name              (int phase, std::string name);

    void set_moduli_precision (double);
    void set_problem (std::string);

    void store_hull (int reg);
    double hull_product (int reg0, int reg1);
 }


#endif // ERPEL_HPP_INCLUDED
