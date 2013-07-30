#include "erpel.hpp"
#include "revision.hpp"
#include "shared.hpp"
#include "sse2_dgemv.hpp"
#include "stencils_iso.hpp"
#include "perf.hpp"
#include "mpisupp.hpp"
#include <pthread.h>
#include <map>
using boost::multi_array;

typedef multi_array <int, 3> index_map_t;
typedef index_map_t::extent_range index_map_range_t;

namespace erpel {

static structure_t structure;

       int itermax = 1000;

static std::vector <matrix> final_stencil;
static std::map <int, std::string> phase_name;
static vector        stencil_temp, stencil_temp2;
static int           num_total_vox;
static int           num_hull_vox, max_extent_vox;
static int           nvox[3];
static multi_array <int, 3> index_map;
static double        prev_res_sq;
static vector        solguess;
static vector        residue;
static vector        lastbase;
static vector        btempo;
static double        cvg_threshold = 1e-4;
static int           int_dof;                 // internal degrees of freedom
static int           accountable_dof;         // dof where this node is accountable
static int           iterctr;
static int           acceptcrit = 3;
static int           max_sequential_cg = 300;  // after 30 steps of CG, re-calc residue
static FILE         *cvg_log = 0;
static FILE         *dist_log = 0;
static int           cvg_log_steps = 300;      // steps between status reports
                                              // (also, steps between cvg tests!)
static double        hull_weight = -1.;       // see sysmat_initialise

static bool          fix_hull_displacements_x = false;  // if false, periodic
static bool          fix_hull_displacements_y = false;  // bc are used.
static bool          fix_hull_displacements_z = false;

// range of voxel layers this note has to manage
static int voxel_plane_begin, voxel_plane_end;
// dito for vertex planes (may differ)
static int vertex_plane_begin, vertex_plane_end;
// vertex plane this node computes the inner products for
static int vertex_accable_end;

static perf_timer_class pt_prod_with_sysmat_e_b ("prod_with_sysmat_enforce_bounds");
static perf_timer_class pt_parallel_inner_prod ("parallel_inner_prod");
static perf_timer_class pt_parallel_norm_inf ("parallel_norm_inf");
static perf_timer_class pt_prod_with_sysmat_noncore ("prod_with_sysmat_noncore");
static perf_timer_class pt_prod_with_sysmat_core_communication ("prod_with_sysmat_core_communication");
static perf_timer_class pt_prod_with_sysmat_core_join_thread ("prod_with_sysmat_core_join_thread");
static perf_timer_class pt_compute_stress_for_slices_bkgd ("compute_stress_for_slices_bkgd");
static perf_timer_class pt_prod_with_sysmat_core_before_fork ("prod_with_sysmat_core_before_fork");
static perf_timer_class pt_prod_with_sysmat_core_fork ("prod_with_sysmat_core_fork");
static perf_timer_class pt_prod_with_sysmat_core_zero ("prod_with_sysmat_core_zero");
static perf_timer_class pt_sysmat_initialise ("sysmat_initialise");

/// perform a parallel inner product on a distributed vector.
/// only the first accountable_dof elements are calculated locally.
/// the rest is assumed to be computed on another node.
/// data is assumed to be synchronised in advance.
static
double parallel_inner_prod (const vector &lhs, const vector &rhs) {

    perf_timer pt (pt_parallel_inner_prod);

    double sdf = serial_inner_prod (lhs, rhs, accountable_dof);
    return all_reduce (*world, sdf, std::plus <double> ());
}

// calculate infinity norm of a distributed vector.
static
double parallel_norm_inf (const vector &in) {

    perf_timer pt (pt_parallel_norm_inf);

    double nrm = serial_norm_inf (in, accountable_dof);
    return all_reduce (*world, nrm, mpi::maximum <double> ());
}

// merge the shared degrees of freedom, i.e. superpose the forces
// calculated locally and those calculated on the neighbour node.
// this only works with support from MPI at the moment, so can't be
// used in the serial code.

static
vector lower_slice, lower_slice_r, upper_slice, upper_slice_r;

enum {
    P_SLICE_DNSYNC = 16,     // MPI message tags
    P_SLICE_UPSYNC = 8
};

static bool message_shown2 = false;

static
void merge_overlapping_dof (vector *out, int flags) {

    if (!message_shown2) {
        fprintf (stderr, "entering merge_overlapping_dof for the first time...\n");
    }

    bool up   = !! (flags & P_SLICE_UPSYNC);
    bool down = !! (flags & P_SLICE_DNSYNC);

    // copy the data we need to send to slices
    // FIXME this could be done more efficiently, using a simple
    // memcpy.  (this would not, however, work once we go from
    // slab decompo to cubes decomposition)
    size_t uctr = 0u, lctr = 0u;
    for (int j = 0; j != nvox[1]+1; ++j)
    for (int k = 0; k != nvox[2]+1; ++k) {

        int lower_x = index_map[vertex_plane_begin][j][k];
        int upper_x = index_map[vertex_plane_end-1][j][k];

        if (upper_x != -1) {
            for (int c = 0; c != 3; ++c)
                upper_slice[uctr++] = (*out)[upper_x + c];
        }

        if (lower_x != -1) {
            for (int c = 0; c != 3; ++c)
                lower_slice[lctr++] = (*out)[lower_x + c];
        }
    }

    // compute which our peers are
    int upper_mate = (world->rank () + 1) % world->size ();
    int lower_mate = (world->rank () + world->size () - 1) % world->size ();

    if (uctr != upper_slice.size ()) {
        if (uctr > upper_slice.size ())
            die ("this should not happen");
        upper_slice.resize (uctr);
    }
    if (lctr != lower_slice.size ()) {
        if (lctr > lower_slice.size ())
            die ("this should not happen");
        lower_slice.resize (lctr);
    }

    /*  this could be faster, however, with MPICH on woody, it isn't.
    mpi::request reqs[4];
    mpi::request *reqs_end = reqs;

    if (down)
        *reqs_end++ = world->isend (lower_mate, P_SLICE_DNSYNC, lower_slice);
    if (up)
        *reqs_end++ = world->irecv (upper_mate, P_SLICE_DNSYNC, upper_slice_r);

    mpi::wait_all (reqs, reqs_end);
    reqs_end = reqs;

    if (!message_shown2) {
        fprintf (stderr, "have received first batch...\n");
    }

    if (up)
        *reqs_end++ = world->isend (upper_mate, P_SLICE_UPSYNC, upper_slice);
    if (down)
        *reqs_end++ = world->irecv (lower_mate, P_SLICE_UPSYNC, lower_slice_r);

    mpi::wait_all (reqs, reqs_end);
    */

#define DB(x)  if (!message_shown2) fprintf(stderr, "%s\n", x)

    if (down && world->rank () % 2) {
        DB("send odd DN");
        world->send (lower_mate, P_SLICE_DNSYNC, lower_slice);
    }
    if (up) {
        DB("recv DN");
        world->recv (upper_mate, P_SLICE_DNSYNC, upper_slice_r);
    }
    if (down && world->rank () % 2 == 0) {
        DB("send even DN");
        world->send (lower_mate, P_SLICE_DNSYNC, lower_slice);
    }

    if (up &&  world->rank () % 2 == 0) {
        DB("send even UP");
        world->send (upper_mate, P_SLICE_UPSYNC, upper_slice);
    }
    if (down) {
        DB("recv UP");
        world->recv (lower_mate, P_SLICE_UPSYNC, lower_slice_r);
    }
    if (up &&  world->rank () % 2) {
        DB("send odd UP");
        world->send (upper_mate, P_SLICE_UPSYNC, upper_slice);
    }
#undef DB

    // after the exchange, merge the data we received.
    // (if required)
    uctr = 0u, lctr = 0u;
    if (up)
        for (int j = 0; j != nvox[1]+1; ++j)
        for (int k = 0; k != nvox[2]+1; ++k) {

            int upper_x = index_map[vertex_plane_end-1][j][k];
            if (upper_x != -1) {
                for (int c = 0; c != 3; ++c)
                    (*out)[upper_x + c] += upper_slice_r[uctr++];
            }
        }

    if (down)
        for (int j = 0; j != nvox[1]+1; ++j)
        for (int k = 0; k != nvox[2]+1; ++k) {

            int lower_x = index_map[vertex_plane_begin][j][k];
            if (lower_x != -1) {
                for (int c = 0; c != 3; ++c)
                    (*out)[lower_x + c] += lower_slice_r[lctr++];
            }
        }

    if (up && uctr != upper_slice_r.size ())
        die ("mismatch in upper packet size: %u %u\n", uctr, (unsigned)upper_slice_r.size ());
    if (down && lctr != lower_slice_r.size ())
        die ("mismatch in lower packet size: %u %u\n", lctr, (unsigned)lower_slice_r.size ());

    if (!message_shown2) {
        fprintf (stderr, "exiting merge_overlapping_dof for the first time...\n");
        message_shown2 = true;
    }
}

// compute the stresses for the given slices
static
void compute_stress_for_slices (vector *out, const vector &in,
                                int slbegin, int slend) {
    assert (slbegin <= slend);
    for (int i = slbegin; i != slend; ++i)
    for (int j = 0; j != nvox[1]; ++j)
    for (int k = 0; k != nvox[2]; ++k) {
        int ph = structure[i][j][k];
        if (ph) {
            // load the vertices into stencil_temp
            int indic[24];
            indic[ 0] = index_map[i][j][k]    ;
            indic[ 1] = index_map[i][j][k] + 1;
            indic[ 2] = index_map[i][j][k] + 2;
            indic[ 3] = index_map[i][j][k+1]    ;
            indic[ 4] = index_map[i][j][k+1] + 1;
            indic[ 5] = index_map[i][j][k+1] + 2;
            indic[ 6] = index_map[i][j+1][k]    ;
            indic[ 7] = index_map[i][j+1][k] + 1;
            indic[ 8] = index_map[i][j+1][k] + 2;
            indic[ 9] = index_map[i][j+1][k+1]    ;
            indic[10] = index_map[i][j+1][k+1] + 1;
            indic[11] = index_map[i][j+1][k+1] + 2;
            indic[12] = index_map[i+1][j][k]    ;
            indic[13] = index_map[i+1][j][k] + 1;
            indic[14] = index_map[i+1][j][k] + 2;
            indic[15] = index_map[i+1][j][k+1]    ;
            indic[16] = index_map[i+1][j][k+1] + 1;
            indic[17] = index_map[i+1][j][k+1] + 2;
            indic[18] = index_map[i+1][j+1][k]    ;
            indic[19] = index_map[i+1][j+1][k] + 1;
            indic[20] = index_map[i+1][j+1][k] + 2;
            indic[21] = index_map[i+1][j+1][k+1]    ;
            indic[22] = index_map[i+1][j+1][k+1] + 1;
            indic[23] = index_map[i+1][j+1][k+1] + 2;
            for (int Z = 0; Z != 24; ++Z)
                stencil_temp[Z] = in[indic[Z]];
            dgemv_24 (&stencil_temp2, final_stencil[ph], stencil_temp);
            for (int Z = 0; Z != 24; ++Z)
                (*out)[indic[Z]] += stencil_temp2[Z];
        }
    }
}

// this mess implements a background worker thread which does most of the
// stress computation while the main thread does MPI comm.
namespace bkgd {

static int have_initialized = 0;
static int exit_flag;
static const vector *in;
static vector *out;
static pthread_mutex_t begin_mutex, end_mutex;
static pthread_t thr_id;

static
void *thread (void *) {
    mpiutil::affinity_report (stderr);

    for (;;) {
        pthread_mutex_lock (&begin_mutex);
        if (exit_flag)
            return 0;
        { perf_timer p (pt_compute_stress_for_slices_bkgd);
        compute_stress_for_slices (out, *in,
            voxel_plane_begin+1, voxel_plane_end-1);
        }
        pthread_mutex_unlock (&end_mutex);
    }
}

static
void create () {
    if (have_initialized)
        return;
    have_initialized = 1;
    int err;
    err = pthread_mutex_init (&begin_mutex, 0);
    if (err) { perror ("pthread_mutex_init"); abort (); }
    err = pthread_mutex_init (&end_mutex, 0);
    if (err) { perror ("pthread_mutex_init"); abort (); }
    pthread_mutex_lock (&begin_mutex);
    pthread_mutex_lock (&end_mutex);
    exit_flag = 0;
    err = pthread_create (&thr_id, 0, &thread, 0);
    if (err) { perror ("pthread_create"); abort (); }
}

#if 0  // commented out to avoid compiler warning 
       // driver code lacks a call to destroy()
static
void destroy () {
    exit_flag = 1;
    // we should hold begin_mutex right now.
    // signal the background thread to exit.
    pthread_mutex_unlock (&begin_mutex);
    void *dummy;
    pthread_join (thr_id, &dummy);
}
#endif

void trigger (vector *out_, const vector &in_) {
    out = out_;
    in = &in_;
    pthread_mutex_unlock (&begin_mutex);
}

static
void join () {
    pthread_mutex_lock (&end_mutex);
}

} // bkgd

static
bool message_shown = false;

static
void prod_with_sysmat_zero (vector *out) {
    perf_timer pt (pt_prod_with_sysmat_core_zero);
    noalias (*out) = zero_vector (int_dof);
}

// compute the product of the displacement vector with the
// system matrix, yielding the stress vector.
// in parallel mode, only a subset is computed.
static
void prod_with_sysmat_core (vector *out, const vector &in, bool enforce_x_pbc) {

    if (!message_shown) {
        fprintf (stderr, "entering prod_with_sysmat_enforce_bounds for the first time...\n");
    }

    prod_with_sysmat_zero (out);

    { perf_timer pt (pt_prod_with_sysmat_core_before_fork);
    // first, do the top & bottom slices.
    // these have to be synced with other nodes, involving
    // communication.  this process can run in the background
    // while we're doing the lion's share of the work.
    compute_stress_for_slices (out, in, voxel_plane_begin, voxel_plane_begin+1);
    compute_stress_for_slices (out, in, voxel_plane_end-1, voxel_plane_end);
    }

    { perf_timer pt (pt_prod_with_sysmat_core_fork);
    // launch worker thread into background
    bkgd::trigger (out, in);
    }

    { perf_timer pt3 (pt_prod_with_sysmat_core_communication);

    // while the main computation is running, sync with other nodes
    int flags = P_SLICE_DNSYNC | P_SLICE_UPSYNC;
    if (world->rank () == 0)
        flags ^= P_SLICE_DNSYNC;
    if (world->rank () == world->size () - 1)
        flags ^= P_SLICE_UPSYNC;
    if (enforce_x_pbc)
        flags = P_SLICE_DNSYNC | P_SLICE_UPSYNC;
    merge_overlapping_dof (out, flags);

    } // end of comm timer

    { perf_timer pt3 (pt_prod_with_sysmat_core_join_thread);
    bkgd::join ();
    }

    if (!message_shown) {
        fprintf (stderr, "exiting prod_with_sysmat_enforce_bounds for the first time...\n");
        message_shown = true;
    }
}

static
void prod_with_sysmat (vector *out, const vector &in) {
    prod_with_sysmat_core (out, in, false);
}

// same as above, but respect the boundary conditions.
// we do this by introducing additional forces as required.
// in fact, we only zero the existing forces of course.
static
void prod_with_sysmat_enforce_bounds (vector *out, const vector &in) {

    perf_timer p (pt_prod_with_sysmat_e_b);

    if (fix_hull_displacements_x) {
        // we zero the stress components on the x-faces anyway -- no
        // need to do communication.
        prod_with_sysmat_core (out, in, false);

        int zero_plane = -1;
        if (vertex_plane_begin == 0)
            zero_plane = vertex_plane_begin;
        else if (vertex_plane_end-1 == nvox[0])
            zero_plane = vertex_plane_end-1;

        if (world->rank () == world->size () - 1 || world->rank () == 0)
            if (zero_plane == -1)
                die ("we should really zero something!");

        if (zero_plane != -1)
        for (int j = 0; j != nvox[1]+1; ++j)
        for (int k = 0; k != nvox[2]+1; ++k) {
            int idx = index_map[zero_plane][j][k];
            if (idx != -1)
                (*out)[idx] = (*out)[idx+1]  = (*out)[idx+2] = 0.;
        }
    } else {
        // periodic bc are handled implicitly by having the topmost
        // and the bottommost (?) node interchange forces.
        prod_with_sysmat_core (out, in, true);
    }

    perf_timer p2 (pt_prod_with_sysmat_noncore);

    // x bounds are done -- handle y now.
    if (fix_hull_displacements_y) {
        for (int j = vertex_plane_begin; j != vertex_plane_end; ++j)
        for (int k = 0; k != nvox[2]+1; ++k) {
            if (index_map[j][0][k] != -1) {
                (*out)[index_map[j][0][k]    ] = 0.;
                (*out)[index_map[j][0][k] + 1] = 0.;
                (*out)[index_map[j][0][k] + 2] = 0.;
                (*out)[index_map[j][nvox[1]][k]    ] = 0.;
                (*out)[index_map[j][nvox[1]][k] + 1] = 0.;
                (*out)[index_map[j][nvox[1]][k] + 2] = 0.;
            }
        }
    } else {
        for (int j = vertex_plane_begin; j != vertex_plane_end; ++j)
        for (int k = 0; k != nvox[2]+1; ++k) {
            if (index_map[j][0][k] != -1) {
                (*out)[index_map[j][0][k]    ] += (*out)[index_map[j][nvox[1]][k]    ];
                (*out)[index_map[j][0][k] + 1] += (*out)[index_map[j][nvox[1]][k] + 1];
                (*out)[index_map[j][0][k] + 2] += (*out)[index_map[j][nvox[1]][k] + 2];
                (*out)[index_map[j][nvox[1]][k]    ] = (*out)[index_map[j][0][k]    ];
                (*out)[index_map[j][nvox[1]][k] + 1] = (*out)[index_map[j][0][k] + 1];
                (*out)[index_map[j][nvox[1]][k] + 2] = (*out)[index_map[j][0][k] + 2];
            }
        }
    }

    // y bounds are done -- handle z now.
    if (fix_hull_displacements_z) {
        for (int j = vertex_plane_begin; j != vertex_plane_end; ++j)
        for (int k = 0; k != nvox[1]+1; ++k) {
            if (index_map[j][k][0] != -1) {
                (*out)[index_map[j][k][0]    ] = 0.;
                (*out)[index_map[j][k][0] + 1] = 0.;
                (*out)[index_map[j][k][0] + 2] = 0.;
                (*out)[index_map[j][k][nvox[2]]    ] = 0.;
                (*out)[index_map[j][k][nvox[2]] + 1] = 0.;
                (*out)[index_map[j][k][nvox[2]] + 2] = 0.;
            }
        }
    } else {
        for (int j = vertex_plane_begin; j != vertex_plane_end; ++j)
        for (int k = 0; k != nvox[1]+1; ++k) {
            if (index_map[j][k][0] != -1) {
                (*out)[index_map[j][k][0]    ] += (*out)[index_map[j][k][nvox[2]]    ];
                (*out)[index_map[j][k][0] + 1] += (*out)[index_map[j][k][nvox[2]] + 1];
                (*out)[index_map[j][k][0] + 2] += (*out)[index_map[j][k][nvox[2]] + 2];
                (*out)[index_map[j][k][nvox[2]]    ] = (*out)[index_map[j][k][0]    ];
                (*out)[index_map[j][k][nvox[2]] + 1] = (*out)[index_map[j][k][0] + 1];
                (*out)[index_map[j][k][nvox[2]] + 2] = (*out)[index_map[j][k][0] + 2];
            }
        }
    }
}

// return the energy stored in the body, in the current state.
static
double energy_density_stored () {
    prod_with_sysmat (&btempo, solguess);
    double ret = parallel_inner_prod (solguess, btempo);
    ret /= num_total_vox;
    return ret;
}

// compute the final stencil from the elastic constants
// lambda is .5 for us, mu = .75
// energy for a compressed cube (hydrostatic mode) is 3!  
// E = 3 lambda + 2mu
// this is why BULKMOD_MODE carries another 1./sqrt(3)
static
void initialise_stencil_for_phase (int phase, double lambda, double mu) {
    // make room
    if ((int)final_stencil.size () <= phase)
        final_stencil.resize (phase+1, matrix ());

    matrix lama_stencil, muma_stencil;
    // load the stencil templates
    lama_stencil.resize (24, 24);
    muma_stencil.resize (24, 24);
    for (int i = 0; i != 24; ++i)
    for (int j = 0; j != 24; ++j) {
        lama_stencil(i,j) = lama_isotropic[24*i+j];
        muma_stencil(i,j) = muma_isotropic[24*i+j];
    }
    // the matrices carry a factor two too much atm.
    lama_stencil *= .5;
    muma_stencil *= .5;
    // linear combination gives final stencil
    final_stencil[phase] = lambda * lama_stencil;     // lambda part
    final_stencil[phase] += (2.*mu) * muma_stencil;   // 2*mu part
    // a few more temps
    stencil_temp.resize (24);
    stencil_temp2.resize (24);
}

// compute the greatest change in displacement w.r.t the initial disp field.
static
mpiutil::range max_distortion_of_single_voxel () {
    mpiutil::range ret;

    for (int i = vertex_plane_begin+1; i != vertex_plane_end; ++i)
    for (int j = 1; j != nvox[1]+1; ++j)
    for (int k = 1; k != nvox[2]+1; ++k) {
        int selfindex;
        if ((selfindex = index_map[i][j][k]) != -1) {

            int otherindex;
            if ((otherindex = index_map[i-1][j][k]) != -1)
                ret.add (1. + solguess[selfindex] - solguess[otherindex]);
            if ((otherindex = index_map[i][j-1][k]) != -1)
                ret.add (1. + solguess[selfindex+1] - solguess[otherindex+1]);
            if ((otherindex = index_map[i][j][k-1]) != -1)
                ret.add (1. + solguess[selfindex+2] - solguess[otherindex+2]);
        }
    }

    ret = all_reduce (*world, ret, mpiutil::range_reduce ());
    return ret;
}

// hydrostatic mode
matrix BULKMOD_MODE;
// cubic shear modes
matrix SHEARMOD_MODE, CUB2ND_MODE;
// Voigt (or rather, Cowin-Mehrabadi) modes
// corresponding to the columns of the CM matrix
matrix XX_MODE, YY_MODE, ZZ_MODE;
matrix XY_MODE, ZX_MODE, YZ_MODE;
const double STRAIN_SCALE = 1;

static
void modes_initialise () {
    BULKMOD_MODE = zero_matrix (3, 3);
    BULKMOD_MODE(0,0) = -STRAIN_SCALE / 3.;
    BULKMOD_MODE(1,1) = -STRAIN_SCALE / 3.;
    BULKMOD_MODE(2,2) = -STRAIN_SCALE / 3.;
    SHEARMOD_MODE = zero_matrix (3, 3);
    SHEARMOD_MODE(0,0) = -STRAIN_SCALE / 2.;
    SHEARMOD_MODE(2,2) = +STRAIN_SCALE / 2.;
    CUB2ND_MODE = zero_matrix (3, 3);
    CUB2ND_MODE(0,2) = -STRAIN_SCALE / 2.;
    CUB2ND_MODE(2,0) = -STRAIN_SCALE / 2.;
    XX_MODE = zero_matrix (3, 3);
    XX_MODE(0,0) = -STRAIN_SCALE;
    YY_MODE = zero_matrix (3, 3);
    YY_MODE(1,1) = -STRAIN_SCALE;
    ZZ_MODE = zero_matrix (3, 3);
    ZZ_MODE(2,2) = -STRAIN_SCALE;
    ZX_MODE = zero_matrix (3, 3);
    ZX_MODE(0,2) = -STRAIN_SCALE / sqrt (2);
    ZX_MODE(2,0) = -STRAIN_SCALE / sqrt (2);
    XY_MODE = zero_matrix (3, 3);
    XY_MODE(0,1) = -STRAIN_SCALE / sqrt (2);
    XY_MODE(1,0) = -STRAIN_SCALE / sqrt (2);
    YZ_MODE = zero_matrix (3, 3);
    YZ_MODE(1,2) = -STRAIN_SCALE / sqrt (2);
    YZ_MODE(2,1) = -STRAIN_SCALE / sqrt (2);
}

// load a subset of the data into a permanent buffer.
// loading the whole structure, even if only on a single node,
// may run in memory problems.
void structure_initialise (structure_reader *reader, int dof_per_vertex) {

    reader->read_header (nvox);
    assert (nvox[0] >= 1);
    assert (nvox[1] >= 1);
    assert (nvox[2] >= 1);

    // set up a few helper variables
    num_total_vox  = nvox[0] * nvox[1] * nvox[2];
    num_hull_vox   = num_total_vox - (nvox[0]-2)*(nvox[1]-2)*(nvox[2]-2);
    max_extent_vox = max3 (nvox[0], nvox[1], nvox[2]);

    // distribute the job among nodes.
    // we make a static and proably suboptimal allocation.
    mpiutil::distribute_work (&voxel_plane_begin, &voxel_plane_end,
                              nvox[0], world->size (), world->rank ());
    vertex_plane_begin = voxel_plane_begin;
    vertex_plane_end = voxel_plane_end + 1;
    vertex_accable_end = voxel_plane_end;
    if (voxel_plane_end == nvox[0])
        vertex_accable_end += 1;

    // if subset is too thin, certain assumptions in the sync code are no longer valid.
    int layers_masternode =
        all_reduce (*world, voxel_plane_end-voxel_plane_begin, mpi::minimum <int> ());
    if (layers_masternode < 2)
        die ("SUBSET IS TOO THIN (TOO MANY CPUS)\n");

    // init sync buffers  -- these are sized down later on in merge_overlapping_dof
    set_up_zero (&upper_slice_r, (nvox[1]+1)*(nvox[2]+1)*dof_per_vertex);
    set_up_zero (&lower_slice_r, (nvox[1]+1)*(nvox[2]+1)*dof_per_vertex);
    set_up_zero (&lower_slice,   (nvox[1]+1)*(nvox[2]+1)*dof_per_vertex);
    set_up_zero (&upper_slice,   (nvox[1]+1)*(nvox[2]+1)*dof_per_vertex);

    // start background thread  FIXME we never get rid of it.
    bkgd::create ();

    structure.resize (boost::extents
        [index_map_range_t (voxel_plane_begin-1, voxel_plane_end+1)]
        [nvox[1]][nvox[2]]);

    if (voxel_plane_begin == 0) {
        for (int i = voxel_plane_begin; i != voxel_plane_end+1; ++i)
            reader->readplane (i, &structure, i);
        reader->readplane (nvox[0]-1, &structure, voxel_plane_begin-1);
    } else if (voxel_plane_end == nvox[0]) {
        reader->readplane (0, &structure, voxel_plane_end);
        for (int i = voxel_plane_begin-1; i != voxel_plane_end; ++i)
            reader->readplane (i, &structure, i);
    } else {
        for (int i = voxel_plane_begin-1; i != voxel_plane_end+1; ++i)
            reader->readplane (i, &structure, i);
    }
    
    // count the total degrees of freedom
    int_dof = 0;
    index_map.resize (boost::extents
        [index_map_range_t (vertex_plane_begin, vertex_plane_end)]
        [nvox[1]+1][nvox[2]+1]);

    for (int i = vertex_plane_begin; i != vertex_plane_end; ++i) {

        for (int j = 0; j != nvox[1]+1; ++j)
        for (int k = 0; k != nvox[2]+1; ++k) {

            // if any adjacent voxel contains matter...
            if (   structure[i]            [j%nvox[1]]            [k%nvox[2]]
                || structure[i]            [j%nvox[1]]            [(k+nvox[2]-1)%nvox[2]]
                || structure[i]            [(j+nvox[1]-1)%nvox[1]][k%nvox[2]]
                || structure[i]            [(j+nvox[1]-1)%nvox[1]][(k+nvox[2]-1)%nvox[2]]
                || structure[(i-1)][j%nvox[1]]            [k%nvox[2]]
                || structure[(i-1)][j%nvox[1]]            [(k+nvox[2]-1)%nvox[2]]
                || structure[(i-1)][(j+nvox[1]-1)%nvox[1]][k%nvox[2]]
                || structure[(i-1)][(j+nvox[1]-1)%nvox[1]][(k+nvox[2]-1)%nvox[2]]) {

                // allocate three degrees of freedom
                index_map[i][j][k] = int_dof;
                int_dof += dof_per_vertex;
            } else {
                // else don't
                index_map[i][j][k] = -1;
            }
        }

        if (i == vertex_accable_end-1)
            accountable_dof = int_dof;
    }
}

static double calc_volfrac_of_phase (int phase) {
    ptrdiff_t cnt = 0;

    for (int i = voxel_plane_begin; i != voxel_plane_end; ++i)
    for (int j = 0; j != nvox[1]; ++j)
    for (int k = 0; k != nvox[2]; ++k)
        cnt += (structure[i][j][k] == phase);

    cnt = mpi::all_reduce (*world, cnt, std::plus <ptrdiff_t> ());

    return 100. * cnt / nvox[0] / nvox[1] / nvox[2];
}

static void report_volfracs () {
    fprintf (stderr, "[report_volfracs] volume fractions of the phases:\n");
    for (int i = 0; i != (int)final_stencil.size (); ++i)
        fprintf (stderr, "      phase %i = %f%%\n", i, calc_volfrac_of_phase (i));
}

// initialise the system matrix (runs after load_structure)
// (actually, it is never stored anywhere, but still...)
// we need to count the variables, and set up the index_map.
// also, we set a sensible initial guess for the solution.
// the initial guess also contains the boundary condition.
void sysmat_initialise (const matrix &mode) {

    perf_timer pt (pt_sysmat_initialise);

    iterctr = 0;

    // this is a little weird, but matrix constants are so ugly in C++.
    modes_initialise ();

    // allocate the temporaries
    residue.resize (int_dof);
    solguess.resize (int_dof);
    lastbase.resize (int_dof);
    btempo.resize (int_dof);
    vector loc = zero_vector (3), loc2 = zero_vector (3);
    hull_weight = 0.;

    // initialise the displacement field
    for (int i = vertex_plane_begin; i != vertex_plane_end; ++i)
    for (int j = 0; j != nvox[1]+1; ++j)
    for (int k = 0; k != nvox[2]+1; ++k) {

        if (index_map[i][j][k] != -1) {

            bool hull_dof = false;
            double r = rand () / (1.+RAND_MAX) * 1e-5;
            if (i == 0 || j == 0 || k == 0 ||
                i == nvox[0] || j == nvox[1] || k == nvox[2]) {

                r = 0.;
                hull_dof = true;
            }

            loc(0) = i;
            loc(1) = j;
            loc(2) = k;
            noalias (loc2) = prod (mode, loc);

            solguess[index_map[i][j][k]    ] = loc2(0) + r;
            solguess[index_map[i][j][k] + 1] = loc2(1) + r;
            solguess[index_map[i][j][k] + 2] = loc2(2) + r;
            if (hull_dof)
                hull_weight += fabs (loc2(0)) + fabs (loc2(1)) + fabs (loc2(2));
        }
    }

    // calculate hull weight
    // the stress on the dataset hull makes enters the moduli at the end,
    // multiplied with the local displacement.
    // here we estimate how much of the hull actually contributes,
    // (this depends on the porosity and geometry of the sample)
    // since later on, for reasons of simplicity, we only have a single real
    // number -- norm_inf(residue) -- to guess how wrong the stress field is.
    hull_weight = all_reduce (*world, hull_weight, std::plus <double> ());

    const double MB = 1024.*1024;
    const int ntempo = 4;

    report_volfracs ();

    fprintf (stderr, "ErPEL version %s has initialised successfully.\n",
             erpel_revision_number ().c_str ());
    fprintf (stderr, "[erpel::sysmat_initialise] memory usage (node %i)\n",
             world->rank ());
    fprintf (stderr, "        %16s: %.6f MB\n",
             "index map", mem_use (index_map) / MB);
    double total = mem_use (index_map) / MB;
    fprintf (stderr, "        %16s: %.6f MB\n",
             "structure", mem_use (structure) / MB);
    total += mem_use (structure) / MB;
    fprintf (stderr, "        %16s: %.6f MB\n",
             "solution", mem_use (solguess) / MB);
    fprintf (stderr, "        %16s: %.6f MB\n",
             "cg temps", mem_use (solguess) / MB * ntempo);
    total += mem_use (solguess)*ntempo / MB;
    fprintf (stderr, "        %16s: %.6f MB\n",
             "sync temps", mem_use (upper_slice_r) / MB * 4);
    total += mem_use (upper_slice_r) / MB * 4;
    fprintf (stderr, "        ------------------------------------------"
                     "---------\n");
    fprintf (stderr, "        %16s: %.6f MB\n",
             "total", total);
    fprintf (stderr, "\n");
    fprintf (stderr, "        %16s: %i\n",
             "local dofs", int_dof);
    fprintf (stderr, "        %16s: %i\n",
             "acc'ables", accountable_dof);
    fprintf (stderr, "        %16s: [%4i; %4i)\n",
             "voxel planes",
             voxel_plane_begin, voxel_plane_end);
    fprintf (stderr, "        %16s: [%4i; %4i)      %i\n",
             "vertex planes",
             vertex_plane_begin, vertex_plane_end, vertex_accable_end);
    fprintf (stderr, "        %16s: %f, .../%i = %f\n",
             "hull weight", hull_weight, num_total_vox, hull_weight/num_total_vox);

    mpiutil::range dist_range = max_distortion_of_single_voxel ();
    fprintf (stderr, "  %22s: %f...%f\n",
             "voxel distortion range", dist_range.lower, dist_range.upper);
}

static
bool cg_single_step (bool force_reinit) {
    bool ret = false;

    if (iterctr % max_sequential_cg == 0 || force_reinit) {
        // (re-)initialise the conjugate gradient solver
        fprintf (stderr, "[erpel::cg_single_step] re-initializing cg solver, iteration = %i\n", iterctr);
        prod_with_sysmat_enforce_bounds (&residue, solguess);
        residue *= -1.;
        lastbase = residue;
        prev_res_sq = parallel_inner_prod (residue, residue);
    }
    prod_with_sysmat_enforce_bounds (&btempo, lastbase);
    double alpha = prev_res_sq / parallel_inner_prod (lastbase, btempo);
    add_to_vector (&solguess, alpha, lastbase);
    add_to_vector (&residue, -alpha, btempo);

    if (iterctr % cvg_log_steps == 0) {
        double last_step_size = fabs (alpha) * parallel_norm_inf (lastbase);
        double res_norm = parallel_norm_inf (residue);
        double dens = energy_density_stored ();
        mpiutil::range dist_range = max_distortion_of_single_voxel ();

        // x is about 3 for a cubic system.
        double x = num_hull_vox * max_extent_vox / num_total_vox;
        x = hull_weight / num_total_vox;
        // we estimate the error in the displacements by last_step_size
        double estimated_error = res_norm * (x + last_step_size);

        if (cvg_log) {
            // this is implicitly done on master node only
            fprintf (cvg_log, "%4i %22.15e %22.15e %22.15e %22.15e\n",
                iterctr, dens, res_norm, last_step_size, estimated_error);
        }

        if (dist_log) {
            // this is implicitly done on master node only
            fprintf (dist_log, "%4i %22.15e %22.15e\n",
                iterctr, dist_range.lower, dist_range.upper);
        }

        // if we're better than our acceptance criterion, quit.
        if (estimated_error < cvg_threshold)
            ret = true;
    }

    double new_res_sq = parallel_inner_prod (residue, residue);
    double beta = new_res_sq / prev_res_sq;
    for (int i = 0; i != int_dof; ++i)
        lastbase(i) = lastbase(i)*beta + residue(i);
    prev_res_sq = new_res_sq;
    ++iterctr;
    return ret;
}

template <typename FIELD>
void produce_field (multi_array <double, 3> *dst, FIELD f, double FILLER = NAN) {

    double average = 0.;

    if_master_node
        dst->resize (boost::extents[nvox[0]][nvox[1]][nvox[2]]);

    multi_array <double, 2> slice;
    slice.resize (boost::extents[nvox[1]][nvox[2]]);

    // we're sending around a whole lot of badly tagged packages,
    // so there better be a barrier
    world->barrier ();
    for (int i = 0; i != nvox[0]; ++i) {

        // produce a slice of the field
        if (i >= voxel_plane_begin && i < voxel_plane_end) {
            average += f (&slice, i, FILLER);
            if_nonmaster_node
                world->send (0, i, slice);
        } else {
            if_master_node
                world->recv (mpi::any_source, i, slice);
        }

        // copy the slice into the final array
        if_master_node
        for (int j = 0; j != nvox[1]; ++j)
        for (int k = 0; k != nvox[2]; ++k)
            (*dst)[i][j][k] = slice[j][k];
    }

    // reduce average
    double average2;
    reduce (*world, average, average2, std::plus <double> (), mpiutil::MASTER_NODE);
    average = average2;

    if_master_node 
        fprintf (stderr, "[produce_field] average = %.6f\n",
                 average/num_total_vox);
}

static
double energy_density_field (multi_array <double, 2> *slice,
                             int i, double FILLER) {
    double ret = 0.;

    for (int j = 0; j != nvox[1]; ++j)
    for (int k = 0; k != nvox[2]; ++k) {
        int ph = structure[i][j][k];
        if (ph) {
            // load the vertices into stencil_temp
            stencil_temp[ 0] = solguess[index_map[i][j][k]    ];
            stencil_temp[ 1] = solguess[index_map[i][j][k] + 1];
            stencil_temp[ 2] = solguess[index_map[i][j][k] + 2];
            stencil_temp[ 3] = solguess[index_map[i][j][k+1]    ];
            stencil_temp[ 4] = solguess[index_map[i][j][k+1] + 1];
            stencil_temp[ 5] = solguess[index_map[i][j][k+1] + 2];
            stencil_temp[ 6] = solguess[index_map[i][j+1][k]    ];
            stencil_temp[ 7] = solguess[index_map[i][j+1][k] + 1];
            stencil_temp[ 8] = solguess[index_map[i][j+1][k] + 2];
            stencil_temp[ 9] = solguess[index_map[i][j+1][k+1]    ];
            stencil_temp[10] = solguess[index_map[i][j+1][k+1] + 1];
            stencil_temp[11] = solguess[index_map[i][j+1][k+1] + 2];
            stencil_temp[12] = solguess[index_map[i+1][j][k]    ];
            stencil_temp[13] = solguess[index_map[i+1][j][k] + 1];
            stencil_temp[14] = solguess[index_map[i+1][j][k] + 2];
            stencil_temp[15] = solguess[index_map[i+1][j][k+1]    ];
            stencil_temp[16] = solguess[index_map[i+1][j][k+1] + 1];
            stencil_temp[17] = solguess[index_map[i+1][j][k+1] + 2];
            stencil_temp[18] = solguess[index_map[i+1][j+1][k]    ];
            stencil_temp[19] = solguess[index_map[i+1][j+1][k] + 1];
            stencil_temp[20] = solguess[index_map[i+1][j+1][k] + 2];
            stencil_temp[21] = solguess[index_map[i+1][j+1][k+1]    ];
            stencil_temp[22] = solguess[index_map[i+1][j+1][k+1] + 1];
            stencil_temp[23] = solguess[index_map[i+1][j+1][k+1] + 2];
            noalias (stencil_temp2) = prod (final_stencil[ph], stencil_temp);
            (*slice)[j][k] = inner_prod (stencil_temp, stencil_temp2);
            ret += (*slice)[j][k];
        } else {
            (*slice)[j][k] = FILLER;
        }
    }

    return ret;
}

template <int component>
double displacement_field (multi_array <double, 2> *slice,
                           int i, double FILLER) {
    double ret = 0.;

    fprintf (stderr, "proc. slice = %i on node %i\n", i, world->rank ());
    for (int j = 0; j != nvox[1]; ++j)
    for (int k = 0; k != nvox[2]; ++k) {
        int x = index_map[i][j][k];
        if (x != -1) {
            // load the vertices into stencil_temp
            (*slice)[j][k] = solguess[x+component];
            ret += (*slice)[j][k];
        } else {
            (*slice)[j][k] = FILLER;
        }
    }

    return ret;
}

    void produce_energy_density_field (multi_array <double, 3> *dst) {
        produce_field (dst, energy_density_field);
    }
    void produce_displacement_field_0 (multi_array <double, 3> *dst) {
        produce_field (dst, &displacement_field <0>);
    }
    void produce_displacement_field_1 (multi_array <double, 3> *dst) {
        produce_field (dst, &displacement_field <1>);
    }
    void produce_displacement_field_2 (multi_array <double, 3> *dst) {
        produce_field (dst, &displacement_field <2>);
    }

// write out the hull field as a single vector.
// * it is easy to reconstruct the effective stiffness tensor using dot products
//   from these vectors.
// * the full result is only available on node zero.
// * hullfiles are not comparable across datasets, or even runs of the code
//   with different parameters (esp. number of cpus)
// * we also record the maximum ignored value, which should be close to zero
//   in converged state.
static
void extract_hull_from (std::vector <double> *full_hull_, double *max_ign_element_,
                        const vector &data) {

    std::vector <double> hull_segment;
    double max_ign_element = 0.;

    // top plane
    int i = vertex_plane_begin;
    if (vertex_plane_begin == 0)
        for (int j = 0; j != nvox[1]+1; ++j)
        for (int k = 0; k != nvox[2]+1; ++k) {
            int x = index_map[i][j][k];
            if (x != -1) {
                for (int c = 0; c != 3; ++c)
                    hull_segment.push_back (data[x+c]);
            }
        }

    // walls
    int imin = vertex_plane_begin;
    int imax = vertex_accable_end;
    if (imin == 0)
        ++imin;
    if (imax == nvox[0]+1)
        --imax;
    for (i = imin; i != imax; ++i) {
        for (int j = 0; j != nvox[1]+1; ++j) {
            int k = 0;
            int x = index_map[i][j][k];
            if (x != -1) {
                for (int c = 0; c != 3; ++c)
                    hull_segment.push_back (data[x+c]);
            }
            k = nvox[2];
            x = index_map[i][j][k];
            if (x != -1) {
                for (int c = 0; c != 3; ++c)
                    hull_segment.push_back (data[x+c]);
            }
        }
        for (int k = 1; k != nvox[2]; ++k) {
            int j = 0;
            int x = index_map[i][j][k];
            if (x != -1) {
                for (int c = 0; c != 3; ++c)
                    hull_segment.push_back (data[x+c]);
            }
            j = nvox[1];
            x = index_map[i][j][k];
            if (x != -1) {
                for (int c = 0; c != 3; ++c)
                    hull_segment.push_back (data[x+c]);
            }
        }
        for (int j = 1; j != nvox[1]; ++j)
        for (int k = 1; k != nvox[2]; ++k) {
            int x = index_map[i][j][k];
            if (x != -1) {
                for (int c = 0; c != 3; ++c)
                    max_ign_element = std::max (
                        max_ign_element,
                        fabs (data[x+c]));
            }
        }
    }

    // bottom plane
    i = nvox[0];
    if (vertex_accable_end == nvox[0]+1)
        for (int j = 0; j != nvox[1]+1; ++j)
        for (int k = 0; k != nvox[2]+1; ++k) {
            int x = index_map[i][j][k];
            if (x != -1) {
                for (int c = 0; c != 3; ++c)
                    hull_segment.push_back (data[x+c]);
            }
        }

    reduce (*world, max_ign_element, *max_ign_element_, mpi::maximum <double> (), mpiutil::MASTER_NODE);
    reduce (*world, hull_segment,    *full_hull_,       mpiutil::concatenate (), mpiutil::MASTER_NODE);
    fprintf (stderr, "[erpel::extract_hull_from] hull size = %i -> %i\n",
             int (hull_segment.size ()), int (full_hull_->size ()));
}

bool cg_loop () {

    // educated guess, more sequential CG is generally more efficient,
    // though errors may accumulate.
    max_sequential_cg = itermax/10;
    int acceptctr = 0;

    for (;;) {
        if (cg_single_step (false)) {
            acceptctr++;
            if (acceptctr >= acceptcrit)
                return true;
            cg_single_step (true);
            continue;
        }
        if (iterctr >= itermax)
            return false;
    }
}

void set_phase_lame_shear_moduli (int phase, double lame, double shear) {
    double bulk = bulk_from_lame_shear (lame, shear);
    double poisson = poisson_from_lame_shear (lame, shear);
    fprintf (stderr, "phase %2i:  "
             "k = %8.5f  g = %8.5f  l = %8.5f  n = %8.5f  # %10s\n",
             phase,
             bulk, shear, lame, poisson,
             phase_name[phase].c_str ());
    initialise_stencil_for_phase (phase, lame, shear);
}

void set_phase_bulk_shear_moduli (int phase, double bulk, double shear) {
    double lame = lame_from_bulk_shear (bulk, shear);
    set_phase_lame_shear_moduli (phase, lame, shear);
}

void set_phase_name (int phase, std::string name) {
    phase_name[phase] = name;
}

void open_convergence_log (std::string filename) {
    if (cvg_log) fclose (cvg_log);
    cvg_log = fopen (filename.c_str (), "wt");
    fprintf (cvg_log, "# steps energy_density norm_inf_residue norm_inf_last_step error_estimate\n");
    // switch to unbuffered. don't write a whole lot of shit into it.
    setbuf (cvg_log, 0);
}

void open_distortion_log (std::string filename) {
    if (dist_log) fclose (dist_log);
    dist_log = fopen (filename.c_str (), "wt");
    if (!dist_log) {
        fprintf (stderr, "cannot open %s\n", filename.c_str ());
    }
    fprintf (dist_log, "# steps min_side max_side\n");
    setbuf (dist_log, 0);
}

void close_convergence_log () {
    if (cvg_log) fclose (cvg_log);
    cvg_log = 0;
}

void close_distortion_log () {
    if (dist_log) fclose (dist_log);
    dist_log = 0;
}

struct hull_info {
    std::vector <double> stress, strain;
    double max_ignored_value;
};

static
std::vector <hull_info> stored_hulls;

void store_hull (int number) {
    if (number <= (int)stored_hulls.size ()) {
        stored_hulls.reserve (10);
        stored_hulls.resize (number+1);
    }

    prod_with_sysmat (&btempo, solguess);
    double valain = 0., valstr = 0.;
    extract_hull_from (&stored_hulls[number].stress, &valstr, btempo);
    extract_hull_from (&stored_hulls[number].strain, &valain, solguess);
    stored_hulls[number].max_ignored_value =
        std::max (valain, valstr);

    fprintf (stderr, "[erpel::store_hull] sizes of hulls:\n");
    for (int i = 0; i != (int)stored_hulls.size (); ++i)
        fprintf (stderr, "   hull[%i] = %i\n", i, (int)stored_hulls[i].stress.size ());
}

double hull_product (int a, int b) {
    return inner_prod (stored_hulls.at (a).stress,
                       stored_hulls.at (b).strain) / num_total_vox;
}

void set_moduli_precision (double mprec) {
    assert (mprec > 0.);
    assert (mprec < 1e-3);
    cvg_threshold = mprec/7;
}

void set_problem (std::string problem) {
    if (problem == "periodic_linear_elasticity")
        fix_hull_displacements_z = fix_hull_displacements_y = fix_hull_displacements_x = false;
    else if (problem == "fixed_displacement_linear_elasticity")
        fix_hull_displacements_z = fix_hull_displacements_y = fix_hull_displacements_x = true;
    else
        die ("[erpel::set_problem] unsupported problem \"%s\"", problem.c_str ());
}

void dump_perf_counters (FILE *fp) {
    perf_timer_class::dump_all (fp);
}

} // namespace erpel
