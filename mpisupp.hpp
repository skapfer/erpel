// a few things we need for MPI
#ifndef MPISUPP_HPP_INCLUDED
#define MPISUPP_HPP_INCLUDED

#include <string>
#include <signal.h>
#include <boost/mpi.hpp>
#include <boost/multi_array.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/binary_object.hpp>

extern boost::mpi::communicator *world;

#define if_master_node if (::world->rank () == ::mpiutil::MASTER_NODE)
#define if_nonmaster_node if (::world->rank () != ::mpiutil::MASTER_NODE)

namespace mpiutil {
    enum { MASTER_NODE = 0 };

    void initialise (int &argc, char **&argv);
    void terminate ();

    void sequentially (), sequentially_end ();

    // FILE is for error messages / status report.
    void affinity_auto_pin (FILE *);
    void affinity_report (FILE *);

    struct concatenate {
        template <typename ELEM>
        std::vector <ELEM> operator() (const std::vector <ELEM> &lhs,
                                       const std::vector <ELEM> &rhs) {
            std::vector <ELEM> ret = lhs;
            ret.insert (ret.end (), rhs.begin (), rhs.end ());
            return ret;
        }
    };

    // redirect stderr to this file. %i is replaced by node number.
    void redirect_stderr (std::string);

    void distribute_work (int *, int *, int ndata, int nnodes, int inode);

    template <typename T>
    void check_identical (const T &loc_stuff) {
        // we don't use broadcast here to avoid having to
        // use operator= (which is inefficient).
        // no really, we do it this way because boost::multi_array
        // has a broken operator=. see
        // http://article.gmane.org/gmane.comp.lib.boost.user/24705
        if_master_node {
            for (int i = 1; i != ::world->size (); ++i)
                ::world->send (i, 91, loc_stuff);
        } else {
            T the_stuff;
            ::world->recv (0, 91, the_stuff);
            if (the_stuff != loc_stuff)
                throw std::runtime_error ("mismatch in mpiutil::check_identical");
        }
    }

    struct range {
        double upper, lower;

        range () {
            upper = -INFINITY;
            lower =  INFINITY;
        }

        void add (double x) {
            lower = std::min (lower, x);
            upper = std::max (upper, x);
        }

        template <typename Archive>
        void serialize (Archive &ar, unsigned) {
            ar & lower & upper;
        }
    };

    struct range_reduce {
        range operator() (const range &lhs, const range &rhs) const {
            range ret;
            ret.lower = std::min (lhs.lower, rhs.lower);
            ret.upper = std::max (lhs.upper, rhs.upper);
            return ret;
        }
    };

    class block_sigs {
    public:
        block_sigs () {
            sigset_t allsigs;
            sigfillset (&allsigs);
            sigprocmask (SIG_SETMASK, &allsigs, &prevmask);
        }
        ~block_sigs () {
            sigprocmask (SIG_SETMASK, &prevmask, NULL);
        }
    private:
        sigset_t prevmask;
    };
}

namespace mpi = boost::mpi;

// serialisation of 2D multi_array's
// written for binary-serialisable datatypes (doubles!)
namespace boost { namespace serialization {
    template <class Archive, typename Elem>
    void save (Archive &ar, const boost::multi_array <Elem, 2> &t, unsigned int) {
        size_t n_elem = t.shape ()[0] * t.shape()[1];
        ar << t.shape ()[0] << t.shape()[1];
        ar << boost::serialization::make_binary_object ((void *)t.data(), n_elem * sizeof (Elem));
    }

    template <class Archive, typename Elem>
    void load (Archive &ar, boost::multi_array <Elem, 2> &t, unsigned int) {
        size_t nx, ny;
        ar >> nx >> ny;
        size_t n_elem = nx * ny;
        t.resize (boost::extents[nx][ny]);
        // hack
        ar >> boost::serialization::make_binary_object (t.data(), n_elem * sizeof (Elem));
    }

    template <class Archive, typename Elem>
    void save (Archive &ar, const boost::multi_array <Elem, 3> &t, unsigned int) {
        size_t n_elem = t.shape ()[0] * t.shape ()[1] * t.shape ()[2];
        ar << t.shape ()[0] << t.shape() [1] << t.shape ()[2];
        ar << boost::serialization::make_binary_object ((void *)t.data(), n_elem * sizeof (Elem));
    }

    template <class Archive, typename Elem>
    void load (Archive &ar, boost::multi_array <Elem, 3> &t, unsigned int) {
        size_t nx, ny, nz;
        ar >> nx >> ny >> nz;
        size_t n_elem = nx * ny * nz;
        t.resize (boost::extents[nx][ny][nz]);
        ar >> boost::serialization::make_binary_object (t.data(), n_elem * sizeof (Elem));
    }
} }

typedef boost::multi_array <double, 2> dummy_boost_multi_array_double_2;
BOOST_SERIALIZATION_SPLIT_FREE (dummy_boost_multi_array_double_2)
typedef boost::multi_array <char, 3> dummy_boost_multi_array_double_3;
BOOST_SERIALIZATION_SPLIT_FREE (dummy_boost_multi_array_double_3)

#endif // MPISUPP_HPP_INCLUDED
