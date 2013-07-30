
#include "config.hpp"
#include "mpisupp.hpp"
#include <stdio.h>
#include <stdarg.h>
#include <numeric>
#include <vector>
#include <assert.h>
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <pthread.h>  // get pthread_setaffinity_np

mpi::communicator *world;
mpi::environment *env;
typedef std::string string;

namespace {
    string cpp_sprintf (const string &msg, ...) {
        va_list va;
        va_start (va, msg);
        int bufsz = vsnprintf (0, 0, msg.c_str (), va);
        va_end (va);
        va_start (va, msg);
        string ret = string (bufsz+1, ' ');
#ifndef NDEBUG
        int bufsz2 =
#endif
        vsnprintf (&ret[0], bufsz+1, msg.c_str (), va);
        assert (bufsz == bufsz2);
        va_end (va);
        assert (ret[ret.size ()-1] == '\0');
        ret.erase (ret.size ()-1);
        return ret;
    }


#if 0   // -- unused, avoid warning
    void dumpargs (int argc, char **argv) {
        for (int i = 0; i != argc; ++i)
            fprintf (stderr, "argv[%i] = %s\n",
                    i, argv[i]);
    }
#endif

    void test_reduce_concat () {
        std::vector <double> locvec (10000);
        std::vector <double> totalvec (world->size () * 10000);
        std::vector <double> recvec;
        for (int i = 0; i != 10000; ++i)
            locvec[i] = world->rank () * 10000 + i;
        for (int i = 0; i != world->size () * 10000; ++i)
            totalvec[i] = i;
        reduce (*world, locvec, recvec, mpiutil::concatenate (), 0);
        if (world->rank () == 0) {
            if (recvec != totalvec) {
                fprintf (stderr, "reduce with mpiutil::concatenate "
                                 "does not work as expected");
                abort ();
            }
        }
    }

    std::string get_hostname () {
        char buffer[256];
        gethostname (buffer, 255);
        buffer[255] = '\0';
        return buffer;
    }
}

namespace mpiutil {

void initialise (int &argc, char **&argv) {
    //dumpargs (argc, argv);
    env = new mpi::environment (argc, argv);
    //dumpargs (argc, argv);
    world = new mpi::communicator ();

    fprintf (stderr, 
             "[mpiutil::initialise] hostname = %s, rank = %i\n",
             get_hostname ().c_str (), world->rank ());

    world->barrier ();
    test_reduce_concat ();
}

void terminate () {
    world->barrier ();
    delete world;
    world = 0;
    delete env;
    env = 0;
}

void sequentially () {
    if (world->rank () > 0) {
        // wait util we're activated
        int dummy;
        world->recv (world->rank () - 1, 0x0815, dummy);
    }
}

void sequentially_end () {
    if (world->rank () + 1 < world->size ()) {
        // activate next
        int dummy = 7;
        world->send (world->rank () + 1, 0x0815, dummy);
    }
}

#ifdef HAVE_LINUX_SETAFFINITY_NP
// CPU_COUNT is not in ancient glibcs. gna.
static int num_cpus_in_cpuset (const cpu_set_t *quux) {
    int ret = 0;
    for (int j = 0; j < CPU_SETSIZE; ++j)
        if (CPU_ISSET (j, quux))
            ++ret;
    return ret;
}

void affinity_auto_pin (FILE *fp) {
    cpu_set_t cpuset;
    const char *func_name = "[mpiutil::affinity_auto_pin]";
    int err = pthread_getaffinity_np (pthread_self (), sizeof (cpu_set_t),
                                      &cpuset);
    if (err) {
        fprintf (fp, "[%s] pthread_getaffinity_np failed %i\n",
                 func_name, err);
        return;
    }

    int cpu_count = num_cpus_in_cpuset (&cpuset);

    if (cpu_count != 1) {
        fprintf (fp, "[mpiutil::affinity_auto_pin] have more than one CPU affinity (%i)\n", cpu_count);

        // check that we run under the conditions this function
        // was designed for
        for (int i = 0; i != cpu_count; ++i)
            if (! (CPU_ISSET (i, &cpuset))) {
                fprintf (fp, "not all CPUs allocated -- "
                         "weird, bailing out\n");
                return;
            }
        if (world->size () % cpu_count) {
            fprintf (fp, "world size is not divisible by "
                     "local CPU count. weird.\n");
            return;
        }

        // hopefully, the ranks are allocated in a linear fashion
        int cpu = world->rank () % cpu_count;
        CPU_ZERO (&cpuset);
        CPU_SET (cpu, &cpuset);
        err = pthread_setaffinity_np (pthread_self (), sizeof (cpu_set_t),
                                      &cpuset);
        if (err) {
            fprintf (fp, "[%s] pthread_setaffinity_np failed %i\n",
                     func_name, err);
            return;
        }
    }
}

void affinity_report (FILE *fp) {
    cpu_set_t cpuset;
    int err = pthread_getaffinity_np (pthread_self (), sizeof (cpu_set_t), &cpuset);
    if (err)
        fprintf (fp, "pthread_getaffinity_np failed %i\n", err);
    fprintf (fp, "[mpiutil::affinity_report] thread affinity mask: ");
    for (int j = 0; j < std::min (32, CPU_SETSIZE); ++j)
        fprintf (fp, "%i", (int)!! (CPU_ISSET (j, &cpuset)));
    fprintf (fp, "\n");
}
#endif // HAVE_LINUX_SETAFFINITY_NP

void redirect_stderr (std::string filename) {
    filename = cpp_sprintf (filename, world->rank ());
    if (!freopen (filename.c_str (), "wt", stderr)) {
        perror ("error redirecting output");
        abort ();
    }
    setbuf (stderr, 0);
}

void distribute_work (int *begi, int *endi, int ndata, int nnodes, int inode) {
    std::vector <int> ndata_per_node (nnodes, ndata/nnodes);
    int missing = ndata - nnodes*ndata_per_node[0];
    for (int i = nnodes-missing; i != nnodes; ++i)
        ++ndata_per_node[i];
    std::vector <int>::iterator beg = ndata_per_node.begin ();
    assert (std::accumulate (beg, beg+nnodes, 0) == ndata);
    *begi = std::accumulate (beg, beg+inode, 0);
    *endi = std::accumulate (beg, beg+inode+1, 0);
}

} // namespace mpiutil
