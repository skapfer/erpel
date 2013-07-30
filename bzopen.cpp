#include "shared.hpp"
#include <stdio.h>
#include <string>
#include <assert.h>
#include <errno.h>
#include <unistd.h>
#include <fcntl.h>

enum {
    BZOPEN_UNREADABLE = -1,
    BZOPEN_UNKNOWN_MAGIC = -2,
    BZOPEN_UNAVAIL = -3
};

static
int zcat (const char *filename, const char *program, int we_write) {
    /* this function duplicates most functionality of popen
     * for fun & profit.  we should really link zlib instead.
     * but it's fun :-)
     */
    errno = 0;
    int file_fd = open (filename, we_write ? O_WRONLY : O_RDONLY);
    if (file_fd < 0) {
        /* user should check errno */
        return BZOPEN_UNREADABLE;
    }
    /* make pipe */
    int pip_fd[2];
    if (pipe (pip_fd) == -1) {
        perror ("pipe");
        return BZOPEN_UNAVAIL;
    }
    int rpip = pip_fd[0];
    int wpip = pip_fd[1];
    switch (fork ()) {
    case -1:
        /* bail out, fork failed */
        close (file_fd);
        perror ("fork");
        return BZOPEN_UNAVAIL;
    case 0:
        /* plumbing */
        if (dup2 (we_write ? rpip : file_fd, 0) == -1) {
            perror ("dup2 (0)");
            abort ();
        }
        if (dup2 (we_write ? file_fd : wpip, 1) == -1) {
            perror ("dup2 (1)");
            abort ();
        }
        close (we_write ? wpip : rpip);
        /* exec zcat */
        execlp (program, program, NULL);
        perror ("execl");
        abort ();
    default:
        close (we_write ? rpip : wpip);
        close (file_fd);
        return we_write ? wpip : rpip;
    }
}

FILE *bzopen (std::string filename, std::string mode) {
    if (ends_with (filename, ".bz2")) {
        int fd;
        if (mode == "r") {
            fd = zcat (filename.c_str (), "bzcat", 0);
            if (fd < 0)
                die ("unable to use bzcat in order to open %s", filename.c_str ());
        } else if (mode == "w") {
            fd = zcat (filename.c_str (), "bzip2", 1);
            if (fd < 0)
                die ("unable to use bzip2 in order to open %s", filename.c_str ());
        } else {
            die ("weird open mode %s", mode.c_str ());
        }
        return fdopen (fd, mode.c_str ());
    } else {
        return fopen (filename.c_str (), mode.c_str ());
    }
}
