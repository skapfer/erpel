#include "shared.hpp"
#include <map>
#include <assert.h>
#include <errno.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>

int main (int, const char **argv) {
    std::map <int, ptrdiff_t> count;
    ptrdiff_t total_sum = 0;
    int nvox[3];
    if (!argv[1]) {
        fprintf (stderr, "missing argument...\n");
        return 1;
    }
    std::string filename = argv[1];
    FILE *fp = bzopen (filename, "r");
    int nsc = fscanf (fp, "%i %i %i", &nvox[0], &nvox[1], &nvox[2]);
    if (nsc != 3)
        die ("[read_bin_file] failed to read dimensions of dataset, got %i,%i,%i", nvox[0], nvox[1], nvox[2]);
    for (int i = 0; i != nvox[0]; ++i)
    for (int j = 0; j != nvox[1]; ++j)
    for (int k = 0; k != nvox[2]; ++k) {
        int c;
        for (;;) {
            c = fgetc (fp);
            switch (c) {
            case ' ':
            case '\n':
            case '\t':
                continue;
            case '0':
            case '1': case '2': case '3': case '4':
            case '5': case '6': case '7': case '8': case '9':
                count[(char)int (c - '0')]++;
                total_sum++;
                goto next_voxel;
            case EOF:
                {
                    std::string error_message = strerror (errno);
                    die ("[read_bin_file] premature end of .bin file at %i,%i,%i\n"
                         "                error code is %s\n", i, j, k,
                         error_message.c_str ());
                }
            default:
                die ("[read_bin_file] invalid character '%c'", (char)c);
            }
        }
    next_voxel:
        continue;
    }

    std::map <int, ptrdiff_t>::const_iterator it;
    for (it = count.begin (); it != count.end (); ++it)
        fprintf (stderr, "%10.3f%%     ", 100. * it->second / total_sum);
    fprintf (stderr, "    # %s\n", filename.c_str ());
    return 0;
}
