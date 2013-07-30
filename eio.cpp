#include "erpel.hpp"
#include <stdio.h>
#include <errno.h>
#include <string.h>

namespace erpel {

    struct bin_reader : public structure_reader {
        bin_reader (std::string filename) {
            zerovox_shown = false;
            pos = 0;
            fp = bzopen (filename, "r");
            if (!fp)
                die ("[bin_reader::bin_reader] opening of %s failed", filename.c_str ());
        }

        virtual ~bin_reader () {
            if (fp) {
                 ignore_final_planes ();
                 ignore_trailing_text ();
                 fclose (fp);
            }
        }

        virtual void read_header (int *nvox) {
            int nsc = fscanf (fp, "%i %i %i", &nvox[2], &nvox[1], &nvox[0]);
            this->nvox[0] = nvox[0];
            this->nvox[1] = nvox[1];
            this->nvox[2] = nvox[2];
            if (nsc != 3)
                die ("[bin_reader::read_header] failed to read dimensions of dataset, got %i,%i,%i",
                     nvox[0], nvox[1], nvox[2]);
            fprintf (stderr, "[bin_reader::read_header] reading .bin, dimensions = %ix%ix%i\n",
                     this->nvox[0], this->nvox[1], this->nvox[2]);
        }

        virtual void readplane (int pln, structure_t *dst, int dst_i) {
            fprintf (stderr, "[bin_reader::readplane] plane %i\n", pln);
            do {
                assert (pln >= pos);
                for (int j = 0; j != nvox[1]; ++j) {
                    for (int k = 0; k != nvox[2]; ++k)
                        (*dst)[dst_i][j][k] = get_voxel ();
                    get_nl ();
                }
            } while (pos++ != pln);
        }

        // read planes and ignore them. this ensures that all nodes read the whole file and
        // bzip2 exits cleanly.
        void ignore_final_planes () {
            while (pos != nvox[0]) {
                for (int j = 0; j != nvox[1]; ++j) {
                    for (int k = 0; k != nvox[2]; ++k)
                        get_voxel ();
                    get_nl ();
                }
                ++pos;
            }
        }

        void ignore_trailing_text () {
            size_t ret = 0u;
            for (;;) {
                if (fgetc (fp) == EOF) {
                    if (ret != 0u) {
                        fprintf (stderr, "[bin_reader::ignore_trailing_text] %u "
                                 "bytes non-structure information in input file\n",
                                 (unsigned)ret);
                    }
                    return;
                } else {
                    ret++;
                }
            }
        }

        phase_index get_voxel () {
            for (;;) {
                int c = fgetc (fp);
                switch (c) {
                case ' ':
                case '\n':
                case '\t':
                    continue;
                case '0':
                    if (!zerovox_shown)
                        fprintf (stderr, "[read_bin_file] WARNING .bin file contains voxels of phase 0!\n"),
                        zerovox_shown = true;
                case '1': case '2': case '3': case '4':
                case '5': case '6': case '7': case '8': case '9':
                    return int (c - '0');
                case EOF:
                    {
                        std::string error_message = strerror (errno);
                        die ("[read_bin_file] premature end of .bin file\n"
                             "                error code is %s\n",
                                              error_message.c_str ());
                    }
                default:
                    die ("[read_bin_file] invalid character '%c'", (char)c);
                }
            }
        }

        void get_nl () {
            for (;;) {
                int c = fgetc (fp);
                switch (c) {
                case ' ':
                case '\t':
                    continue;
                case '\n':
                    return;
                default:
                    fprintf (stderr, "[read_bin_file] WARNING scanline in .bin is not terminated by a linefeed. probably there's something fishy.\n");
                    return;
                }
            }
        }

        bool zerovox_shown;
        FILE *fp;
        int nvox[3];
        int pos;
    };

    structure_reader *open_bin_file (std::string filename) {
        return new bin_reader (filename);
    }

    structure_reader::~structure_reader () {
    }
}
