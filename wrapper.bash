#!/bin/bash -l
# wrapper script containing the ugly set-up bits for various systems.
ERPEL=$0.bin

woody () {
    mpirun -mpich -- "$ERPEL" "$@"
}

tinyblue () {
    mpirun -mpich -npernode 8 -pin "0 1 2 3 4 5 6 7" -- "$ERPEL" "$@"
}

lima () {
    /apps/rrze/bin/mpirun_rrze-intelmpd -intelmpd -pin "0,12 1,13 2,14 3,15 4,16 5,17 6,18 7,19 8,20 9,21 10,22 11,23" -- "$ERPEL" "$@"
}

theo1adhoc () {
    ${MPIRUN:-mpirun -np 4} -- "$ERPEL" "$@"
}

case "$(hostname)" in
klein|noether|ries|duerer|van-gogh|da-vinci)
    mpiwrap=theo1adhoc
    echo '[Erpel] Ad-hoc configuration for Theo1 systems' >&2
    # module system is broken in noninteractive shells.
    . /etc/mod_system/tcl/init/bash
    module load boost
    ;;
# RRZE systems
tb*)
    mpiwrap=tinyblue
    echo '[Erpel] Ignore the following warnings.'
    module load mpich/p4-intel64 intel64/11.1.064 2>&1
    module load acml-intel64/4.1.0_mp 2>&1
    ;;
w*)
    mpiwrap=woody
    echo '[Erpel] Ignore the following warnings.'
    module load mpich/p4-intel64 intel64/11.1.064 2>&1
    module load acml-intel64/4.1.0_mp 2>&1
    ;;
l*)
    mpiwrap=lima
    . /apps/rrze/etc/use-rrze-modules.sh
    module load boost/1.44.0-intel11.1-impi4.0 intel64/11.1.073 2>&1
    ;;
# FIXME add cfg for townsend
*)
    echo "[Erpel] WARNING no configuration for this host" >&2
    exit 1
esac
echo -n '[Erpel] '
module list
( while sleep 5; do cat /proc/loadavg ; done ) >loadavg &
echo '[Erpel] Executing mpirun.'
$mpiwrap "$@"
kill %1
