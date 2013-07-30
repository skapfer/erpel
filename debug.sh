#!/bin/bash
echo running Erpel under debugger...
gdb ./Erpel <<EOF
run $@
bt



kill
quit
EOF
