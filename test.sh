#!/bin/bash
FLAG="$1"
if [[ ! "$FLAG" ]]
then
	FLAG="all"
fi

set -o errexit
rm -f MD
#./format.sh
#http://www.unknownroad.com/rtfm/gdbtut/gdbsegfault.html
make "$FLAG"
./MD
