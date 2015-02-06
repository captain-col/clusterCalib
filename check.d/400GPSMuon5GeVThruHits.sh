#! /bin/sh

# Enable error trapping.
set -e
set -x

INPUT=elecsim-200GPSMuon5GeVThruDigits.root
if [ ! -f $INPUT ]; then
    echo "Missing file: " ${INPUT}
    echo MISSING INPUT
    exit
fi

OUTPUT=calib-400GPSMuon5GeVThruHits.root
if [ -f $OUTPUT ]; then
    rm $OUTPUT
fi

CLUSTERCALIB.exe -o $OUTPUT $INPUT

