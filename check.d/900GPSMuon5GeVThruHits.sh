#!/bin/bash

INPUT=calib-400GPSMuon5GeVThruHits.root

# Check that the pmt hit container was created
if ! (dump-event.exe ${INPUT} | grep 'THitSelection.*pmt'); then
    echo FAILED Missing pmt hit container in ${INPUT}
    exit 1
fi

# Check the number of PMTs created
if [ $(dump-event.exe -O hits ${INPUT} | grep 'TDataHit.*MC-P-' | wc -l) -lt 5 ]; then
    echo FAILED Missing pmt hits in ${INPUT}
    exit 1
fi

# Check that the drift digit container was created
if ! (dump-event.exe ${INPUT} | grep 'THitSelection.*drift'); then
    echo FAILED Missing drift digit container in ${INPUT}
    exit 1
fi

# Check the number of X wires created
if [ $(dump-event.exe -O hits ${INPUT} | grep 'TDataHit.*MC-W-X-' | wc -l) -lt 10 ]; then
    echo FAILED Missing X wire hits in ${INPUT}
    exit 1
fi

# Check the number of U wires created
if [ $(dump-event.exe -O hits ${INPUT} | grep 'TDataHit.*MC-W-U-' | wc -l) -lt 10 ]; then
    echo FAILED Missing U wire hits in ${INPUT}
    exit 1
fi

# Check the number of V wires created
if [ $(dump-event.exe -O hits ${INPUT} | grep 'TDataHit.*MC-W-V-' | wc -l) -lt 10 ]; then
    echo FAILED Missing V wire hits in ${INPUT}
    exit 1
fi


echo SUCCESS
