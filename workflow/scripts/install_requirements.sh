#!/bin/bash

if [ -f "requirements_met.flag" ]; then
    echo "Requirements are already met. Skipping execution."
    exit 0
fi

if [ ! -d "PhaBOX" ]; then
    git clone https://github.com/KennthShang/PhaBOX.git
fi
cd PhaBOX
python  -m pip install .
cd ..
rm -rf PhaBOX

wget https://github.com/KennthShang/PhaBOX/releases/download/v2/phabox_db_v2.zip
unzip phabox_db_v2.zip > /dev/null

touch requirements_met.flag
