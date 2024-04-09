#!/bin/bash

if [ -f "requirements_met.flag" ]; then
    echo "Requirements are already met. Skipping execution."
    exit 0
fi

if [ ! -d "PHABOX" ]; then
    git clone https://github.com/KennthShang/PhaBOX.git
fi

pip install gdown
gdown  --id 1hjACPsIOqqcS5emGaduYvYrCzrIpt2_9

gdown  --id 1E94ii3Q0O8ZBm7UsyDT_n06YekNtfV20

unzip phagesuite_database.zip  > /dev/null
unzip phagesuite_parameters.zip  > /dev/null


if [ -n "$CONDA_PREFIX" ]; then
    cp PhaBOX/blastxml_to_tabular.py "$CONDA_PREFIX/bin/"
    chmod 777 "$CONDA_PREFIX/bin/blastxml_to_tabular.py"
fi

touch requirements_met.flag
