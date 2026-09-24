#!/usr/bin/env bash
set -e
cd /mnt/d/T2D_lncRNA_Standalone
python3 -m venv venv
source venv/bin/activate
pip install --quiet --upgrade pip
pip install --quiet -r requirements.txt
echo "VENV_READY"
python3 -c "import pandas, numpy, scipy, statsmodels, matplotlib, requests; print('versions:', pandas.__version__, numpy.__version__, scipy.__version__, statsmodels.__version__, matplotlib.__version__, requests.__version__)"
