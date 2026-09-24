#!/usr/bin/env bash
set -e
cd /mnt/d/T2D_lncRNA_Standalone
python3 -m venv venv
source venv/bin/activate
pip install --quiet --upgrade pip
pip install --quiet pandas numpy scipy scikit-learn statsmodels joblib matplotlib seaborn pyyaml requests openpyxl
echo "VENV_READY"
python3 -c "import pandas, sklearn, statsmodels; print('versions:', pandas.__version__, sklearn.__version__, statsmodels.__version__)"
