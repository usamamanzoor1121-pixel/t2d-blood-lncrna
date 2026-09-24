#!/usr/bin/env bash
set -e
export JAVA_HOME=$HOME/jdk17
export PATH=$JAVA_HOME/bin:$PATH
cd /mnt/d/T2D_lncRNA_Standalone
source venv/bin/activate
export MLFLOW_ALLOW_FILE_STORE=true
./nextflow run main.nf
