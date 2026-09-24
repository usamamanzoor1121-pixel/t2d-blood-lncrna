#!/usr/bin/env bash
set -e
export JAVA_HOME=$HOME/jdk17
export PATH=$JAVA_HOME/bin:$PATH
cd /mnt/d/T2D_lncRNA_Standalone
java -version
curl -s https://get.nextflow.io | bash
chmod +x nextflow
./nextflow -version
