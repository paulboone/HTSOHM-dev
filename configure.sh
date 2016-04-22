#!/bin/bash

export HTSOHM_DIR=${PWD}          # specifies HTSOHM directory
export RASPA_DIR=${HOME}/RASPA    # specifies RASPA directory


#writing environment variables to .bashrc
DEST=~/.bashrc
echo $'\n# HTSOHM directories' >> $DEST
echo "export HTSOHM_DIR=${HTSOHM_DIR}" >> $DEST
echo "export RASPA_DIR=${RASPA_DIR}" >> $DEST
echo "export SRC_DIR=\${HTSOHM_DIR}/bin" >> $DEST
echo "export FF_DIR=\${RASPA_DIR}/share/raspa/forcefield" >> $DEST
echo "export MAT_DIR=\${RASPA_DIR}/share/raspa/structures/cif" >> $DEST
echo "# " >> ~/.bashrc

#create local database/table...
python $HTSOHM_DIR/bin/runDB_declarative.py
