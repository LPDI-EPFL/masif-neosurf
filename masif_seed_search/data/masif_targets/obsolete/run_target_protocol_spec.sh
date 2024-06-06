#!/bin/bash
# For specific file
# PDBFilePath PDB_chain Small_molecule_ID mol2file
./data_prepare_one_spec.sh $1 $2 $3 $4
./predict_site.sh $2
./color_site.sh $2
./compute_descriptors.sh $2
echo "Creating running directory targets/$2 "
cp -r targets/template/ targets/$2
## UNCOMMENT THESE LINES IF YOU WANT D-AMINO ACIDS
# Now compute descriptors and predictions for the mirror image target (d-amino acids)
#./predict_site.sh d$1
#./color_site.sh d$1
#./compute_descriptors.sh d$1
echo "Creating running directory for the mirror image (d-amino acids): targets/d$2 "
#cp -r targets/template/ targets/d$1

