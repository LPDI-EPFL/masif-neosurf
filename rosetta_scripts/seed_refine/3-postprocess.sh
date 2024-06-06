#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Error: Please provide either 'S' (Sheet seeds) or 'H' (Helical seeds), and the ligand identifier as an argument."
    echo "Example: ./3-postprocess.sh S BB2"
    exit 1
fi

if [ "$1" == "S" ]; then
    echo "Executing post-processing for sheet-based seeds."
    echo "--- Calculating percentage of sheet-based contacts ---"
    python3 ./utils/sheet/get_dssp_sheet.py
    echo "--- Calculating similarity score ---"
    python3 ./utils/sheet/get_similarity_sheet.py
    echo "--- Calculating atom contact counts with molecule ---"
    python3 ./utils/sheet/get_contacts_sheet.py $2 True
elif [ "$1" == "H" ]; then
    echo "Executing post-processing for helical seeds."
    echo "--- Calculating similarity score ---"
    python3 ./utils/helix/get_similarity_helix.py
    echo "--- Calculating atom contact counts with molecule ---"
    python3 ./utils/helix/get_contacts_helix.py $2
else
    echo "Error: Please provide either 'S' (Sheet seeds) or 'H' (Helical seeds) as an argument."
    echo "Example: ./3-postprocess.sh S BB2"
    exit 1
fi

    echo "--- Run finished ---"
