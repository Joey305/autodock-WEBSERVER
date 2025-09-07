#!/bin/bash
# 1_Bcount.sh - count files inside any directory with "_PDBQT_" in its name

echo "📂 Counting files in directories with '_PDBQT_'..."
echo ""

# Loop over matching directories
for d in *_PDBQT_*; do
    if [ -d "$d" ]; then
        count=$(find "$d" -type f | wc -l)
        echo "$d : $count files"
    fi
done

