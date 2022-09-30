#!/bin/bash
INPUT="$@"
if [ -z "$INPUT" ]
then
	echo "----"
    echo "The script runs PHANOTATE will all combinations of STOP codons."
    echo "Usage: run_phanotate.sh <output_directory> <input_file(s)>"
    echo "Example: ./run_phanotate.sh genbank/ test-data/*"
    echo "----"
else
    OUTDIR="$1"
    for INFILE in ${@:2}
    do
        # get the file name without the extension
        FILENAME=$(basename "$INFILE")
        FILENAME="${FILENAME%.*}"
        echo " - working on $INFILE"
        echo "TAG, TGA, TAA" #all three get assigned as stop codons 
        phanotate.py "$INFILE" -e TAG,TAA,TGA -f genbank -o "$OUTDIR/$FILENAME"-TAG-TGA-TAA.gbk 
        echo "TAG, TGA"
        phanotate.py "$INFILE" -e TAG,TGA -f genbank -o "$OUTDIR/$FILENAME"-TAG-TGA.gbk
        echo "TGA, TAA"
        phanotate.py "$INFILE" -e TGA,TAA -f genbank -o "$OUTDIR/$FILENAME"-TGA-TAA.gbk
        echo "TAG, TAA"
        phanotate.py "$INFILE" -e TAG,TAA -f genbank -o "$OUTDIR/$FILENAME"-TAG-TAA.gbk
        echo "TAG"
        phanotate.py "$INFILE" -e TAG -f genbank -o "$OUTDIR/$FILENAME"-TAG.gbk
        echo "TGA"
        phanotate.py "$INFILE" -e TGA -f genbank -o "$OUTDIR/$FILENAME"-TGA.gbk
        echo "TAA"
        phanotate.py "$INFILE" -e TAA -f genbank -o "$OUTDIR/$FILENAME"-TAA.gbk
    done
fi

