#!/bin/sh

./run.py THetA \
    /media/d/data/NGS/projects/subclone/ucsc_benchmark/1954/getGCMap/getGCMapHCC1954.mix1.n20t80.bam.bicseq.gc.txt \
    /media/d/data/NGS/projects/subclone/ucsc_benchmark/1954/snp/HCC1954.n20t80.withCounts \
    /media/d/data/NGS/projects/subclone/ucsc_benchmark/1954/snp/HCC1954.NORMAL.30x.withCounts \
    --gc_correction_method visual\
    --baseline_selection_method visual\



