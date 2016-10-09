import sys
import time
import pickle as pkl
from multiprocessing import Pool

import numpy as np
import pysam

from mixclone import constants
from mixclone.preprocess.data import Data
from mixclone.preprocess.io import PairedCountsIterator, PairedPileupIterator
from mixclone.preprocess.utils import *

def run_preprocess(args):
    time_start = time.time()

    if "MixClone" == args.IOFormat:
        converter = BamToDataConverter(
            args.normal_bam,
            args.tumor_bam,
            args.reference_genome,
            args.input_filename_base,
            args.segments_bed,
            min_depth=args.min_depth,
            min_bqual=args.min_base_qual,
            min_mqual=args.min_map_qual,
            process_num = args.process_num
        )
    elif "THetA" == args.IOFormat:
        converter = BICseqSNPToDataConverter(
            args.BICseq_bed,
            args.tumor_SNP,
            args.normal_SNP,
        )

    converter.convert()

    time_end = time.time()

    print 'Run time: {0:.2f} seconds'.format(time_end - time_start)
    sys.stdout.flush()
