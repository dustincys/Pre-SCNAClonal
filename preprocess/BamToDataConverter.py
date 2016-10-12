'''
# =============================================================================
#      FileName: BamToDataConverter.py
#          Desc:
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2016-10-08 10:54:48
#       History:
# =============================================================================
'''

import sys
import pickle as pkl
from multiprocessing import Pool
import numpy as np
import pysam

from GCBASELIN.preprocess.data import Data
from GCBASELIN.preprocess.iofun import PairedCountsIterator, PairedPileupIterator
from GCBASELIN.preprocess.utils import *


class BamToDataConverter:

    def __init__(self, normal_bam_filename, tumor_bam_filename,
                 reference_genome_filename, input_filename_base, segments_bed,
                 min_depth=20, min_bqual=10, min_mqual=10, process_num=1):
        self.normal_bam_filename = normal_bam_filename
        self.tumor_bam_filename = tumor_bam_filename
        self.reference_genome_filename = reference_genome_filename
        self.input_filename_base = input_filename_base
        self.segments_bed = segments_bed

        self.min_depth = min_depth
        self.min_bqual = min_bqual
        self.min_mqual = min_mqual
        self.process_num = process_num

        self.data = Data()

    def convert(self, methods):

        gc_correction_method, baseline_selection_method = methods

        self._load_segments()

        self._get_counts()

        self._get_LOH_frac()

        self._get_APM_frac()

        data_file_name = self.input_filename_base + '.MixClone.input.pkl'
        outfile = open(data_file_name, 'wb')
        pkl.dump(self.data, outfile, protocol=2)

        outfile.close()

    def _load_segmentsn(self):
        """
        :returns: TODO

        """
        normal_bam = pysam.Samfile(self.normal_bam_filename, 'rb')
        tumor_bam = pysam.Samfile(self.tumor_bam_filename, 'rb')

        print 'Loading normalized segments by {0}...'.format(self.segments_bed)
        sys.stdout.flush()
        self.data.load_segmentsn(normal_bam, tumor_bam, self.segments_bed)

        normal_bam.close()
        tumor_bam.close()

    def _load_segments(self):
        normal_bam = pysam.Samfile(self.normal_bam_filename, 'rb')
        tumor_bam = pysam.Samfile(self.tumor_bam_filename, 'rb')

        print 'Loading segments by {0}...'.format(self.segments_bed)
        sys.stdout.flush()
        self.data.load_segments(normal_bam, tumor_bam, self.segments_bed)

        normal_bam.close()
        tumor_bam.close()

    def _get_counts(self):
        seg_num = self.data.seg_num
        process_num = self.process_num

        if process_num > seg_num:
            process_num = seg_num

        pool = Pool(processes=process_num)

        args_list = []

        for j in range(0, seg_num):
            seg_name = self.data.segments[j].name
            chrom_name = self.data.segments[j].chrom_name
            chrom_idx = self.data.segments[j].chrom_idx
            start = self.data.segments[j].start
            end = self.data.segments[j].end

            args_tuple = (
                seg_name,
                chrom_name,
                chrom_idx,
                start,
                end,
                self.normal_bam_filename,
                self.tumor_bam_filename,
                self.reference_genome_filename,
                self.min_depth,
                self.min_bqual,
                self.min_mqual)

            args_list.append(args_tuple)

        counts_tuple_list = pool.map(process_by_segment, args_list)

        for j in range(0, seg_num):
            paired_counts_j, BAF_counts_j = counts_tuple_list[j]

            self.data.segments[j].paired_counts = paired_counts_j
            self.data.segments[j].BAF_counts = BAF_counts_j

    def _get_LOH_frac(self):
        self.data.get_LOH_frac()

    def _get_APM_frac(self):
        self.data.get_APM_frac()
#===============================================================================
# Function
#===============================================================================


def process_by_segment(args_tuple):
    seg_name, chrom_name, chrom_idx, start, end, normal_bam_filename, tumor_bam_filename, \
        reference_genome_filename, min_depth, min_bqual, min_mqual = args_tuple

    print 'Preprocessing segment {0}...'.format(seg_name)
    sys.stdout.flush()

    normal_bam = pysam.Samfile(normal_bam_filename, 'rb')
    tumor_bam = pysam.Samfile(tumor_bam_filename, 'rb')
    ref_genome_fasta = pysam.Fastafile(reference_genome_filename)

    normal_pileup_iter = normal_bam.pileup(chrom_name, start, end)
    tumor_pileup_iter = tumor_bam.pileup(chrom_name, start, end)

    paired_pileup_iter = PairedPileupIterator(
        normal_pileup_iter, tumor_pileup_iter, start, end)
    paired_counts_iter = PairedCountsIterator(
        paired_pileup_iter,
        ref_genome_fasta,
        chrom_name,
        chrom_idx,
        min_depth,
        min_bqual,
     min_mqual)

    paired_counts_j, BAF_counts_j = iterator_to_counts(paired_counts_iter)
    counts_tuple_j = (paired_counts_j, BAF_counts_j)

    normal_bam.close()
    tumor_bam.close()
    ref_genome_fasta.close()

    return counts_tuple_j


def iterator_to_counts(paired_counts_iter):
    buff = 100000

    paired_counts_j = np.array([[], [], [], [], [], []], dtype=int).transpose()
    BAF_counts_j = np.zeros((100, 100))
    buff_counts = []
    i = 0

    for counts in paired_counts_iter:
        buff_counts.append(counts)
        i = i + 1

        if i < buff:
            continue

        buff_counts = np.array(buff_counts)

        if buff_counts.shape[0] != 0:
            BAF_counts_buff = get_BAF_counts(buff_counts)
            BAF_counts_j += BAF_counts_buff

        buff_counts_filtered = normal_heterozygous_filter(buff_counts)

        if buff_counts_filtered.shape[0] != 0:
            paired_counts_j = np.vstack((paired_counts_j, buff_counts_filtered))

        buff_counts = []
        i = 0

    buff_counts = np.array(buff_counts)

    if buff_counts.shape[0] != 0:
        BAF_counts_buff = get_BAF_counts(buff_counts)
        BAF_counts_j += BAF_counts_buff

    buff_counts_filtered = normal_heterozygous_filter(buff_counts)

    if buff_counts_filtered.shape[0] != 0:
        paired_counts_j = np.vstack((paired_counts_j, buff_counts_filtered))

    return (paired_counts_j, BAF_counts_j)
