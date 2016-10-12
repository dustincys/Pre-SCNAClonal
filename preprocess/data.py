'''
Created on 2013-08-13

@author: Yi Li

pyloh.preprocess.data

================================================================================

Modified on 2014-04-09

@author: Yi Li
'''
import sys
import numpy as np
import scipy.cluster.hierarchy as hcluster
from collections import Counter

#from GCBASELINE import constants
#from GCBASELINE.preprocess.utils import *


class Segment:

    def __init__(self):
        self.name = ""
        self.chrom_idx = -1
        self.chrom_name = ""
        self.start = -1
        self.end = -1
        self.normal_reads_num = -1
        self.tumor_reads_num = -1
        self.sites_num = 0
        self.LOH_frac = 0.0
        self.LOH_status = 'NONE'
        self.APM_frac = 0.0
        self.APM_status = 'NONE'
        self.baseline_label = 'FALSE'
        self.log2_ratio = 0.0
        self.paired_counts = None
        self.BAF_counts = None
        self.copy_number = -1
        self.allele_type = 'NONE'
        self.subclone_prev = -1
        self.subclone_cluster = 'NONE'

        self.gc = -1


#class Data:
#
#    def __init__(self):
#        self.seg_num = 0
#        self.Lambda_S = -1
#        self.segments = []
#
#    def reload_segmentsn(self, normal_bam, tumor_bam, bed_file_normalized_name):
#        """reload the gc corrected tumor counts
#
#        :normal_bam: TODO
#        :tumor_bam: TODO
#        :bed_file_normalized_name: TODO
#        :returns: TODO
#
#        """
#
#        chrom_idx_list = constants.CHROM_IDX_LIST
#        chrom_start = constants.CHROM_START
#
#        sam_SQ = normal_bam.header['SQ']
#        sam_chrom_format = get_chrom_format(map(lambda x: x['SN'], sam_SQ))
#        chrom_lens, chrom_idxs = get_chrom_lens_idxs(chrom_idx_list, sam_SQ)
#        bed_chroms, bed_starts, bed_ends, tumor_reads, normal_reads =\
#            BEDnParser(bed_file_normalized_name)
#        bed_num = len(bed_chroms)
#
#        segment_index = 0
#        for i in range(0, bed_num):
#            chrom_idx = chrom_name_to_idx(bed_chroms[i])
#            chrom_name = chrom_idx_to_name(chrom_idx, sam_chrom_format)
#            seg_name = get_segment_name(chrom_name, bed_starts[i], bed_ends[i])
#
#            if chrom_idx not in chrom_idx_list:
#                print 'Chromsome {0} not found, segment {1} excluded...'.format(
#                    bed_chroms[i], seg_name)
#                sys.stdout.flush()
#                continue
#
#            tumor_reads_num = tumor_reads[i]
#
#            if self.segments[segment_index].name == seg_name and \
#                    self.segments[segment_index].chrom_idx == chrom_idx and\
#                    self.segments[segment_index].chrom_name == chrom_name and\
#                    self.segments[segment_index].start == start and\
#                    self.segments[segment_index].end == end:
#                self.segments[segment_index].tumor_reads_num = tumor_reads_num
#            else:
#                print "segments reload corrected tumor reads error"
#            segment_index = segment_index + 1
#
#    def load_segmentsn(self, bed_file_normalized_name):
#        """
#
#        :bed_file_normalized_name: TODO
#        :returns: TODO
#
#        """
#
#        bed_chroms, bed_starts, bed_ends, tumor_reads, normal_reads, gcs =\
#            BEDnParser(bed_file_normalized_name)
#        get_chrom_format(bed_chroms)
#        bed_num = len(bed_chroms)
#
#        for i in range(0, bed_num):
#            chrom_idx = chrom_name_to_idx(bed_chroms[i])
#            seg_name = get_segment_name(
#                bed_chroms[i], bed_starts[i], bed_ends[i])
#
#            normal_reads_num = normal_reads[i]
#            tumor_reads_num = tumor_reads[i]
#
#            segment_i = Segment()
#            segment_i.name = seg_name
#            segment_i.chrom_idx = chrom_idx
#            segment_i.chrom_name = bed_chroms[i]
#            segment_i.start = bed_starts[i]
#            segment_i.end = bed_ends[i]
#            segment_i.normal_reads_num = normal_reads_num
#            segment_i.tumor_reads_num = tumor_reads_num
#            segment_i.log2_ratio = np.log2(1.0*tumor_reads_num/normal_reads_num)
#
#            segment_i.gc = gcs[i]
#
#            self.segments.append(segment_i)
#            self.seg_num += 1
#
#    def load_segments(self, normal_bam, tumor_bam, bed_file_name):
#        chrom_idx_list = constants.CHROM_IDX_LIST
#        chrom_start = constants.CHROM_START
#
#        sam_SQ = normal_bam.header['SQ']
#        sam_chrom_format = get_chrom_format(map(lambda x: x['SN'], sam_SQ))
#        chrom_lens, chrom_idxs = get_chrom_lens_idxs(chrom_idx_list, sam_SQ)
#
#        bed_chroms, bed_starts, bed_ends, gcs = BEDParser(bed_file_name)
#        get_chrom_format(bed_chroms)
#        bed_num = len(bed_chroms)
#
#        for i in range(0, bed_num):
#            chrom_idx = chrom_name_to_idx(bed_chroms[i])
#            chrom_name = chrom_idx_to_name(chrom_idx, sam_chrom_format)
#            seg_name = get_segment_name(chrom_name, bed_starts[i], bed_ends[i])
#
#            if chrom_idx not in chrom_idx_list:
#                print 'Chromsome {0} not found, segment {1} excluded...'.format(
#                    bed_chroms[i], seg_name)
#                sys.stdout.flush()
#                continue
#
#            chrom_lst_idx = chrom_idxs.index(chrom_idx)
#
#            if bed_starts[i] < chrom_start or bed_ends[
#                    i] > chrom_lens[chrom_lst_idx]:
#                print 'Out of range chromsome {0}, segment {1} excluded...'.\
#                    format(bed_chroms[i], seg_name)
#                sys.stdout.flush()
#                continue
#
#            normal_reads_num = normal_bam.count(
#                chrom_name, bed_starts[i], bed_ends[i])
#            tumor_reads_num = tumor_bam.count(
#                chrom_name, bed_starts[i], bed_ends[i])
#
#            segment_i = Segment()
#            segment_i.name = seg_name
#            segment_i.chrom_idx = chrom_idx
#            segment_i.chrom_name = chrom_name
#            segment_i.start = bed_starts[i]
#            segment_i.end = bed_ends[i]
#            segment_i.normal_reads_num = normal_reads_num
#            segment_i.tumor_reads_num = tumor_reads_num
#            segment_i.log2_ratio = np.log2(1.0 *
#                                           tumor_reads_num/normal_reads_num)
#
#            self.segments.append(segment_i)
#            self.seg_num += 1
#
#    def load_counts_fromSNP(self, tumorData, normalData):
#        """Load the SNP data into each segment,
#        the format is paired_counts and BAF_counts
#
#        :tumorData: TODO
#        :normalData: TODO
#        :returns: TODO
#
#        """
#        for j in range(0, len(self.segments)):
#            normalData_temp, tumorData_temp = get_row_by_segment(
#                tumorData, normalData, self.segments[j])
#            paired_counts_temp = get_paired_counts(
#                normalData_temp, tumorData_temp)
#            self.segments[j].paired_counts = paired_counts_temp
#
#    def get_LOH_frac(self):
#        for j in range(0, self.seg_num):
#            self.segments[j].LOH_frac = get_LOH_frac(
#                self.segments[j].paired_counts)
#
#    def get_APM_frac(self):
#        """
#        :returns: TODO
#
#        """
#        for j in range(0, self.seg_num):
#            self.segments[j].APM_frac = get_APM_frac_MAXMIN(
#                self.segments[j].paired_counts)
#
#    def get_LOH_status(self, baseline_thred):
#        LOH_num = 0
#        for j in range(0, self.seg_num):
#            self.segments[j].LOH_status = get_LOH_status(
#                self.segments[j].LOH_frac, baseline_thred)
#            if self.segments[j].LOH_status == "TRUE":
#                LOH_num = LOH_num + 1
#
#        print "LOH_num/seg_num = {0}/{1}".format(LOH_num, self.seg_num)
#
#    def get_APM_status(self, baseline_thred_APM):
#        APM_num = 0
#        for j in range(0, self.seg_num):
#            self.segments[j].APM_status = get_APM_status(
#                self.segments[j].APM_frac, baseline_thred_APM)
#            if self.segments[j].APM_status == "TRUE":
#                APM_num = APM_num + 1
#
#        print "APM_num/seg_num = {0}/{1}".format(APM_num, self.seg_num)
#
#    def compute_Lambda_S(self, max_copynumber, subclone_num):
#        """ compute the Lambda S, through hierarchy clustering
#        :returns: TODO
#
#        """
#
#        thresh = constants.HC_THRESH
#
#        reads_depth_ratio_log = []
#        reads_depth_ratio = []
#        for j in range(0, self.seg_num):
#            if self.segments[j].APM_status == 'TRUE':
#                ratio = self.segments[
#                    j].tumor_reads_num*1.0/self.segments[j].normal_reads_num
#                reads_depth_ratio_log.append(np.log(ratio))
#                reads_depth_ratio.append(ratio)
#
#        reads_depth_ratio = np.array(reads_depth_ratio)
#        reads_depth_ratio_log = np.array(reads_depth_ratio_log)
#        if reads_depth_ratio_log.shape[0] == 0:
#            print 'Error: APM position found, existing...'
#            sys.exit(1)
#
#        reads_depth_ratio_log = reads_depth_ratio_log.reshape(
#            reads_depth_ratio_log.shape[0], 1)
#        y = np.ones(reads_depth_ratio_log.shape)
#        reads_depth_ratio_log = np.hstack((reads_depth_ratio_log, y))
#        clusters = hcluster.fclusterdata(
#            reads_depth_ratio_log, thresh, criterion="distance")
#        mccs = Counter(clusters).most_common(max_copynumber * subclone_num)
#
#        rdr_min = float('Inf')
#        cluster_min = -1
#        for i in range(0, len(mccs)):
#            cluster_temp = mccs[i][0]
#            print "cluster temp : {}".format(cluster_temp)
#            rdr_temp = reads_depth_ratio[clusters == cluster_temp].mean()
#            print "rdr_temp"
#            print rdr_temp
#            if rdr_min > rdr_temp:
#                rdr_min = rdr_temp
#                cluster_min = cluster_temp
#
#        print mccs
#        print "log Lambda_S: {}".format(np.log(rdr_min))
#
#        cluster_flag = (clusters == cluster_min)
#        baseline_num = 0
#        rdr_i = 0
#        for j in range(0, self.seg_num):
#            if self.segments[j].APM_status == 'TRUE':
#                if cluster_flag[rdr_i]:
#                    if self.segments[j].LOH_status == 'FALSE':
#                        self.segments[j].baseline_label == 'TRUE'
#                        baseline_num = baseline_num + 1
#                    else:
#                        self.segments[j].baseline_label == 'FALSE'
#                else:
#                    self.segments[j].baseline_label == 'FALSE'
#
#                rdr_i = rdr_i + 1
#            else:
#                self.segments[j].baseline_label == 'FALSE'
#
#        print "baseline_num: {}".format(baseline_num)
#
#        if baseline_num == 0:
#            print 'Error: No diploid segments found, existing...'
#            sys.exit(1)
#
#        self.Lambda_S = rdr_min
#
#    def compute_normal_boundary(self):
#        """ compute the Lambda S, through hierarchy clustering
#        :returns: TODO
#
#        """
#        baselines = filter(lambda seg: seg.baseline_label == 'TRUE',
#                           self.segments)
#        sum_r_n = sum(map(lambda seg: seg.normal_reads_num, self.segments))
#        sum_r_t = sum(map(lambda seg: seg.tumor_reads_num, self.segments))
#
#        upper_bound = -float('inf')
#        lower_bound = float('inf')
#
#        for seg in baselines:
#            ratio = (float(seg.tumor_reads_num) / float(sum_r_t)
#                     ) / (float(seg.normal_reads_num) / float(sum_r_n))
#            if upper_bound < ratio:
#                upper_bound = ratio
#            if lower_bound > ratio:
#                lower_bound = ratio
#
#        return upper_bound, lower_bound
