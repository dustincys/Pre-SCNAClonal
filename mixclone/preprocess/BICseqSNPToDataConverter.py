'''
# =============================================================================
#      FileName: BICseqSNPToDataConverter.py
#          Desc:
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2016-10-11 10:06:07
#       History:
# =============================================================================
'''

import sys

import pickle as pkl

import numpy as np
from scipy.stats import gaussian_kde
from scipy.signal import argrelextrema

from plotGC import GCStripePlot

from data import Data
from utils import read_snp_file
from mcmc import MCMCLM


class THetA_Converter:

    """Docstring for BICseqSNPToDataConverter. """

    def __init__(self, input_filename_base, BICseq_bed_fileName,
                 BICseq_bed_fileName_corrected, tumor_SNP_fileName,
                 normal_SNP_fileName, pkl_path,
                 max_copynumber, subclone_num, sampleNumber,
                 baseline_thred_LOH, baseline_thred_APM,
                 baseline_thred_SNPDENSITY,
                 gamma, process_num):
        """
            BICseq_bed_fileName: bicseq file name
        """

        self.input_filename_base = input_filename_base
        self.BICseq_bed_fileName = BICseq_bed_fileName
        self.BICseq_bed_fileName_corrected = BICseq_bed_fileName_corrected
        self.tumor_SNP_fileName = tumor_SNP_fileName
        self.normal_SNP_fileName = normal_SNP_fileName
        self.pkl_path = pkl_path

        self.max_copynumber = max_copynumber
        self.subclone_num = subclone_num
        self.sampleNumber = sampleNumber

        self.baseline_thred_LOH = baseline_thred_LOH
        self.baseline_thred_APM = baseline_thred_APM
        self.baseline_thred_SNPDENSITY = baseline_thred_SNPDENSITY
        self.gamma = gamma
        self.process_num = process_num

        self.data = Data()

    def convert(self, method, pkl_flag=False):
        """convert data to the THetA style
        :returns: BICseq bed file, gc corrected, and relevant parameters

        """
        print "pkl_flag"
        print pkl_flag
        if pkl_flag and self.pkl_path != "":
            infile = open(self.pkl_path, 'rb')
            self.data = pkl.load(infile)
            infile.close()
        else:
            self._load_segments()
            print "THetA converter converting"

            if "auto" == method:
                self._MCMC_gccorrection()
            elif "visual" == method:
                self._visual_gccorrection()

            self._load_SNP()

        self._baseline_selection()

        self._output()

        self.data.outSNV(self.input_filename_base + '.snv.txt')

        data_file_name = self.input_filename_base + '.THetA.input.pkl'
        outfile = open(data_file_name, 'wb')
        pkl.dump(self.data, outfile, protocol=2)

        outfile.close()
        pass

    def _baseline_selection(self):
        """
        :returns: TODO

        """
        print "LOH = "
        print self.baseline_thred_LOH
        print "APM = "
        print self.baseline_thred_APM
        self.data.get_LOH_frac_SNP()
        self.data.get_LOH_status(self.baseline_thred_LOH, True)
        self.data.get_APM_frac_SNP()
        self.data.get_APM_status_SNP(self.baseline_thred_APM,
                                     self.baseline_thred_SNPDENSITY)
        self.data.compute_Lambda_S(self.max_copynumber, self.subclone_num, True)

    def _MCMC_gccorrection(self):
        """
        The interception is irrelevant for correction, set as median
        MCMCLM only returns the m and c, then correct the data here
        """
        mcmclm = MCMCLM(self.data, 0, self.subclone_num, self.max_copynumber)
        m, c = mcmclm.run()
        print "MCMC slope = {}".format(m)
#       self.data.pr = mcmclm.getPeakRange(m)
#       print "MCMC peak range = {}".format(self.data.pr)
        self._correct(m, c)

    def _correct(self, slope, intercept):

        x = np.array(map(lambda seg: seg.gc, self.data.segments))
        y = np.array(map(lambda seg: np.log(seg.tumor_reads_num + 1) -
                         np.log(seg.normal_reads_num + 1),
                         self.data.segments))

        K = np.percentile(y, 50)
        A = slope * x + intercept
        y_corrected = y - A + K

        for i in range(len(y_corrected)):
            self.data.segments[i].tumor_reads_num = np.exp(
                y_corrected[i] +
                np.log(self.data.segments[i].normal_reads_num + 1)
            ) - 1

        print "gc corrected, with slope = {0}, intercept = {1}".\
            format(slope, intercept)

    def _visual_gccorrection(self):
        gsp = GCStripePlot(self.data.segments, self.sampleNumber)
        print "total number: {}".format(self.data.seg_num)

#       Sampling then linear regression, poor performance
#       gsp.sampleln([i * 1000 for i in range(1,9)], 100)
        gsp.plot()

# todo   trimed x, y position
        x, y, m, c = gsp.output()
#       self.data.pr = (max(y) - min(y)) / self.max_copynumber

        print "x, y, m, c"
        print x, y, m, c

        self._correct(m, c)

    def getPeakRange(self):
        """
        :returns: TODO

        """
        baseLine = np.log(self.data.Lambda_S)

        y_corrected = np.array(map(lambda seg: np.log(seg.tumor_reads_num + 1) -
                                   np.log(seg.normal_reads_num + 1),
                                   self.data.segments))
        y_density = gaussian_kde(y_corrected)
        y_down = min(y_corrected)
        y_up = max(y_corrected)
        y_xs = np.linspace(y_down, y_up,
                           1000*self.subclone_num*self.max_copynumber)
        y_ys = y_density(y_xs)
        peaks = argrelextrema(y_ys, np.greater)[0]
        valleys = argrelextrema(y_ys, np.less)[0]

        basePeak = peaks[np.argmin(np.abs(y_ys[peaks] - baseLine))]

        valley_u = len(y_ys)
        valley_l = -1
        for v in valleys:
            if v < basePeak:
                if v > valley_l:
                    valley_l = v
            if v > basePeak:
                if v < valley_u:
                    valley_u = v

        hillside_u = y_xs[valley_u] - y_xs[basePeak]
        hillside_l = y_xs[basePeak] - y_xs[valley_l]

        pr = 2 * min(hillside_l, hillside_u)

        return pr

    def _output(self):
        """Output the parameter for THetA

        The Upper and Lower Boundaries for normal heuristic
        The GC corrected interval_count_file
        """
        self.data.pr = self.getPeakRange()
        upper_bound, lower_bound = self.data.compute_normal_boundary(
            self.data.pr)
        print "upper_bound = {0}\n lower_bound = {1}".format(upper_bound,
                                                             lower_bound)
        sys.stdout.flush()

        interval_count_file = open(self.BICseq_bed_fileName_corrected, 'w')
#       interval_count_file.write(
#           "chrm\tstart\tend\ttumorCount\tnormalCount\tgc\n")

        for i in range(len(self.data.segments)):
            chrm_i = self.data.segments[i].chrom_name
            start_i = self.data.segments[i].start
            end_i = self.data.segments[i].end
            tumorCount_i = int(float(self.data.segments[i].tumor_reads_num))
            normalCount_i = int(float(self.data.segments[i].normal_reads_num))
            gc_i = self.data.segments[i].gc
            interval_count_file.write(
                "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(chrm_i,
                                                        start_i,
                                                        end_i,
                                                        tumorCount_i,
                                                        normalCount_i,
                                                        gc_i))

        interval_count_file.close()

        print "GC corrected interval file generated!"
        sys.stdout.flush()

    def _load_segments(self):
        print 'Loading normalized segments by {0}...'.\
            format(self.BICseq_bed_fileName)
        self.data.load_segmentsn(self.BICseq_bed_fileName)

    def _load_SNP(self):
        """ load snp from file
        :returns: TODO

        """
        tumorData = read_snp_file(self.tumor_SNP_fileName)
        normalData = read_snp_file(self.normal_SNP_fileName)

        # generate the paired_counts and BAF_counts from snp
        self.data.load_counts_fromSNP(tumorData, normalData, self.gamma,
                                      self.process_num)

        pass
