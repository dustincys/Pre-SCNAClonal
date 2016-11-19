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

from plotGC import GCStripePlot
from plotBaseline import BaselinePlot

from data import Data
from utils import read_snp_file
from mcmc import MCMCLM


class THetA_Converter:

    """Docstring for BICseqSNPToDataConverter. """

    def __init__(self, BICseq_bed_fileName, tumor_SNP_fileName,
                 normal_SNP_fileName, seg_length, max_copynumber,
                 subclone_num, sampleNumber, lm_lowerbound, lm_upperbound,
                 delta):
        """
            BICseq_bed_fileName: bicseq file name
        """

        self.BICseq_bed_fileName = BICseq_bed_fileName
        self.tumor_SNP_fileName = tumor_SNP_fileName
        self.normal_SNP_fileName = normal_SNP_fileName

        self.seg_length = seg_length
        self.max_copynumber = max_copynumber
        self.subclone_num = subclone_num
        self.sampleNumber = sampleNumber
        self.lm_lowerbound = lm_lowerbound
        self.lm_upperbound = lm_upperbound
        self.delta = delta

        self.data = Data()

    def convert(self, methods):
        """convert data to the THetA style
        :returns: BICseq bed file, gc corrected, and relevant parameters

        """
        self._load_segments()

        print "THetA converter converting"
        gc_correction_method, baseline_selection_method = methods

        if "auto" == gc_correction_method:
            print "auto gc correction"
            self._MCMC_gccorrection()
        elif "visual" == gc_correction_method:
            print "visual gc correction"
            self._visual_gccorrection()

        self._load_SNP()
        self._baseline_selection()

# baseline visual selection, not meaningful, if the baseline segments are not
# obviously located
#        if "visual" == gc_correction_method:
#            self._visual_baseline_selection()
        self._output()
        pass

    def _baseline_selection(self):
        """
        :returns: TODO

        """
        self.data.get_LOH_frac()
        self.data.get_LOH_status(self.baseline_thred_LOH, True)
        self.data.get_APM_frac()
        self.data.get_APM_status(self.baseline_thred_APM)
        self.data.compute_Lambda_S(self.max_copynumber, self.subclone_num, True)

    def _visual_baseline_selection(self):
        """
        :returns: TODO

        """

        blp = BaselinePlot(self.data, self.max_copynumber, self.subclone_num)
        blp.plot()

    def _MCMC_gccorrection(self):
        """
        The interception is irrelevant for correction, set as median
        MCMCLM only returns the m and c, then correct the data here
        """
        mcmclm = MCMCLM(self.data, 0, self.subclone_num, self.max_copynumber)
        m, c = mcmclm.run()
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

#todo   trimed x, y position
        x, y, m, c = gsp.output()

        print "x, y, m, c"
        print x, y, m, c

        self._correct(m,c)


    def _output(self):
        """Output the parameter for THetA

        The Upper and Lower Boundaries for normal heuristic
        """
        upper_bound, lower_bound = self.data.compute_normal_boundary()
        print "upper_bound = {0}\n lower_bound = {1}".format(upper_bound,
                                                             lower_bound)
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
        self.data.load_counts_fromSNP(tumorData, normalData)

        pass
