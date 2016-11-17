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
        gc_correction_method, baseline_selection_method = methods

        print "THetA converter converting"
#        if "auto" == gc_correction_method:
#            print "auto gc correction"
#            self._gccorrection()

        self._load_segments()

        if "visual" == gc_correction_method:
            print "visual gc correction"
            self._visual_gccorrection()
            self._MCMC_gccorrection()

        # self._load_SNP()
        # self._baseline_selection()

        # if "visual" == gc_correction_method:
        #    self._visual_baseline_selection()
        # self._output()
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
        """
        mcmclm = MCMCLM(self.data, 0, self.subclone_num, self.max_copynumber)
        mcmclm.run()

    def _visual_gccorrection(self):
        gsp = GCStripePlot(self.data.segments, self.sampleNumber)
        print "total number: {}".format(self.data.seg_num)
#        gsp.sampleln([i * 1000 for i in range(1,9)], 100)
        gsp.plot()

        x, y, m, c = gsp.output()
        print "x, y, m, c"
        print x, y, m, c

        return x, y, m, c

    def _gccorrection(self):
        gc_corrected_BICseq_bed_fileName = self.BICseq_bed_fileName + ".temp"
        args = (self.BICseq_bed_fileName,
                gc_corrected_BICseq_bed_fileName,
                self.seg_length,
                self.max_copynumber,
                self.subclone_num,
                self.lm_lowerbound,
                self.lm_upperbound,
                self.delta)
        gccorrect(args)
        self.BICseq_bed_fileName = gc_corrected_BICseq_bed_fileName

    def _output(self):
        """Output the parameter for THetA

        The Upper and Lower Boundaries for normal heuristic
        """
        upper_bound, lower_bound = self.data.compute_normal_boundary()
        print "upper_bound = {0}\n lower_bound = {1}".format(upper_bound,
                                                             lower_bound)
        sys.stdout.flush()

    def _load_segments(self):
        """
        :returns: TODO

        """

        print 'Loading normalized segments by {0}...'.format(self.BICseq_bed_fileName)
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
