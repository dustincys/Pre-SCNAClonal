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


class BICseqSNPToDataConverter:

    """Docstring for BICseqSNPToDataConverter. """

    def __init__(self, BICseq_bed_fileName, tumor_SNP_fileName,
                 normal_SNP_fileName, gc_corrected_BICseq_bed_fileName,
                 seg_length, max_copy_number, subclone_num, lm_lowerbound,
                 lm_upperbound, delta
                 ):
        """
            BICseq_bed_fileName: bicseq file name
        """

        self.BICseq_bed_fileName = BICseq_bed_fileName
        self.tumor_SNP_fileName = tumor_SNP_fileName
        self.normal_SNP_fileName = normal_SNP_fileName

        self.gc_corrected_BICseq_bed_fileName = gc_corrected_BICseq_bed_fileName
        self.seg_length = seg_length
        self.max_copy_number = max_copy_number
        self.subclone_num = subclone_num
        self.lm_lowerbound = lm_lowerbound
        self.lm_upperbound = lm_upperbound
        self.delta = delta

        self.data = Data()

    def convert(self):
        """convert data to the THetA style
        :returns: BICseq bed file, gc corrected, and relevant parameters

        """
        self._load_segments()
        self._load_SNP()
        self._get_Baseline()

        self._output()
        pass

    def _output(self):
        """Output the parameter for THetA

        The Upper and Lower Boundaries for normal heuristic
        """
        self.data.get_LOH_frac()
        self.data.get_LOH_status(self.baseline_thred)
        self.data.get_APM_frac()
        self.data.get_APM_status(self.baseline_thred_APM)
        self.data.compute_Lambda_S_hc_APM(
            self.max_copynumber, self.subclone_num)

        upper_bound, lower_bound = self.data.compute_normal_boundary()
        print "upper_bound = {0}\n lower_bound = {1}".format(upper_bound,
                                                             lower_bound)
        sys.stdout.flush()

    def _load_segments(self):
        """
        :returns: TODO

        """
        args = (self.BICseq_bed_fileName,
                self.gc_corrected_BICseq_bed_fileName,
                self.seg_length,
                self.max_copy_number,
                self.subclone_num,
                self.lm_lowerbound,
                self.lm_upperbound,
                self.delta)

        gccorrect(args)

        self.BICseq_bed_fileName = gc_corrected_BICseq_bed_fileName

        print 'Loading normalized segments by {0}...'.format(self.segments_bed)
        self.data.load_segmentsn(self.segments_bed)

    def _load_SNP(self):
        """ load snp from file
        :returns: TODO

        """
        tumorData = read_snp_file(tumorfile)
        normalData = read_snp_file(normalfile)

        # generate the paired_counts and BAF_counts from snp
        self.data.load_counts_fromSNP(tumorData, normalData)

        pass
