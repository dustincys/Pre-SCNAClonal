'''
# =============================================================================
#      FileName: BICseqSNPToDataConverter.py
#          Desc:
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2016-10-07 18:21:55
#       History:
# =============================================================================
'''


class BICseqSNPToDataConverter:

    """Docstring for BICseqSNPToDataConverter. """

    def __init__( self, BICseq_bed_fileName, tumor_SNP_fileName,
                 normal_SNP_fileName):
        """
            BICseq_bed_fileName: bicseq file name
        """

        self.BICseq_bed_fileName = BICseq_bed_fileName
        self.tumor_SNP_fileName = tumor_SNP_fileName
        self.normal_SNP_fileName = normal_SNP_fileName

        self.data = Data()

    def convert(self):
        """convert data to the THetA style
        :returns: BICseq bed file, gc corrected, and relevant parameters

        """
        self._load_segments()
        self._load_SNP()
        self._get_BAF()
        self._get_LOH()
        self._get_APM()
        self._get_Baseline()

        self._output()
        pass

    def _load_segments(self):
        """
        :returns: TODO

        """
        args = (self.BICseq_bed_fileName, gc_corrected_BICseq_bed_fileName, seg_length, max_copy_number, subclone_num, lm_lowerbound, lm_upperbound, delta)

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
        # paired_counts_j = np.array([[], [], [], [], [], []], dtype=int).transpose()
        self.data.load_counts_fromSNP(tumorData, normalData)

        pass


