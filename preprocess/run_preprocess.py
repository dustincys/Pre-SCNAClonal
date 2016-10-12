'''
# =============================================================================
#      FileName: run_preprocess.py
#          Desc: run_preprocess
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2016-10-11 10:57:17
#       History: Yi Li
# =============================================================================
'''

import time


def run_preprocess(args):
    '''
    args.gc_correction_method: manual, auto
    args.baseline_selection_method: manual, auto
    '''

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
            process_num=args.process_num
        )
    elif "THetA" == args.IOFormat:
        converter = BICseqSNPToDataConverter(
            args.BICseq_bed,
            args.tumor_SNP,
            args.normal_SNP,
        )

    methods = (args.gc_correction_method, args.baseline_selection_method)
    converter.convert(methods)

    time_end = time.time()

    print 'Run time: {0:.2f} seconds'.format(time_end - time_start)
    sys.stdout.flush()
