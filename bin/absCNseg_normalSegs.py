#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: absCNseg_normalSegs.py
#          Desc: generate the normlized segments for absCNseg
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2016-12-20 15:02:04
#       History:
# =============================================================================
'''

import numpy as np


class ABSCNSEGIN(object):

    """get absCNseg input"""

    def __init__(self):
        pass

    def getSegFn(self, segments_file_path, segFn_file_path):
        with open(segFn_file_path, "w") as outFile:
            outFile.write(
                "chrom\tloc.start\tloc.end\teff.seg.len\tnormalized.ratio\n")
            with open(segments_file_path) as inFile:
                for line in inFile:
                    line = line.strip()
                    if line == "":
                        continue
                    listLine = line.split("\t")
                    if listLine[0] == "chrom":
                        continue
                    chrom = listLine[0]
                    start = listLine[1]
                    end = listLine[2]
                    segLen = str(int(end) - int(start))
                    tumorCount = float(listLine[3])
                    normalCount = float(listLine[4])
                    ratio = str((tumorCount + 1.0) / (normalCount + 1.0))

                    outFile.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(chrom,
                                                                     start,
                                                                     end,
                                                                     segLen,
                                                                     ratio))

    def getSegFnCorrected(self,
                          segments_file_path, segFn_file_path, baseline=0):
        with open(segFn_file_path, "w") as outFile:
            outFile.write(
                "chrom\tloc.start\tloc.end\teff.seg.len\tnormalized.ratio\n")
            with open(segments_file_path) as inFile:
                for line in inFile:
                    line = line.strip()
                    if line == "":
                        continue
                    listLine = line.split("\t")
                    if listLine[0] == "ID":
                        continue
                    chrom = listLine[1]
                    start = listLine[2]
                    end = listLine[3]
                    segLen = str(int(end) - int(start))
                    tumorCount = float(listLine[4])
                    normalCount = float(listLine[5])
                    ratio = str(tumorCount + 1 / normalCount + 1)
                    if baseline != 0:
                        ratio = np.log(tumorCount / normalCount) - baseline
                        ratio = str(np.exp(ratio))

                    outFile.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(chrom,
                                                                     start,
                                                                     end,
                                                                     segLen,
                                                                     ratio))

def main():

    absIO = ABSCNSEGIN()

    idx = ["n5t95", "n20t80", "n40t60", "n60t40", "n80t20", "n95t5"]
    baseline = [-0.356, -0.113, -0.188, -0.126, -0.089, -0.040]

    for item, bl in zip(idx, baseline):
        originalSegPath = "/data/yschu/projects/subclone/data/BICseq/ucsc_benchmark/1954/getGCMap/getGCMapHCC1954.mix1.{}.bam.bicseq.gc.txt".format(item)
        correcstedSegPath = "/data/yschu/projects/subclone_GCBASELINE/pipelineTest/MixClone/corrected/getGCMapHCC1954.mix1.{}.bam.bicseq.gc.txt.50e4.bed.gccorrected".format(item)
        outOriginalFn = "/data/yschu/projects/subclone_GCBASELINE/pipelineTest/absCNseq/original/{}.original.fn".format(item)
        outCorrectedFn = "/data/yschu/projects/subclone_GCBASELINE/pipelineTest/absCNseq/result/{}.corrected.fn".format(item)
        absIO.getSegFn(originalSegPath, outOriginalFn)
        absIO.getSegFnCorrected(correcstedSegPath, outCorrectedFn, bl)
        print "finished {}!".format(item)


if __name__=="__main__":
    main()
