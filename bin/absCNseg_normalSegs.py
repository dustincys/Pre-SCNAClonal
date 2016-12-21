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
                    ratio = str(tumorCount / normalCount)

                    outFile.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(chrom,
                                                                     start,
                                                                     end,
                                                                     segLen,
                                                                     ratio))

    def getSegFnCorrected(self,
                          segments_file_path, segFn_file_path, baseline=-1):
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
                    ratio = str(tumorCount / normalCount)
                    if baseline != -1:
                        ratio = np.log(tumorCount / normalCount) - baseline
                        ratio = str(np.exp(ratio))

                    outFile.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(chrom,
                                                                     start,
                                                                     end,
                                                                     segLen,
                                                                     ratio))

