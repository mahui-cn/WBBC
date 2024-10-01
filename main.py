# -*- coding: utf-8 -*-

import sys, os
import argparse
import time
from multiprocessing import cpu_count
import wbbc

sys.path.append(os.path.dirname(__file__))


def arguments():
    # argument parser
    parser = argparse.ArgumentParser()

    # TSV raw data file
    parser.add_argument(
        "-tf",
        "--tsvFile",
        nargs="?",
        required=False,
        help="file name of TSV raw data",
    )

    # TSV raw data path
    parser.add_argument(
        "-tp",
        "--tsvPath",
        nargs="?",
        required=False,
        help="path of TSV raw data",
    )

    # filename of allele frequency without extension
    parser.add_argument(
        "-af",
        "--alleleFrqFile",
        nargs="?",
        default="wbbc_{}".format(time.strftime("%Y-%m-%d", time.localtime())),
        help="file name of allele frequency without extension",
    )

    # decimal digits of alleles frequency
    parser.add_argument(
        "-ad",
        "--afDigits",
        nargs="?",
        type=int,
        default=6,
        help="decimal digits of alleles frequency",
    )

    # threshold of standard deviation of alleles frequency
    parser.add_argument(
        "-sd",
        "--stdev",
        nargs="?",
        type=float,
        default=0.03,
        help="threshold of standard deviation of alleles frequency",
    )

    # filename of regions of High Linkage Disequilibrium
    parser.add_argument(
        "-hld",
        "--highLD",
        nargs="?",
        default="",
        help="filename of regions of High Linkage Disequilibrium. Refer to: https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)",
    )

    # threshold of standard deviation of alleles frequency
    parser.add_argument(
        "-th",
        "--thread",
        nargs="?",
        type=int,
        default=cpu_count() if cpu_count() > 0 else 4,
        help="threads count for concurrency",
    )

    return parser.parse_args()


# 基于基因芯片的TSV数据，从WBBC的VCF中筛选出符合条件的突变和频率数据，生成alleles和frequency文件作为祖源计算器的参考数据集
def main():
    # get arguments
    args = arguments()

    # TSV文件集合
    tsvFiles = []

    # 指定TSV文件
    if args.tsvFile != None:
        tsvFiles.extend(args.tsvFile.strip().split())

    # 指定TSV目录
    if args.tsvPath != None:
        tsvPath = args.tsvPath.strip()
        for file in os.listdir(tsvPath):
            if os.path.isfile(os.path.join(tsvPath, file)):
                tsvFiles.append(os.path.join(tsvPath, file))

    if len(tsvFiles) > 0:
        isTsvExist = True
        for file in tsvFiles:
            if not os.access(file, os.F_OK):
                print("TSV file {} is not available.".format(file))
                isTsvExist = False
    else:
        print("Please set one more TSV files or path.")
        isTsvExist = False

    if isTsvExist:
        wbbc.make_allele_frq(
            tsvFiles,
            args.highLD,
            args.alleleFrqFile,
            args.afDigits,
            args.stdev,
            args.thread,
        )


if __name__ == "__main__":
    main()
