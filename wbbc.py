# -*- coding: utf-8 -*-

import os
import statistics
from concurrent.futures import ThreadPoolExecutor, as_completed

VCF_FILENAME = "wbbc_vcf/WBBC.chr{}.GRCh37_PhaseI.vcf"


# https://wbbc.westlake.edu.cn/
# 根据西湖中国样本库，生成祖源模型文件
# The VCF is annotated with rsIDs from dbSNP151, and the following INFO fields:
# AC:Allele count in called genotypes in WBBC
# AF:Allele frequency in called genotypes in WBBC
# AN:Total number of alleles in called genotypes in WBBC
# NS:Total number of samples in called genotypes in WBBC
# North_AF:Allele frequency in North Han Chinese
# North_AN:Total number of alleles in North Han Chinese
# Central_AF:Allele frequency in Central Han Chinese
# Central_AN:Total number of alleles in Central Han Chinese
# South_AF:Allele frequency in South Han Chinese
# South_AN:Total number of alleles in South Han Chinese
# Lingnan_AF:Allele frequency in Lingnan Han Chinese
# Lingnan_AN:Total number of alleles in Lingnan Han Chinese
# DP:Raw read depth
# VQSLOD:Variant Recalibration Score from GATK
def make_allele_frq(
    tsvFiles,
    highLDFileName="",
    modelPath=".",
    alleleFrqFile="wbbc",
    afDigits=6,
    stdDevThreshold=0.03,
    maxWorkers=4,
):
    try:
        # 从TSV文件获取RSID集合
        rsid_set = get_rsid(tsvFiles, highLDFileName)

        # SNP总计数
        snp_total_count = 0

        with open(
            "{}/{}.alleles".format(modelPath, alleleFrqFile), "w", encoding="utf-8"
        ) as alleleFile:
            with open(
                "{}/{}.F".format(modelPath, alleleFrqFile), "w", encoding="utf-8"
            ) as frqFile:
                # 多线程遍历所有VCF文件
                with ThreadPoolExecutor(max_workers=maxWorkers) as t:
                    task_list = []
                    for chr in range(1, 23):
                        print("\tTask {} is launching...".format(chr))
                        task_list.append(
                            t.submit(
                                match_snp_from_vcf,
                                rsid_set,
                                VCF_FILENAME.format(chr),
                                afDigits,
                                stdDevThreshold,
                            )
                        )

                    # 异步等待线程，每个VCF文件处理结果
                    for task in as_completed(task_list):
                        allele_list, frq_list = task.result()

                        if len(allele_list) == len(frq_list):
                            snp_total_count += len(allele_list)
                            # 结果写入alleles文件
                            alleleFile.writelines(allele_list)
                            alleleFile.flush()
                            # 结果写入frequency文件
                            frqFile.writelines(frq_list)
                            frqFile.flush()
                            print(
                                "\tOne task success, {:,} SNPs have been saved in {}.alleles and {}.F respectively".format(
                                    len(allele_list), alleleFrqFile, alleleFrqFile
                                )
                            )
                        else:
                            print(
                                "\tOne task failed and not saved to file. The rows of alleles is {}, while the rows of frequency is {}".format(
                                    len(allele_list), len(frq_list)
                                )
                            )

        print(
            "All {} tasks finished. Totally {:,} SNPs have been saved in {}.alleles and {}.F respectively.".format(
                len(task_list), snp_total_count, alleleFrqFile, alleleFrqFile
            )
        )

    except FileNotFoundError:
        print("文件不存在！")
    except PermissionError:
        print("无权限访问文件！")
    except Exception as e:
        print(f"发生未知错误：{e}")


# 在指定的VCF文件中筛选符合条件的SNP
def match_snp_from_vcf(rsid_set, vcfFileName, afDigits, stdDevThreshold):
    # 当前VCF文件中筛选出的alleles和frequency列表
    allele_list = []
    frq_list = []

    # 遍历指定的VCF文件
    with open(
        vcfFileName,
        "r",
        encoding="utf-8",
    ) as vcfFile:
        print("Processing WBBC vcf file: " + vcfFileName)
        # 构造vcf字典，以RSID作为键，REF，ALT和INFO作为键值，只需要SNP，不包括Indel，生成vcf字典集
        vcfDict = {}
        for vcfLine in vcfFile.readlines():
            if len(vcfLine) > 0 and not vcfLine.startswith("#"):
                vcfLineList = vcfLine.split("\t")
                if (
                    len(vcfLineList) == 8
                    and vcfLineList[2] != "."
                    and len(vcfLineList[3]) == 1
                    and len(vcfLineList[4]) == 1
                ):
                    vcfDict[vcfLineList[2]] = {
                        "ref": vcfLineList[3],
                        "alt": vcfLineList[4],
                        "info": vcfLineList[7],
                    }

        # 遍历rsid模板集合，在vcf字典集查询对应的rsid，如有则加入
        for rsid in rsid_set:
            if rsid in vcfDict:
                infoList = vcfDict[rsid]["info"].split(";")
                if len(infoList) == 15:
                    af = float(infoList[1].split("=")[1])
                    if af != 0 and af != 1:
                        # 解析每个人群的alt allele频率
                        north_af = float(infoList[4].split("=")[1])
                        central_af = float(infoList[6].split("=")[1])
                        south_af = float(infoList[8].split("=")[1])
                        lingnan_af = float(infoList[10].split("=")[1])

                        # 只统计人群allele频率标准差大于阈值的SNP
                        if (
                            statistics.pstdev(
                                [
                                    north_af,
                                    central_af,
                                    south_af,
                                    lingnan_af,
                                ]
                            )
                            < stdDevThreshold
                        ):
                            continue

                        # 构造allele数据行，较大的是major突变, 较小的是minor突变
                        major_allele = (
                            vcfDict[rsid]["ref"] if af < 0.5 else vcfDict[rsid]["alt"]
                        )

                        minor_allele = (
                            vcfDict[rsid]["ref"] if af >= 0.5 else vcfDict[rsid]["alt"]
                        )

                        alleleLine = "{} {} {}\n".format(
                            rsid,
                            minor_allele,
                            major_allele,
                        )
                        allele_list.append(alleleLine)
                        # print("\tCollecting allele data: " + alleleLine)

                        # 构造frequency数据行，获取每个祖源成分的major突变频率
                        frqLine = "{} {} {} {}\n".format(
                            round(
                                (north_af if af >= 0.5 else 1 - north_af),
                                afDigits,
                            ),
                            round(
                                (central_af if af >= 0.5 else 1 - central_af),
                                afDigits,
                            ),
                            round(
                                (south_af if af >= 0.5 else 1 - south_af),
                                afDigits,
                            ),
                            round(
                                (lingnan_af if af >= 0.5 else 1 - lingnan_af),
                                afDigits,
                            ),
                        )
                        frq_list.append(frqLine)
                        # print("\tCollecting frequency data: " + frqLine)

    return allele_list, frq_list


# 从TSV文件获取RSID集合
def get_rsid(tsvFiles, highLDFileName):
    # TSV数据源的RSID集合
    rsid_set = set()

    # 获取高连锁不平衡数据
    hld_list = (
        get_high_ld(highLDFileName)
        if highLDFileName != "" and highLDFileName != None
        else []
    )
    rsid_in_hld_count = 0

    # 常染色体集合
    chr_set = {str(chr) for chr in range(1, 23)}

    if tsvFiles != None:
        for tsvFileName in tsvFiles:
            print("Collecting RSID in TSV file {0}...".format(tsvFileName))
            with open(tsvFileName, "r", encoding="utf-8") as tsvFile:
                for tsvLine in tsvFile.readlines():
                    if len(tsvLine) > 0 and not tsvLine.startswith(
                        ("#", "\n", "\t", '"')
                    ):
                        tsvLineList = tsvLine.split("\t")
                        if (
                            len(tsvLineList) >= 2
                            and tsvLineList[0] not in rsid_set
                            and tsvLineList[1] in chr_set
                        ):
                            # 去除高连锁不平衡的SNP
                            inHLD = False
                            for i in range(len(hld_list)):
                                if (
                                    tsvLineList[1] == hld_list[i]["chr"]
                                    and int(tsvLineList[2]) >= hld_list[i]["start_pos"]
                                    and int(tsvLineList[2]) <= hld_list[i]["end_pos"]
                                ):
                                    rsid_in_hld_count += 1
                                    inHLD = True
                                    break
                            if not inHLD:
                                rsid_set.add(tsvLineList[0])

    if len(rsid_set) > 0:
        print(
            "Finish collecting {:,} RSID in {}\nApart from {:,} RSID in high LD region.".format(
                len(rsid_set), " ".join(tsvFiles), rsid_in_hld_count
            )
        )
    else:
        raise Exception("There is no available SNP in TSV files.")

    return rsid_set


# 获取连锁不平衡数据high linkage disequilibrium
# https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)
def get_high_ld(ldFileName):
    if not os.access(ldFileName, os.F_OK):
        raise Exception(
            "High Linkage Disequilibrium file '{}' is not available.".format(ldFileName)
        )

    ld_list = []
    with open(
        ldFileName,
        "r",
        encoding="utf-8",
    ) as hdFile:
        for hdLine in hdFile.readlines():
            if len(hdLine) > 0 and not hdLine.startswith("#"):
                hdLineList = hdLine.split("\t")
                if len(hdLineList) == 3:
                    ld_list.append(
                        {
                            "chr": hdLineList[0],
                            "start_pos": int(hdLineList[1]),
                            "end_pos": int(hdLineList[2]),
                        }
                    )

    return ld_list
