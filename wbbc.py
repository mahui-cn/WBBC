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
    tsv_files: list,
    high_ld_filename: str = "",
    model_path: str = ".",
    allele_frq_file: str = "wbbc",
    af_digits: int = 6,
    std_dev_threshold: float = 0.03,
    max_workers: int = 4,
) -> None:
    try:
        # 从TSV文件获取RSID集合
        rsid_set = get_rsid(tsv_files, high_ld_filename)

        # SNP总计数
        snp_total_count = 0

        with open(
            "{}/{}.alleles".format(model_path, allele_frq_file), "w", encoding="utf-8"
        ) as allele_file:
            with open(
                "{}/{}.F".format(model_path, allele_frq_file), "w", encoding="utf-8"
            ) as frq_file:
                # 多线程遍历所有VCF文件
                with ThreadPoolExecutor(max_workers=max_workers) as t:
                    task_list = []
                    for chr in range(1, 23):
                        print("\tTask {} is launching...".format(chr))
                        task_list.append(
                            t.submit(
                                match_snp_from_vcf,
                                rsid_set,
                                VCF_FILENAME.format(chr),
                                af_digits,
                                std_dev_threshold,
                            )
                        )

                    # 异步等待线程，每个VCF文件处理结果
                    for task in as_completed(task_list):
                        allele_list, frq_list = task.result()

                        if len(allele_list) == len(frq_list):
                            snp_total_count += len(allele_list)
                            # 结果写入alleles文件
                            allele_file.writelines(allele_list)
                            allele_file.flush()
                            # 结果写入frequency文件
                            frq_file.writelines(frq_list)
                            frq_file.flush()
                            print(
                                "\tOne task success, {:,} SNPs have been saved in {}.alleles and {}.F respectively".format(
                                    len(allele_list), allele_frq_file, allele_frq_file
                                )
                            )
                        else:
                            print(
                                "\tOne task failed and not saved to file due to not matched amount between alleles and frequency. The rows of alleles is {}, while the rows of frequency is {}".format(
                                    len(allele_list), len(frq_list)
                                )
                            )

        print(
            "All {} tasks finished. Totally {:,} SNPs have been saved in {}.alleles and {}.F respectively.".format(
                len(task_list), snp_total_count, allele_frq_file, allele_frq_file
            )
        )

    except FileNotFoundError:
        print("文件不存在！")
    except PermissionError:
        print("无权限访问文件！")
    except Exception as e:
        print(f"发生未知错误：{e}")


# 在指定的VCF文件中筛选符合条件的SNP
def match_snp_from_vcf(
    rsid_set: set, vcf_filename: str, af_digits: int, std_dev_threshold: float
):
    # 当前VCF文件中筛选出的alleles和frequency列表
    allele_list = []
    frq_list = []

    # 遍历指定的VCF文件
    with open(
        vcf_filename,
        "r",
        encoding="utf-8",
    ) as vcf_file:
        print("\nProcessing WBBC vcf file: " + vcf_filename)
        # 构造vcf字典，以RSID作为键，REF，ALT和INFO作为键值，只需要SNP，不包括Indel，生成vcf字典集
        vcf_dict = {}
        for vcf_line in vcf_file.readlines():
            if len(vcf_line) > 0 and not vcf_line.startswith("#"):
                vcf_line_list = vcf_line.split("\t")
                if (
                    len(vcf_line_list) == 8
                    and vcf_line_list[2] != "."
                    and len(vcf_line_list[3]) == 1
                    and len(vcf_line_list[4]) == 1
                ):
                    vcf_dict[vcf_line_list[2]] = {
                        "ref": vcf_line_list[3],
                        "alt": vcf_line_list[4],
                        "info": vcf_line_list[7],
                    }

        # 遍历rsid模板集合，在vcf字典集查询对应的rsid，如有则加入
        for rsid in rsid_set:
            if rsid in vcf_dict:
                info_list = vcf_dict[rsid]["info"].split(";")
                if len(info_list) == 15:
                    af = float(info_list[1].split("=")[1])
                    if af != 0 and af != 1:
                        # 解析每个人群的alt allele频率
                        north_af = float(info_list[4].split("=")[1])
                        central_af = float(info_list[6].split("=")[1])
                        south_af = float(info_list[8].split("=")[1])
                        lingnan_af = float(info_list[10].split("=")[1])

                        # 过滤allele频率标准差大于阈值的SNP
                        if (
                            statistics.pstdev(
                                [
                                    north_af,
                                    central_af,
                                    south_af,
                                    lingnan_af,
                                ]
                            )
                            < std_dev_threshold
                        ):
                            continue

                        # 构造allele数据行，A1是ref allele, A2是alt allele
                        allele_line = "{} {} {}\n".format(
                            rsid,
                            vcf_dict[rsid]["ref"],
                            vcf_dict[rsid]["alt"],
                        )
                        allele_list.append(allele_line)
                        # print("\tCollecting allele data: " + allele_line)

                        # frequency数据行，每个祖源成分的ref allele频率
                        frq_line = "{} {} {} {}\n".format(
                            round(
                                1 - north_af,
                                af_digits,
                            ),
                            round(
                                1 - central_af,
                                af_digits,
                            ),
                            round(
                                1 - south_af,
                                af_digits,
                            ),
                            round(
                                1 - lingnan_af,
                                af_digits,
                            ),
                        )
                        frq_list.append(frq_line)
                        # print("\tCollecting frequency data: " + frq_line)

    return allele_list, frq_list


# 从TSV文件获取RSID集合
def get_rsid(tsv_files: list, high_ld_filename: str) -> set:
    # TSV数据源的RSID集合
    rsid_set = set()

    # 获取高连锁不平衡数据
    hld_list = (
        get_high_ld(high_ld_filename)
        if high_ld_filename != "" and high_ld_filename != None
        else []
    )
    rsid_in_hld_count = 0

    # 常染色体集合
    chr_set = {str(chr) for chr in range(1, 23)}

    if tsv_files != None:
        for tsv_filename in tsv_files:
            print("Collecting RSID in TSV file {0}...".format(tsv_filename))
            with open(tsv_filename, "r", encoding="utf-8") as tsv_file:
                for tsv_line in tsv_file.readlines():
                    if len(tsv_line) > 0 and not tsv_line.startswith(
                        ("#", "\n", "\t", '"')
                    ):
                        tsv_line_list = tsv_line.split("\t")
                        if (
                            len(tsv_line_list) >= 2
                            and tsv_line_list[0] not in rsid_set
                            and tsv_line_list[1] in chr_set
                        ):
                            # 去除连锁不平衡区域的SNP
                            inHLD = False
                            for i in range(len(hld_list)):
                                if (
                                    tsv_line_list[1] == hld_list[i]["chr"]
                                    and int(tsv_line_list[2])
                                    >= hld_list[i]["start_pos"]
                                    and int(tsv_line_list[2]) <= hld_list[i]["end_pos"]
                                ):
                                    rsid_in_hld_count += 1
                                    inHLD = True
                                    break
                            if not inHLD:
                                rsid_set.add(tsv_line_list[0])

    if len(rsid_set) > 0:
        print(
            "Finish collecting {:,} RSID in {}\nApart from {:,} RSID in high LD region.".format(
                len(rsid_set), " ".join(tsv_files), rsid_in_hld_count
            )
        )
    else:
        raise Exception("There is no available SNP in TSV files.")

    return rsid_set


# 获取连锁不平衡数据high linkage disequilibrium
# https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)
def get_high_ld(ld_filename: str) -> list:
    if not os.access(ld_filename, os.F_OK):
        raise Exception(
            "High Linkage Disequilibrium file '{}' is not available.".format(
                ld_filename
            )
        )

    ld_list = []
    with open(
        ld_filename,
        "r",
        encoding="utf-8",
    ) as hd_file:
        for hd_line in hd_file.readlines():
            if len(hd_line) > 0 and not hd_line.startswith("#"):
                hd_line_list = hd_line.split("\t")
                if len(hd_line_list) == 3:
                    ld_list.append(
                        {
                            "chr": hd_line_list[0],
                            "start_pos": int(hd_line_list[1]),
                            "end_pos": int(hd_line_list[2]),
                        }
                    )

    return ld_list
