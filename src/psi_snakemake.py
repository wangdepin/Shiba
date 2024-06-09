import warnings
warnings.simplefilter('ignore')
import argparse

def str2bool(v):
    if v == "True":
        return True
    elif v == "False":
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def get_args():
    ## Get arguments from command line

    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = "PSI calculation for alternative splicing events"
    )

    parser.add_argument("junctions", type = str, help = "A bed file of Junction read counts generated by bam2junc.sh")
    parser.add_argument("event", type = str, help = "Directory that contains text files of alternative splicing events generated by gtf2event.py")
    parser.add_argument("output", type = str, help = "Directory for output files")
    parser.add_argument("-p", "--num-process", type = int, help = "Number of processors to use", default = 1)
    parser.add_argument("-g", "--group", type = str, help = "Group information for detecting differential events")
    parser.add_argument("-f", "--fdr", type = float, help = "FDR for detecting differential events", default = 0.05)
    parser.add_argument("-d", "--psi", type = float, help = "Threshold of delta PSI for detecting differential events", default = 0.1)
    parser.add_argument("-r", "--reference", type = str, help = "Reference group for detecting differential events")
    parser.add_argument("-a", "--alternative", type = str, help = "Alternative group for detecting differential events")
    parser.add_argument("-m", "--minimum-reads", type = int, help = "Minumum value of total reads for each junction for detecting differential events", default = 10)
    parser.add_argument("-i", "--individual-psi", help = "Print PSI for individual samples to output files", type = str2bool, nargs = "?", const = True, default = False)
    parser.add_argument("-t", "--ttest", help = "Perform Welch's t-test between reference and alternative group", type = str2bool, nargs = "?", const = True, default = False)
    parser.add_argument("--onlypsi", help = "Just calculate PSI for each sample, not perform statistical tests", type = str2bool, nargs = "?", const = True, default = False)
    parser.add_argument("--onlypsi-group", help = "Just calculate PSI for each group, not perform statistical tests (Overrides --onlypsi when used together)", type = str2bool, nargs = "?", const = True, default = False)
    parser.add_argument("--excel", help = "Make result files in excel format", type = str2bool, nargs = "?", const = True, default = False)

    args = parser.parse_args()

    return(args)


def main():
    ## Main

    args = get_args()
    junction_path = args.junctions
    event_path = args.event
    output_path = args.output
    group_path = args.group
    num_process = args.num_process
    FDR = args.fdr
    dPSI = args.psi
    reference = args.reference
    alternative = args.alternative
    minimum_reads = args.minimum_reads
    individual_psi = args.individual_psi
    ttest_bool = args.ttest
    onlypsi = args.onlypsi
    onlypsi_group = args.onlypsi_group
    excel = args.excel

    # start_time = time.time()  # Time

    # Load modules
    import sys
    import os
    current_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
    sys.path.append(parent_dir)
    import time
    import pandas as pd
    from lib import shibalib

    # Read events
    print("Load event files....", file = sys.stdout)
    event_df_dict = shibalib.read_events(event_path)

    # Read junction count
    print("Load junction files....", file = sys.stdout)
    junc_df = shibalib.read_junctions(junction_path)

    # Dictionary for junction read counts from each sample
    junc_dict_all = shibalib.junc_dict(junc_df)

    # Get sample list
    sample_list = shibalib.make_sample_list(junc_df)
    # print("Samples: " + str(sample_list), file = sys.stdout)

    # Sum read count
    if onlypsi_group or group_path:

        print("Sum read count for each group....", file = sys.stdout)
        # Read group information
        group_df = shibalib.read_group(group_path)
        # Set groups
        group_list = shibalib.set_group(group_df, onlypsi_group, reference, alternative)
        print("Groups: " + str(group_list), file = sys.stdout)
        # Get sample list for each group
        if onlypsi_group == False:
            sample_list = shibalib.sample_in_group_list(group_df, group_list)
        # print("Samples: " + str(sample_list), file = sys.stdout)
        # Sum read count for each group
        junc_group_df = shibalib.sum_reads(onlypsi_group, junc_df, group_df, group_list)
        # print("Sum read count for each group" + str(junc_group_df), file = sys.stdout)
        # Dictionary for junction read counts from each group
        junc_dict_group = shibalib.junc_dict(junc_group_df)

    # Set of junctions
    junc_set = shibalib.make_junc_set(junc_df)

    # print("Load files: " + str(time.time() - start_time), file = sys.stdout)

    ################################################################################
    # Skipped exon
    # start_time = time.time()  # Time
    print("PSI for each skipped exon....", file = sys.stdout)

    event_for_analysis_df = shibalib.event_for_analysis_se(event_df_dict["SE"], junc_set)

    # Make PSI matrix
    if onlypsi_group:

        psi_table_group_df = shibalib.make_psi_table_group(
            group_list, event_for_analysis_df, junc_dict_group, shibalib.se,
            shibalib.col_se, num_process, minimum_reads
        )

    elif onlypsi:

        psi_table_sample_df = shibalib.make_psi_table_sample(
            sample_list, event_for_analysis_df, junc_dict_all, shibalib.se,
            shibalib.col_se, num_process, minimum_reads
        )

    else:

        psi_table_group_df = shibalib.make_psi_table_group(
            group_list, event_for_analysis_df, junc_dict_group, shibalib.se,
            shibalib.col_se, num_process, minimum_reads
        )

        psi_table_sample_df = shibalib.make_psi_table_sample(
            sample_list, event_for_analysis_df, junc_dict_all, shibalib.se,
            shibalib.col_se, num_process, minimum_reads
        )

    # PSI matrix for each group
    if onlypsi_group:

        SE_nodiff_group_df, output_mtx_SE_group_df = shibalib.make_psi_mtx(psi_table_group_df)

    # PSI matrix for each sample
    elif onlypsi:

        SE_nodiff_sample_df, output_mtx_SE_sample_df = shibalib.make_psi_mtx(psi_table_sample_df)

    # Differential analysis and PSI matrix
    else:

        SE_diff_df = shibalib.diff_event(
            event_for_analysis_df, psi_table_group_df, junc_dict_all,
            group_df, group_list, sample_list,
            shibalib.diff_se, shibalib.se_ind, num_process,
            FDR, dPSI, individual_psi, ttest_bool
        )

        SE_nodiff_sample_df, output_mtx_SE_sample_df = shibalib.make_psi_mtx(psi_table_sample_df)
        SE_nodiff_group_df, output_mtx_SE_group_df = shibalib.make_psi_mtx(psi_table_group_df)

    ################################################################################
    # Alternative five ss
    print("PSI for each alternative five prime ss....", file = sys.stdout)

    event_for_analysis_df = shibalib.event_for_analysis_five_three(event_df_dict["FIVE"], junc_set)

    # Make PSI matrix
    if onlypsi_group:

        psi_table_group_df = shibalib.make_psi_table_group(
            group_list, event_for_analysis_df, junc_dict_group, shibalib.five,
            shibalib.col_five_three, num_process, minimum_reads
        )

    elif onlypsi:

        psi_table_sample_df = shibalib.make_psi_table_sample(
            sample_list, event_for_analysis_df, junc_dict_all, shibalib.five,
            shibalib.col_five_three, num_process, minimum_reads
        )

    else:

        psi_table_group_df = shibalib.make_psi_table_group(
            group_list, event_for_analysis_df, junc_dict_group, shibalib.five,
            shibalib.col_five_three, num_process, minimum_reads
        )

        psi_table_sample_df = shibalib.make_psi_table_sample(
            sample_list, event_for_analysis_df, junc_dict_all, shibalib.five,
            shibalib.col_five_three, num_process, minimum_reads
        )

    # PSI matrix for each group
    if onlypsi_group:

        FIVE_nodiff_group_df, output_mtx_FIVE_group_df = shibalib.make_psi_mtx(psi_table_group_df)

    # PSI matrix for each sample
    elif onlypsi:

        FIVE_nodiff_sample_df, output_mtx_FIVE_sample_df = shibalib.make_psi_mtx(psi_table_sample_df)

    # Differential analysis and PSI matrix
    else:

        FIVE_diff_df = shibalib.diff_event(
            event_for_analysis_df, psi_table_group_df, junc_dict_all,
            group_df, group_list, sample_list,
            shibalib.diff_five_three, shibalib.five_ind, num_process,
            FDR, dPSI, individual_psi, ttest_bool
        )

        FIVE_nodiff_sample_df, output_mtx_FIVE_sample_df = shibalib.make_psi_mtx(psi_table_sample_df)
        FIVE_nodiff_group_df, output_mtx_FIVE_group_df = shibalib.make_psi_mtx(psi_table_group_df)

    ################################################################################
    # Alternative three ss
    print("PSI for each alternative three prime ss....", file = sys.stdout)

    event_for_analysis_df = shibalib.event_for_analysis_five_three(event_df_dict["THREE"], junc_set)

    # Make PSI matrix
    if onlypsi_group:

        psi_table_group_df = shibalib.make_psi_table_group(
            group_list, event_for_analysis_df, junc_dict_group, shibalib.three,
            shibalib.col_five_three, num_process, minimum_reads
        )

    elif onlypsi:

        psi_table_sample_df = shibalib.make_psi_table_sample(
            sample_list, event_for_analysis_df, junc_dict_all, shibalib.three,
            shibalib.col_five_three, num_process, minimum_reads
        )

    else:

        psi_table_group_df = shibalib.make_psi_table_group(
            group_list, event_for_analysis_df, junc_dict_group, shibalib.three,
            shibalib.col_five_three, num_process, minimum_reads
        )

        psi_table_sample_df = shibalib.make_psi_table_sample(
            sample_list, event_for_analysis_df, junc_dict_all, shibalib.three,
            shibalib.col_five_three, num_process, minimum_reads
        )

    # PSI matrix for each group
    if onlypsi_group:

        THREE_nodiff_group_df, output_mtx_THREE_group_df = shibalib.make_psi_mtx(psi_table_group_df)

    # PSI matrix for each sample
    elif onlypsi:

        THREE_nodiff_sample_df, output_mtx_THREE_sample_df = shibalib.make_psi_mtx(psi_table_sample_df)

    # Differential analysis and PSI matrix
    else:

        THREE_diff_df = shibalib.diff_event(
            event_for_analysis_df, psi_table_group_df, junc_dict_all,
            group_df, group_list, sample_list,
            shibalib.diff_five_three, shibalib.three_ind, num_process,
            FDR, dPSI, individual_psi, ttest_bool
        )

        THREE_nodiff_sample_df, output_mtx_THREE_sample_df = shibalib.make_psi_mtx(psi_table_sample_df)
        THREE_nodiff_group_df, output_mtx_THREE_group_df = shibalib.make_psi_mtx(psi_table_group_df)

    ################################################################################
    # Mutually exclusive exons
    print("PSI for each mutually exclusive exons....", file = sys.stdout)

    event_for_analysis_df = shibalib.event_for_analysis_mxe(event_df_dict["MXE"], junc_set)

    # Make PSI matrix
    if onlypsi_group:

        psi_table_group_df = shibalib.make_psi_table_group(
            group_list, event_for_analysis_df, junc_dict_group, shibalib.mxe,
            shibalib.col_mxe, num_process, minimum_reads
        )

    elif onlypsi:

        psi_table_sample_df = shibalib.make_psi_table_sample(
            sample_list, event_for_analysis_df, junc_dict_all, shibalib.mxe,
            shibalib.col_mxe, num_process, minimum_reads
        )

    else:

        psi_table_group_df = shibalib.make_psi_table_group(
            group_list, event_for_analysis_df, junc_dict_group, shibalib.mxe,
            shibalib.col_mxe, num_process, minimum_reads
        )

        psi_table_sample_df = shibalib.make_psi_table_sample(
            sample_list, event_for_analysis_df, junc_dict_all, shibalib.mxe,
            shibalib.col_mxe, num_process, minimum_reads
        )

    # PSI matrix for each group
    if onlypsi_group:

        MXE_nodiff_group_df, output_mtx_MXE_group_df = shibalib.make_psi_mtx(psi_table_group_df)

    # PSI matrix for each sample
    elif onlypsi:

        MXE_nodiff_sample_df, output_mtx_MXE_sample_df = shibalib.make_psi_mtx(psi_table_sample_df)

    # Differential analysis and PSI matrix
    else:

        MXE_diff_df = shibalib.diff_event(
            event_for_analysis_df, psi_table_group_df, junc_dict_all,
            group_df, group_list, sample_list,
            shibalib.diff_mxe, shibalib.mxe_ind, num_process,
            FDR, dPSI, individual_psi, ttest_bool
        )

        MXE_nodiff_sample_df, output_mtx_MXE_sample_df = shibalib.make_psi_mtx(psi_table_sample_df)
        MXE_nodiff_group_df, output_mtx_MXE_group_df = shibalib.make_psi_mtx(psi_table_group_df)

    ################################################################################
    # Retained intron
    print("PSI for each retained intron....", file = sys.stdout)

    event_for_analysis_df = shibalib.event_for_analysis_ri(event_df_dict["RI"], junc_set)

    # Make PSI matrix
    if onlypsi_group:

        psi_table_group_df = shibalib.make_psi_table_group(
            group_list, event_for_analysis_df, junc_dict_group, shibalib.ri,
            shibalib.col_ri, num_process, minimum_reads
        )

    elif onlypsi:

        psi_table_sample_df = shibalib.make_psi_table_sample(
            sample_list, event_for_analysis_df, junc_dict_all, shibalib.ri,
            shibalib.col_ri, num_process, minimum_reads
        )

    else:

        psi_table_group_df = shibalib.make_psi_table_group(
            group_list, event_for_analysis_df, junc_dict_group, shibalib.ri,
            shibalib.col_ri, num_process, minimum_reads
        )

        psi_table_sample_df = shibalib.make_psi_table_sample(
            sample_list, event_for_analysis_df, junc_dict_all, shibalib.ri,
            shibalib.col_ri, num_process, minimum_reads
        )

    # PSI matrix for each group
    if onlypsi_group:

        RI_nodiff_group_df, output_mtx_RI_group_df = shibalib.make_psi_mtx(psi_table_group_df)

    # PSI matrix for each sample
    elif onlypsi:

        RI_nodiff_sample_df, output_mtx_RI_sample_df = shibalib.make_psi_mtx(psi_table_sample_df)

    # Differential analysis and PSI matrix
    else:

        RI_diff_df = shibalib.diff_event(
            event_for_analysis_df, psi_table_group_df, junc_dict_all,
            group_df, group_list, sample_list,
            shibalib.diff_ri, shibalib.ri_ind, num_process,
            FDR, dPSI, individual_psi, ttest_bool
        )

        RI_nodiff_sample_df, output_mtx_RI_sample_df = shibalib.make_psi_mtx(psi_table_sample_df)
        RI_nodiff_group_df, output_mtx_RI_group_df = shibalib.make_psi_mtx(psi_table_group_df)

    ################################################################################
    # Make directory
    os.makedirs(output_path, exist_ok = True)
    ################################################################################
    # Save PSI matrix
    if onlypsi_group:

        simple_psi_group_df = pd.concat(
            [output_mtx_SE_group_df, output_mtx_FIVE_group_df, output_mtx_THREE_group_df, output_mtx_MXE_group_df, output_mtx_RI_group_df]
        )

        simple_psi_group_df.to_csv(

            output_path + "/PSI_matrix_group.txt",
            sep = "\t",
            index = False,

        )

        del simple_psi_group_df

    elif onlypsi:

        simple_psi_sample_df = pd.concat(
            [output_mtx_SE_sample_df, output_mtx_FIVE_sample_df, output_mtx_THREE_sample_df, output_mtx_MXE_sample_df, output_mtx_RI_sample_df]
        )

        simple_psi_sample_df.to_csv(

            output_path + "/PSI_matrix_sample.txt",
            sep = "\t",
            index = False,

        )

        del simple_psi_sample_df

    else:

        simple_psi_group_df = pd.concat(
            [output_mtx_SE_group_df, output_mtx_FIVE_group_df, output_mtx_THREE_group_df, output_mtx_MXE_group_df, output_mtx_RI_group_df]
        )
        simple_psi_sample_df = pd.concat(
            [output_mtx_SE_sample_df, output_mtx_FIVE_sample_df, output_mtx_THREE_sample_df, output_mtx_MXE_sample_df, output_mtx_RI_sample_df]
        )

        simple_psi_group_df.to_csv(

            output_path + "/PSI_matrix_group.txt",
            sep = "\t",
            index = False,

        )

        simple_psi_sample_df.to_csv(

            output_path + "/PSI_matrix_sample.txt",
            sep = "\t",
            index = False,

        )

        del simple_psi_group_df
        del simple_psi_sample_df

    ################################################################################
    # Save summary file (up-regulated and down-regulated exons)
    if group_path and (onlypsi == False) and (onlypsi_group == False):

        event_counter_SE = shibalib.EventCounter(SE_diff_df, dPSI)
        event_counts_SE = event_counter_SE.count_all_events()
        event_counter_FIVE = shibalib.EventCounter(FIVE_diff_df, dPSI)
        event_counts_FIVE = event_counter_FIVE.count_all_events()
        event_counter_THREE = shibalib.EventCounter(THREE_diff_df, dPSI)
        event_counts_THREE = event_counter_THREE.count_all_events()
        event_counter_MXE = shibalib.EventCounter(MXE_diff_df, dPSI)
        event_counts_MXE = event_counter_MXE.count_all_events()
        event_counter_RI = shibalib.EventCounter(RI_diff_df, dPSI)
        event_counts_RI = event_counter_RI.count_all_events()

        summary_l = [
            ["SE", "Up", "annotated", event_counts_SE["up_annotated_num"]],
            ["SE", "Down", "annotated", event_counts_SE["down_annotated_num"]],
            ["SE", "Up", "unannotated", event_counts_SE["up_unannotated_num"]],
            ["SE", "Down", "unannotated", event_counts_SE["down_unannotated_num"]],
            ["FIVE", "Up", "annotated", event_counts_FIVE["up_annotated_num"]],
            ["FIVE", "Down", "annotated", event_counts_FIVE["down_annotated_num"]],
            ["FIVE", "Up", "unannotated", event_counts_FIVE["up_unannotated_num"]],
            ["FIVE", "Down", "unannotated", event_counts_FIVE["down_unannotated_num"]],
            ["THREE", "Up", "annotated", event_counts_THREE["up_annotated_num"]],
            ["THREE", "Down", "annotated", event_counts_THREE["down_annotated_num"]],
            ["THREE", "Up", "unannotated", event_counts_THREE["up_unannotated_num"]],
            ["THREE", "Down", "unannotated", event_counts_THREE["down_unannotated_num"]],
            ["MXE", "Up", "annotated", event_counts_MXE["up_annotated_num"]],
            ["MXE", "Down", "annotated", event_counts_MXE["down_annotated_num"]],
            ["MXE", "Up", "unannotated", event_counts_MXE["up_unannotated_num"]],
            ["MXE", "Down", "unannotated", event_counts_MXE["down_unannotated_num"]],
            ["RI", "Up", "annotated", event_counts_RI["up_annotated_num"]],
            ["RI", "Down", "annotated", event_counts_RI["down_annotated_num"]],
            ["RI", "Up", "unannotated", event_counts_RI["up_unannotated_num"]],
            ["RI", "Down", "unannotated", event_counts_RI["down_unannotated_num"]]
        ]

        summary_df = pd.DataFrame(

            summary_l,
            columns = ["AS", "Direction", "Label", "Number"]

        )

        summary_df.to_csv(

            output_path + "/summary.txt",
            sep = "\t",
            index = False,

        )

    ###############################################################################
    # Save results to files
    if onlypsi_group:

        SE_df = SE_nodiff_group_df
        FIVE_df = FIVE_nodiff_group_df
        THREE_df = THREE_nodiff_group_df
        MXE_df = MXE_nodiff_group_df
        RI_df = RI_nodiff_group_df

    elif onlypsi:

        SE_df = SE_nodiff_sample_df
        FIVE_df = FIVE_nodiff_sample_df
        THREE_df = THREE_nodiff_sample_df
        MXE_df = MXE_nodiff_sample_df
        RI_df = RI_nodiff_sample_df

    else:

        SE_df = SE_diff_df
        FIVE_df = FIVE_diff_df
        THREE_df = THREE_diff_df
        MXE_df = MXE_diff_df
        RI_df = RI_diff_df

    SE_df.to_csv(

        output_path + "/PSI_SE.txt",
        sep = "\t",
        index = False

    )

    FIVE_df.to_csv(

        output_path + "/PSI_FIVE.txt",
        sep = "\t",
        index = False

    )

    THREE_df.to_csv(

        output_path + "/PSI_THREE.txt",
        sep = "\t",
        index = False

    )

    MXE_df.to_csv(

        output_path + "/PSI_MXE.txt",
        sep = "\t",
        index = False

    )

    RI_df.to_csv(

        output_path + "/PSI_RI.txt",
        sep = "\t",
        index = False

    )

    ################################################################################
    # Excel file
    if excel:

        print("Export to an excel file....", file = sys.stdout)
        shibalib.save_excel(output_path, SE_df, FIVE_df, THREE_df, MXE_df, RI_df)

    ################################################################################

    print("Output file: " + output_path, file = sys.stdout)


if __name__ == '__main__':

    main()

