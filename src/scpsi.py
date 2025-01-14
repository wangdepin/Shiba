import warnings
warnings.simplefilter('ignore')
import argparse
import logging
import sys
import os
import pandas as pd
from lib import shibalib

# Configure logging
logger = logging.getLogger(__name__)

def get_args():
    ## Get arguments from command line

    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = "PSI calculation for alternative splicing events in scRNA-seq data"
    )
    parser.add_argument("junctions", type = str, help = "A bed file of Junction read counts generated by bam2junc.sh")
    parser.add_argument("event", type = str, help = "Directory that contains text files of alternative splicing events generated by gtf2event.py")
    parser.add_argument("output", type = str, help = "Directory for output files")
    parser.add_argument("-p", "--num-process", type = int, help = "Number of processors to use", default = 1)
    parser.add_argument("-f", "--fdr", type = float, help = "FDR for detecting differential events", default = 0.05)
    parser.add_argument("-d", "--psi", type = float, help = "Threshold of delta PSI for detecting differential events", default = 0.1)
    parser.add_argument("-r", "--reference", type = str, help = "Reference group for detecting differential events")
    parser.add_argument("-a", "--alternative", type = str, help = "Alternative group for detecting differential events")
    parser.add_argument("-m", "--minimum-reads", type = int, help = "Minumum value of total reads for each junction for detecting differential events", default = 10)
    parser.add_argument("--onlypsi", help = "Just calculate PSI for each sample, not perform statistical tests", action = 'store_true')
    parser.add_argument("--excel", help = "Make result files in excel format", action = 'store_true')
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    args = parser.parse_args()
    return(args)

def main():
    ## Main
    args = get_args()

    # Set up logging
    logging.basicConfig(
        format="[%(asctime)s] %(levelname)7s %(message)s",
        level=logging.DEBUG if args.verbose else logging.INFO,
    )
    logger.info("Starting PSI calculation")
    logger.debug(args)

    # Parse arguments
    paths = {
        "junction": args.junctions,
        "event": args.event,
        "output": args.output,
    }
    params = {
        "num_process": args.num_process,
        "FDR": args.fdr,
        "dPSI": args.psi,
        "reference": args.reference,
        "alternative": args.alternative,
        "minimum_reads": args.minimum_reads,
        "onlypsi": args.onlypsi,
        "excel": args.excel,
    }

    # Load event and junction data
    logger.info("Loading event and junction files...")
    event_df_dict = shibalib.read_events(paths["event"])
    junc_df = shibalib.read_junctions(paths["junction"])
    junc_dict_all = shibalib.junc_dict(junc_df)
    sample_list = shibalib.make_sample_list(junc_df)
    junc_set = shibalib.make_junc_set(junc_df)

    # Group handling
    group_data = {}
    if not params["onlypsi"]:
        logger.info("Processing group information...")
        if params["reference"] in sample_list and params["alternative"] in sample_list:
            logger.info(f"{params['reference']} vs. {params['alternative']}")
            group_list = [params["reference"], params["alternative"]]
            junc_group_df = junc_df[["chr", "start", "end", "ID"] + group_list]
            junc_dict_group = shibalib.junc_dict(junc_group_df)
            group_data = {"group_list": group_list, "junc_dict_group": junc_dict_group, "group_df": junc_group_df, "sample_list_diff": sample_list}
        else:
            logger.error(f"Error: {params['reference']} or {params['alternative']} is not in the sample list")
            sys.exit(1)

    # Define event processing
    def process_event(event_type, event_func, func, col_func, diff_func=None, index_func=None):
        logger.info(f"Processing {event_type} events...")
        logger.debug(f"Event type: {event_type}")
        logger.debug(f"Event function: {event_func.__name__}")
        logger.debug(f"PSI function: {func.__name__}")
        logger.debug(f"Column function: {col_func.__name__}")
        logger.debug(f"Differential function: {diff_func.__name__ if diff_func else None}")
        logger.debug(f"Individual PSI function: {index_func.__name__ if index_func else None}")
        event_for_analysis_df = event_func(event_df_dict[event_type], junc_set)

        # Generate PSI tables
        psi_table_group_df, psi_table_sample_df = None, None
        if not params["onlypsi"]:
            psi_table_group_df = shibalib.make_psi_table_group(group_data["group_list"], event_for_analysis_df, group_data["junc_dict_group"], func, col_func, params["num_process"], params["minimum_reads"])
        else:
            psi_table_sample_df = shibalib.make_psi_table_sample(sample_list, event_for_analysis_df, junc_dict_all, func, col_func, params["num_process"], params["minimum_reads"])

        # Perform differential analysis
        diff_df = None
        nodiff_sample_df = None
        output_mtx_sample_df = None
        if not params["onlypsi"]:
            diff_df = shibalib.diff_event(
                event_for_analysis_df, psi_table_group_df, junc_dict_all, False,
                [params["reference"], params["alternative"]], sample_list,
                diff_func, index_func, params["num_process"], params["FDR"], params["dPSI"], False, False
            )
        else:
            # Generate PSI matrices
            nodiff_sample_df, output_mtx_sample_df = (shibalib.make_psi_mtx(psi_table_sample_df) if psi_table_sample_df is not None else (None, None))

        return {"nodiff_sample": nodiff_sample_df, "output_mtx_sample": output_mtx_sample_df, "diff": diff_df}

    # Define events and functions
    event_definitions = {
        "SE": (shibalib.event_for_analysis_se, shibalib.se, shibalib.col_se, shibalib.diff_se, shibalib.se_ind),
        "FIVE": (shibalib.event_for_analysis_five_three_afe_ale, shibalib.five_three_afe_ale, shibalib.col_five_three_afe_ale, shibalib.diff_five_three_afe_ale, shibalib.five_three_afe_ale_ind),
        "THREE": (shibalib.event_for_analysis_five_three_afe_ale, shibalib.five_three_afe_ale, shibalib.col_five_three_afe_ale, shibalib.diff_five_three_afe_ale, shibalib.five_three_afe_ale_ind),
        "MXE": (shibalib.event_for_analysis_mxe, shibalib.mxe, shibalib.col_mxe, shibalib.diff_mxe, shibalib.mxe_ind),
        "MSE": (shibalib.event_for_analysis_mse, shibalib.mse, shibalib.col_mse, shibalib.diff_mse, shibalib.mse_ind),
        "AFE": (shibalib.event_for_analysis_five_three_afe_ale, shibalib.five_three_afe_ale, shibalib.col_five_three_afe_ale, shibalib.diff_five_three_afe_ale, shibalib.five_three_afe_ale_ind),
        "ALE": (shibalib.event_for_analysis_five_three_afe_ale, shibalib.five_three_afe_ale, shibalib.col_five_three_afe_ale, shibalib.diff_five_three_afe_ale, shibalib.five_three_afe_ale_ind),
    }

    # Process each event
    event_results = {event: process_event(event, *functions) for event, functions in event_definitions.items()}

    # Save PSI matrices
    logger.info("Saving PSI matrices...")
    os.makedirs(paths["output"], exist_ok=True)

    if params["onlypsi"]:
        simple_psi_sample_df = pd.concat([result["output_mtx_sample"] for result in event_results.values()])
        simple_psi_sample_df.to_csv(f"{paths['output']}/PSI_matrix_sample.txt", sep="\t", index=False)
        for event, result in event_results.items():
            if result["nodiff_sample"] is not None:
                result["nodiff_sample"].to_csv(f"{paths['output']}/PSI_{event}.txt", sep="\t", index=False)
    else:
        for event, result in event_results.items():
            if result["diff"] is not None:
                result["diff"].to_csv(f"{paths['output']}/PSI_{event}.txt", sep="\t", index=False)

    # Save summary file
    logger.info("Saving summary file...")
    if not params["onlypsi"]:
        # Define the event types
        event_types = ["SE", "FIVE", "THREE", "MXE", "MSE", "AFE", "ALE"]
        # Collect event counts
        summary_l = []
        for event in event_types:
            logger.debug(f"Counting events for {event}...")
            event_counter = shibalib.EventCounter(event_results[event]["diff"], params["dPSI"])
            event_counts = event_counter.count_all_events()
            summary_l.extend([
                [event, "Up", "annotated", event_counts["up_annotated_num"]],
                [event, "Down", "annotated", event_counts["down_annotated_num"]],
                [event, "Up", "unannotated", event_counts["up_unannotated_num"]],
                [event, "Down", "unannotated", event_counts["down_unannotated_num"]],
            ])
        # Create and save the summary DataFrame
        summary_df = pd.DataFrame(
            summary_l,
            columns=["AS", "Direction", "Label", "Number"]
        )
        summary_df.to_csv(
            os.path.join(paths["output"], "summary.txt"),
            sep="\t",
            index=False,
        )

    # Optionally save to Excel
    if params["excel"]:
        logger.info("Exporting results to Excel...")
        results_to_save = [result["nodiff_sample"] if params["onlypsi"] else result["diff"] for result in event_results.values()]
        shibalib.save_excel_sc(paths["output"], *results_to_save)

    logger.info("All processes completed.")

if __name__ == '__main__':

    main()
