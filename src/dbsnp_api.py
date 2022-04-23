"""Utilities for retrieving data from dbSNP."""

import json
import time

import requests

DB_SNP = "snp"

# dbsnp efetch max refsnp results:
MAX_DBSNP_QUERIES = 15

# dbsnp API min time between requests.
SLEEP_SECONDS = 3

SUBCOL_DELIM = ";"

# Use this study for reporting population mean allele frequency (MAF).
PREFERED_AF_STUDY = "dbGaP_PopFreq"


def get_mafs_for_refsnps(snp_ids):
    """Given a dbSNP ID, fetches allele frequency data about the SNP from dbSNP.
    Returns a dict of dicts, where outer key is SNP ID, inner dict key is allele, and value is MAF for that allele for
    that SNP. dbSNP API has some rate limiting features that make it a bit tedious and slow to work with...
    """
    start = 0
    stop = len(snp_ids)
    data = {}
    # A batch may retrieve invalid JSON. If this happens, store the IDs and retry them afterwards one by one,
    # skipping any in the batch that have invalid JSON individually.
    retry_ids = []
    failed_id_count = 0
    print(f"Retrieving data from dbSNP...")
    while start < stop:
        cutoff = start + MAX_DBSNP_QUERIES
        cutoff = min(stop, cutoff)
        batch_ids = snp_ids[start:cutoff]

        print_progress_bar_(start + 1, stop)

        time.sleep(SLEEP_SECONDS)
        try:
            data_batch = get_mafs_for_refsnps_internal_(batch_ids)
            # Pythonic way to merge dicts (3.5+)
            data = {**data, **data_batch}
        except json.decoder.JSONDecodeError:
            retry_ids = retry_ids + batch_ids

        start += MAX_DBSNP_QUERIES

    print_progress_bar_(stop, stop)
    if len(retry_ids) > 0:
        print(f"Retrying {len(retry_ids)} SNPs...")

    for retry_id in retry_ids:
        try:
            data_batch = get_mafs_for_refsnps_internal_([retry_id])
            # Pythonic way to merge dicts (3.5+)
            data = {**data, **data_batch}
        except json.decoder.JSONDecodeError:
            # Just accept that we won't have this SNP.
            failed_id_count += 1

    if failed_id_count > 0:
        print(f"Unable to parse response for {failed_id_count} SNPs")
    return data


def get_mafs_for_refsnps_internal_(snp_ids):
    """Makes API call and parses response. Wrapped to avoid exceeding API request/response size limitations."""
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": DB_SNP,
        "id": snp_ids,
        "rettype": "json",
        "retmode": "text",
    }

    maf_dicts = {}
    request = requests.get(url=url, params=params)
    parseable_json = (
        "[{" + request.text[1:].replace('{"refsnp_id":', ',{"refsnp_id":') + "]"
    )
    response = json.loads(parseable_json)
    for snp_response in response:
        if "primary_snapshot_data" not in snp_response:
            continue

        snp_id = "rs" + snp_response["refsnp_id"]
        allele_to_maf = {}
        allele_annotations = snp_response["primary_snapshot_data"]["allele_annotations"]
        for allele_annotation in allele_annotations:
            frequencies = allele_annotation["frequency"]
            pop_freq_entries = [
                entry
                for entry in frequencies
                if entry["study_name"] == PREFERED_AF_STUDY
            ]
            if len(pop_freq_entries) == 0:
                continue
            pop_freq_entry = pop_freq_entries[0]
            pop_maf = pop_freq_entry["allele_count"] / pop_freq_entry["total_count"]
            allele = pop_freq_entry["observation"]["inserted_sequence"]
            allele_to_maf[allele] = pop_maf
        # Sometimes all alleles are reported, even with 0.0 values.
        # Drop those, assuming they are truly 0 and should be ignored.
        allele_to_maf = {
            key: value for key, value in allele_to_maf.items() if value > 0.0
        }
        maf_dicts[snp_id] = allele_to_maf

    return maf_dicts


# From https://stackoverflow.com/questions/3173320/text-progress-bar-in-terminal-with-block-characters
def print_progress_bar_(
    iteration,
    total,
    prefix="",
    suffix="",
    decimals=1,
    length=100,
    fill="█",
    printEnd="\r",
):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + "-" * (length - filled_length)
    print(f"\r{prefix} |{bar}| {percent}% {suffix}", end=printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()
