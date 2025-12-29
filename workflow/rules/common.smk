# ============================================
# Common Functions - Sample Management
# ============================================

import pandas as pd
import sys

samples_df = pd.read_csv("config/samples.tsv", sep="\t")


def get_active_samples():
    """
    Get list of samples to process with validation.

    Returns:
        list: Valid sample names

    Raises:
        ValueError: If requested samples don't exist in samples.tsv
    """
    requested = config.get("samples_to_run", None)
    available = set(samples_df["sample"].tolist())

    if requested is None:
        return list(available)
    else:
        requested_set = set(requested)
        invalid = requested_set - available

        if invalid:
            print(f"\n   ERROR: Invalid samples in config.yaml", file=sys.stderr)
            print(f"   Requested: {requested_set}", file=sys.stderr)
            print(f"   Available: {available}", file=sys.stderr)
            print(f"   Invalid: {invalid}", file=sys.stderr)
            sys.exit(1)

        return list(requested_set)


def get_process_type(wildcards):
    """
    Get how this sample should be processed (WXS or WGS).

    Returns:
        str: "WXS" or "WGS"
    """
    sample_row = samples_df[samples_df["sample"] == wildcards.sample].iloc[0]
    return sample_row["process_as"]


def get_microsatellite_bed(wildcards):
    """
    Get the correct microsatellite BED file for this sample.

    WXS samples: exonic microsatellites only
    WGS samples: all microsatellites
    """
    process_type = get_process_type(wildcards)

    if process_type == "WXS":
        return "resources/microsatellites/wxs_microsatellites.bed"
    else:
        return "resources/microsatellites/wgs_microsatellites.bed"


def get_required_ms_beds():
    """
    Figure out which microsatellite BED files we need to create.

    Only creates BEDs for process types actually used.
    """
    active_samples = get_active_samples()
    process_types = samples_df[samples_df["sample"].isin(active_samples)][
        "process_as"
    ].unique()

    beds = []
    if "WXS" in process_types:
        beds.append("resources/microsatellites/wxs_microsatellites.bed")
    if "WGS" in process_types:
        beds.append("resources/microsatellites/wgs_microsatellites.bed")

    return beds


# Print summary
print("=" * 60)
print("Workflow Configuration Summary")
print("=" * 60)
print(f"Samples to process: {get_active_samples()}")
print(f"Total: {len(get_active_samples())} samples")
print(f"Chromosome: {config['ref'].get('chromosome', 'all')}")
print("=" * 60)
