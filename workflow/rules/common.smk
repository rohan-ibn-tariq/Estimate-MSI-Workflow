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

# Print summary
print("=" * 60)
print("Workflow Configuration Summary")
print("=" * 60)
print(f"Data source: {config['data_source']}")
print(f"Samples to process: {get_active_samples()}")
print(f"Total: {len(get_active_samples())} samples")
print(f"Chromosome: {config['ref'].get('chromosome', 'all')}")
print("=" * 60)
