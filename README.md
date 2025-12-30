# Estimate-MSI-Workflow
Estimate-MSI-Workflow is a complete pipeline that processes raw FASTQ files, calls indels, computes posterior variant probabilities using Varlociraptor, and derives a final MSI score for each sample using Varlociraptor's estimation msi subcommand.

## Setup (Temporary - Development Phase)
This workflow currently uses a custom varlociraptor build with MSI estimation support.
Once the feature is merged, this step(#1) won't be needed.

### 1. Create Local Config
Copy the template and add your paths:

```bash
cp config/config.local.yaml.template config/config.local.yaml
nano config/config.local.yaml
```

Edit with your cluster-specific paths to varlociraptor.

### 2. Run Workflow
```bash
snakemake -n --cores 12  # Dry run
snakemake --cores 12     # Real run
```
