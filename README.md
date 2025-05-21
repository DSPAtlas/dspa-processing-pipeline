# Processing LiP-MS and LiP-Quant Data for Dynaprot

## Step-by-Step Guide

1. *Prepare the YAML Configuration File*
Create a configuration file based on the provided params_template.yaml.
This file should include:

  - The path to the Spectronaut output file

  - Associated metadata required for database integration

2. *Run the Appropriate R Script from the Command Line*
For LiP/TrP Experiments:
```bash

Rscript LiPTrPpipeline.R path/to/your_params.yaml
```
For LiP-Quant Experiments with Dose-Response:
```bash
Rscript LiPQuantPipeline.R path/to/your_params.yaml
```

3. *Upload to Dynaprot*
After successful processing, use the generated output folder to upload results to the Dynaprot database on the virtual machine (VM).