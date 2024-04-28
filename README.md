# **eor_fisher_pipeline**

This repository holds some semi-automatic scripts for simulating a large suites of reionization lightcones, with parameters varying around some fiducial values over a number of realizations. The jobs can then be submitted to a slurm managed cluster for computation. This is useful for running Fisher analysis.

## Dependencies
- [21cmfast](https://github.com/21cmfast)
- [hiimtool](https://github.com/zhaotingchen/hiimtool) and dependencies therein

## Usage
1. Change the parameter values in [config file](config_template.py)
2. Run `python eor_pipeline.sh`. A `run.sh` file should appear (see `run_pipeline.sh` for example)
3. Run `sbatch run.sh`
4. If any jobs timed out or failed for no reason, you can run `python tidy_missing_parts.py` and then `bash tidy.sh` to finish any missing parts.