# Manuscript Figure Code

This repository contains the MATLAB and Python code used to generate the figures in the accompanying manuscript. It also includes the summary spreadsheets and example trained recurrent neural network (RNN) checkpoints used for the model analyses.

> **Manuscript:** *Title to be added*

## Repository structure

| Directory | Description |
| --- | --- |
| `analysis/` | MATLAB code for behavioral and neural-data analyses (Figures 1–4). |
| `mTDR_analysis/` | MATLAB code for model-based targeted dimensionality reduction (mTDR) analyses (Figures 5–6). |
| `RNN/` | Python code, summary spreadsheets, and trained-network checkpoints for the RNN analyses (Figure 7 and Supplementary Figures S10–S11). |

## Figure-to-code guide

### Main figures

| Figure panel(s) | Script(s) | Related files | Description |
| --- | --- | --- | --- |
| Figure 1C–E | `analysis/Psychometric_all.m` | Experimental behavioral data (not included) | Psychometric analysis and plotting. |
| Figure 2A–D | `analysis/CPturn_temporally.m` | Experimental neural data (not included) | Time-resolved choice-probability analysis. |
| Figure 2E–F | `analysis/CP_corr_2M.m` | Experimental neural data (not included) | Choice-probability comparison across two monkeys. |
| Figure 2G–H | `analysis/neural_proportion.m` | Experimental neural data (not included) | Proportions of neurons encoding perceptual choice, saccadic choice, or both. |
| Figure 3A | `analysis/signal_noise_corr.m` | Experimental neural data (not included) | Signal- and noise-correlation analysis. |
| Figure 3B | `analysis/noise_corr_contrast.m` | Experimental neural data (not included) | Comparison of noise correlations across conditions. |
| Figure 4A–F | `analysis/popu_CPcompare_new.m` | Experimental neural data (not included) | Population choice-probability analysis. |
| Figure 5B–E | `mTDR_analysis/brainarea_compare_new.m`<br>`mTDR_analysis/compare_starttime_p34_same_area.m`<br>`mTDR_analysis/time_test_new.m` | Processed mTDR results or experimental neural data (not included) | Comparison of temporal dynamics across brain areas. |
| Figure 6A, D | `mTDR_analysis/plot_angle.m` | Processed mTDR results (not included) | Angle analysis and plotting. |
| Figure 6B, C, E, F | `mTDR_analysis/axis_compare.m` | Processed mTDR results (not included) | Comparison of coding axes. |
| Figure 7A | `RNN/inputoutput.py` | Checkpoints in `RNN/checkpoint/` and `RNN/checkpoint_128_256/` | RNN input/output and connectivity analyses. |
| Figure 7B | `RNN/checkneuron2.py` | Checkpoints in `RNN/checkpoint/` and `RNN/checkpoint_128_256/` | Example-neuron and population analyses. |
| Figure 7C | `RNN/batch_comparemodule_temporal.py` | `RNN/choicetemporal_128_256.xlsx` | Temporal comparison of model choice signals. |
| Figure 7D | `RNN/batch_comparemodule_temporal.py` | `RNN/sactemporal_128_256.xlsx` | Temporal comparison of model saccadic signals. |
| Figure 7E | `RNN/batch_comparemodule_temporal.py` | `RNN/sactemporal_128_256_cut.xlsx` | Temporal comparison using the cut-network condition. |

### Supplementary figures

| Figure panel(s) | Script(s) | Related files | Description |
| --- | --- | --- | --- |
| Figure S10A | `RNN/inputoutput.py` | RNN checkpoints | RNN input/output analysis. |
| Figure S10B–D | `RNN/checkneuron2.py` | RNN checkpoints | Additional example-neuron and population analyses. |
| Figure S11A | `RNN/batch_compare_behavior.py`<br>`RNN/check_batch_behavior.py` | `RNN/behavior_128_256.xlsx`<br>`RNN/behavior_128_256_cut.xlsx` | Behavioral comparison of intact and cut-network conditions. |
| Figure S11B | `RNN/batch_comparemodule_temporal.py` | `RNN/choicetemporal_128_256_cut.xlsx` | Choice-signal comparison using the cut-network condition. |

## RNN checkpoints

The repository includes example trained networks used by the plotting and analysis scripts:

| Path | Description |
| --- | --- |
| `RNN/checkpoint_128_256/` | Example networks with colored-target input to the second hidden layer. |
| `RNN/checkpoint/color2h1.t7` | Example network with colored-target input to the first hidden layer. |
| `RNN/checkpoint/` | Additional trained-network checkpoints used by the RNN analyses. |

## Software requirements

The experimental-data and mTDR analyses require MATLAB. The RNN analyses require Python and use the following third-party packages:

- NumPy
- pandas
- Matplotlib
- SciPy
- seaborn
- PyTorch
- six

Exact software versions used for the manuscript should be added before archival release.

## Data and paths

The raw experimental data are not included in this repository. Several MATLAB scripts currently contain local Windows paths such as `Z:\Data\...`. Before running these scripts, replace those paths with the location of the corresponding data on your system.

Add the manuscript's data-access statement or public dataset URL here before creating the archival release.

## Running the analyses

1. Clone or download this repository.
2. Open MATLAB in `analysis/` or `mTDR_analysis/`, or run the Python scripts from `RNN/` so that local helper modules can be found.
3. Update the data paths in the relevant scripts.
4. Use the tables above to identify the script and associated file for each figure panel.

The scripts reflect the analysis workflow used for the manuscript and may require manuscript-specific intermediate data that are not distributed here.

## Citation

A DOI and full citation will be added after the first archived release is published through Zenodo.

## License

License information will be added before the archival release.
