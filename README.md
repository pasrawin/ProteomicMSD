# TheModifiedNMF
Peak identification and quantification in proteomic mass spectrogram using modified non-negative matrix factorization

## Overview
Here, we developed a novel approach based on non-negative matrix factorization, in which the isotopic distribution profiles of peptide ions, the learned noise subspace, and the protein-peptide hierarchy for group sparsity constraint are incorporated into the original NMF to identify and quantify peaks in proteomic mass spectrogram.

## Requirements
Our modified NMF (mNMF) was written in Python 2.7 and tested on Window systems.
#### Dependencies
* numpy, pandas, scipy, collections, sklearn, and pyteomics

#### The inputs are: 
1. A .txt file of LC/MS experiment, converted by ProteoWizard MSconvert (required)
    * **Example:** Example_athousandions_msconvert.txt
2. A .xlsx of theoretical proteins and their sequences (required)
    * **Required column names:** Protein, Sequence
    * **Example:** Example_theoreticalproteins.xlsx
3. A .xlsx of peptide retention time calibration (optional)
    * **Required column names:** Protein, Sequence, PepRtimeStart, PepRtimePeak, PepRtimeEnd
    * **Example:** Example_RTcalibration.xlsx
#### The outputs are: 
1. A .pkl file of LC/MS experiment
2. A .xlsx file of NMF result
    * **Worksheets:** 1.Noisemean, 2.Peakresult, 3.XIC
3. Three .npz files of V, W, and H matrices

## Installation
git clone https://github.com/pasrawin/TheModifiedNMF.git
## Running
1. Prepare a directory containing your input files
2. Define your input files and parameters in ```mNMF00_handler.py```
    * The mass spectrograms from different instruments and proteomic experiments provide different features. In order to obtain an optimal result from mNMF, we strongly suggest that the parameters should be carefully set for each observed mass spectrogram. 
    * The *m/z* range, retention time range, smoothing window and shift, bins, *in silico* digestion criteria, number of best peaks reported etc. can be modified easily by replacing default values here.
3. Run ```python mNMF00_handler.py```
    * Yor command prompt will show mNMF process from 1 to 9

## Benchmark Datasets
For benchmarking in the paper, there are 2 standard protein mixtures
* A separation of a thousand ions: 1,057 peptide ions are generated from 476 unique peptides of 4 proteins
* A separation of nine thousand ions: 9,319 peptide ions are generated from 4,045 unique peptides of 48 proteins

## Support
If you have any questions about mNMF, please contact Pasrawin Taechawattananant (pasrawin@gmail.com)

## Citation


## License


## Acknowledgement
