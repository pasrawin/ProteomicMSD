# Non-negative matrix factorization for proteomic mass spectrogram
Peak identification and quantification in proteomic mass spectrogram using modified non-negative matrix factorization

## Overview
Here, we developed a novel approach based on non-negative matrix factorization, in which the isotopic distribution profiles of peptide ions, the learned noise subspace, and the protein-peptide hierarchy for group sparsity constraint are incorporated into the original NMF to identify and quantify peaks in proteomic mass spectrogram.

## Requirements
Our proteomic NMF (protNMF) was written in Python 2.7 and tested on Window systems.
#### Dependencies
* division from future, numpy, pandas, scipy, math, collections, itertools, sklearn, gc, and pyteomics

#### The inputs are: 
1. An .mzML or a .txt file of LC/MS experiment converted by ProteoWizard MSconvert, or a .pkl file of protNMF output (required)
    * **Example:** Example_mNMF_standard4proteins.txt
    * **Example:** Example_mNMF_standard4proteins.pkl (download [here](https://repository.jpostdb.org/preview/4730115445e7091ad1182a). Currently, available for reviewers only 4330)
       * After running protNMF, .pkl will be obtained. You may reuse this in following protNMF run to save the computational time.
2. A .FASTA or an .xlsx of theoretical proteins and their sequences (required)
    * **Required column names:** Protein, Sequence
    * **Example:** Example_theoreticalproteins.fasta or Example_theoreticalproteins.xlsx (available here)
3. An .xlsx of peptide retention time calibration (optional)
    * **Required column names:** Protein, Sequence, PepRtimeStart, PepRtimePeak, PepRtimeEnd
    * **Example:** Example_RTcalibration.xlsx (available here)
#### The outputs are: 
1. A .pkl file of LC/MS experiment
2. An .xlsx file of NMF result
    * **Worksheets:** 1.XIC, 2.Peakresult
3. Three .npz files of V, W, and H matrices

## Installation
```git clone https://github.com/pasrawin/ProteomicNMF.git```
## Running
1. Prepare a directory containing your input files
2. Define your input file names and parameters in ```mNMF00_handler.py``` 
    * The mass spectrograms from different instruments and proteomic experiments provide different features. In order to obtain an optimal result from protNMF, we strongly suggest that the parameters should be carefully set for each observed mass spectrogram. 
    * The *m/z* range, retention time range, smoothing window and shift, bins, *in silico* digestion criteria, number of best peaks reported etc. can be modified easily by replacing default values here.
3. Run ```python mNMF00_handler.py```
    * Yor command prompt will show protNMF process from 1 to 10

## Benchmark Datasets
For benchmarking in the paper, there are 2 standard protein mixtures
* A separation of a thousand ions: 1,057 peptide ions are generated from 476 unique peptides of 4 proteins
* A separation of nine thousand ions: 9,319 peptide ions are generated from 4,045 unique peptides of 48 proteins

The MS raw data were deposited at the ProteomeXchange Consortium via jPOST partner repository with identifier JPST000663 (currently, available for reviewers only).

## Support
If you have any questions about protNMF, please contact Pasrawin Taechawattananant (pasrawin@gmail.com), Kazuyoshi Yoshii (yoshii@kuis.kyoto-u.ac.jp), or Yasushi Ishihama (yishiham@pharm.kyoto-u.ac.jp)

