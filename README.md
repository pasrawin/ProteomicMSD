# Non-negative matrix decomposition for proteomic mass spectrogram
Peak identification and quantification in proteomic mass spectrogram using modified non-negative matrix decomposition

## Overview
Here, we developed a novel approach based on non-negative matrix decomposition, in which the isotopic distribution profiles of peptide ions, the learned noise subspace, and the protein-peptide hierarchy for group sparsity constraint are incorporated into the original matrix decomposition to identify and quantify peaks in proteomic mass spectrogram.

## Requirements
Our Proteomic Mass Spectrogram Decomposition (protMSD) was written in Python 2.7 and tested on Window systems.
#### Dependencies
* division from future, numpy, pandas, scipy, math, collections, itertools, sklearn, gc, and pyteomics

#### The inputs are: 
1. An .mzML or a .txt file of LC/MS experiment converted by ProteoWizard MSconvert, or a .pkl file of protMSD output (required)
    * **Example:** Example_mNMF_standard4proteins.txt
    * **Example:** Example_mNMF_standard4proteins.pkl 
       * Currently, available for reviewers only. Please use the access key in Supplementary Information.
       * After running protMSD, .pkl will be obtained. You may reuse this in following protMSD run to save the computational time.
2. A .FASTA or an .xlsx of theoretical proteins and their sequences (required)
    * **Required column names:** Protein, Sequence
    * **Example:** Example_theoreticalproteins.fasta or Example_theoreticalproteins.xlsx (available here)
3. An .xlsx of peptide retention time calibration (optional)
    * **Required column names:** Protein, Sequence, PepRtimeStart, PepRtimePeak, PepRtimeEnd
    * **Example:** Example_RTcalibration.xlsx (available here)
#### The outputs are: 
1. A .pkl file of LC/MS experiment
2. An .xlsx file of protMSD result
    * **Worksheets:** 1.XIC, 2.Peakresult
3. Three .npz files of V, W, and H matrices

## Installation
```git clone https://github.com/pasrawin/ProteomicMSD.git```
## Running
1. Prepare a directory containing your input files
2. Define your input file names and parameters in ```mNMF00_handler.py``` 
    * The mass spectrograms from different instruments and proteomic experiments provide different features. In order to obtain an optimal result from protMSD, we strongly suggest that the parameters should be carefully set for each observed mass spectrogram. 
    * The *m/z* range, retention time range, smoothing window and shift, bins, *in silico* digestion criteria, number of best peaks reported etc. can be modified easily by replacing default values here.
3. Run ```python mNMF00_handler.py```
    * Yor command prompt will show protMSD process from 1 to 10

## Benchmark Datasets
For benchmarking in the paper, there are 2 standard protein mixtures
* A separation of a thousand ions: 1,057 peptide ions are generated from 476 unique peptides of 4 proteins
* A separation of nine thousand ions: 9,319 peptide ions are generated from 4,045 unique peptides of 48 proteins

The MS raw data were deposited at the ProteomeXchange Consortium via jPOST partner repository with identifier [JPST000663](https://repository.jpostdb.org/preview/4730115445e7091ad1182a) (Currently, available for reviewers only. Please use the access key in Supplementary Information).

## Support
If you have any questions about protMSD, please contact Pasrawin Taechawattananant (pasrawin@gmail.com), Kazuyoshi Yoshii (yoshii@kuis.kyoto-u.ac.jp), or Yasushi Ishihama (yishiham@pharm.kyoto-u.ac.jp)

