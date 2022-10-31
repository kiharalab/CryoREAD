# cryo_READ

<a href="https://github.com/marktext/marktext/releases/latest">
   <img src="https://img.shields.io/badge/cryo_READ-v1.0.0-green">
   <img src="https://img.shields.io/badge/platform-Linux%20%7C%20Mac%20-green">
   <img src="https://img.shields.io/badge/Language-python3-green">
   <img src="https://img.shields.io/badge/dependencies-tested-green">
   <img src="https://img.shields.io/badge/licence-GNU-green">
</a>  

Cryo_READ is a computational tool using deep learning to automatically build full DNA/RNA atomic structure from cryo-EM map.  

Copyright (C) 2022 Xiao Wang, Genki Terashi, Daisuke Kihara, and Purdue University. 

License: GPL v3. (If you are interested in a different license, for example, for commercial use, please contact us.) 

Contact: Daisuke Kihara (dkihara@purdue.edu), Xiao Wang (wang3702@purdue.edu)

## Citation:

Xiao Wang, Genki Terashi & Daisuke Kihara. Cryo-READ: DNA/RNA dE novo Atomic structure moDeling in cryo-EM maps with deep learning. bioArxiv, 2022.
```
@article{wang2021emap2secplus,   
  title={Cryo-READ: DNA/RNA dE Novo Atomic Structure MoDeling in cryo-EM Maps with Deep Learning},   
  author={Xiao Wang,  Genki Terashi, and Daisuke Kihara},    
  journal={bioArxiv},    
  year={2022}    
}   
```

## Project website: 
## Online Platform:  (Run easily and freely)

## Introduction

## Overall Protocol 
1) Structure Detection by deep neural network Cryo-READ networks;   
2) Tracing backbone according to detections;   
3) Fragment-based nucleotide assignment;  
4) Full atomic structure modeling.   

<p align="center">
  <img src="figures/framework.jpg" alt="cryo-READ framework" width="90%">
</p> 


## Pre-required software
### Required 
Python 3 : https://www.python.org/downloads/     
Phenix: https://phenix-online.org/documentation/install-setup-run.html   
### Optional
Pymol (for map visualization): https://pymol.org/2/    
Chimera (for map visualization): https://www.cgl.ucsf.edu/chimera/download.html  

## Installation  
### 1. [`Install git`](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) 
### 2. Clone the repository in your computer 
```
git clone git@github.com:kiharalab/cryo_READ.git && cd cryo_READ
```

### 3. Build dependencies.   
You have two options to install dependency on your computer:
#### 3.1 Install with pip and python(Ver 3.6.9).
##### 3.1.1[`install pip`](https://pip.pypa.io/en/stable/installing/).
##### 3.1.2  Install dependency in command line.
```
pip3 install -r requirements.txt --user
```
If you encounter any errors, you can install each library one by one:
```
pip3 install biopython
pip3 install numpy
pip3 install numba
pip3 install scipy
pip3 install ortools
pip3 install sklearn
pip3 install mrcfile
pip3 install torch==1.6.0
```

#### 3.2 Install with anaconda
##### 3.2.1 [`install conda`](https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html). 
##### 3.2.2 Install dependency in command line
```
conda create -n cryo_READ python=3.8.10
conda activate cryo_READ
pip install -r requirements.txt 
```
Each time when you want to run this software, simply activate the environment by
```
conda activate cryo_READ
conda deactivate(If you want to exit) 
```

#### 4 Verify the pre-installed software
To verify phenix is correctly installed for final refinement step, please run
```
phenix.real_space_refine -h
```
If it can print out the help information of this function, then the refinemnt step of our program can be supported.
**If not, please always remove --refine command line in all the command, then cryoREAD should output structure without refinement.**


## Usage
```
usage: main.py [-h] [-F F] [-M M] [-P P] --mode MODE [--contour CONTOUR] [--stride STRIDE] [--box_size BOX_SIZE] [--gpu GPU] [--batch_size BATCH_SIZE] [-f F] [-m M]
               [-g G] [-k K] [-R R] [--rule_soft RULE_SOFT] [--frag_size FRAG_SIZE] [--frag_stride FRAG_STRIDE] [--top_select TOP_SELECT] [--resolution RESOLUTION]
               [--num_workers NUM_WORKERS] [--prediction_only PREDICTION_ONLY] [--no_seqinfo NO_SEQINFO]

optional arguments:
  -h, --help            show this help message and exit
  -F F                  Input map file path. (str)
  -M M                  Pre-trained model path. (str) Default value: "best_model"
  -P P                  Optional fasta sequence file path. (str)
  --mode MODE           Control Mode for program: 0: cryo_READ structure modeling. Required parameter. (Integer), Default value: 0
  --contour CONTOUR     Contour level for input map, suggested 0.5*[author_contour]. (Float), Default value: 0.0
  --stride STRIDE       Stride for scanning of deep learning model. (Integer), Default value: 16.
  --box_size BOX_SIZE   Input box size for deep learning model. (Integer), Default value: 64
  --gpu GPU             Specify the gpu we will use. (str), Default value: None.
  --batch_size BATCH_SIZE
                        Batch size for inference of network. (Integer), Default value: 8.
  -f F                  Filter for representative points, for LDPs, removing points' normalized density<=-f (Float), Default value: 0.05
  -m M                  After meanshifting merge points distance<[float]. (Float), Default value: 2.0.
  -g G                  Bandwidth of the Gaussian filter, (Float), Default value: 3.0.
  -k K                  Always keep edges where d<k parameter. (Float), Default value: 0.5
  -R R                  Maximum length of local edges. (Float), Default value: 10.0.
  --rule_soft RULE_SOFT
                        Use strict/soft rules to assemble collected fragments in DP step. (Integer), Default value: 0 (strict rules)
  --frag_size FRAG_SIZE
                        Fragment size for sequence split.(Integer), Default value: 20
  --frag_stride FRAG_STRIDE
                        Frag stride step. (Integer), Default value: 2
  --top_select TOP_SELECT
                        Select top fragment candidate here. (Integer), Default value: 20
  --resolution RESOLUTION
                        resolution of maps, used for final structure refinement. (Float), Default value: 2.5
  --num_workers NUM_WORKERS
                        number of workers to fetch data for GPU inference. (Integer), Default value: 4
  --prediction_only PREDICTION_ONLY
                        Optional input. Only run the deep learning prediction step. (True/False) Default value: False
  --no_seqinfo NO_SEQINFO
                        Optional input. Build structures when no sequence information is available. (True/False) Default value: False
```

### 1. Only Make Structure Information Predictions by cryo-READ.
```
python3 main.py --mode=0 -F=[Map_Path] -M=[Model_Path] --contour=[half_contour_level] --gpu=[GPU_ID] --batch_size=[batch_size] --prediction_only 
```
[Map_Path] is the path of the experimental cryo-EM map, [Model_Path] is the path of our pre-trained deep learning model, [half_contour_level] is 0.5* contour_level (suggested by author) to remove outside regions to save processing time, [GPU_ID] specifies the gpu used for inference, [batch_size] is the number of examples per batch in the inference (we used 8 with a 24GB GPU). 



The predicted probability maps are saved in [Predict_Result/(map_name)/2nd_stage_detection] with mrc format. It will include 8 mrc files corresponding to 8 different classes.

Example Command:
```
python3 main.py --mode=0 -F=example/21051.mrc -M=best_model --contour=0.3 --gpu=0 --batch_size=8 --prediction_only
```

### 2. Build atomic structure without sequence information
```
python3 main.py --mode=0 -F=[Map_Path] -M=[Model_Path] --contour=[half_contour_level] --gpu=[GPU_ID] --batch_size=[batch_size] --resolution=[Map_Resolution] --no_seqinfo --refine
```
[Map_Path] is the path of the experimental cryo-EM map, [Model_Path] is the path of our pre-trained deep learning model,  [half_contour_level] is 0.5* contour_level (suggested by author) to remove outside regions to save processing time, [GPU_ID] specifies the gpu used for inference, [batch_size] is the number of examples per batch in the inference (we used 8 with a 24GB GPU), [Map_Resolution] is the resolution of the deposited maps.

"--refine" should be removed if you can not successfully install Phenix correctly.

The automatically build atomic structure is saved in [Predict_Result/(map-name)/graph_atomic_modeling/Output_Structure_noseq/Final_Refined*.pdb] in pdb format.

Example Command:
```
python3 main.py --mode=0 -F=example/21051.mrc -M=best_model --contour=0.3 --gpu=0 --batch_size=8 --resolution=3.7 --no_seqinfo --refine
```


### 3. Build atomic structure with sequence information
```
python3 main.py --mode=0 -F=[Map_Path] -M=[Model_Path] -P=[Fasta_Path] --contour=[half_contour_level] --gpu=[GPU_ID] --batch_size=[batch_size] --rule_soft=[assignment_rule] --resolution=[Map_Resolution] --refine
```
[Map_Path] is the path of the experimental cryo-EM map, [Model_Path] is the path of our pre-trained deep learning model, [Fasta_Path] is the path of the input fasta file about sequence information, [half_contour_level] is 0.5* contour_level (suggested by author) to remove outside regions to save processing time, [GPU_ID] specifies the gpu used for inference, [batch_size] is the number of examples per batch in the inference (we used 8 with a 24GB GPU), [rule_soft] specifies the assignment rule, default is 0 to use the strict assignment assembling rule, [Map_Resolution] is the resolution of the deposited maps.

"--refine" should be removed if you can not successfully install Phenix correctly.

Example Command:
```
python3 main.py --mode=0 -F=example/21051.mrc -M=best_model -P=example/21051.fasta --contour=0.3 --gpu=0 --batch_size=8 --rule_soft=0 --resolution=3.7  --refine
```
The automatically build atomic structure is saved in [Predict_Result/(map-name)/graph_atomic_modeling/Output_Structure/Final_Refined*.pdb] in pdb format.



## Example
### Input File
Cryo-EM map with mrc format. 
(Optional) Sequence information with fasta format.

### Output File 
1 *.mrc: a mrc file saved our detected probabilites by our deep learning model.    
2 *.pdb: a PDB file that stores the atomic DNA/RNA structure by our method.

