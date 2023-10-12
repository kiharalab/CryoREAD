# CryoREAD

<a href="https://github.com/marktext/marktext/releases/latest">
   <img src="https://img.shields.io/badge/CryoREAD-v1.0.0-green">
   <img src="https://img.shields.io/badge/platform-Linux%20%7C%20Mac%20-green">
   <img src="https://img.shields.io/badge/Language-python3-green">
   <img src="https://img.shields.io/badge/dependencies-tested-green">
   <img src="https://img.shields.io/badge/licence-GNU-green">
</a>  

Cryo_READ is a computational tool using deep learning to automatically build full DNA/RNA atomic structure from cryo-EM map.  

Copyright (C) 2022 Xiao Wang, Genki Terashi, Daisuke Kihara, and Purdue University. 

License: GPL v3. (If you are interested in a different license, for example, for commercial use, please contact us.) 

Contact: Daisuke Kihara (dkihara@purdue.edu)

For technical problems or questions, please reach to Xiao Wang (wang3702@purdue.edu).

## Citation:

Xiao Wang, Genki Terashi & Daisuke Kihara. De novo structure modeling for nucleic acids in cryo-EM maps using deep learning. Nature Methods, 2023.
[https://www.nature.com/articles/s41592-023-02032-5](https://www.nature.com/articles/s41592-023-02032-5)
```
@article{wang2023CryoREAD,   
  title={De novo structure modeling for nucleic acids in cryo-EM maps using deep learning},   
  author={Xiao Wang, Genki Terashi, and Daisuke Kihara},    
  journal={Nature Methods},    
  year={2023}    
}   
```


# Online Platform:

## Server(Recommended): https://em.kiharalab.org/algorithm/CryoREAD
<details>
We have three publicly available platforms, which basically offer the same functionality.
Input: cryo-EM map+sequence file (optional). Output: modeled structure. The input and output are the same across all platforms.
</details>

### Google Colab: https://bit.ly/CryoREAD
<details> 
   Step-by-step instructions are available. For free user, colab has 4-hour running time limit and may not work for large structure(>=1000 nucleotides).
</details>

### Local installation with source code at Github
<details>
Full code is available here and it is easier for user to modify to develop their own tools.
<br>It provides two additional supports:
<br>1. Detection Output: This option outputs probability values of detected phosphate, sugar, base, and base types, computed by deep learning, in the map, for users reference.
<br>2. Refinement pipeline: structures from other source can be refined in the specified EM map.
</details>

### Project website: https://kiharalab.org/emsuites
### Detailed pipeline instructions can be found https://kiharalab.org/emsuites/cryoread.php

## Introduction
<details>
   <summary>Cryo_READ is a computational tool using deep learning to automatically build full DNA/RNA atomic structure from cryo-EM map.  </summary>
DNA and RNA play fundamental roles in various cellular processes, where the three-dimensional (3D) structure provides critical information to understand molecular mechanisms of their functions.  Although an increasing number of structures of nucleic acids and their complexes with proteins are determined by cryogenic electron microscopy (cryo-EM), structure modeling for DNA and RNA is still often challenging particularly when the map is determined at sub-atomic resolution. Moreover, computational methods are sparse for nucleic acid structure modeling.

Here, we developed a deep learning-based fully automated de novo DNA/RNA atomic structure modeling method, CryoREAD. CryoREAD identifies phosphate, sugar, and base positions in a cryo-EM map using deep learning, which are traced and modeled into a 3D structure. When tested on cryo-EM maps determined at 2.0 to 5.0 Å resolution, CryoREAD built substantially accurate models than existing methods. We have further applied the method on cryo-EM maps of biomolecular complexes in SARS-CoV-2.
</details>

## Overall Protocol 
<details>
1) Structure Detection by deep neural network CryoREAD networks;   
2) Tracing backbone according to detections;   
3) Fragment-based nucleotide assignment;  
4) Full atomic structure modeling.   


<p align="center">
  <img src="https://user-images.githubusercontent.com/50850224/199084130-34b35a89-3c0c-4647-b693-82fbcc10c820.jpg" alt="CryoREAD framework" width="70%">
</p>
</details>

## Installation
<details>

### System Requirements
CPU: >=8 cores <br>
Memory (RAM): >=50Gb. For maps with more than 3,000 nucleotides, memory space should be higher than 200GB if the sequence is provided. <br>
GPU: any GPU supports CUDA with at least 12GB memory. <br>
GPU is required for CryoREAD and no CPU version is available for CryoREAD since it is too slow.

## Pre-required software
### Required 
Python 3 : https://www.python.org/downloads/     
Phenix: https://phenix-online.org/documentation/install-setup-run.html   
Coot: https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/
### Optional
Pymol (for map visualization): https://pymol.org/2/    
Chimera (for map visualization): https://www.cgl.ucsf.edu/chimera/download.html  

## Installation  
### 1. [`Install git`](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) 
### 2. Clone the repository in your computer 
```
git clone  https://github.com/kiharalab/CryoREAD.git && cd CryoREAD
```

### 3. Build dependencies.   
You have two options to install dependency on your computer:
#### 3.2 Install with anaconda (Recommended)
##### 3.2.1 [`install anaconda`](https://www.anaconda.com/download). 
##### 3.2.2 Install dependency in command line
Make sure you are in the CryoREAD directory and then run 
```
conda env create -f environment.yml
```
Each time when you want to run this software, simply activate the environment by
```
conda activate CryoREAD
conda deactivate(If you want to exit) 
```


#### 3.2 Install with pip and python (Not Suggested).
##### 3.2.1[`install pip`](https://pip.pypa.io/en/stable/installing/).
##### 3.2.2  Install dependency in command line.
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
pip3 install mrcfile
pip3 install torch==1.6.0
```



#### 4 Verify the pre-installed software
To verify phenix is correctly installed for final refinement step, please run
```
phenix.real_space_refine -h
```
To veryify coot is correctly installed for final refinement step, please run
```commandline
coot
```
If it can print out the help information of this function, then the refinemnt step of our program can be supported.
**If not, please always remove --refine command line in all the commands, then CryoREAD should output structure without refinement.**

</details>





# Usage

<details>
<summary>Command Parameters</summary>
   
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
                        Batch size for inference of network. (Integer), Default value: 4.
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
  --thread THREAD
                        Use multiple threads for fragment-based sequence assignment,default:1 (multi-threading is disabled)
```

</details>




<details>
   <summary> DNA/RNA structure modeling without FASTA sequence</summary>
   
### Build atomic structure without sequence information
```
python3 main.py --mode=0 -F=[Map_Path] -M=[Model_Path] --contour=[half_contour_level] --gpu=[GPU_ID] --batch_size=[batch_size] --resolution=[Map_Resolution] --no_seqinfo --refine
```
[Map_Path] is the path of the experimental cryo-EM map, 
<br>[Model_Path] is the path of our pre-trained deep learning model, 
<br>[half_contour_level] is 0.5* contour_level (suggested by author) to remove outside regions to save processing time, 
<br>[GPU_ID] specifies the gpu used for inference,
<br> [batch_size] is the number of examples per batch in the inference (we used 8 with a 24GB GPU), 
<br>[Map_Resolution] is the resolution of the deposited maps.

"--refine" should be removed if you can not successfully install Phenix/coot correctly, which may result in nucleotides that do not satisfy some geometry and chemical constraints.

The automatically build atomic structure is saved in [Predict_Result/(map-name)/Output/Refine_cycle[k].pdb] in pdb format, here default k is 3. However, it may fail if your dependencies are not properly installed, then you may only find Refine_cycle1.pdb or Refine_cycle2.pdb.

#### Example Command:
```
python3 main.py --mode=0 -F=example/21051.mrc -M=best_model --contour=0.3 --gpu=0 --batch_size=4 --resolution=3.7 --no_seqinfo --refine
```

</details>

<details>

<summary> DNA/RNA structure modeling with FASTA sequence</summary>

### Build atomic structure with sequence information
```
python3 main.py --mode=0 -F=[Map_Path] -M=[Model_Path] -P=[Fasta_Path] --contour=[half_contour_level] --gpu=[GPU_ID] --batch_size=[batch_size] --rule_soft=[assignment_rule] --resolution=[Map_Resolution] --refine --thread=[num_threads]
```
[Map_Path] is the path of the experimental cryo-EM map, 
<br>[Model_Path] is the path of our pre-trained deep learning model, 
<br>[Fasta_Path] is the path of the input fasta file about sequence information, 
<br>[half_contour_level] is 0.5* contour_level (suggested by author) to remove outside regions to save processing time, 
<br>[GPU_ID] specifies the gpu used for inference, 
<br>[batch_size] is the number of examples per batch in the inference (we used 8 with a 24GB GPU), 
<br>[rule_soft] specifies the assignment rule, default is 0 to use the strict assignment assembling rule, 
<br>[Map_Resolution] is the resolution of the deposited maps.
<br>[num_thread] specifies the number of CPUs used for fragment-based sequence assignment.

"--refine" should be removed if you can not successfully install Phenix/coot correctly,which may result in nucleotides that do not satisfy some geometry and chemical constraints.



#### Example Command:
```
python3 main.py --mode=0 -F=example/21051.mrc -M=best_model -P=example/21051.fasta --contour=0.3 --gpu=0 --batch_size=4 --rule_soft=0 --resolution=3.7  --refine --thread 4 
```
The automatically build atomic structure is saved in [Predict_Result/(map-name)/Output/Refine_cycle[k].pdb] in pdb format, here default k is 3. However, it may fail if your dependencies are not properly installed, then you may only find Refine_cycle1.pdb or Refine_cycle2.pdb. Modeled structures without considering sequence information are also saved as [Predict_Result/(map-name)/Output/CryoREAD_noseq.pdb] (without refinement). Meanwhile, structures only considering the sequence information without connecting gap regions are saved in [Predict_Result/(map-name)/Output/CryoREAD_seqonly.pdb] (without refinement) for reference.
<br>Please adjust --thread based on your available CPU numbers (more is better).

</details>

<details>

<summary>CryoREAD detection (if you only want to check detection by deep learning)</summary>

### Structure Information Predictions by CryoREAD.
```
python3 main.py --mode=0 -F=[Map_Path] -M=[Model_Path] --contour=[half_contour_level] --gpu=[GPU_ID] --batch_size=[batch_size] --prediction_only 
```
[Map_Path] is the path of the experimental cryo-EM map, 
<br>[Model_Path] is the path of our pre-trained deep learning model, 
<br>[half_contour_level] is 0.5* contour_level (suggested by author) to remove outside regions to save processing time, 
<br>[GPU_ID] specifies the gpu used for inference, 
<br>[batch_size] is the number of examples per batch in the inference (we used 8 with a 24GB GPU). 

The predicted probability maps are saved in [Predict_Result/(map_name)/2nd_stage_detection] with mrc format. It will include 8 mrc files corresponding to 8 different classes.

#### Example Command:
```
python3 main.py --mode=0 -F=example/21051.mrc -M=best_model --contour=0.3 --gpu=0 --batch_size=4 --prediction_only
```
</details>

<details>

<summary>Structure Refinement Pipeline in CryoREAD (if you only want to refine your structure)</summary>

### Structure refinement
The full refinement pipeline involving Phenix and coot is also available for refinement-only purposes. 
```
python3 main.py --mode=1 -F=[input_structure_pdb] -M=[input_map_path] -P=[output_dir]
```
This refinement pipeline can work for any given structure (not limited to DNA/RNA) and a corresponding map.
<br>[input_structure_pdb] is the path of the input structure in pdb format, 
<br>[input_map_path] corresponds to the input map path. 
<br>[output_dir]: the directory you specify to save the outputs during refinement process. The final output Refine_cycle3.pdb will be generated in this directory.

#### Example Command:
```
python3 main.py --mode=1 -F=example/6v5b_drna.pdb -M=example/21051.mrc -P=refine_test
```
This will refine the input structure according to density and output the refined structure in [refine_test] directory. 

</details>

## Example

<details>

### Input File
Cryo-EM map with mrc format. 
(Optional) Sequence information with fasta format.
Our example input can be found [here](https://github.com/kiharalab/CryoREAD/tree/main/example)

### Output File 
1 *.mrc: a mrc file saved our detected probabilites by our deep learning model.    
2 *.pdb: a PDB file that stores the atomic DNA/RNA structure by our method.
Our example output can be found [here](https://kiharalab.org/emsuites/cryoread/output_21051.tar.gz). All the intermediate results are also kept here. 
</details>
