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

Contact: Daisuke Kihara (dkihara@purdue.edu)

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
Python 3 : https://www.python.org/downloads/    
Phenix: https://phenix-online.org/documentation/install-setup-run.html
Pymol(for visualization): https://pymol.org/2/  

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



