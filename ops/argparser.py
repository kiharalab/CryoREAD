#
# Copyright (C) 2020 Xiao Wang
# Email:xiaowang20140001@gmail.com wang3702@purdue.edu
#

import argparse

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-F',type=str, help='Input map file path. (str)')
    parser.add_argument('-M', type=str,  default="best_model",help='Pre-trained model path.  (str) Default value: "best_model"')
    parser.add_argument("-P",type=str,help="Optional fasta sequence file path. (str) ")
    parser.add_argument("--output",type=str,help="Output directory")
    parser.add_argument('--mode',type=int,required=True,help='Control Mode for program: 0: cryo_READ structure modeling. Required parameter. (Integer), Default value: 0')
    parser.add_argument("--contour",type=float,default=0,help="Contour level for input map, suggested 0.5*[author_contour]. (Float), Default value: 0.0")
    parser.add_argument("--stride",type=int,default=16,help="Stride for scanning of deep learning model. (Integer), Default value: 16.")
    parser.add_argument("--box_size",type=int,default=64,help="Input box size for deep learning model. (Integer), Default value: 64")
    parser.add_argument("--gpu",type=str,default=None,help="Specify the gpu we will use. (str), Default value: None.")
    parser.add_argument('--batch_size', type=int, default=8, help='Batch size for inference of network. (Integer), Default value: 8.')
    parser.add_argument('-f', type=float, default=0.05, help="Filter for representative points, for LDPs, "
                                                             "removing points' normalized density<=-f (Float), Default value: 0.05")  #
    parser.add_argument('-m', type=float, default=2.0, help="After meanshifting merge points distance<[float]. (Float), Default value: 2.0. ")  # merge distance
    parser.add_argument('-g', type=float, default=3.0, help="Bandwidth of the Gaussian filter, (Float), Default value: 3.0.")  # gaussian filter
    parser.add_argument('-k', type=float, default=0.5, help="Always keep edges where d<k parameter. (Float), Default value: 0.5 ")
    parser.add_argument('-R', type=float, default=10.0, help="Maximum length of local edges. (Float), Default value: 10.0. ")
    parser.add_argument("--rule_soft",type=int,default=0,help="Use strict/soft rules to assemble collected fragments in DP step. (Integer), Default value: 0 (strict rules)")
    parser.add_argument("--frag_size",type=int,default=20,help="Fragment size for sequence split.(Integer), Default value: 20")
    parser.add_argument("--frag_stride",type=int,default=2, help="Frag stride step. (Integer), Default value: 2")
    parser.add_argument("--top_select",type=int,default=20,help="Select top fragment candidate here. (Integer), Default value: 20")
    parser.add_argument("--resolution",type=float,default=2.5,help="resolution of maps, used for final structure refinement. (Float), Default value: 2.5")
    parser.add_argument("--num_workers",type=int,default=4,help="number of workers to fetch data for GPU inference. (Integer), Default value: 4")
    parser.add_argument('--prediction_only',  action='store_true', help="Optional input. Only run the deep learning prediction step. (True/False) Default value: False")
    parser.add_argument('--no_seqinfo',  action='store_true', help="Optional input. Build structures when no sequence information is available. (True/False) Default value: False")
    parser.add_argument("--refine",action="store_true",help="Optional Input. Use Phenix to do the last step refinement or not (Suggested to set as True).")
    parser.add_argument("--colab",action="store_true",help="Optional Input. Used to specify the phenix path for colab settings")
    args = parser.parse_args()
    params = vars(args)
    return params
