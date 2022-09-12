#
# Copyright (C) 2020 Xiao Wang
# Email:xiaowang20140001@gmail.com wang3702@purdue.edu
#

from collections import defaultdict
import parser
import argparse

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-F',type=str, help='Input map file path')
    parser.add_argument('-M', type=str,  help='pre-trained model path')
    parser.add_argument("-P",type=str,help="optional fasta sequence file path")
    parser.add_argument('--mode',type=int,required=True,help='control mode')
    parser.add_argument("--contour",type=float,default=0,help="contour level for input map, suggested 0.5*[author_contour]")
    parser.add_argument("--stride",type=int,default=16,help="stride for scanning of deep learning model")
    parser.add_argument("--box_size",type=int,default=64,help="input box size for deep learning model")
    parser.add_argument("--gpu",type=str,default=None,help="specify the gpu we will use")
    parser.add_argument('--batch_size', type=int, default=8, help='batch size for inference')
    parser.add_argument('-f', type=float, default=0.05, help="Filter for representative points, for LDPs, "
                                                             "removing points' normalized density<=-f")  #
    parser.add_argument('-m', type=float, default=1.0, help="after meanshifting merge points distance<[float]")  # merge distance
    parser.add_argument('-g', type=float, default=3.0, help="Bandwidth of the Gaussian filter")  # gaussian filter
    parser.add_argument('-k', type=float, default=0.5, help="keep edges where d<k parameter")
    parser.add_argument('-R', type=float, default=10.0, help="radius of local MST")
    parser.add_argument("--rule_soft",type=int,default=0,help="Use strict/soft rules to assemble collected fragments in DP step")
    parser.add_argument("--resolution",type=float,default=2.5,help="resolution of maps, used for final structure refinement")
    parser.add_argument("--num_workers",type=int,default=4,help="number of workers to fetch data")
    args = parser.parse_args()
    params = vars(args)
    return params
