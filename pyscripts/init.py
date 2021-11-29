import argparse
import pandas as pd
import subprocess
import csv
import itertools

parser = argparse.ArgumentParser(description='name')
parser.add_argument('-v1', dest='unfiltered', type=str, help="name of unfilterd merged vcf file")
parser.add_argument('-v2', dest='filtered', type=str, help="name of filtered merged vcf file")
parser.add_argument('-b1', dest='bam', type=str, help="name of bam file")
parser.add_argument('-b2', dest='bed', type=str, help="name of bed file")
parser.add_argument('-o1', dest='name', type=str, help="name of output")
parser.add_argument('-e1', dest='candidates', type=str, help="name of candidates")
parser.add_argument('-e2', dest='voi', type=str, help="name of variants of interest")
parser.add_argument('-f1', dest='fasta', type=str, help="name of fasta file")
args = parser.parse_args()