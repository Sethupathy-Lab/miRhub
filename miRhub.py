#!/home/pr46_0001/shared/bin/python3

import glob
import os
import argparse
import sys
from miRhub.miRhub_job import mirHubJob
from miRhub.xls_write import write_miRhub_results

usage = '''
miRhub determines enrichment of predicted miR binding sites within
a list of genes. 

For more information on miRhub, refer to:

Baran-Gale J, Fannin EE, Kurtz CL, Sethupathy P (2013) Beta Cell 
5' Shifted isomiRs Are Candidate Regulatory Hubs in Type 2 Diabetes. 
PLOS ONE 8(9): e73240. 

Supplementary information, section 'Candidate miRNA regulatory hub 
identification pipeline'. NOTE: Protein regulatory hub is depreciated.
'''


def run_miRhub(arg):
    '''
    Build perl miRhub job
    '''
    sheets = []
    report = {}
    for geneList in arg.gene_lists:
        print(geneList)
        for c in arg.cons: 
            job = mirHubJob()
            job.add_DEG(geneList)
            job.add_proj(arg.project)
            job.add_spec_files(arg.species)
            job.add_cons(c)
            job.add_iter(arg.iter)
            job.build_output_name()
            job.build_mirhub_job()
            job.build_summary_job()
            sheets.append(job.outname)
            report[job.outname] = job
    print(sheets)
    write_miRhub_results(report)
    

def main(arg):
    run_miRhub(arg)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
             description=usage,
             formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-c', '--cons',
                        action='store',
                        default='012',
                        help='Required miR conservation across included species (default = 012)')
    parser.add_argument('-p', '--project',
                        action='store',
                        default='miRhubProject',
                        help='Name for the project to be affixed to output files (default = miRhubProject)')
    parser.add_argument('-i', '--iter',
                        action='store',
                        default=1000,
                        type=int,
                        help='Number of iterations to run (default = 1000)')
    parser.add_argument('-sp', '--species',
                        action='store',
                        dest='species',
                        default='mouse',
                        help='Indicate which species miRhub is being run for; should be mouse, rat or human (default = mouse)')
    parser.add_argument('gene_lists',
                        nargs='+',
                        help='Lists of differentially expressed genes from DESeq')
    main(parser.parse_args())
