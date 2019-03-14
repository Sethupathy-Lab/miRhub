#!/usr/bin/python3

import csv
import glob
import os
import sys
import xlwt

def write_miRhub_results(job_di):
    wb = xlwt.Workbook()

    ws = wb.add_sheet('miRhub_run_info')
    for i, h in enumerate(['sheet', 'cons', 'iters']):
        ws.write(0, i, h) 
    for rowx, filename in enumerate(sorted(job_di)):
        job = job_di[filename]
        sheet = '{}_{}'.format(job.DEG, job.cons)
        for colx, val in enumerate([sheet, job.cons, job.iter]):
            ws.write(rowx + 1, colx, val)



    for fi in sorted(job_di):
        job = job_di[fi]
        filename = '{}.out2.txt'.format(job.outname)
        sheet_name = '{}_{}'.format(os.path.splitext(os.path.basename(job.DEG))[0][:28], job.cons)
        ws = wb.add_sheet(sheet_name)
        with open(filename) as f:
            for rowx, row in enumerate(f):
                for colx, value in enumerate(row.split('\t')):
                    ws.write(rowx, colx, value)
    wb.save("{}.xls".format(job.project))

