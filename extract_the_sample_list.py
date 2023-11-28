#!/bin/python3
#Program:
#   Extract the sample list
'''
    @Autor: Jiang XH
    @Date:2023/11/10
    @First release: 2023/11/11
'''

import os
import glob
import re
import json
import argparse


class Sample_list():

    def __init__(self,path):

        self.path = path # the abs path of sample sequencing data (expected clean data)

    def out_sample_json(self):
    
        Reads1 = glob.glob(self.path+'/*/*R1.fq.gz')
        Reads2 = glob.glob(self.path+'/*/*R2.fq.gz')
        
        All_files = {} # all files container
        Samples = [re.split(r'[./]',name)[-4] for name in Reads1] # get all samplenames default_format:{xx}/{xx}.R1.fq.gz
        
        for sample in Samples:
            R1 = lambda Reads1: sample in Reads1 # abs path of Reads1
            R2 = lambda Reads2: sample in Reads2 # abs path of Reads2
            All_files[sample] = {}
            All_files[sample]['R1'] = sorted(filter(R1,Reads1)) 
            All_files[sample]['R2'] = sorted(filter(R2,Reads2))
        
        sample_json = json.dumps(All_files,indent=4,sort_keys=True) # generate the jsom format of sample 
        open('All_sample.json','w').writelines(sample_json)
        
if __name__ == '__main__':

    '''
        The main
                '''

    parser = argparse.ArgumentParser(description='Extract the sample list from the raw data path')
    parser.add_argument('--path','-p',required=True,help='The clean sequencing data path, abs path is recommended')
    args = parser.parse_args()

    extract_sf = Sample_list(args.path)
    
    if not os.path.exists(args.path):
        print('The fastq folder is not exist, Please checking your data path !')
        exit()
    
    extract_sf.out_sample_json() # sample json

