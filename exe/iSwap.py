#!/usr/bin/python3
#filter of data
from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import Counter
import os, sys, argparse, time, multiprocessing, math, itertools, datetime, editdistance
from multiprocessing import Queue,Process
from random import randint
import seaborn as sns
from itertools import combinations


### DEFAULT VALUES AS MODULE ##############################################################################
# Print on screen
verbose = True
print_time = True  # eval if verbose is True
# About output
output_matrix_sep = '\t'
output_matrixes_encoding = 'utf-8'
output_na_rep = ''  # how to WRITE missing data
# About input
input_matrixes_sep = '\t'
input_matrixes_encoding = 'utf-8'
input_matrixes_prefix = None  # None or a list paired with input_paths

use_attributes = None  # None or list of int 0-based (do relabelling)
old_attibute_sep = '_'  # Exploited in relabelling (use_attributes is not None) but also by add_prefix!
new_attribute_sep = '_'  # It makes sense only if use_attributes is not None (do relabelling)
g = '\033[92m'
e = '\033[0m'
w = '\033[93m'
r = '\033[91m'
###########################################################################################################

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

header = """
+------------------------------------------------------+
         *** {g}ISwap{e} {g}I{e}ntegration {g}S{e}ite {g}Swap{e}  ***
+------------------------------------------------------+
 Author:    Adriano De Marino
 Date:      May 2018
 Contact:   demarino.adriano@hsr.it,calabria.andrea@hsr.it,
 Revision:  1.2 Origin
 Revision:  1.3 Modified October 2018 by Adriano De Marino
 Revision:  1.4 Modified March 2019 by Adriano De Marino
 Revision:  1.5 Modified June 2019 by Adriano De Marino
+------------------------------------------------------+
""".format(g=g,e=e)

description = """ 

Description:
This program clean the dataset from Switch of the LTR & LC and partial Wabbly of Shearsite

"""
print(header)

usage_example = """
 {g}-f{e}   'path to the *_contamination.tsv' to process
 {g}-j{e}   number of jobs exp. ( 5 means that the total work will be splitted in 5 part )   
 {g}-o{e}   name of the output file, only name no extension [Default 'untitled']
 {g}-h{e}   print help

Examples of usage: {w}nohup python ISwap.py -f /absolute_path/filename_contamination.csv.gz -j 4 -o /absolute_path/project_name {e}


""".format(g=g,e=e,w=w)

parser = argparse.ArgumentParser(usage = usage_example, epilog = "[ HSR-TIGET - Vector Integration Core - Bioinformatics ] \n", description = description)
parser.add_argument('-f', dest="contamina_file", help="*_contamination.tsv. No default option.", action="store", required=True)
parser.add_argument("-j", dest='jobs', default=2, action="store", required=False, help="set the cpu to use for the Parallelism, default 2")
parser.add_argument('-o', dest="out", help="output name", action="store", required=True, default='untitled')

args = parser.parse_args()

contaminationFile = args.contamina_file
jobs = int(args.jobs)
output = args.out

start_time = time.time()


import re
def humanSorted(l):
    def tryint(s):
        try:
            return int(s)
        except:
            return s
    def alphanum_key(s):
        return [ tryint(c) for c in re.split('([0-9]+)', s) ]
    def sort_nicely(l):
        return sorted(l, key=alphanum_key)
    return sort_nicely(l)

def verbosePrint(x, verbose=verbose, print_time=print_time):
    '''
    Purpose: print( immediately if verbose is True, with
             date and time if print(_time is True too.
    
    IN: 
    x - anything supporting __str__ method
    
    OUT:
    x on stdout, immediately flushed
    
    NOTE:
    uncommon/complex patterns of "\n" may be not properly
    handled in case of 'print_time'. 
    '''
    if verbose:
        if print_time:
            from time import localtime, strftime
            y = str(x)
            nnl = y.count("\n", 0, int(len(x)/2+1))
            y = y.replace("\n", "", nnl)
            nl_str = "".join(['\n']*nnl)
            print(nl_str+"[{time}] ".format(time=strftime("%Y-%m-%d %H:%M:%S", localtime())), y)
        else:
            print( x)
        sys.stdout.flush()

def check_input_paths (*paths):
    '''
    Purpose: check input matrix file path(s);
             exit and explain why if something is wrong.
    
    IN: 
    *paths - complete path(s) of matrix file(s)
    
    OUT:
    0
    '''
    verbosePrint("\n[CHECK IN PATH(S)]")
    for p in paths:
        # normalize path
        p = str(p)
        try:
            p = os.path.normpath(p)
        except( Exception, err_message):
            print( "\n[ERROR] input path='{p}' is not formatted as valid path!".format(p=str(p)))
            print( "os.path.normpath returned: ", err_message)
            sys.exit("\n[QUIT]\n")
        p = os.path.expanduser(p)
        p = os.path.expandvars(p)
        # check if path is absolute, read permission and if path points to an existing file
        if not os.path.isabs(p):
            print( "\n[ERROR] input path must be an absolute path! Your input path='{p}'".format(p=str(p)))
            sys.exit("\n[QUIT]\n")
        if not os.access(os.path.dirname(p), os.R_OK):
            print( "\n[ERROR] You have not read permission in folder='{f}'".format(f=str(os.path.dirname(p))))
            sys.exit("\n[QUIT]\n")
        if not os.path.isfile(p):
            print( "\n[ERROR] input path must point to an existing file! Your input path='{p}'".format(p=str(p)))
            sys.exit("\n[QUIT]\n")
    verbosePrint("...OK!")
    verbosePrint("[DONE]")
    return 0

def buildOUTDIR(ground_dir, *subfolders):
    '''
    Purpose: check (and create) an absolute folder path where WRITE FILES
    
    IN: 
    ground_dir - absolute folder path (checked)
    *subfolders - optional, further subfolders to concatenate (checked every time)
    
    OUT:
    OUTDIR - absolute folder path where write, permissions are checked
    '''
    # Final path variable
    OUTDIR = None  # <- os.path.normpath(ground_dir) + join *subfolders
    # Try OUTDIR=normpath(OUTDIR) and formal check
    try:
        OUTDIR = os.path.normpath(ground_dir)
    except(Exception, err_message):
        print( "\n[ERROR] ground_dir='{ground_dir}' is not formatted as valid path!".format(ground_dir=str(ground_dir)))
        print( "os.path.normpath returned: ", err_message)
        sys.exit("\n[QUIT]\n")
    if os.path.isfile(OUTDIR):
        print( "\n[ERROR] ground_dir='{OUTDIR}' must be a folder path, not a file!".format(OUTDIR=str(OUTDIR)))
        sys.exit("\n[QUIT]\n")
    if not os.path.isabs(OUTDIR):
        print( "\n[ERROR] ground_dir='{OUTDIR}' must be an absolute path!".format(OUTDIR=str(OUTDIR)))
        sys.exit("\n[QUIT]\n")
    # Check ground_dir: try to create it if not exists
    if not os.path.exists(OUTDIR):
        try:
            os.makedirs(OUTDIR)
        except (Exception, err_message):
            print( "\n[ERROR] ground_dir='{OUTDIR}' is not a valid path or you don't have permissions to create it!".format(OUTDIR=str(OUTDIR)))
            print( "os.makedirs returned: ", err_message)
            sys.exit("\n[QUIT]\n")
        verbosePrint("> ground_dir created: {OUTDIR}".format(OUTDIR=str(OUTDIR)))
    else:
        verbosePrint("> ground_dir found: {OUTDIR}".format(OUTDIR=str(OUTDIR)))
    # Check ground_dir: write permissions (needed in case you don't give any *subfolders and ground_dir already exist)
    if not subfolders:
        if not os.access(OUTDIR, os.W_OK):
            print( "\n[ERROR] You don't have write permissions in ground_dir='{OUTDIR}'".format(OUTDIR=str(OUTDIR)))
            sys.exit("\n[QUIT]\n")
    # Loop over *subfolders
    for sf in subfolders:
        # Try OUTDIR=join(OUTDIR, sf) and formal check
        try:
            OUTDIR = os.path.normpath(os.path.join(OUTDIR, sf))
        except( Exception, err_message):
            print( "\n[ERROR] Cannot join OUTDIR='{OUTDIR}' and SUBFOLDER='{SUBFOLDER}' as a valid path!".format(OUTDIR=str(OUTDIR), SUBFOLDER=str(sf)))
            print( "os.path.normpath(os.path.join(...)) returned: ", err_message)
            sys.exit("\n[QUIT]\n")
        if os.path.isfile(OUTDIR):
            print( "\n[ERROR] *subfolders args must be folder(s), not file(s)! Input subfolders: {subfolders}.".format(subfolders=str(subfolders)))
            sys.exit("\n[QUIT]\n")
        # Check: try to create OUTDIR if not exists
        if not os.path.exists(OUTDIR):
            try:
                os.makedirs(OUTDIR)
            except (Exception, err_message):
                print( "\n[ERROR] OUTDIR='{OUTDIR}' + SUBFOLDER='{SUBFOLDER}' is not a valid path or you don't have permissions to create it!".format(OUTDIR=str(OUTDIR), SUBFOLDER=str(sf)))
                print( "os.makedirs returned: ", err_message)
                sys.exit("\n[QUIT]\n")
            verbosePrint("> subfolder created: {OUTDIR}".format(OUTDIR=str(OUTDIR)))
        else:
            verbosePrint("> subfolder found: {OUTDIR}".format(OUTDIR=str(OUTDIR)))
    # Check: write permissions (needed in case you specify *subfolders that already exists)
    if not os.access(OUTDIR, os.W_OK):
        print( "\n[ERROR] You don't have write permissions in OUTDIR='{OUTDIR}'".format(OUTDIR=str(OUTDIR)))
        sys.exit("\n[QUIT]\n")
    # Return OUTDIR
    verbosePrint("> OUTDIR: {OUTDIR}".format(OUTDIR=str(OUTDIR)))
    verbosePrint("> write permission: OK.")
    return OUTDIR
    
def build_outpath (out_path):
    '''
    Purpose:
    1) check if out_path is ok, normpath, expanduser, expandvars.
    2) exploit buildOUTDIR to create desired folder's tree till the last directory.
    
    IN: 
    out_path - complete abs path of the file you want to write
    
    OUT:
    out_path - complete abs path of the file you want to write
    (checked, normed, expanded and existing till the last directory)
    '''
    verbosePrint("\n[CHECK AND CREATE OUT PATH]")
    out_path = str(out_path)
    try:
        out_path = os.path.normpath(out_path)
    except (Exception, err_message):
        print( "\n[ERROR] out_path='{out_path}' is not formatted as valid path!".format(out_path=str(out_path)))
        print( "os.path.normpath returned: ", err_message)
        sys.exit("\n[QUIT]\n")
    out_path = os.path.expanduser(out_path)
    out_path = os.path.expandvars(out_path)
    if not os.path.isabs(out_path):
        print( "\n[ERROR] out_path must be an absolute path! Your out_path='{out_path}'".format(out_path=str(out_path)))
        sys.exit("\n[QUIT]\n")
    items = out_path.split(os.sep)
    buildOUTDIR(items[0]+os.sep, *tuple(items[1:-1]))
    verbosePrint("[DONE]")
    return str(out_path)

def groups(inputc):

    def LC(x):
        key = 'FAKE'
        if key in x:
            return x.split('_')[4].split('.')[1]
        else:
            return 'LC'+x.split('_')[2].split('LC')[1]

    def LTR(x):
        key = 'FAKE'
        if key in x:
            return x.split('_')[4].split('.')[0]
        else:
            return x.split('_')[2].replace(".","").split('LC')[0]

    any_df = pd.read_csv(inputc, sep='\t')

    if 'Unnamed: 0' in any_df.columns:
        any_df.drop('Unnamed: 0', axis=1, inplace=True)
    
    any_df['LTR'] = any_df.association_ID.apply(LTR)
    any_df['LC'] = any_df.association_ID.apply(LC)

    any_df.shearsite = any_df['shearsite'].astype(np.int32)
    any_df.seq_count = any_df['seq_count'].astype(np.int32)

    gb = any_df.groupby('genomic_coordinates') 
    gb_list = []

    for _chr_ in gb:
        gb_list.append(_chr_[1])
    return gb_list,any_df

def concatenat(lis):
    return pd.concat(lis,ignore_index=True,sort=False)

# split a list into evenly sized chunks
def chunks(l, n):
    return [l[i:i+int(n)] for i in range(0, len(l), int(n))]

def split_list(alist, wanted_parts=1):
    length = len(alist)
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] 
             for i in range(wanted_parts) ]

def paradox(classSize,span=800):
    '''
    Calculate the prob to observe same shearsite in a IS
    '''
    prob = 1.0
    total = span
    for i in range(classSize):
        prob = prob * (total-i)/total
    return 1-prob
    
def swap_LC(job_id, data_slice, return_dict, events, trash, threshold):

    for gb in data_slice:
        if gb.seq_count.nunique() == 1 and gb.seq_count.unique()[0] == 1 and gb.association_ID.nunique() > 1:
            trash[gb.iloc[0].genomic_coordinates] = gb
            continue
        else:
            #new implementation for CEM samples TO DO....important...because in this way is underestimate
            test = gb[['association_ID','seq_count']].copy().groupby('association_ID').sum()
            inizio = test[test.seq_count > 0].to_dict()['seq_count']
            sos = gb.copy().reset_index(drop=True)
            x = sos.groupby(['shearsite','LTR'], sort=False)
            ss1 = x.sum().dropna(0)
            nome1 = x.apply(lambda subf: subf['association_ID'][subf['seq_count'].idxmax()])
            LC2 = x.apply(lambda subf: subf['LC'][subf['seq_count'].idxmax()])
            randomBC1 = x.apply(lambda subf: subf['randomBC'][subf['seq_count'].idxmax()])
            swapLC = pd.concat([nome1,LC2,randomBC1,ss1],1).reset_index()#.drop('LTR',axis=1)
            swapLC['genomic_coordinates'] = gb.iloc[0].genomic_coordinates
            swapLC.columns = ['shearsite','LTR','association_ID','LC','randomBC', 'seq_count','genomic_coordinates']
            swapLC = swapLC[['association_ID', 'genomic_coordinates','shearsite', 'randomBC', 'seq_count','LTR','LC']]
            swapLC.columns = swapLC.columns.get_level_values(0)
            test2 = swapLC[['association_ID','seq_count']].copy().groupby('association_ID').sum()
            fine = test2[test2.seq_count > 0].to_dict()['seq_count']
            contaminating = set(inizio.keys()) - set(fine.keys())
            reads = np.sum([inizio[x] for x in contaminating])
            source_tmp = pd.Series(inizio)-pd.Series(fine)
            source_T = source_tmp[source_tmp < 0].to_dict()
            try:
                umi_removed = gb.randomBC.nunique()-swapLC.randomBC.nunique()
            except Exception as e:
                umi_removed = 0
            #in questo caso entreranno nell loop soltanto i campioni scartati e dunque che si trovano nel contaminating 
            if all(inizio[x]*threshold <= max(fine.values()) for x in contaminating):
                events[gb.iloc[0].genomic_coordinates] = [reads,umi_removed,{ x : gb[gb.association_ID == x].seq_count.values for x in contaminating},source_T,gb.iloc[0].genomic_coordinates]
                return_dict[gb.iloc[0].genomic_coordinates] = swapLC
            else:
                trash[gb.iloc[0].genomic_coordinates] = swapLC
                continue

def cluster(data, maxgap):
    '''Arrange data into groups where successive elements
       differ by no more than *maxgap*

        >>> cluster([1, 6, 9, 100, 102, 105, 109, 134, 139], maxgap=10)
        [[1, 6, 9], [100, 102, 105, 109], [134, 139]]

        >>> cluster([1, 6, 9, 99, 100, 102, 105, 134, 139, 141], maxgap=10)
        [[1, 6, 9], [99, 100, 102, 105], [134, 139, 141]]
    '''
    data.sort()
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - np.median(groups[-1])) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups

def wabbly(guru,output):
    # output = []
    for i in guru:
        #Se IS ha tutti 1 come seq_count = BUTTA 
        if i.seq_count.nunique() == 1 and i.seq_count.unique()[0] == 1:
            continue #non salvare questa integrazione
        elif paradox(i.shearsite.count()) < 0.5:
            i.association_ID = i.association_ID.astype(object)
            #per ogni campione controllo lo wobbling dello shearsite
            for s1,s1_df in i.groupby('association_ID'):
                #cluster dello wabbling
                gruppos = cluster(s1_df.copy().shearsite.values,7)
                '''
                per ogni gruppo unisco poiche si trovano molto vicini tra di loro ed essendo pochi per essere
                entrati in questo if allora significa che e' altamente improbabile che questa cosa sia successa per caso
                ma e' dovuta ad un errore di wabbling
                '''
                for grb in gruppos:
                    datafram_small = s1_df[s1_df.shearsite.isin(grb)]
                    x = datafram_small.groupby(['association_ID'])
                    ss1 = x.sum().dropna(0).drop('shearsite',axis=1)
                    ssh = x.apply(lambda subf: subf['shearsite'][subf['seq_count'].idxmax()])
                    LC2 = x.apply(lambda subf: subf['LC'][subf['seq_count'].idxmax()])
                    LTR1 = x.apply(lambda subf: subf['LTR'][subf['seq_count'].idxmax()])
                    randomBC1 = x.apply(lambda subf: subf['randomBC'][subf['seq_count'].idxmax()])
                    wabbling = pd.concat([ssh,LC2,LTR1,randomBC1,ss1],1).reset_index()
                    wabbling['genomic_coordinates'] = datafram_small.iloc[0].genomic_coordinates
                    wabbling.columns = ['association_ID','shearsite','LC','LTR','randomBC', 'seq_count','genomic_coordinates']
                    wabbling = wabbling[['association_ID', 'genomic_coordinates','shearsite', 'randomBC', 'seq_count','LTR','LC']]
                    output.append(wabbling)
        elif paradox(i.shearsite.count()) > 0.5:
            i.association_ID = i.association_ID.astype(object)
            #per ogni campione controllo lo wobbling dello shearsite
            for s1,s1_df in i.groupby('association_ID'):
                #cluster dello wabbling
                gruppos = cluster(s1_df.copy().shearsite.values,2)
                '''
                per ogni gruppo unisco poiche si trovano molto vicini tra di loro ed essendo pochi per essere
                entrati in questo if allora significa che e' altamente improbabile che questa cosa sia successa per caso
                ma e' dovuta ad un errore di wabbling
                '''
                for grb in gruppos:
                    datafram_small = s1_df[s1_df.shearsite.isin(grb)]
                    x = datafram_small.groupby(['association_ID'])
                    ss1 = x.sum().dropna(0).drop('shearsite',axis=1)
                    ssh = x.apply(lambda subf: subf['shearsite'][subf['seq_count'].idxmax()])
                    LC2 = x.apply(lambda subf: subf['LC'][subf['seq_count'].idxmax()])
                    LTR1 = x.apply(lambda subf: subf['LTR'][subf['seq_count'].idxmax()])
                    randomBC1 = x.apply(lambda subf: subf['randomBC'][subf['seq_count'].idxmax()])
                    wabbling = pd.concat([ssh,LC2,LTR1,randomBC1,ss1],1).reset_index()
                    wabbling['genomic_coordinates'] = datafram_small.iloc[0].genomic_coordinates
                    wabbling.columns = ['association_ID','shearsite','LC','LTR','randomBC', 'seq_count','genomic_coordinates']
                    wabbling = wabbling[['association_ID', 'genomic_coordinates','shearsite', 'randomBC', 'seq_count','LTR','LC']]
                    output.append(wabbling)

def swap_LTR(job_id, data_slice, return_dict, events):
    CEM_LIST = ['chr8_8866486_+','chr11_64537168_-','chr17_2032352_-','chr17_47732339_-','chr2_24546571_-','chr2_73762398_-','chr16_28497498_-']
    for gb in data_slice:
        test = gb[['association_ID','seq_count']].copy().groupby('association_ID').sum()
        inizio = test[test.seq_count > 0].to_dict()['seq_count']
        sos = gb.copy()
        y = sos.groupby(['randomBC','LC'])
        ss2 = y.sum().drop('shearsite',axis=1).dropna(0)
        shearsite1 = y.apply(lambda subf: subf['shearsite'][subf['seq_count'].idxmax()])
        nome2 = y.apply(lambda subf: subf['association_ID'][subf['seq_count'].idxmax()])
        LTR1 = y.apply(lambda subf: subf['LTR'][subf['seq_count'].idxmax()])
        swapLTR = pd.concat([nome2,LTR1,shearsite1,ss2],1).reset_index()
        swapLTR['genomic_coordinates'] = gb.iloc[0].genomic_coordinates
        swapLTR.columns = ['randomBC','LC','association_ID','LTR','shearsite', 'seq_count','genomic_coordinates']
        swapLTR = swapLTR[['association_ID', 'genomic_coordinates','shearsite', 'randomBC', 'seq_count','LTR','LC']]
        
        # if ((len(swapLTR[swapLTR.genomic_coordinates.isin(CEM_LIST)]) == 0) and (len(swapLTR) < 500)):
        #     # Parte di edit distance
        #     umi_edi = swapLTR.randomBC.values
        #     for comb in combinations(umi_edi,2):
        #         if editdistance.eval(comb[0],comb[1]) <= 4:
        #             _new_ = swapLTR[swapLTR.randomBC.isin(list(comb))]
        #             swapLTR_edit = (_new_[_new_.seq_count == _new_.seq_count.max() ]).reset_index(drop=True)
        #             swapLTR_edit.loc[0,'seq_count'] = _new_.seq_count.sum()
        #             swapLTR = swapLTR[~swapLTR.randomBC.isin(list(comb))]
        #             swapLTR = pd.concat([swapLTR_edit,swapLTR],sort=False)
        #             swapLTR = swapLTR.dropna()
        
        swapLTR.columns = swapLTR.columns.get_level_values(0)
        test2 = swapLTR[['association_ID','seq_count']].copy().groupby('association_ID').sum()
        fine = test2[test2.seq_count > 0].to_dict()['seq_count']
        contaminating = set(inizio.keys()) - set(fine.keys())
        reads = np.sum([inizio[x] for x in contaminating])
        source_tmp = pd.Series(inizio)-pd.Series(fine)
        source_T = source_tmp[source_tmp < 0].to_dict()
        try:
            umi_removed = gb.randomBC.nunique()-swapLTR.randomBC.nunique()
        except Exception as e:
            umi_removed = 0
        events[gb.iloc[0].genomic_coordinates] = [reads,umi_removed,{ x : gb[gb.association_ID == x].seq_count.values for x in contaminating},source_T,gb.iloc[0].genomic_coordinates]
        return_dict[gb.iloc[0].genomic_coordinates] = swapLTR



def swappy_LC(data, job_number, threshold=3.2):
    total = len(data)
    chunk_size = total / job_number
    slice = chunks(data, chunk_size)
    jobs = []
    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    events = manager.dict()
    trash = manager.dict()

    for i, s in enumerate(slice):
        j = multiprocessing.Process(target=swap_LC, args=(i, s, return_dict, events, trash,threshold))
        jobs.append(j)
    for j in jobs:
        j.start()
    for p in jobs:
        p.join()

    return return_dict.values(),events.values(),trash.values()

def swappy_LTR(data, job_number):
    total = len(data)
    chunk_size = total / job_number
    slice = chunks(data, chunk_size)
    jobs = []
    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    events = manager.dict()

    for i, s in enumerate(slice):
        j = multiprocessing.Process(target=swap_LTR, args=(i, s, return_dict, events))
        jobs.append(j)
    for j in jobs:
        j.start()
    for p in jobs:
        p.join()

    return return_dict.values(),events.values()

def swappy_wabbly(data, job_number):
    total = len(data)
    chunk_size = total / job_number
    slice = chunks(data, chunk_size)
    jobs = []
    manager = multiprocessing.Manager()
    output = manager.list()

    for i, s in enumerate(slice):
        j = multiprocessing.Process(target=wabbly, args=(s,output))
        jobs.append(j)
    for j in jobs:
        j.start()
    for p in jobs:
        p.join()

    return output

def stats(df,final,genom_grouped,stats1,stats2,out_path):
    import warnings
    warnings.filterwarnings("ignore") 
    # define the name of the directory to be created
    path = out_path

    # sys.stdout = open('{path}/log_file.txt'.format(path=path),'wt')

    # try:  
    #     os.mkdir(path)
        
    # except OSError:  
    #     verbosePrint("Creation of the directory %s failed because exist" % path)
    # else:  
    #     verbosePrint("Successfully created the directory %s " % path)

    log_file = open("{path}/log_file.txt".format(path=path), "a")

    Contaminated_LTR = []
    pollutant_LTR = []

    IS_TOTAL_with_swap_LTR = []
    read_swapped_total_LTR = 0
    mean_read_swapped_total_LTR = []

    Contaminated_LC = []
    pollutant_LC = []

    IS_TOTAL_with_swap_LC = []
    read_swapped_total_LC = 0
    mean_read_swapped_total_LC = []

    IS_Contaminated_LTR = {} 
    IS_Contaminated_LC =  {}

    umi_removed = {}
    maxs1 = {}
    maxs2 = {}
    c1 = 0
    for n1,arr1 in enumerate(stats1):
        if arr1[0] > c1:
            c1 = arr1[0]
            maxs1 = {arr1[4]:c1}
        if arr1[0] != 0:
            read_swapped_total_LTR += arr1[0]  
            mean_read_swapped_total_LTR.extend([k for row in arr1[2].values() for k in row])
            IS_TOTAL_with_swap_LTR.append(arr1[4])
            Contaminated_LTR.extend([k for k in arr1[2].keys()])
            pollutant_LTR.extend([k for k in arr1[3].keys()])
            IS_Contaminated_LTR[arr1[4]] = arr1[0]
        elif arr1[1] != 0:
            umi_removed[arr1[4]] = arr1[1]
        else:
            continue
    if len(maxs1) == 0:
        maxs1['NA'] = 0
    c2 = 0
    for n2,arr2 in enumerate(stats2):
        if arr2[0] > c2:
            c2 = arr2[0]
            maxs2 = {arr2[4]:c2}            
        if arr2[0] != 0:
            read_swapped_total_LC += arr2[0]   
            mean_read_swapped_total_LC.extend([k for row in arr2[2].values() for k in row])
            IS_TOTAL_with_swap_LC.append(arr2[4])
            Contaminated_LC.extend([k for k in arr2[2].keys()])
            pollutant_LC.extend([k for k in arr2[3].keys()])
            IS_Contaminated_LC[arr2[4]] = arr2[0]
        elif arr2[1] != 0:
            umi_removed[arr2[4]] = arr2[1]
        else:
            continue
    if len(maxs2) == 0:
        maxs2['NA'] = 0

    IS_Contaminated_LTR = pd.Series(IS_Contaminated_LTR).to_frame()
    IS_Contaminated_LC = pd.Series(IS_Contaminated_LC).to_frame()
    IS_Contaminated_LTR['type_swap'] = 'LTR'
    IS_Contaminated_LC['type_swap'] = 'LC'
    IS_Contaminated = pd.concat([IS_Contaminated_LTR,IS_Contaminated_LC]).reset_index().rename({'index':'IS_Contaminated',0 : '#reads'},axis=1)
    if len(IS_Contaminated) != 0:
        IS_Contaminated.to_csv('{path}/{x}_IS_contaminated.tsv'.format(path=path,x=contaminationFile.split('.')[0].split('/')[-1].replace('_contamination','')),index=False,sep='\t')
    
    contaminated_tmp1_ = pd.Series(Contaminated_LTR).value_counts().to_frame().reset_index().rename({'index':'Contaminated_sample',0 : '#IS_Contaminated'},axis=1)#.to_csv('path/samples_contaminated_LTR_swap.tsv',index=False,sep='\t')
    contaminated_tmp2_ = pd.Series(Contaminated_LC).value_counts().to_frame().reset_index().rename({'index':'Contaminated_sample',0 : '#IS_Contaminated'},axis=1)#.to_csv('path/samples_contaminated_LC_swap.tsv',index=False,sep='\t')
    contaminated_tmp1_['type_swap'] = 'LTR'
    contaminated_tmp2_['type_swap'] = 'LC'

    #campioni contaminati
    contaminated_ = pd.concat([contaminated_tmp1_,contaminated_tmp2_])
    if len(contaminated_) != 0:
        contaminated_.to_csv('{path}/{x}_samples_contaminated.tsv'.format(path=path,x=contaminationFile.split('.')[0].split('/')[-1].replace('_contamination','')),index=False,sep='\t')

    pollutant_tmp1_ = pd.Series(pollutant_LTR).value_counts().to_frame().reset_index().rename({'index':'Pollutant_sample',0 : '#IS_Pollutant'},axis=1)
    pollutant_tmp2_ = pd.Series(pollutant_LC).value_counts().to_frame().reset_index().rename({'index':'Pollutant_sample',0 : '#IS_Pollutant'},axis=1)
    pollutant_tmp1_['type_swap'] = 'LTR'
    pollutant_tmp2_['type_swap'] = 'LC'

    #campioni che hanno contaminato
    pollutant_ = pd.concat([pollutant_tmp1_,pollutant_tmp2_])
    if len(pollutant_) != 0:
        pollutant_.to_csv('{path}/{x}_samples_pollutant.tsv'.format(path=path,x=contaminationFile.split('.')[0].split('/')[-1].replace('_contamination','')),index=False,sep='\t')

    import pylab
    #numero di integrazioni prima e dopo
    IS_initial = len(genom_grouped)
    CEM_LIST = ['chr8_8866486_+','chr11_64537168_-','chr17_2032352_-','chr17_47732339_-','chr2_24546571_-','chr2_73762398_-','chr16_28497498_-']
    
    try:
        o_cem = df[(df.association_ID.apply(lambda x: 'CEM' in x)) & ~(df.genomic_coordinates.isin(CEM_LIST))]
        n_cem = final[(final.association_ID.apply(lambda x: 'CEM' in x)) & ~(final.genomic_coordinates.isin(CEM_LIST))]
        cem_is_o = o_cem.genomic_coordinates.nunique()
        cem_is_n = n_cem.genomic_coordinates.nunique()
        reduction = round((1-(cem_is_n/float(cem_is_o)))*100,2)
        

        ccc = df[df.association_ID.apply(lambda x: 'CEM' in x)]
        TP = ccc[ccc.genomic_coordinates.isin(CEM_LIST)].seq_count.sum()
        FP = ccc[~ccc.genomic_coordinates.isin(CEM_LIST)].seq_count.sum()
        FN = df[(df.genomic_coordinates.isin(CEM_LIST)) & (df.association_ID.apply(lambda x: 'CEM' not in x))].seq_count.sum()
        precision = TP/float((TP+FP))
        recall = TP/float((TP+FN))
        f_measure = 2*(precision*recall)/float((precision+recall))
        f_measure = round(f_measure,8)

        aaa = final[final.association_ID.apply(lambda x: 'CEM' in x)]
        TP2 = aaa[aaa.genomic_coordinates.isin(CEM_LIST)].seq_count.sum()
        FP2 = aaa[~aaa.genomic_coordinates.isin(CEM_LIST)].seq_count.sum()
        FN2 = final[(final.genomic_coordinates.isin(CEM_LIST)) & (final.association_ID.apply(lambda x: 'CEM' not in x))].seq_count.sum()
        precision2 = TP2/float((TP2+FP2))
        recall2 = TP2/float((TP2+FN2))
        f_measure2 = 2*(precision2*recall2)/float((precision2+recall2))
        f_measure2 = round(f_measure2,8)
    except Exception as e:
        o_cem = 'NA'
        n_cem = 'NA'
        cem_is_o = 'NA'
        cem_is_n = 'NA'
        reduction = 'NA'
        reduction = 'NA'
        TP ='NA'
        FP ='NA'
        FN ='NA'
        precision = 'NA'
        recall = 'NA'
        f_measure = 'NA'
        TP2 = 'NA'
        FP2 = 'NA'
        FN2 = 'NA'
        precision2 = 'NA'
        recall2 = 'NA'
        f_measure2 = 'NA'

    IS_final = final.genomic_coordinates.nunique()
    n_IS_removed = IS_initial-IS_final
    percentage_of_IS_removed = (1-(IS_final/float(IS_initial)))*100
    percentage_of_IS_contaminated = (len(IS_Contaminated)/float(IS_final))*100
    IS_LTR = len(set(IS_Contaminated_LTR.index))
    IS_LC = len(set(IS_Contaminated_LC.index))
    shared_IS_LTR_LC = len(set(IS_Contaminated_LTR.index).intersection(IS_Contaminated_LC.index))
    #numero di sample contaminated
    total_samples = final.association_ID.nunique()
    n_polluntant_samples = pollutant_.Pollutant_sample.nunique()
    n_Contaminated_samples = contaminated_.Contaminated_sample.nunique()
    percentage_of_samples_contaminated = (n_Contaminated_samples/float(total_samples))*100
    try:
        more_polluntant_sample = pollutant_[pollutant_['#IS_Pollutant'] == pollutant_['#IS_Pollutant'].max()].Pollutant_sample.tolist()[0]
        more_contaminated_sample = contaminated_[contaminated_['#IS_Contaminated'] == contaminated_['#IS_Contaminated'].max()].Contaminated_sample.tolist()[0]
    except Exception as e:
        more_polluntant_sample = 'NA'
        more_contaminated_sample = 'NA'
    #numero di reads
    reads_totali = final.seq_count.sum()
    percentage_of_reads_contaminated = round((read_swapped_total_LTR+read_swapped_total_LC)/float(reads_totali)*100,4)
    type_of_swap = ['LTR', 'LC']    
    try:
        x_pos = np.arange(len(type_of_swap))
        CTEs = [np.mean(mean_read_swapped_total_LTR),np.mean(mean_read_swapped_total_LC)]
        error = [np.std(mean_read_swapped_total_LTR),np.std(mean_read_swapped_total_LC)]
        fig, ax = plt.subplots()
        ax.bar(x_pos, CTEs, yerr=error, align='center', alpha=0.5, ecolor='black', capsize=10)
        #First Figure
        ax.set_ylabel('Mean Sequence count in Switch of index')
        ax.set_xticks(x_pos)
        ax.set_xticklabels(type_of_swap)
        ax.set_title('Sequence count average in index Switching in pool {x}'.format(x=contaminationFile.split('.')[0].split('/')[-1].replace('_contamination','')), fontsize=9)
        ax.yaxis.grid(True)
        ax.xaxis.grid(True)
        plt.tight_layout()
        ax.get_figure().savefig('{path}/{x}_bar_plot_with_error_bars_average_seq_count.png'.format(path=path,x=contaminationFile.split('.')[0].split('/')[-1].replace('_contamination','')),format='png', dpi=500)
    except Exception as e:
        pass
    #Second Figure LTR
    try:
        sq_times = pd.Series(mean_read_swapped_total_LTR).value_counts().sort_index().to_frame().rename({0:'Number of times found'},axis=1)
        sq_times['Percentage %'] = sq_times['Number of times found']/sq_times['Number of times found'].sum()*100
        fig1 = plt.figure()
        ax2 = fig1.add_subplot(111)
        ax2 = sq_times.plot(kind='barh', legend=False)
        plt.xscale('log')
        ax2.set_title('Sequence count found with LTR index Switching in pool {x}'.format(x=contaminationFile.split('.')[0].split('/')[-1].replace('_contamination','')), fontsize=8,y=1.13)
        ax2.set_xlabel('Value (Log10)')
        ax2.set_ylabel('Seq_count number')
        # ax2_t = ax2.twiny()
        # ax2_t.set_xlabel('Value (linear)')
        # plt.xticks([0,1,10,100,1000])
        ax2.yaxis.grid(True)
        ax2.xaxis.grid(True)
        plt.tight_layout()
        plt.legend()
        ax2.get_figure().savefig('{path}/{x}_bar_plot_seq_count_found_LTR_SWITCH.png'.format(path=path,x=contaminationFile.split('.')[0].split('/')[-1].replace('_contamination','')),format='png', dpi=500)
        plt.cla()
    except Exception as e:
        pass
    #third Figure LC 
    try:
        sq_times2 = pd.Series(mean_read_swapped_total_LC).value_counts().sort_index().to_frame().rename({0:'Number of times found'},axis=1)
        sq_times2['Percentage %'] = sq_times2['Number of times found']/sq_times2['Number of times found'].sum()*100
        fig2 = plt.figure()
        ax3 = fig2.add_subplot(111)
        ax3 = sq_times2.plot(kind='barh', legend=False)
        plt.xscale('log')
        ax3.set_title('Sequence count found with LC index Switching in pool {x}'.format(x=contaminationFile.split('.')[0].split('/')[-1].replace('_contamination','')), fontsize=8,y=1.13)
        ax3.set_xlabel('Value (Log10)')
        ax3.set_ylabel('Seq_count number')
        ax3.yaxis.grid(True)
        ax3.xaxis.grid(True)
        plt.tight_layout()
        plt.legend()
        ax3.get_figure().savefig('{path}/{x}_bar_plot_seq_count_found_LC_SWITCH.png'.format(path=path,x=contaminationFile.split('.')[0].split('/')[-1].replace('_contamination','')),format='png', dpi=500)
        plt.cla()
    except Exception as e:
        pass
    #numero di shearsite
    try:
        shs_initial = df.groupby('genomic_coordinates').shearsite.size().describe()
        shs_final = final.groupby('genomic_coordinates').shearsite.size().describe()
        file1 = df.groupby('genomic_coordinates').shearsite.size()
        file2 = final.groupby('genomic_coordinates').shearsite.size()
        file3 = file1.to_frame().join(file2.to_frame(),lsuffix='_before', rsuffix='_after').sort_values('shearsite_before',ascending=False)
        file3.to_csv('{path}/{x}_cells_contaminated_per_IS.tsv'.format(path=path,x=contaminationFile.split('.')[0].split('/')[-1].replace('_contamination','')),index=True,sep='\t')
    except Exception as e:
        pass
    try:
        fig3 = plt.figure()
        ax4 = fig3.add_subplot(111)
        shs1 = df.shearsite.copy()
        shs1.name = 'Shearsite before filter'
        ax4 = sns.kdeplot(shs1, shade=True,color="r")
        shs2 = final.shearsite.copy()
        shs2.name = 'Shearsite after filter'
        ax4 = sns.kdeplot(shs2, shade=True,color="g")
        ax4.set_xlabel('Shearsite type')
        ax4.set_ylabel('Distribution')
        ax4.set_title('Distribution of shearsite in pool {x}'.format(x=contaminationFile.split('.')[0].split('/')[-1].replace('_contamination','')), fontsize=9)
        ax4.yaxis.grid(True)
        ax4.xaxis.grid(True)
        ax4.get_figure().savefig('{path}/{x}_Distribution_shearsite.png'.format(path=path,x=contaminationFile.split('.')[0].split('/')[-1].replace('_contamination','')),format='png', dpi=500)
    except Exception as e:
        pass
    #numero di cellule/UMI
    umi_cell = pd.Series(umi_removed)

   
    ########### Display ####################################
    print('1. Information on ISs:',file=log_file)
    print('>>> Number of initial ISs Found:\t'+'{x}'.format(x=IS_initial),file=log_file)
    print('>>> Number of Final ISs Found:\t\t'+'{x}'.format(x=IS_final),file=log_file)
    print('>>> Number of ISs filtered out:\t\t'+'{x}'.format(x=n_IS_removed),file=log_file)
    print('>>> .%. of ISs filtered out:\t\t'+'{p:4.2f} %'.format(p=percentage_of_IS_removed),file=log_file)
    print('>>> Number of ISs Contaminated:\t\t'+'{x}'.format(x=len(IS_Contaminated))+' out of {y}'.format(y=IS_final),file=log_file)
    print('>>> %. of ISs Contaminated:\t\t'+'{p:4.2f} %'.format(p=percentage_of_IS_contaminated),file=log_file)
    print('>>> ISs more contaminated for LTR, LC swaps respectively:\t'+'{x}'.format(x=(list(maxs1.keys()),list(maxs2.keys()))),file=log_file)
    print('>>> Number of shared IS btw swap LTR and LC (Venny-DiaGram):\t{y1}  '.format(y1=IS_LTR)+'{x}'.format(x=(shared_IS_LTR_LC))+'  {y2}'.format(y2=IS_LC),file=log_file)
    if cem_is_o != 0:
        print('>>> Number of ISs Contaminated in CEM samples before:\t\t'+'{x}'.format(x=cem_is_o),file=log_file)
        print('>>> Number of ISs Contaminated in CEM samples After:\t\t'+'{x}'.format(x=cem_is_n),file=log_file)
        print('>>> Reduction of ISs Contaminated in CEM samples:\t\t'+'{x}'.format(x=reduction),file=log_file)
        try:
            n_cem.to_csv('{path}/{x}_IS_remained_inCEM.tsv'.format(path=path,x=contaminationFile.split('.')[0].split('/')[-1].replace('_contamination','')),index=False,sep='\t')
        except Exception as e:
            pass
    ########### Display ####################################
    print('2. Information on Samples:',file=log_file)
    print('>>> Number of Contaminated samples:\t\t'+'{x}'.format(x=n_Contaminated_samples)+' out of {y}'.format(y=df.association_ID.nunique()),file=log_file)
    print('>>> Number of Pollutant samples:\t\t'+'{x}'.format(x=n_polluntant_samples)+' out of {y}'.format(y=df.association_ID.nunique()),file=log_file)
    print('>>> %. of Contaminated samples:\t\t\t'+'{p:4.2f} %'.format(p=percentage_of_samples_contaminated),file=log_file)
    print('>>> Sample more contaminated:\t\t'+'{x}'.format(x=(more_contaminated_sample)),file=log_file)
    print('>>> Sample more pollutant:\t\t'+'{x}'.format(x=(more_polluntant_sample)),file=log_file)
    ########### Display ####################################
    print('3. Information on Reads:',file=log_file)
    print('>>> Number of Contaminated Reads:\t\t'+'{x}'.format(x=read_swapped_total_LTR+read_swapped_total_LC)+' out of {y}'.format(y=reads_totali),file=log_file)
    print('>>> %. of Contaminated Reads:\t\t\t'+'{p:6.4f} %'.format(p=percentage_of_reads_contaminated),file=log_file)
    print('>>> Contaminated Reads for LTR:\t\t\t'+'{x}'.format(x=read_swapped_total_LTR),file=log_file)
    try:
        print('>>> Mean contaminated Reads for LTR:\t\t'+'{x}'.format(x=round(np.mean(mean_read_swapped_total_LTR),2)),file=log_file)
        print('>>> Median contaminated Reads for LTR:\t\t'+'{x}'.format(x=round(np.median(mean_read_swapped_total_LTR),2)),file=log_file)
        print('>>> standard dev contaminated Reads for LTR:\t'+'{x}'.format(x=round(np.std(mean_read_swapped_total_LTR),2)),file=log_file)
        print('>>> max() contaminated Reads for LTR:\t\t'+'{x}'.format(x=np.max(mean_read_swapped_total_LTR)),file=log_file)
    except:
        print('>>> Mean contaminated Reads for LTR:\t\t'+'{x}'.format(x='NA'),file=log_file)
        print('>>> Median contaminated Reads for LTR:\t\t'+'{x}'.format(x='NA'),file=log_file)
        print('>>> standard dev contaminated Reads for LTR:\t'+'{x}'.format(x='NA'),file=log_file)        
        print('>>> max() contaminated Reads for LTR:\t\t'+'{x}'.format(x='NA'),file=log_file)
    print('>>> 1-2-3 seq_count contaminated for LTR found:\t'+'{x} %'.format(x=round(sq_times.loc[0:3,'Percentage %'].sum(),2)),file=log_file)
    print('>>> Contaminated Reads for LC:\t\t\t'+'{x}'.format(x=read_swapped_total_LC),file=log_file)
    try:
        print('>>> Mean contaminated Reads for LC:\t\t'+'{x}'.format(x=round(np.mean(mean_read_swapped_total_LC),2)),file=log_file)
        print('>>> Median contaminated Reads for LC:\t\t'+'{x}'.format(x=round(np.median(mean_read_swapped_total_LC),2)),file=log_file)
        print('>>> standard dev contaminated Reads for LC:\t'+'{x}'.format(x=round(np.std(mean_read_swapped_total_LC),2)),file=log_file)
        print('>>> max() contaminated Reads for LC:\t\t'+'{x}'.format(x=np.max(mean_read_swapped_total_LC)),file=log_file)
    except:
        print('>>> Mean contaminated Reads for LC:\t\t'+'{x}'.format(x='NA'),file=log_file)
        print('>>> Median contaminated Reads for LC:\t\t'+'{x}'.format(x='NA'),file=log_file)
        print('>>> standard dev contaminated Reads for LC:\t'+'{x}'.format(x='NA'),file=log_file)
        print('>>> max() contaminated Reads for LC:\t\t'+'{x}'.format(x='NA'),file=log_file)
    print('>>> 1-2-3 seq_count contaminated for LC found:\t'+'{x} %'.format(x=round(sq_times2.loc[0:3,'Percentage %'].sum(),2)),file=log_file)
    ########### Display ####################################
    print('4. Information on Shearsite/Cells:',file=log_file)
    print('4.1 Before Filter:',file=log_file)
    print('>>> Mean Cells per IS before filter:\t\t'+'{x}'.format(x=shs_initial['mean'].round(2)),file=log_file)
    print('>>> standard dev Cells per IS before filter:\t'+'{x}'.format(x=shs_initial['std'].round(2)),file=log_file)
    print('>>> Median Cells per IS before filter:\t\t'+'{x}'.format(x=shs_initial['50%']),file=log_file)
    print('>>> max() Cells per IS before filter:\t\t'+'{x}'.format(x=shs_initial['max']),file=log_file)
    print('4.2 AFTER Filter:',file=log_file)
    print('>>> Mean Cells per IS AFTER filter:\t\t'+'{x}'.format(x=shs_final['mean'].round(2)),file=log_file)
    print('>>> standard dev Cells per IS AFTER filter:\t'+'{x}'.format(x=shs_final['std'].round(2)),file=log_file)
    print('>>> Median Cells per IS AFTER filter:\t\t'+'{x}'.format(x=shs_final['50%']),file=log_file)
    print('>>> max() Cells per IS AFTER filter:\t\t'+'{x}'.format(x=shs_final['max']),file=log_file)
    ########### Display ####################################
    print('4. Information on UMI/Cells:',file=log_file)
    print('>>> UMI/Genomes removed Total:\t\t\t'+'{x}'.format(x=umi_cell.sum()),file=log_file)
    
    # run_time = str(datetime.timedelta(seconds=(round(time.time() - start_time ))))+' hh:mm:ss'
    # output_name = '{x}_filtered_contamination.csv.gz'.format(x=output)
    # stats = { 'nCpu' : jobs, 'line old file' : len(a), 'line new file' : len(b),'number of line removed' : diff, 'Percentage of line removed' : percentage, 'UMI removed' : spread_umi, 'swap_event_LTR' : swap_event_LTR, 'swap_event_LC' : swap_event_LC ,'Shearsite removed' : shs_obsolete, 'Run time' : run_time, 'output name' : output_name }
    # (pd.Series(stats, name = 'stats').to_frame()).to_csv('log_file_{x}'.format(x=output), sep='\t', index=True)
        

if __name__ == "__main__":
    done = '[' + ' Done ' + ']'
    end = '[' + ' END ' + ']'
    check_input_paths (contaminationFile)
    buildOUTDIR(output)
    out_path = build_outpath (output)
    verbosePrint('\n> reading input file...')
    genom_grouped,df = groups(contaminationFile)
    verbosePrint(done)

    ######## SWAP LTR #################################
    ###################################################
    verbosePrint('> SWAP LTR...')
    step1,stats1 = swappy_LTR(genom_grouped, jobs)
    verbosePrint(done)
    ######## SWAP LC ##################################
    ###################################################
    verbosePrint('> SWAP LC...')
    step3,stats2,trash = swappy_LC(step1, jobs)
    verbosePrint(done)
    ######## Wabbling #################################
    ###################################################
    # print('> Wabbling', end='')
    # step3 = swappy_wabbly(step2, jobs) 
    # print('\t\t[' + ' Done ' + ']')

    ##### ASSEMBLING
    verbosePrint('> concatenation...')
    final = concatenat(step3)
    verbosePrint(done)
    
    verbosePrint('\n> writing log_file of this operation - ')
    stats(df,final,genom_grouped,stats1,stats2,out_path)
    verbosePrint(done)
    # print('\t\t[' + ' Done ' + ']')
    verbosePrint('\n> writing results in [ {out_path} ]'.format(out_path=out_path))
    (final.drop(['LTR','LC'],axis=1)).to_csv('{out_path}/{x}_swap_removed.csv.gz'.format(out_path=out_path,x=contaminationFile.split('.')[0].split('/')[-1].replace('_contamination','')), sep='\t', index=False, compression='gzip')
    os.system('python /home/adriano/skip_project/create_matrix_adriano/_matrix_MAIN.py a,b,{out_path},{out_path}/{filename}_swap_removed.csv.gz --dataset_ID {filename}'.format(out_path=out_path,filename=contaminationFile.split('.')[0].split('/')[-1].replace('_contamination','')))
    verbosePrint(done)
    print()
    verbosePrint(end)

    run_time = str(datetime.timedelta(seconds=(round(time.time() - start_time ))))+' hh:mm:ss'

    # print(run_time,file=open("{out_path}/log_file.txt".format(out_path=out_path)), "a")
    verbosePrint('Run Time = '+run_time)


