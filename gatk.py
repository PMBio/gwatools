import scipy as sp
import gzip as gz
import pdb
import os
import fnmatch

def gt_lookup(gt):
    
    if gt == '0/0':
        return 0
    if gt in ['0/1', '1/0']:
        return 1
    if gt == '1/1':
        return 2
    if gt == './.':
        return -1
    return -2

def sumAll(aldic):
    ids = []
    pos = []
    snps = []
    for subdic in aldic:
        ids.extend(subdic['ids'])
        pos.extend(subdic['pos'])
        snps.extend(subdic['snps'])
    ids = sp.array(ids)
    pos = sp.array(pos)
    snps = sp.array(snps)
    
  #  posIdx = sp.lexsort((pos[:,1],pos[:,0]))
  #  
  #  ids = ids[posIdx]
  #  pos = pos[posIdx]
  #  snps = snps[posIdx]
    return {'ids':ids, 'pos':pos, 'snps':snps}
    
def parseDir(baseDir):
    allDic = []
    for root, dirs, files in os.walk(baseDir):
        files = [f for f in fnmatch.filter(files, '*.vcf')]
        for f in sorted(files, key = lambda x: (int(x.split('.')[2]) if x.split('.')[2].isdigit() else x.split('.')[2], int(x.split('.')[3]), int(x.split('.')[4]))):
            print f
            subDic = readVCF(os.path.join(root,f))
            allDic.append(subDic)
    return sumAll(allDic)
            
def readVCF(filename, delimiter = '\t', missing = './.', valid_exome = None, valid_gt = None, valid_sample = None):

    if filename[-2:] == 'gz':
        infile = gz.open(filename, 'rb')
    else:
        infile = open(filename, 'r')
    
    data = sp.loadtxt(infile, dtype = 'str', delimiter = delimiter, comments = '##') 

    ### vcf file was empty
    if len(data.shape) < 2:
        return dict()

    ### data[0, :] is line of headers
    ids = sp.array(data[0, 9:], dtype = 'str')

    ### replace X, Y chrms with PLINK convention (x-> 23, Y->24, pseudoA->25, MT->26)
    #data[data[:, 0] == 'X', 0] = '23'
    #data[data[:, 0] == 'Y', 0] = '24'

    ### remove all lines with more than 2 alleles
    ### TODO make max allele count filter configurable
    data = data[[i for i in xrange(data.shape[0]) if not ',' in data[i, 4]], :]

    pos = data[1:,0:2]#sp.array(data[1:,0:2], dtype = 'float')
    qual = sp.array(data[1:, 5], dtype = 'float')
    low_qual = sp.array(data[1:, 4] == 'LowQual', dtype = 'bool')

    data = data[:, 9:]
    if not valid_exome == None:
        v_idx = sp.where(sp.in1d(ids, valid_exome))[0]
        ids = ids[v_idx]
        data = data[:, v_idx]

    ids_shortened = False
    if not valid_gt == None:
        ids = sp.array([x.rsplit('-', 4)[0] for x in ids], dtype = 'str')
        v_idx = sp.where(sp.in1d(ids, valid_gt))[0]
        ids = ids[v_idx]
        data = data[:, v_idx]
        ids_shortened = True

    if not valid_sample == None:
        if not ids_shortened:
            ids = sp.array([x.rsplit('-', 4)[0] for x in ids], dtype = 'str')
        v_idx = sp.where(sp.in1d(ids, valid_sample))[0]
        ids = ids[v_idx]
        data = data[:, v_idx]

    snps = sp.zeros((data.shape[0] - 1, data.shape[1]), dtype = 'float')
    snps[:] = sp.nan
    isnan = data[1:,:] == missing
    rm_idx = [];
    for i in xrange(snps.shape[0]):
        tmp_idx = sp.where(~isnan[i, :])
        if tmp_idx[0].shape[0] > 0:
            #snps[i, tmp_idx] = map(lambda x: gt_lookup(x[0].split('\t')[0].split(':')[0]), data[i+1, tmp_idx])
            snps[i, tmp_idx] = map(lambda x: gt_lookup(x.split('\t')[0].split(':')[0]), data[i+1, tmp_idx][0].tolist())
                         
    infile.close()

    return {'ids':ids, 'pos':pos, 'snps':snps, 'qual':qual, 'low_qual':low_qual}

