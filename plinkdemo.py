import plink
import os
plinkpath = "/ebi/research/stegle/projects/1000GenomesRNASeq/data/genotypes/plink"

bfile = os.path.join(plinkpath,'maf1_chr22')

print "the plink file reader returns a structure containing the SNPs and meta information."
if 0:
    print "load all of chr 22:"
    s1 = plink.readBED(bfile)
    print "done.\n"
if 1:
    print "load SNPs with index (starting with 0) 20000-20020:"
    s2 = plink.readBED(bfile, start=20000,nSNPs=20)
    print "done.\n"
if 1:
    print "load all SNPs on chr 22 between basepair positions 18,000,000 to 18,100,000\nnote that the second number refers to genomic distance (cM) and is not specified in our files."
    s3 = plink.readBED(bfile,startpos=[22,0,18000000],endpos = [22,0,18100000])
    print "done.\n"
    print 'keys of the genotype struct loaded:'
    print s3.keys()
    print '\nsize of the SNP matrix:'
    print s3['snps'].shape
