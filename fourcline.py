import fileinput
import os
import pdb
import re
import seeq
import sys
import subprocess
import tempfile
from collections import defaultdict
from gzopen import gzopen
from itertools import izip

TOMAPfname = sys.argv[1] + '_2map'
#pdb.set_trace()
with gzopen(sys.argv[1]) as f, open(TOMAPfname,'w') as g:
   for lineno,line in enumerate(f):
      # Is a fastq keep only sequence
      if lineno % 4 != 1: continue
      # Exact search of NlaIII
      brcd = line.rstrip().split('CATG')[0]
      if len(brcd) == len(line.rstrip()): continue
      seq  = line.rstrip().split('CATG')[1]   
      # Cut if there is a MlucI site
      dna  = seq.split('AATT')[0]
      # Write fasta to map it 
      if not 10 < len(brcd) < 22 : continue
      if not 5 < len(dna): continue 
      g.write('>%s\n%s\n' % (brcd,dna))


# Map the sequences
      
outfname = sys.argv[1]  
INDEX = '/mnt/shared/seq/gem/dm3R5/dm3R5_pT2_unmasked.gem'


# System call to `gem-mapper` passing the desired arguments.
subprocess.call([
   'gem-mapper',
   '-I', INDEX ,
   '-i', TOMAPfname,
   '-o', outfname,
   '-m3',
   '--unique-mapping',
   '-T4'])

# Open the mapped file and keep barcode position if uniquely mapped
experiment = open(sys.argv[1] + '_experiment.tsv','w')
BCDICT = defaultdict(int)
MAPPED = sys.argv[1] + '.map'
previousln = ''
with open(MAPPED) as g:
   for lineno,line in enumerate(g):
      #pdb.set_trace()
      items = line.rstrip().split()
      if items[-1] == '-': continue
      # last column contains chr:strand:pos:matches
      bcd = items[0]
      loc = ':'.join(items[-1].split(':')[:3])
      if len(loc.split(':')) < 3:
         print bcd
      BCDICT[bcd + ':' + loc] += 1

   for k in BCDICT:
      items = k.split(':')
      if len(items) != 4:
         print 'Error in ', k
      bcd,chrom,strand,pos = k.split(':')
      count = BCDICT[k]
      out = (bcd,count,chrom,strand,pos)
      experiment.write('%s\t%s\t%s\t%s\t%s\n' % out)

# Append original barcodes with experiment barcode to starcode

originals  = sys.argv[2] # bcd count
experiment = sys.argv[1] + '_experiment.tsv'
allbcdsfn  = sys.argv[1] + '_allbcd.tsv'

with open(allbcdsfn, 'w') as outf, \
     open(originals, 'r') as ori, \
     open(experiment, 'r') as ex:
   try:
      for lineno,line in enumerate(ori):
         looping = 'Originals'
         bcd,count = line.rstrip().split()[:2] 
         outf.write(bcd + '\t' + count + '\n')
   except ValueError:
      pdb.set_trace()
      print 'The error is in the Originals file'
      print line
      
   try:
      for lineno,line in enumerate(ex):
         looping = 'barcodes'
         bcd,count = line.rstrip().split()[:2] 
         outf.write(bcd + '\t' + count + '\n')
   except ValueError:
      print 'The error is in the experiment file'
      print line
      
      
# Call starcode on the files
starcoded  = sys.argv[1] + '_starcoded.txt'
subprocess.call([
   'starcode',
   '-t4',
   '-d2',
   '--print-clusters',
   '-i' ,
   allbcdsfn,
   '-o',
   starcoded])

# Change barcodes for it's canonical and print final file

CANON = defaultdict(str)
with open(starcoded) as f:
   for line in f:
      items = line.split()
      canonical = items[0]
      alternate = list(items[2].split(','))
      for bcd in alternate:
         CANON[bcd] = canonical
      

INF  = sys.argv[1] + '_experiment.tsv'         
FINAL = sys.argv[1] + '_final.tsv'

with open(INF,'r') as f, open(FINAL,'w') as g:
   for lineno,line in enumerate(f):
      items  = line.split()
      bcd    = items[0]
      #count  = int(items[1]) - 10^6
      rest   = items[1:]
      newbcd = CANON[bcd]
      if newbcd == '': continue
      g.write(newbcd + ' ' + ' '.join(rest) + '\n')
      
      

         
      
   

         
