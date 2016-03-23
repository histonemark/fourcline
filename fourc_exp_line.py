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

# This script processes the gDNA and cDNA of barcoded insertions per clone. All
# the replicates should be passed in order gdna1,gdna2,cdna1,cdna2

bcd_dict = defaultdict(int)
processed = []
spikessed = []
def call_starcode_fastq_file(fastq):
   #pdb.set_trace()
   MIN_BRCD = 15
   MAX_BRCD = 25
      
   brcd_outfname = fname + '_barcodes.tsv'
   spk_outfname  = fname + '_spikes.tsv'
   
   GFP   = seeq.compile('CATGCTAGTTGTGGTTTGTCCAAACT', 4)
   SPIKE = seeq.compile('CATGATTACCCTGTTATC', 2)
   barcode_tempf = tempfile.NamedTemporaryFile(delete=False)
   spike_tempf   = tempfile.NamedTemporaryFile(delete=False)
   
   with gzopen(fastq) as f:
      outf = None
      for lineno,line in enumerate(f):
         if lineno % 4 != 1: continue
         hit = GFP.match(line)
         if hit is not None:
            outf = barcode_tempf
         else:
            hit = SPIKE.match(line)
            if hit is not None:
               outf = spike_tempf
            else:
               continue
         pos = hit.matchlist[0][0]
         if MIN_BRCD <= pos <= MAX_BRCD:
            outf.write(line[:pos] + '\n')
   barcode_tempf.close()
   spike_tempf.close()

   subprocess.call([
      'starcode',
      '-t4',
      '-i', barcode_tempf.name,
      '-o', brcd_outfname,])
      
   subprocess.call([
      'starcode',
      '-t4',
      '-i', spike_tempf.name,
      '-o', spk_outfname,])

   # Delete temporary files.
   os.unlink(barcode_tempf.name)
   os.unlink(spike_tempf.name)

   # Save the names of the files processsed
   #processed.append([brcd_outfname,spk_outfname])
   processed.append(brcd_outfname)
   spikessed.append(spk_outfname)
   #pdb.set_trace()
   return

countsd = defaultdict(lambda: defaultdict(int))
spikesd = defaultdict(lambda: defaultdict(int))
names = []

def gather_counts():
   if len(processed) == 2:
      names = ['gDNA','cDNA']
   elif len(processed) == 4:
      names = ['gDNA1','gDNA2','cDNA1','cDNA2']
   else:
      print 'You Can only process 2 or 4 files(replicates)'

   for idx,procf in enumerate(processed):
      with open(procf,'r') as f:
         for line in f:
            bcd,count = line.rstrip().split()
            countsd[bcd][names[idx]] = count
   #pdb.set_trace()
   
   exp_fname = processed[0].split('.')[0] + '_expression.tsv'
      
   with open(exp_fname, 'w') as g:
      g.write('barcode\tgDNA1\tgDNA2\tcDNA1\tcDNA2\n')
      lnames = len(names)
      for k in countsd:
         if lnames == 2:
            gdna = countsd[k]['gDNA']
            cdna = countsd[k]['cDNA']
            g.write('%s\t%s\t%s\n' % (k,gdna,cdna))
         else:
            gdna1 = countsd[k]['gDNA1']
            gdna2 = countsd[k]['gDNA2']
            cdna1 = countsd[k]['cDNA1']
            cdna2 = countsd[k]['cDNA2']
            g.write('%s\t%s\t%s\t%s\t%s\n' % \
                    (k,gdna1,gdna2,cdna1,cdna2))
      
      # We add the spikes
      #pdb.set_trace()
      for fnum,fname in enumerate(spikessed):
         with open(fname) as f:
            for line in f:
               columns = ['0' for i in names]
               bcd,count = line.split()
               columns[fnum] = str(count)
               out = '\t'.join(columns)
               g.write('spike\t%s\n' % out)
               columns = ['0' for i in names]


for fname in sys.argv[1:]:
   call_starcode_fastq_file(fname)

gather_counts()