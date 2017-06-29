#!/usr/bin/env python
import sys
def main(filename,outfilename):
  f=open(filename)
  g=open(outfilename,'w')
  g2=open('zp_transcripts.txt','w')
  line=f.readline()
  for line in f: 
    line=line.rstrip()
    parts=line.split('\t')
    if 'vitellogenin' in parts[1].lower() or 'vitellogenin' in parts[2].lower():
      idx=parts[0]
      g.write(idx+'\t'+parts[1].split('[')[0]+'\n')
    if 'zona pellucida' in parts[1].lower() or 'zona pellucida' in parts[2].lower():
      idx=parts[0]
      g2.write(idx+'\t'+parts[1].split('[')[0]+'\n')
  f.close()
  g.close()
  g2.close()

if __name__=='__main__':
  main(sys.argv[1],sys.argv[2])