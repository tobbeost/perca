#!/usr/bin/env python
import sys
def main(transcriptfile,annotationfile):
  alltranscripts={}
  f=open(transcriptfile)
  line=f.readline()
  for line in f:
    line=line.rstrip()
    line=line.split('\t')
    idx=line[0]
    if idx not in alltranscripts:
      alltranscripts[idx]=[]
  f.close()
  f=open(annotationfile)
  for line in f:
    if not line.startswith('#'):
      line=line.rstrip()
      parts=line.split('\t')
      idx=parts[0]
      target=parts[8].lstrip('Target=').split()[0]
      if target not in alltranscripts[idx]:
        alltranscripts[idx].append(target)
  f.close()
  g=open('annotation_transcripts.txt','w')
  for transcript in alltranscripts:
    if len(alltranscripts[transcript])==0:
      g.write(transcript+'\tNA\n')
    else:
      g.write(transcript+'\t'+';'.join(alltranscripts[transcript])+'\n')
  g.close()
    
if __name__=='__main__':
  main(sys.argv[1],sys.argv[2])