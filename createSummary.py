#!/usr/bin/env python
def readAnnot(transcriptsfile,annotationfile):
  alltranscripts={}
  f=open(transcriptsfile)
  line=f.readline()
  for line in f:
    line=line.rstrip()
    line=line.split('\t')
    idx=line[0]
    if idx not in alltranscripts:
      alltranscripts[idx]={}
      alltranscripts[idx]['annot1']='NA'
      alltranscripts[idx]['annot2']='NA'
      alltranscripts[idx]['annot3']='NA'
      alltranscripts[idx]['gobp']=[]
      alltranscripts[idx]['gomf']=[]
      alltranscripts[idx]['gocc']=[]
      
  f.close()
  f=open(annotationfile)
  line=f.readline()
  line=line.rstrip()
  header=line.split('\t')
  transcriptindex=header.index('TranscriptName')
  a1=header.index('DescriptionSP')
  a1id=header.index('HSPNameSP')
  a1org=header.index('OSNameSP')
  a2=header.index('DescriptionUf')
  a2org=header.index('Taxonomy')
  gobp=header.index('BPId')
  gomf=header.index('MFId')
  gocc=header.index('CCId')
  pfam=header.index('CDDesc')
  
  for line in f:
    if not line.startswith('#'):
      line=line.rstrip()
      parts=line.split('\t')
      idx=parts[transcriptindex]
      if parts[a1] not in '-':
        alltranscripts[idx]['annot1']=parts[a1]+' ['+parts[a1id]+', '+parts[a1org]+']'
      if parts[a2] not in '-':
        alltranscripts[idx]['annot2']=parts[a2]+' ['+parts[a2org]+']'
      if parts[pfam] not in '-':
        alltranscripts[idx]['annot3']=parts[pfam]
      go=parts[gobp]
      if go not in '-':
        gos=go.split(']---[')
        alltranscripts[idx]['gobp']=gos
      go=parts[gomf]
      if go not in '-':
        gos=go.split(']---[')
        alltranscripts[idx]['gomf']=gos
      go=parts[gocc]
      if go not in '-':
        gos=go.split(']---[')
        alltranscripts[idx]['gocc']=gos
      
        
  f.close()
  return(alltranscripts)
  
def main(inputfile,outputfile):
  abundancematrix='wholetranscript/west_coast/abundance_matrix_westcoast.txt'
  annotationfile='/storage/tobiaso/perca/02_annotation/Perca_fluviatilis_west-coast/output/west_coast_good-transcripts_CD-HIT_97_TWO_SAMPLES_uniref_2016_08_filt_ann_out.txt'
  annot=readAnnot(abundancematrix,annotationfile)
  print annot['P3503_115_TRINITY_DN20600_c0_g2_i1']
  f=open(inputfile)
  g=open(outputfile,'w')
  for line in f:
    line=line.rstrip()
    parts=line.split('\t')
    idx=
    
if __name__=='__main__':
  import sys
  inputfile=sys.argv[1]
  outputfile=sys.argv[2]
  main(inputfile,outputfile)