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
  outputfilename='/storage/tobiaso/perca/'+transcriptfile.replace('abundance_matrix','annotation_transcripts2')
  #gobpfilename=transcriptfile.replace('abundance_matrix','go_bp_')
  #gomffilename=transcriptfile.replace('abundance_matrix','go_mf_')
  #goccfilename=transcriptfile.replace('abundance_matrix','go_cc_')
  g1=open(outputfilename,'w')
  #g2=open(gobpfilename,'w')
  #g3=open(gomffilename,'w')
  #g4=open(goccfilename,'w')
  g1.write('transcript_id\tAnnotation_uniprot\tAnnotation_fish\PFAM_COG_CDD\n')
  numtot=len(alltranscripts)
  numannot=0
  numannotpfam=0
  for transcript in alltranscripts:
    g1.write(transcript+'\t'+alltranscripts[transcript]['annot1']+'\t'+alltranscripts[transcript]['annot2']+'\t'+alltranscripts[transcript]['annot3']+'\n')
    if not (alltranscripts[transcript]['annot1']=='NA' and alltranscripts[transcript]['annot2']=='NA'):
      numannot+=1
    if not (alltranscripts[transcript]['annot1']=='NA' and alltranscripts[transcript]['annot2']=='NA' and alltranscripts[transcript]['annot3']=='NA'):
      numannotpfam+=1
  
#    if len(alltranscripts[transcript]['gobp'])>0:
#      for go in alltranscripts[transcript]['gobp']:
#        g2.write(transcript+'\t'+go+'\n')
    #if len(alltranscripts[transcript]['gomf'])>0:
    #  for go in alltranscripts[transcript]['gomf']:
    #    g3.write(transcript+'\t'+go+'\n')
    #if len(alltranscripts[transcript]['gocc'])>0:
    #  for go in alltranscripts[transcript]['gocc']:
    #    g4.write(transcript+'\t'+go+'\n')
  g1.close()
  #g2.close()
  #g3.close()
  #g4.close()
  print "total number of transcripts: ",numtot
  print "annotated transcripts: ",numannot
  print "annotated transcripts including pfam annotation: ",numannotpfam

if __name__=='__main__':
  abundancematrix='abundance_matrix_eastcoast.txt'
  annotationfile='/storage/tobiaso/perca/02_annotation/Perca_fluviatilis_east-coast/output/east_coast_good-transcripts_CD-HIT_97_TWO_SAMPLES_uniref_2016_08_filt_ann_out.txt'
  #numberofannotations=1
  main(abundancematrix,annotationfile)