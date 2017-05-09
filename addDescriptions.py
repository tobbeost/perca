#!/usr/bin/env python 
import sys
def main(GSEAfile,allGoterms):
  f=open(allGoterms)
  golist={}
  for line in f:
    line=line.rstrip()
    line=line.split('\t')
    goid=line[0]
    ontology=line[1]
    goterm=line[2]
    if goid not in golist:
      golist[goid]={}
      golist[goid]['ont']=ontology
      golist[goid]['term']=goterm
    
  f.close()
  f=open(GSEAfile)
  newfilename=GSEAfile.rstrip('.txt')+'4.txt'
  g=open(newfilename,'w')
  descr=f.readline()
  descr=descr.split('\t')
  newdescr=[descr[0],'Ontology','Description']+descr[1:]
  g.write('\t'.join(newdescr)+'\n')
  for line in f:
    line=line.rstrip()
    line=line.split('\t')
    goid=line[0]
    if goid in golist:
      ont=golist[goid]['ont']
      term=golist[goid]['term']
    else:
      ont='NA'
      term='NA'
    newline=[line[0],ont,term]+line[1:]
    g.write('\t'.join(newline)+'\n')
  f.close()
  g.close()
  
if __name__=='__main__':
  main(sys.argv[1],sys.argv[2])