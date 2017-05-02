#!/usr/bin/env python
def main(filelist):

    import numpy as np
    print filelist
    #find unique samples
    samplelist=[]
    domains={}
    for filename in filelist:
        sample=filename.split('/')[-1].split('.')[0]
        if sample not in samplelist:
           samplelist.append(sample)
    print samplelist
    for filename in filelist:
        f=open(filename)
        sample=filename.split('/')[-1].split('.')[0]
        for line in f:
            line=line.rstrip()
            tigrfam=line.split('\t')[0].split(' ')[0]
            count=int(line.split('\t')[2])
            if tigrfam not in domains:
               domains[tigrfam]={}
               for i in samplelist:
                   domains[tigrfam][i]=0
            domains[tigrfam][sample] += count
        f.close()
    g=open('abundance_matrix3.txt','w')
    g.write('domain')
    for sample in samplelist:
        g.write('\t'+sample)
    g.write('\n')
    for domain in domains.keys():
        if np.sum(domains[domain].values())>0:
          g.write(domain)
          for sample in samplelist:
            g.write('\t'+str(domains[domain][sample]))
          g.write('\n')
    g.close()
    
                                  

if __name__=='__main__':
   import sys
   main(sys.argv[1:])