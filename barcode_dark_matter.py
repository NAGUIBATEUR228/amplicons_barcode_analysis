import os
import pandas as pd 
import numpy as np
import re
from datetime import datetime

u1 = "GATGTCCACGAGGTCTCT"
u2 = "CGTACGCTGCAGGTCGAC"

def rc(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = ""
    for nt in seq:
        if nt in complement:
            reverse_complement = complement[nt] + reverse_complement
    return reverse_complement

u={'u1':u1,'u2':u2,'u1_rc':rc(u1),'u2_rc':rc(u2)}

print('directory with script: C:\\Users\\zokmi\\Desktop\\study\\coursework\\')
directory=input('enter directory with data: ')
path=f'C:\\Users\\zokmi\\Desktop\\study\\coursework\\{directory}\\'
print(str(datetime.now())+' '+path)
ref=pd.read_csv(f'{path}..\\reference.txt')

dmp=path+'..\\dark_matter'
if not os.path.exists(dmp):os.mkdir(dmp)
if not os.path.exists(dmp+'\\out'):os.mkdir(dmp+'\\out')
f=open(f'{dmp}\\query_u.fasta',"w")
for i,j in u.items():
    f.write(f'>{i}\n{j}\n')
f.close()

os.system(f'makeblastdb -in dark_matter/query_u.fasta -dbtype nucl -out dark_matter/uref')

dirs=list()
for i in os.listdir(path):#list of directory and file names
    if os.path.isdir(os.path.join(path, i)):
        dirs.append(i)
print(str(datetime.now())+' '+str(dirs))#directories with fastq files

dirs=['1Pet_L001-ds.c61c642938334ad581ddd020b0d2602c']

exp=[]
totalreadcount=[]
readcountnb=[]
readcountb=[]
uniqueb=[]
matched=[]
notmatched=[]
#same in percent
for i in dirs:
    p=path+i+'\\artem'
    content = os.listdir(p)
    files = []
    for file in content:
        if os.path.isfile(os.path.join(p, file)) and file.endswith('not_barcoded_raw_qual_count.csv'):
            files.append(file)
    #parsing not_matched files in folders
    for j in files:
        name=re.split(r'_not_barcoded_raw_qual_count\.csv',j)[0]
        exp.append(name)
        nb=pd.read_csv(f'{p}\\{j}')#not_barcoded table
        f=open(f'{dmp}\\{name}.fasta',"w")#making fasta with BLAST queries
        for k in nb['seq']:
            f.write(f'>{k}-unknownseqname\n{k}\n')
        f.close()
        print(str(datetime.now())+' '+f'{dmp}\\{name}.fasta')
        
        #BLAST command, searching not_matched sequences in database ref. output file end with blastout.txt. e-value 10, alignment initiating word size 6, search on the same strand, finds only the best hit. outfile contains different information about alignment in tsv format.
        #os.system(f'blastn -query dark_matter/{name}.fasta -db dark_matter/uref -out dark_matter/out/{name}_blastout.txt -evalue 0.001 -word_size 6 -strand plus')
        
        os.system(f'blastn -query dark_matter/{name}.fasta -db dark_matter/uref -out dark_matter/out/{name}_blastout.txt -evalue 0.001 -word_size 6 -strand plus -max_target_seqs 2 -outfmt \"6 qacc qlen sacc slen length nident evalue qstart qend sstart send\"')

        #subtracting -unknownseqname in outfile.
        with open(f'{dmp}\\out\\{name}_blastout.txt','r') as f:
            old_data=f.read()
        new_data = 'qacc\tqlen\tsacc\tslen\tlength\tnident\tevalue\tqstart\tqend\tsstart\tsend\n'+old_data.replace('-unknownseqname\t','\t')
        with open(f'{dmp}\\out\\{name}_blastout.txt','w') as f:
            f.write(new_data)

        #tables:output_dm_count (matched barcodes),not_matched_dm (nm barcodes)
        #join three tables: blast_barcodes, reference, not_barcoded. make dict with orfs; filter where nan in orf

        check=pd.read_csv(f'{dmp}\\out\\{name}_blastout.txt',sep='\t')
        #chooses the best alignment of a subject if two by max length
        check=check[(check.groupby(['qacc','sacc'])['length'].transform('max'))==check['length']]
        #filtering meaningful data and extracting barcodes
        def primers(col):
          u=col.unique()
          c1=len(u)==2
          c2=('u1' in u) and ('u2' in u)
          c3=('u1_rc' in u) and ('u2_rc' in u)
          return c1 and (c2 or c3)
        check=check[check.groupby(['qacc'])['sacc'].transform(primers)]
        check['ind']=np.where(check['sacc']=='u1',check['qend'],
                              np.where(check['sacc']=='u2',check['qstart']-1,
                                       np.where(check['sacc']=='u2_rc',check['qend'],
                                                np.where(check['sacc']=='u1_rc',check['qstart']-1,
                                                    None))))
        u=check['sacc'].unique().tolist()
        w=pd.pivot(data=check, values='ind', index='qacc',columns='sacc').reset_index()
        v=np.vectorize(lambda x,y,z:x[y:z])
        vecrc=np.vectorize(rc)
        w['barcode']=w['qacc']
        if 'u1' in u: w['barcode']=v(w['qacc'],w['u1'],w['u2'])
        if 'u1_rc' in u: w['barcode']=np.where(w['barcode']==w['qacc'],vecrc(v(w['qacc'],w['u2_rc'],w['u1_rc'])),w['barcode'])
        #adding barcodes to table with not_barcoded reads data and adding information on barcodes from reference table
        full=pd.merge(w,nb,left_on='qacc',right_on='seq',how='outer')
        full=pd.merge(full,ref,left_on='barcode',right_on='UPTAG_seqs',how='left')
        #orfs={}-??
        #либо пройти словариками либо эту жуткую аггрегацію использовать и если чо вернуться къ фулл

print(str(datetime.now())+' '+'DONE')