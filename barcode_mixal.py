import os
import numpy as np
import pandas as pd 
import re
from datetime import datetime
from sys import argv
if len(argv)>1:
    need_to_BLAST=bool(int(argv[1]))
else: 
    need_to_BLAST=True

print('directory with script: C:\\Users\\zokmi\\Desktop\\study\\coursework\\')
directory=input('enter directory with data: ')
path=f'C:\\Users\\zokmi\\Desktop\\study\\coursework\\{directory}\\'
print(str(datetime.now())+' '+path)
#import reference barcode table, made by R scirpt
ref=pd.read_csv(f'{path}..\\reference.txt')
print(str(datetime.now()))
print(ref.head(5))

#function to find best matches of barcode
def get_min(mass):
    res=[228,[]]
    for i in mass:
        if i[1]<res[0]:
            res[0]=i[1]
            res[1]=[i[0]]
            continue
        if i[1]==res[0]:
            res[1].append(i[0])
    return res

#function globally aligns a barcode with a trie and returns a list of pairs [row in ref table with matched ORF,distance]
def go(node,pat,pp=0,dist=0,maxdist=2):#node of trie, pattern, position in pattern, current distance, maximum possible distance.
    res=[]
    if dist>maxdist: return res
    if pp==len(pat):
        if node.end==None: return res 
        else: 
            res.append([node.end,dist])
            for k in node.tr.keys():
                res.extend(go(node.tr[k],pat,pp,dist+1))
            return res
    c=pat[pp]
    if c in node.tr.keys():
        res.extend(go(node.tr[c],pat,pp+1,dist))
        for k in node.tr.keys():
            if k != c: 
                res.extend(go(node.tr[k],pat,pp,dist+1))
                res.extend(go(node.tr[k],pat,pp+1,dist+1))
    else:
        res.extend(go(node,pat,pp+1,dist+1))
        for k in node.tr.keys():
            res.extend(go(node.tr[k],pat,pp,dist+1))
            res.extend(go(node.tr[k],pat,pp+1,dist+1))
    return res

dirs=list()
for i in os.listdir(path):#list of directory and file names
    if os.path.isdir(os.path.join(path, i)):
        dirs.append(i)
print(str(datetime.now())+' '+str(dirs))#directories with fastq files
if need_to_BLAST:    
    #making fasta file for BLAST database
    f=open('./ref_db.fasta',"w")
    for i in range(len(ref['UPTAG_seqs'])):
        i1=ref['UPTAG_seqs'][i]
        i2=i1+'-refseqname'
        f.write(f'>{i2}\n{i1}\n')
    f.close()
    #BLAST command, making nucleotide database 'ref'
    os.system('makeblastdb -in ./ref_db.fasta -dbtype nucl -out ref')

summixal = pd.DataFrame({
    'exp': [], 'all_count': [], 'm_count': [], 'percentage': [], 'qual': []
})
summixaldm = pd.DataFrame({
    'exp': [], 'all_count': [], 'm_count': [], 'percentage': [], 'qual': []
})

barcodes=pd.Series(dtype=object)
for i in dirs:
    p=path+i+'\\artem'
    content = os.listdir(p)
    files = []
    for file in content:
        if os.path.isfile(os.path.join(p, file)) and (file.endswith('not_matched.csv') or file.endswith('not_matched_dm.csv')):
            files.append(file)
    #parsing not_matched files in folders
    for j in files:
        nb=pd.read_csv(f'{p}\\{j}')#not_matched table
        #print(len(nb[~nb['barcode'].isin(barcodes)]))
        barcodes=pd.Series(pd.concat([barcodes,nb['barcode']],ignore_index=True).unique())

print(len(barcodes))
if need_to_BLAST:
    f=open(f'{path}\\all.fasta',"w")#making fasta with BLAST queries
    for k in barcodes:
        f.write(f'>{k}-unknownseqname\n{k}\n')
    f.close()
    os.system(f'cd {path}')
    print(str(datetime.now())+' '+f'{p}\\all.fasta')
    #BLAST command, searching not_matched sequences in database ref. output file end with blastout3.txt. e-value 10, alignment initiating word size 6, search on the same strand, finds only the best hit. outfile contains different information about alignment in tsv format.
    os.system(f'blastn -query {path}\\all.fasta -db ref -out {path}\\all_blastout.txt -evalue 10 -word_size 6 -strand plus -max_target_seqs 1 -max_hsps 1 -outfmt \"6 qacc qlen sacc slen length nident evalue qstart qend sstart send\"')#-num_threads 8
    #subtracting -unknownseqname in outfile.
    with open(f'{path}\\all_blastout.txt','r') as f:
        old_data=f.read()
    new_data = 'qacc\tqlen\tsacc\tslen\tlength\tnident\tevalue\tqstart\tqend\tsstart\tsend\n'+old_data.replace('-unknownseqname\t','\t')
    with open(f'{path}\\all_blastout.txt','w') as f:
        f.write(new_data)
    with open(f'{path}\\all_blastout.txt','r') as f:
        old_data=f.read()
    new_data = old_data.replace('-refseqname\t','\t')
    with open(f'{path}\\all_blastout.txt','w') as f:
        f.write(new_data)

check=pd.read_csv(f'{path}\\all_blastout.txt',sep='\t')

#filtering distance between sequences
check['dist']=(check['slen']+check['qlen']-check['length']-check['nident'])
check=check[check['dist']<=2]
print(str(datetime.now())+' BLASTed')
print(str(datetime.now())+' adding reference and not_matched data')

barcode=pd.DataFrame({'barcode':barcodes})
full=pd.merge(barcode,check,left_on='barcode',right_on='qacc',how='left')
full=pd.merge(full,ref,left_on='sacc',right_on='UPTAG_seqs',how='left')
full=full[~pd.isna(full['qacc'])].reset_index()[['barcode','Confirmed_deletion','UPTAG_seqs','UPTAG_notes']]

print(str(datetime.now())+' aligning')
toal=barcodes[~barcodes.isin(check['qacc'])]
al=pd.DataFrame({'barcode':[],'Confirmed_deletion':[],'UPTAG_seqs':[],'UPTAG_notes':[],'dist':[]})

print(str(datetime.now())+' building trie')

#node of trie with transition dictionary and row number of barcode, that ends in this node, in ref table
class Node:
    def __init__(self,end=None,num=None):
        self.tr=dict()#потомки по буквамъ
        self.end=end

#adding barcodes from reference table
root=Node()
txt=toal.tolist()
for j in range(len(txt)):
    v=root
    pat=txt[j]
    if not isinstance(pat, str):continue
    for i in range(len(pat)):
        if not (pat[i] in v.tr.keys()):
            u=Node()
            v.tr[pat[i]]=u
        v=v.tr[pat[i]]
    v.end=j

print('toal',len(toal))

for elem in ref.iterrows():
      bcd=elem[1]['UPTAG_seqs']
      if not isinstance(bcd, str):continue
      mass=go(root,bcd)#[[row,dist],...]
      if len(mass)==0: continue
      cd=elem[1]['Confirmed_deletion']
      seq=elem[1]['UPTAG_seqs']
      note=elem[1]['UPTAG_notes']
      for sm in mass:
        bcd=txt[sm[0]]
        al.loc[ len(al.index)] = [
            bcd,
            cd,
            seq,
            note,
            sm[1]
            ]
print(str(datetime.now())+' aligned')
print('al',len(al))
al=al.drop_duplicates()
#or uncomment this to get only the best hits
#al=al[al.groupby('barcode')['dist'].transform('min')==al['dist']]
print('al',len(al))
for i in dirs:
    p=path+i+'\\artem'
    content = os.listdir(p)
    files = []
    for file in content:
        if os.path.isfile(os.path.join(p, file)) and file.endswith('not_matched.csv'):
            files.append(file)
    #parsing not_matched files in folders
    for j in files:
        nb=pd.read_csv(f'{p}\\{j}')#not_matched table
        name=re.split(r'_not_matched\.csv',j)[0]
        f=pd.merge(full,nb,left_on='barcode',right_on='barcode',how='right')
        nm=f[pd.isna(f['Confirmed_deletion'])].reset_index()[['barcode','n','count','qual']]
        m=f[~pd.isna(f['Confirmed_deletion'])].reset_index()[['barcode',    'Confirmed_deletion',   'UPTAG_seqs', 'UPTAG_notes',    'n',    'count',    'qual']]
        nm=pd.merge(al,nm,left_on='barcode',right_on='barcode',how='right')
        nm=nm[~pd.isna(nm['Confirmed_deletion'])].reset_index()[['barcode', 'Confirmed_deletion',   'UPTAG_seqs', 'UPTAG_notes',    'n',    'count',    'qual']]
        maxnm=nm.groupby('barcode')['count'].transform('max').values 
        lennm=nm.groupby('barcode')['count'].transform('count').values
        nm['count']=maxnm/lennm
        m=pd.concat([m,nm],ignore_index=True)
        m = m.groupby('Confirmed_deletion').agg(
                barcode= ('barcode',lambda x: '|'.join(x[~pd.isna(x)].unique())),
                n=('n', 'sum'),
                count=('count', 'sum'),
                notes=('UPTAG_notes', lambda x: '|'.join(x[~pd.isna(x)].unique())),
                original_barcode=('UPTAG_seqs',lambda x: '|'.join(x[~pd.isna(x)].unique()))).reset_index()

        m=m.sort_values (by = ['count'], ascending = [ False ])
        print(str(datetime.now())+' '+f'{p}\\{name}_mixaled.csv')
        m.to_csv (f'{p}\\{name}_mixaled.csv', index= False )

        summixal.loc[ len(summixal.index )] = [
            name,
            round(nb['count'].sum()),
            round(m['count'].sum()),
            round(m['count'].sum())/round(nb['count'].sum()),
            nb['qual'].mean().item()
            ]

print(str(datetime.now()))
print(str(summixal['percentage'].tolist())+' '+'DONE')
summixal.to_csv (f'{path}summixal.csv', index= False )


for i in dirs:
    p=path+i+'\\artem'
    content = os.listdir(p)
    files = []
    for file in content:
        if os.path.isfile(os.path.join(p, file)) and file.endswith('not_matched_dm.csv'):
            files.append(file)
    #parsing not_matched files in folders
    for j in files:
        nb=pd.read_csv(f'{p}\\{j}')#not_matched table
        name=re.split(r'_not_matched_dm\.csv',j)[0]
        #merging blastout3,not_matched and reference table in order to obtain necessary information
        f=pd.merge(full,nb,left_on='barcode',right_on='barcode',how='right')
        nm=f[pd.isna(f['Confirmed_deletion'])].reset_index()[['barcode','n','count','qual']]
        m=f[~pd.isna(f['Confirmed_deletion'])].reset_index()[['barcode',    'Confirmed_deletion',   'UPTAG_seqs', 'UPTAG_notes',    'n',    'count',    'qual']]
        nm=pd.merge(al,nm,left_on='barcode',right_on='barcode',how='right')
        nm=nm[~pd.isna(nm['Confirmed_deletion'])].reset_index()[['barcode', 'Confirmed_deletion',   'UPTAG_seqs', 'UPTAG_notes',    'n',    'count',    'qual']]
        maxnm=nm.groupby('barcode')['count'].transform('max').values 
        lennm=nm.groupby('barcode')['count'].transform('count').values
        nm['count']=maxnm/lennm
        m=pd.concat([m,nm],ignore_index=True)
        m = m.groupby('Confirmed_deletion').agg(
                barcode= ('barcode',lambda x: '|'.join(x[~pd.isna(x)].unique())),
                n=('n', 'sum'),
                count=('count', 'sum'),
                notes=('UPTAG_notes', lambda x: '|'.join(x[~pd.isna(x)].unique())),
                original_barcode=('UPTAG_seqs',lambda x: '|'.join(x[~pd.isna(x)].unique()))).reset_index()

        m=m.sort_values (by = ['count'], ascending = [ False ])
        print(str(datetime.now())+' '+f'{p}\\{name}_mixaled_dm.csv')
        m.to_csv (f'{p}\\{name}_mixaled_dm.csv', index= False )

        summixaldm.loc[ len(summixaldm.index )] = [
            name,
            round(nb['count'].sum()),
            round(m['count'].sum()),
            round(m['count'].sum())/round(nb['count'].sum()),
            nb['qual'].mean().item()
            ]

print(str(datetime.now()))
print(str(summixaldm['percentage'].tolist())+' '+'DONE')
summixaldm.to_csv (f'{path}summixaldm.csv', index= False )
