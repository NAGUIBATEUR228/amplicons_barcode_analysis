import os
import numpy as np
import pandas as pd 
import re
from datetime import datetime

print('directory with script: C:\\Users\\zokmi\\Desktop\\study\\coursework\\')
directory=input('enter directory with data: ')
path=f'C:\\Users\\zokmi\\Desktop\\study\\coursework\\{directory}\\'
print(str(datetime.now())+' '+path)
#import reference barcode table, made by R scirpt
ref=pd.read_csv(f'{path}..\\reference.txt')
print(ref.head(10))
dirs=list()
for i in os.listdir(path):#list of directory and file names
    if os.path.isdir(os.path.join(path, i)):
        dirs.append(i)
print(str(datetime.now())+' '+str(dirs))#directories with fastq files
#making fasta file for BLAST database
f=open('./ref_db.fasta',"w")
for i in range(len(ref['UPTAG_seqs'])):
    i1=ref['UPTAG_seqs'][i]
    i2=ref['Confirmed_deletion'][i]
    f.write(f'>{i2}\n{i1}\n')
f.close()
#BLAST command, making nucleotide database 'ref'
os.system('makeblastdb -in ./ref_db.fasta -dbtype nucl -out ref')

#lists for summary table
distar=[]#list of number of reads, which are too different from any reference barcode
distcount=[]#their percentage
allcount=[]#total read count
exp=[]#name of an experiment
qual=[]
for i in dirs:
    p=path+i+'\\artem'
    content = os.listdir(p)
    files = []
    for file in content:
        if os.path.isfile(os.path.join(p, file)) and file.endswith('not_matched_dm.csv'):
            files.append(file)
    #parsing not_matched files in folders
    for j in files:
        name=re.split(r'\.csv',j)[0]
        nb=pd.read_csv(f'{p}\\{j}')#not_matched table
        f=open(f'{p}\\{name}.fasta',"w")#making fasta with BLAST queries
        for k in nb['barcode']:
            f.write(f'>{k}-unknownseqname\n{k}\n')
        f.close()
        os.system(f'cd {path}')
        print(str(datetime.now())+' '+f'{p}\\{name}.fasta')
        #BLAST command, searching not_matched sequences in database ref. output file end with blastout.txt. e-value 10, alignment initiating word size 6, search on the same strand, finds only the best hit. outfile contains different information about alignment in tsv format.
        os.system(f'blastn -query {p}\\{name}.fasta -db ref -out {p}\\{name}_blastout.txt -evalue 10 -word_size 6 -strand plus -max_target_seqs 1 -outfmt \"6 qacc qlen sacc slen length nident evalue qstart qend sstart send\"')
        #subtracting -unknownseqname in outfile.
        with open(f'{p}\\{name}_blastout.txt','r') as f:
            old_data=f.read()
        new_data = 'qacc\tqlen\tsacc\tslen\tlength\tnident\tevalue\tqstart\tqend\tsstart\tsend\n'+old_data.replace('-unknownseqname\t','\t')
        with open(f'{p}\\{name}_blastout.txt','w') as f:
            f.write(new_data)
        #debug for checking if there are exact matches between not_matched and reference barcodes.
        check=pd.read_csv(f'{p}\\{name}_blastout.txt',sep='\t')
        check['norm']=(check['length']>check['nident']) | (check['qlen']!=check['length']) | (check['slen']!=check['length']) | ((check['qstart']-check['qend'])*(check['sstart']-check['send'])>=0)
        abnorm=check[check['norm']==False]
        if not abnorm.empty:
            print(str(datetime.now())) 
            print(abnorm)
        else:
            print(str(datetime.now())+'\n'+'OK')
        
        #filtering distance between sequences
        check['dist']=(check['slen']+check['qlen']-check['length']-check['nident'])
        check=check[check['dist']<=3]

        name=re.split(r'_not_matched_dm\.csv',j)[0]
        exp.append(name)
        print(str(datetime.now())+' '+f'{p}\\{name}_blasted_dm.csv')
        f=open(f'{p}\\{name}_blasted_dm.csv','w')
        f.write('Confirmed_deletion,barcode,n,count,notes,original_barcode')

        #merging blastout,not_matched and reference table in order to obtain necessary information
        full=pd.merge(check,nb,left_on='qacc',right_on='barcode',how='outer')
        full=pd.merge(full,ref,left_on='sacc',right_on='Confirmed_deletion',how='left')
        orfs={}
        toofar=0
        #making hash table with orfs for final result file
        for elem in full.iterrows():
          if pd.isna(elem[1]['qacc']): toofar+=elem[1]['count']; continue
          k=elem[1]["sacc"]#orf
          if k not in orfs.keys():
            orfs[k]=["",0,0,"",elem[1]['UPTAG_seqs']]
          bc=elem[1]["qacc"]#barcode
          ev="("+str(elem[1]["evalue"])+")"
          orfs[k][1]+=elem[1]["n"]
          orfs[k][2]+=elem[1]["count"]
          orfs[k][0]+=bc+"|"
          orfs[k][3]+=ev+"-"

        distcount.append(toofar)
        allcount.append(nb['count'].sum().item())
        qual.append(round(nb['qual'].mean().item(),3))
        distar.append(round(distcount[-1]/allcount[-1],3))
        for elem in orfs.keys():
          f.write(f"\n{elem},{orfs[elem][0][:-1]},{orfs[elem][1]},{orfs[elem][2]},{orfs[elem][3][:-1]},{orfs[elem][4]}")
        f.close()
print(str(datetime.now())+' '+'DONE')
print(str(datetime.now())+' '+str(distar))

f=open(f'{path}sumblastdm.csv','w')
f.write('exp,all_count,toofar_count,percentage,qual\n')
for i in range(len(exp)):
    f.write(f'{exp[i]},{allcount[i]},{distcount[i]},{distar[i]},{qual[i]}\n')
f.close()