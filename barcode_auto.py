import os
import numpy as np
import pandas as pd 
import re
from datetime import datetime

#script inputs are reference.txt with reference barcode table and directory with folders, denoting experiments and containing tables with unmatched barcodes.
print('directory with script: C:\\Users\\zokmi\\Desktop\\study\\coursework\\')
directory=input('enter directory with data: ')
path=f'C:\\Users\\zokmi\\Desktop\\study\\coursework\\{directory}\\'
print(str(datetime.now())+' '+path)
ref=pd.read_csv(f'{path}..\\reference.txt')
txt=list(ref['UPTAG_seqs'])

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

print(str(datetime.now())+' building trie')

#node of trie with transition dictionary and row number of barcode, that ends in this node, in ref table
class Node:
    def __init__(self,end=None,num=None):
        self.tr=dict()#потомки по буквамъ
        self.end=end

#adding barcodes from reference table
root=Node()
for j in range(len(txt)):
    v=root
    pat=txt[j]
    for i in range(len(pat)):
        if not (pat[i] in v.tr.keys()):
            u=Node()
            v.tr[pat[i]]=u
        v=v.tr[pat[i]]
    v.end=j

#getting list of directories
dirs=list()
for i in os.listdir(path):
    if os.path.isdir(os.path.join(path, i)):
        dirs.append(i)
print(str(datetime.now())+' '+str(dirs))

#these are for summary table in the end
distar=[]
distcount=[]
allcount=[]
exp=[]
qual=[]
for i in dirs:
    p=path+i+'\\artem'
    content = os.listdir(p)
    files = []
    for file in content:
        if os.path.isfile(os.path.join(p, file)) and file.endswith('not_matched.csv'):
            files.append(file)
    for j in files:
        nb=pd.read_csv(f'{p}\\{j}')#table of not matched barcodes
        name=re.split(r'_not_matched\.csv',j)[0]
        exp.append(name)
        print(str(datetime.now())+' '+f'{p}\\{name}_aligned.csv')
        f=open(f'{p}\\{name}_aligned.csv','w')#res table
        f.write('Confirmed_deletion,barcode,n,count,notes,original_barcode')
        
        toofar=0
        orfs={}
        for elem in nb.iterrows():
            bcd=elem[1]['barcode']
            if not isinstance(bcd, str):toofar+=elem[1]['count']; continue
            mass=go(root,bcd)#[[row,dist],...]
            if len(mass)==0: toofar+=elem[1]['count']; continue
            subm=get_min(mass)#getting the best hits
            dist=subm[0]
            mass=subm[1]
            l=len(mass)#to distribute counts between hits
            for sm in mass:
                k=ref.iloc[sm]['Confirmed_deletion']#orf
                if k not in orfs.keys(): orfs[k]=["",0,0,str(dist)+'|'+ref.iloc[sm]['UPTAG_notes'],ref.iloc[sm]['UPTAG_seqs']]
                orfs[k][1]+=elem[1]['n']/l
                orfs[k][2]+=elem[1]['count']/l
                orfs[k][0]+=bcd+"|"
                #orfs[k][4]
        
        #writing tables
        distcount.append(toofar)
        allcount.append(nb['count'].sum().item())
        qual.append(round(nb['qual'].mean().item(),3))
        distar.append(round(distcount[-1]/allcount[-1],3))
        for elem in orfs.keys():
          f.write(f"\n{elem},{orfs[elem][0][:-1]},{orfs[elem][1]},{orfs[elem][2]},{orfs[elem][3]},{orfs[elem][4]}")
        f.close()

print(str(datetime.now())+' '+str(distar)+' DONE')

f=open(f'{path}sumalign.csv','w')
f.write('exp,all_count,toofar_count,percentage,qual\n')
for i in range(len(exp)):
    f.write(f'{exp[i]},{allcount[i]},{distcount[i]},{distar[i]},{qual[i]}\n')
f.close()