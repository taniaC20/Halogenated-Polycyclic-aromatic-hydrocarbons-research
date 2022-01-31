import numpy as np
import pandas as pd
from pandas import DataFrame
import os
import json
import argparse

parser = argparse.ArgumentParser(description='Reading index of acbap.')
parser.add_argument("-i",'--index',type=int, help='acbap index')
parser.add_argument("-p",'--practice',type=bool, help='Is Practice?', default=False)
args = parser.parse_args()

print(f"Command-line arguments: {args}")

aromatic_index = args.index
is_practice = args.practice

if aromatic_index == None:
  raise Exception("aromatic INDEX IS NONE...")

print(f"aromatic INDEX = {aromatic_index}")


root_path         = "/home/tnchevez/aromatic_batch"
if is_practice:
    root_path         = "."
output_folder     = f"{root_path}/output_2"
aromatic_filepath = f'{root_path}/Aromatic_GeneralStructures.csv'
output_filepath   = f'{output_folder}/Aromatic-output_{aromatic_index}.csv'

try:
  if not os.path.exists(output_folder):
      os.makedirs(output_folder)
except Exception as e:
  print("ERROR CREATING FOLDER: ", e)
  
def readData(aromatic_filepath):
  aromatic_df   = pd.read_csv(aromatic_filepath)
  return aromatic_df

def split(word):
        return [char for char in word]
    
def permutation(nCl=0,nF=0,nBr=0,numH=0):
    elementsList=[]
    numHalo=nCl+nF+nBr
    for i in range(0,numH-numHalo): elementsList.append('H')
    if nCl > 0:
        for x in range(0,nCl): elementsList.append('C')
    if nF > 0:
        for x in range(0,nF): elementsList.append('F')
    if nBr > 0:
        for x in range(0,nBr): elementsList.append('B')
    perms = [[]]
    for n in elementsList:
        new_perm = []
        for perm in perms:
            for i in range(len(perm) + 1):
                new_perm.append(perm[:i] + [n] + perm[i:])
                # handle duplication
                if i < len(perm) and perm[i] == n: break
        perms = new_perm
    return perms
def duplicates(seq,item):
    start_at = -1
    locs = []
    while True:
        try:
            loc = seq.index(item,start_at+1)
        except ValueError:
            break
        else:
            locs.append(loc)
            start_at = loc
    return locs

def replaceSmileFromPermutation(positions_1=[],positions_2=[],positions_3=[],maxH=0, halo_1="",halo_2="", halo_3="",smile="",cont=0):
    cont=0
    smile=split(smile)
    for index in range(0,len(smile)):
            if smile[index] != 'H':
                continue
            else:
                cont+=1
                # print("contH",cont)
                
                if positions_1:
                    if (cont-1) in positions_1:
                        smile[index]=halo_1
                        # print("smile in function",smile)
                if positions_2:
                    if (cont-1) in positions_2:
                        smile[index]=halo_2
                        # print("smile in function",smile)
                if positions_3:
                    if (cont-1) in positions_3:
                        smile[index]=halo_3
                        # print("smile in function",smile)
                
                    
    smile="".join(smile)
    return smile

def generateGroupOfOne(group="",smile="",numH=0):
    allPossibleSmiles=[]
    elementToLook=''
    halo=""

    for h in range(1,numH+1):
        # print(h)
        if group=="Cl":
            perms=permutation(nCl=h,numH=numH)
            elementToLook='C'
            halo="Cl"
        if group=="F":
            perms=permutation(nF=h,numH=numH)
            elementToLook='F'
            halo="F"
        if group=="Br":
            perms=permutation(nBr=h,numH=numH)
            elementToLook='B'
            halo="Br"
        
        for p in perms:
                # print("ENTERED PERMS")
                elementInPerms=''.join(p)
                # print("elementInPerms",elementInPerms)
                positions=duplicates(elementInPerms,elementToLook)
                # print("positions",positions)
                
                smile_Final=replaceSmileFromPermutation(positions_1=positions,maxH=numH,halo_1=halo,smile=smile)
                # print("smile",smile_Final)
                allPossibleSmiles.append(smile_Final)
                
    return allPossibleSmiles

def generateGroupOfTwo(group="",smile="",numH=0):
    allPossibleSmiles=[]
    elementToLook_1=''
    elementToLook_2=''
    halo_1=""
    halo_2=""
    
    for h1 in range(1,numH+2):
        
        for h2 in range(1,(numH+1)-h1):
            # print("Halo 1, Halo 2",h1,h2)
            if group =="ClBr" or group =="BrCl":
                perms=permutation(nCl=h1,nBr=h2)
                elementToLook_1='C'
                elementToLook_2='B'
                halo_1="Cl" 
                halo_2="Br"
            if group =="ClF" or group =="FCl":
                perms=permutation(nCl=h1,nF=h2)
                elementToLook_1='C'
                elementToLook_2='F'
                halo_1="Cl" 
                halo_2="F"
            if group == "FBr" or group == "BrF":
                perms=permutation(nBr=h1,nF=h2)
                elementToLook_1='B'
                elementToLook_2='F'
                halo_1="Br" 
                halo_2="F"

            for p in perms:
                    # print("ENTERED PERMS")
                    elementInPerms=''.join(p)
                    # print("elementInPerms",elementInPerms)
                    positions_1=duplicates(elementInPerms,elementToLook_1)
                    positions_2=duplicates(elementInPerms,elementToLook_2)
                    
                    # print("positions 1, positions 2",positions_1,positions_2)

                    smile_Final=replaceSmileFromPermutation(positions_1=positions_1,positions_2=positions_2,maxH=numH,halo_1=halo_1,halo_2=halo_2,smile=smile)
                    # print("smile",smile_Final)
                    allPossibleSmiles.append(smile_Final)
                
    return allPossibleSmiles
    
def generateGroupOfThree(group="",smile="",numH=0):
    allPossibleSmiles=[]
    elementToLook_1=''
    elementToLook_2=''
    elementToLook_3=''
    halo_1=""
    halo_2=""
    halo_3=""
    
    for h1 in range(1,numH+2):
        
        for h2 in range(1,(numH+1)-h1):
            
            for h3 in range(1,(numH+1)-(h1+h2)):
            
                # print("Halo 1, Halo 2, Halo 3:",h1,h2,h3)
                perms=permutation(nCl=h1,nBr=h2,nF=h3)
                elementToLook_1='C'
                elementToLook_2='B'
                elementToLook_3='F'
                halo_1="Cl" 
                halo_2="Br"
                halo_3="F"

                for p in perms:
                        # print("ENTERED PERMS")
                        elementInPerms=''.join(p)
                        # print("elementInPerms",elementInPerms)
                        positions_1=duplicates(elementInPerms,elementToLook_1)
                        positions_2=duplicates(elementInPerms,elementToLook_2)
                        positions_3=duplicates(elementInPerms,elementToLook_3)

                        # print("positions 1, positions 2, positions 3:",positions_1,positions_2,positions_3)

                        smile_Final=replaceSmileFromPermutation(positions_1=positions_1,positions_2=positions_2,positions_3=positions_3,maxH=numH,halo_1=halo_1,halo_2=halo_2,halo_3=halo_3,smile=smile)
                        # print("smile",smile_Final)
                        allPossibleSmiles.append(smile_Final)
                
    return allPossibleSmiles
    
    
def generateAllPossible(group,smile,numH):
    allPossibleSmiles=[]
    
    if group =="Cl" :
        allPossibleSmiles=generateGroupOfOne(group,smile,numH)
    elif group =="Br":
        allPossibleSmiles=generateGroupOfOne(group,smile,numH)
    elif group =="F":
        allPossibleSmiles=generateGroupOfOne(group,smile,numH)
    elif group =="ClBr" or group =="BrCl":
        allPossibleSmiles=generateGroupOfTwo(group,smile,numH)
    elif group =="ClF" or group =="FCl":
        allPossibleSmiles=generateGroupOfTwo(group,smile,numH)
    elif group == "FBr" or group == "BrF":
        allPossibleSmiles=generateGroupOfTwo(group,smile,numH)
    elif group == "BrClF":
        allPossibleSmiles=generateGroupOfThree(group,smile,numH)
    else:
        print("Group not recognized")
        return

    return allPossibleSmiles 

def allPosgenerateAllPossible_AllGroupssibleSmiles():
    aromatic_df = readData(aromatic_filepath)
    aromatic_arr = list(aromatic_df.iterrows()) 
    
    if aromatic_index >= len(aromatic_arr):
        raise Exception(f"\n\nTHE INDEX PROVIDED ({aromatic_index}) IS OUT OF RANGE FOR A ARRAY OF SIZE: {len(aromatic_arr)} WITH INDEXES: 0->{len(aromatic_arr)-1}\n\n")
  
    _, aromatic = aromatic_arr[aromatic_index]
    
    numH=aromatic["H"]
    smile=aromatic["General Smiles"]
    
    all_Cl=generateAllPossible("Cl",smile,numH)
    all_Br=generateAllPossible("Br",smile,numH)
    all_F=generateAllPossible("F",smile,numH)
    all_FBr=generateAllPossible("BrF",smile,numH)
    all_BrCl=generateAllPossible("BrCl",smile,numH)
    all_ClF=generateAllPossible("ClF",smile,numH)
    all_BrClF=generateAllPossible("BrClF",smile,numH)
    
    all_Possible=all_Cl+all_Br+all_F+ all_FBr+all_BrCl+all_ClF+ all_BrClF
    
    results_df = DataFrame (all_Possible,columns=['Smile'])
    results_df.to_csv(output_filepath, index=False)
    
if __name__ == "__main__":
  allPosgenerateAllPossible_AllGroupssibleSmiles()