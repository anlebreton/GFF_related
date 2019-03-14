#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
#author          :Annie Lebreton
#date            :20190314
#version         :1
#python_version  :3.6
"""

import argparse

def parseoptions( ):
    """ Docstring 
    .... """

    parser = argparse.ArgumentParser( description="Initialy done for adapt JGI gff annotation without gene and mRNA to genometools requierements " )
    parser.add_argument( '-g',  '--gtf',default="test.gff",  help="[input gff without gene nor mRNA]")
    
    global ARGS         # Update the global ARGS variable 
    ARGS = parser.parse_args()

##----- INPUT ---------
# scaffold_1      JGI     exon    446     1141    .       +       .       name "CE1_279"; transcriptId 137
# scaffold_1      JGI     CDS     642     1073    .       +       0       name "CE1_279"; proteinId 2; exonNumber 1
# scaffold_1      JGI     start_codon     642     644     .       +       0       name "CE1_279"
# scaffold_1      JGI     stop_codon      1071    1073    .       +       0       name "CE1_279"

##----- OUTPUT --------
# scaffold_1      JGI     gene    446     1141    .       +       .       ID=CE1_279, transcriptId=137
# scaffold_1      JGI     mRNA    446     1141    .       +       .       Parent=CE1_279, ID=CE1_279.m; transcriptId=137
# scaffold_1      JGI     exon    446     1141    .       +       .       Parent=CE1_279.m, transcriptId=137
# scaffold_1      JGI     CDS     642     1073    .       +       0       Parent=CE1_279.m, proteinId=2; exonNumber=1
# scaffold_1      JGI     start_codon     642     644     .       +       0       Parent=CE1_279.m
# scaffold_1      JGI     stop_codon      1071    1073    .       +       0       Parent=CE1_279.m



def getTxt(lS,txt,sep):
    txt=str(txt)
    out="no_"+str(txt)
    l=lS
#    print txt
    try:
        d=l.find(txt)
        f=l.find(sep, (d+len(txt)))
        if f== -1:
            out=str(l[d+len(txt):])
        else:
            out=str(l[d+len(txt):f])

        return out
    except:
       return out

def getIdent(lS):
    tmp=getTxt(lS[8],"name",";")
    tmp=tmp.replace('"','')
    ident=tmp.strip()
    return ident

def createDicoID():
    dicoID={}
    with open(ARGS.gtf,'r') as ef :
        lines=ef.readlines()

        for line in lines:
            if line.startswith('##FASTA'): break
            if line.startswith('#'): continue
            lS=line.split("\t")
            ident=getIdent(lS)
            #print(ident)
            if ident not in dicoID:
                dicoID[ident]=[min(int(lS[3]),int(lS[4])), max(int(lS[3]),int(lS[4]))]
            else:
                start=min(dicoID[ident][0],int(lS[3]),int(lS[4])) 
                stop =max(dicoID[ident][1],int(lS[3]),int(lS[4]))
                dicoID[ident]=[start,stop]
    return dicoID




def main():
    parseoptions( )    
    dicoID=createDicoID()
    ## dicoID[name/ID]=(start,stop)
    print("##gff-version 3")
    with open(ARGS.gtf,'r') as ef :
        lines=ef.readlines()
        previousIdent=""
        ##keep comments and FASTA at the end of the gff file            
        fasta=False
        for line in lines:
            line=line.strip()
            if line.startswith('##FASTA') or fasta==True:
                print(line)
                fasta=True
                continue
            if line.startswith('#'):
                print(line)
                continue

            ##process other lines
            lS=line.split('\t')
            ident=getIdent(lS)
            if ident != previousIdent:
                previousIdent = ident
                m="ID="+ident
                m2="Parent=" + ident + ";ID=" + ident + ".m"
                print("\t".join([lS[0], lS[1],"gene", str(dicoID[ident][0]), str(dicoID[ident][1]), lS[5],lS[6],".", m]))
                print("\t".join([lS[0], lS[1],"mRNA", str(dicoID[ident][0]), str(dicoID[ident][1]), lS[5],lS[6],".", m2]))

            line=line.replace('"','')
            tmp="name "+ident
            tmp2="Parent=" + ident + ".m"
            line=line.replace(tmp, tmp2)
            line=line.replace("; ",";") 	
            line=line.replace(' ','=')
            print(line)


main()        
