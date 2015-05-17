#!/usr/bin/python

def about():
    from time import strftime
    __time__   = strftime("%Y-%m-%d %H:%M:%S")
    __author__ = "By taushif @ JNU : taushifkhan@gmail.com"
    __path__   = "Created From /home/taushif/PyCodes/ProLegoMods/pdbParse.py"
    aboutString = "#%s\n#%s\n#On %s\n"%(__path__,__author__,__time__)
    return aboutString

import re

def chainInfo(pdbfl):
    '''
    pdbfl : opened pdb file in .realines()
    '''
    for k in pdbfl:
        if re.match(r'^COMPND',k) and re.search(r'CHAIN',k):
            chain = k.split(":")[1].rstrip(";").split(",")
            return chain
        else:
            continue

class chainParse:
    def __init__(self):
        self.pdbid = ''
        self.chain = ''
        self.atomLines = []
        self.resDic = {}
        
    def parsePDB(self,pdbFile,chain):
        for i in pdbFile:
            i = i.rstrip()
            if re.match(r'^HEADER',i):
                self.pdbid = i[62:]
            if re.match(r'^ATOM',i) and i[21] == chain and re.match(r'[N,O,C,S]',i[77]):
                self.atomLines.append(i)
                self.resDic[int(i[22:26].strip())] = i[17:20].strip()
        
    def getatom(self,resId):
        atomCoord = {}
        for i in self.atomLines:
            if int(i[22:26].strip()) == resId:
                atomCoord[i[12:16].strip()] = (float(i[30:38]),float(i[38:46]),float(i[46:54]))
        return atomCoord

    def writePDB(self,outfl,pdb_detail):
        wfl = open(outfl,'w')
        authstring = about()
        wfl.write("%s%s\n"%(pdb_detail,authstring)) # writes header
        # writing chain atom lines from pdb

        for l in self.atomLines:
            wfl.write("%s\n"%l)

        wfl.write("END\n")
        wfl.close()
        return 1


if __name__ == '__main__':
    import os
    import argparse
    parse = argparse.ArgumentParser()
    parse.add_argument('-p',action='store',dest='pdbfl',\
        help='pdbid with path')
    parse.add_argument('-c',action='store',dest='chain',\
        help='chain of the pdb file')
    presult = parse.parse_args()
    
    path  = os.getcwd().rstrip("/")

    try:
        pdbfl = open(presult.pdbfl,'r').readlines()
    except:
        pdbfl = path+"/"+presult.pdbfl
        pdbfl = open(pdbfl,'r').readlines()

    pdbid = presult.pdbfl[-8:-4]
    chain = presult.chain

    outfl = path+"/"+pdbid+"_"+chain+".pdb"

    print pdbid,outfl

    q = chainParse()
    q.parsePDB(pdbfl,chain)
    pdbresidues = q.resDic

    print pdbresidues
    pdb_detail = "#pdbId: %s; Chain: %s; Residues %d\n"%(pdbid,chain,len(pdbresidues.keys()))
    # atC = q.getatom(30)
    # print atC
    wstat = q.writePDB(outfl,pdb_detail)
    print "Write status %d"%wstat
