#!/usr/bin/python
import re,sys
import numpy as Np

def calDist(atm1,atm2):
	return(Np.sum((atm1-atm2)**2))

def getAtomLine(l):
	atomDetail = {}
	atomDetail['atom_num']=int(l[5:11])
	atomDetail['atom_name'] = l[12:16].strip()
	atomDetail['res_name'] = l[17:20].strip()
	atomDetail['res_num'] = int(l[22:26])
	cord = (float(l[30:38]),float(l[38:46]),float(l[46:54]))
	atomDetail['cord'] = Np.array(cord)
	return atomDetail
	

def rco(pdbfl,chain):
    pdbFile = open(pdbfl,'r').readlines()
    print "processing %s"%pdbfl
    chain = chain
    atomList = []
    for i in pdbFile:
        i = i.rstrip()
        if re.match(r'^ATOM',i) and i[21] == chain:
            atomList.append(getAtomLine(i))
        if re.match(r'^TER',i):
            break
    #print atomList
    # Checking id contact or not
    length = atomList[len(atomList)-1]['res_num'] - atomList[0]['res_num']
    # Mapping pdb residue ids from 0 to length-1
    resDictionary = {}
    for k in atomList:
        resDictionary[k['res_num']] = len(resDictionary.keys()) -1
    

    contacts = 0
    sequenceDist = 0
    residueContact = Np.zeros((length+1,length+1))
    #import ipdb; ipdb.set_trace();

    for atom1 in atomList:
        for atom2 in atomList:
            seqDist = atom1['res_num']-atom2['res_num']
            if seqDist > 0:	# calculating non-adjecent residues
                dist = calDist(atom1['cord'],atom2['cord'])
                if dist < 36:	# 6*6 d comare with square dist
                    contacts +=1
                    sequenceDist += seqDist     
                    residueContact[resDictionary[atom1['res_num']],resDictionary[atom2['res_num']]]=1

    print "Analyzing contact"

    longContactRes = Np.count_nonzero(residueContact)
    #import ipdb; ipdb.set_trace();
    print longContactRes

    rco_fl = sequenceDist/float(length*contacts)
    absCO_fl = rco_fl*length
    lco_fl = longContactRes/float(150)

    #print contacts,sequenceDist, rco
    print length
    return rco_fl,absCO_fl,lco_fl
    

if __name__ == '__main__':
    import argparse
    parse = argparse.ArgumentParser()
    parse.add_argument('-p',action='store',dest='pdbFl',help='pdb file')
    parse.add_argument('-c',dest='chain',action='store',help='chain')
    presult = parse.parse_args()

    print "pdb fl:"+presult.pdbFl
    print "chain :"+presult.chain

    rco_fl,absCo_fl,lco_fl = rco(presult.pdbFl,presult.chain)
    print rco_fl,absCo_fl,lco_fl
