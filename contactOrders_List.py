#!/usr/bin/python
import re,sys
import numpy as Np

def calDist(atm1,atm2):
	return(Np.sum((atm1-atm2)**2))

def chainInfo(pdbfl):
    '''
    pdbfl : opened pdb file in .realines()
    '''
    pdbfl = open(pdbfl,'r').readlines()

    for k in pdbfl:
        if re.match(r'^COMPND',k) and re.search(r'CHAIN',k):
            chain = k.split(":")[1].strip().strip(";").split(",")
            return chain
        else:
            continue


def getAtomLine(l):
	atomDetail = {}
	atomDetail['atom_num']=int(l[5:11])
	atomDetail['atom_name'] = l[12:16].strip()
	atomDetail['res_name'] = l[17:20].strip()
	atomDetail['res_num'] = int(l[22:26])
	cord = (float(l[30:38]),float(l[38:46]),float(l[46:54]))
	atomDetail['cord'] = Np.array(cord)
	return atomDetail

def contactParameter(contactMatrix):
    """
    S_max,S_min,S_mid
    1.RCO: relative contact order,
    2.ACO: absolute contact order,
    3.SRO: short Range order,(2-4)
    4.MRO: Medium Range Order,(5-8)
    5.LRO: Long range Order,(>12)
    6.TCD: Total Contact Distance, ()
    7.N_alpha:
    8.Effective Length (L,L^(1/2),L^(2/3),L^(3/5))

    From contact matrix get indices(i,j) whose values are not zero. i-j is the sequence separion
    between two contacting residues.

    """	
    from collections import Counter
    contactRes_list = []
    length = Np.shape(contactMatrix)[0]
    print length
    for i in range(length-1):
        for j in range(i+1,length):
            if contactMatrix[j,i]:
                #print i,j
                contactRes_list.append(abs(j-i))
            else:
                continue
    contactRes_Dic = Counter(contactRes_list)
    #print contactRes_list,contactRes_Dic
    totalContact = Np.count_nonzero(contactMatrix)
    print "Total Contact: %d"%totalContact

    sro = []
    lro = []
    mrc_percent = []
    src_percent = []
    lrc_percent = []

    tcd = []
    local_CO = []
    nonlc_CO = []

    for k in contactRes_Dic.keys():
        if k >=2 and k <=8:
            local_CO.append(contactRes_Dic[k]*k)
            sro.append(contactRes_Dic[k])

            if k <= 4:
                src_percent.append(contactRes_Dic[k])
            else:
                mrc_percent.append(contactRes_Dic[k])
        elif k >= 12:
            nonlc_CO.append(contactRes_Dic[k]*k)
            lro.append(contactRes_Dic[k])
        elif k >= 3:
            tcd.append(contactRes_Dic[k]*k)
    #index = 0,1,2,3,4,5,6,

    if length and totalContact:

        sro_val = sum(sro)/float(length)
        lro_val = sum(lro)/float(length)

        src_per = (sum(src_percent)/float(totalContact))*100
        mrc_per = (sum(mrc_percent)/float(totalContact))*100
        lrc_per = (sum(lro)/float(totalContact))*100

        localCO_val = sum(local_CO)/float(totalContact)
        nonloCO_val = sum(nonlc_CO)/float(totalContact)

        tcd_val = (sum(tcd)/float(length*length))

        if totalContact:
            return totalContact,sro_val,lro_val,src_per,mrc_per,lrc_per,localCO_val,nonloCO_val,tcd_val
        else:
            return 0,0,0,0,0,0,0,0,0
    else:
        return 0,0,0,0,0,0,0,0,0

def gen_contactMatrix(pdbfl,chain):
    try :
        pdbFile = open(pdbfl,'r').readlines()
    except:
        return (0,0,0,0,0,0,0,0,0,0,0,0,0)

    print "processing %s"%pdbfl
    #chain = chain
    atomList = []
    for i in pdbFile:
        i = i.rstrip()
        if re.match(r'^ATOM',i) and i[21] == chain:
            atomList.append(getAtomLine(i))
        if re.match(r'^TER',i):
            break

    if len(atomList):
        length = atomList[len(atomList)-1]['res_num'] - atomList[0]['res_num']
    else:
        return 0,0,0,0

    #print atomList
    # Checking id contact or not
    resDictionary = {}
    for k in atomList:
        resDictionary[k['res_num']] = len(resDictionary.keys()) -1

    contacts = 0
    sequenceDist = 0
    residueContact = Np.zeros((length+1,length+1))

    for atom1 in atomList:
        for atom2 in atomList:
            seqDist = atom1['res_num']-atom2['res_num']
            if seqDist > 2:	# calculating non-adjecent residues
                dist = calDist(atom1['cord'],atom2['cord'])
                if dist < 36:	# 6*6 d comare with square dist
                    contacts +=1
                    sequenceDist += seqDist
                    residueContact[resDictionary[atom1['res_num']],resDictionary[atom2['res_num']]]=1

    print "Analyzing contact"
    #import ipdb;ipdb.set_trace()

    if length and contacts:
        totalRes_Contact,sro,lro,src_per,mrc_per,lrc_per,lCO,nloCO,tcd = contactParameter(residueContact)
        #print sro,lro,src_per,mrc_per,lrc_per,lCO,nloCO,tcd
        rco_fl = sequenceDist/float(length*contacts)
        absCO_fl = rco_fl*length
        #print contacts,sequenceDist, rco
        #print length 
        #longContactRes = Np.count_nonzero(residueContact)/float(length)
        # import ipdb; ipdb.set_trace()
        return (length,contacts,totalRes_Contact,rco_fl,absCO_fl,sro,lro,src_per,mrc_per,lrc_per,lCO,nloCO,tcd)
    else:
        return (0,0,0,0,0,0,0,0,0,0,0,0,0)


if __name__ == '__main__':
    """
    contact order parameters (relative and absolute) for pdbids
    Length
    Atom Contacts
    Residue contacts
    Relative contact order, Absolute contact order
    """
    import argparse
    parse = argparse.ArgumentParser()
    parse.add_argument('-p',action='store',dest='pdbFl',help='pdb file and chain')
    parse.add_argument('-d',dest='pdbDir',action='store',help='PDB directory As in Only Beta: /home/taushif/ThesisWork/allPDbs/Data/PDB_files/')
    parse.add_argument('-o',dest='outfl',action='store',help='output file name')
    presult = parse.parse_args()

    print "pdb file with pid and chain:"+presult.pdbFl
    
    pfl = open(presult.pdbFl,'r').readlines()
    ofl = open(presult.outfl,'w')
    elog = open("error.log",'w')
    ofl.write("p.ID\tChain\tlength\tAtom_Con\tRes_Cont\tRCO\tACO\tSRO\tLRO\tSR_per\tMR_per\tLR_per\tLCO\tNlCO\tTCD\n")
    pdbDir = presult.pdbDir
    #pdbDir = "/home/taushif/ThesisWork/Codes_Data/alphaWork/Data2_7June/PDB/"
    #pdbDir = "/home/taushif/PyCodes/ContactOrders/pdbs/"
    #pdbDir = "/home/taushif/PyCodes/ContactOrders/kineticDataset/pdbs/"
    chainCount_suc = 0
    chainCount_fail = 0

    for k in pfl:
        k = k.strip().split()
        #dbfl = pdbDir+"pdb"+k[0]+".ent"
        pdbfl = pdbDir+k[0]+".pdb"
        #chain = chainInfo(pdbfl)
        #print chain
        chain = k[1]
    #    for c in chain:
        res = gen_contactMatrix(pdbfl,chain)
        #protl,tAContact,tRContact,rco_fl,absCo_fl,sro,lro,src_p,mrc_p,lrc_p,lCO,nloCO,tcd 
        if res[0]:
            print "%d\t%d\t%d\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.6f\t"\
                %(res[0],res[1],res[2],res[3],res[4],res[5],res[6],res[7],res[8],res[9],res[10],res[11],res[12])
            ofl.write("%s\t%s\t%d\t%d\t%d\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.6f\n"\
                %(k[0],chain,res[0],res[1],res[2],res[3],res[4],res[5],res[6],res[7],res[8],res[9],res[10],res[11],res[12]))
            chainCount_suc += 1
        else:
            print "problem in protein %s"%k
            chainCount_fail += 1
            #elog.write("%s\t%s\n"%(k[0],k[1]))
            elog.write("%s\n"%(k))
            continue
        print "\nDONE : [%d/%d] FAIL: [%d/%d]\n"%(chainCount_suc,len(pfl),chainCount_fail,len(pfl))
        #import ipdb; ipdb.set_trace()
    ofl.close()
    elog.close()
