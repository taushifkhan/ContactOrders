# ContactOrders

Availability of experimental data enabled statistical inference of topological 
descriptors to its folding kinetics. Folding and unfolding rates have been correlated 
to short and long range contacts. Short contacts increases the propensity of 
local structural motifs like helix, beta hairpins, turns, whereas, long range 
interactions contributes the formation of beta structures and complex topologies. 
However, protein length has shown as positive correlation with logarithm of rate of 
folding, over the years many topological descriptors have been proposed to increase 
the predictability of protein folding/unfolding kinetics. 


contact order parameters (relative and absolute) for pdbids
for a defined Minimum residue distance to maximum residue distance cutoff

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

#Usage

--------------Analyzing multiple pdb files (input as list line terminated and tab separated)
    $ python contactOrders_List.py -h
    usage: contactOrders_List.py [-h] [-p PDBFL] [-d PDBDIR] [-o OUTFL]

    optional arguments:
      -h, --help  show this help message and exit
      -p PDBFL    pdb file and chain
      -d PDBDIR   PDB directory As in Only Beta:
                  ./allPDbs/Data/PDB_files/
      -o OUTFL    output file name

    ----------- For individual pdb file--------------

    python contactOrders.py -h
    usage: contactOrders.py [-h] [-p PDBFL] [-c CHAIN]
    
    optional arguments:
      -h, --help  show this help message and exit
      -p PDBFL    pdb file
      -c CHAIN    chain
    

#Reference

Multiple research article have been considered to analyze above topological descriptors,
some are listed as folowing,
Zou, Ozkon (2011, Proteins), Selveraja and Gorima (2004-2009 papers), Zhou and Zohu (2004),
Ouyang and Liang (2002), Ivanakov and Baker (2003) 
+
