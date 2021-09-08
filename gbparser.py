import Bio.SeqIO as SeqIO
from Bio.SeqFeature import FeatureLocation
import pandas as pd
import re
import os
import numpy as np
scardict ={'CGAC' : 'I_scar', 'CTCC' : 'J_scar', 'TTGC' : 'K_scar', 'TAGT' : 'L_scar', 'GTCG' : 'M_scar', 'AAGC' : 'N_scar', 'CATT' : 'S_scar', 'AGTA' : 'alpha_scar', 'AGCG' : 'beta_scar', 'GGCA' : 'gamma_scar', 'ACTA' : 'delta_scar', 'GGGA' : 'epsilon_scar', 'TGCC' : 'zeta_scar', 'GCAA' : 'eta_scar', 'CGCT' : 'F_scar', 'GCTT' : 'E_scar', 'AGGT' : 'D_scar', 'AATG' : 'C_scar', 'TACG' : 'B_scar', 'GGAG' : 'A_scar', 'ATTG' : 'Y_scar', 'TGTC' : 'X_scar', 'TCTG' : 'V_scar', 'GGGC' : 'U_scar', 'GTAA' : 'R_scar', 'GAGT' : 'Q_scar', 'CCTA': 'P_scar', 'ATGC' : 'O_scar', 'ATAG' : 'G_scar', 'TACT' : 'H_scar'}

pdict = eval(open("~/Level_0_CDSs__RBSs/Level_0_CDSs__RBSs/dict.txt").read())
strdirectory = "~/gbfiles/type/Level_2_Destination_Vectors/"
directory = os.fsencode(strdirectory)
filelist = {}
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    
    filelist[filename] = 1
    if filename.endswith(".gb"):
        print(filename)
        scars = False
        SequenceDictionary = []
        handle = open(os.path.join(strdirectory, filename))
        record = SeqIO.read(handle, "gb")
        features = record.features
        for f in features:
                q = f.qualifiers
                loc = f.location
                start = loc.start
                end = loc.end
                s = loc.strand
                if str(f.type).lower() == "rbs & cds":
                        
                        types = q['label'][0].split("_")
                        if types[1] == "PhIF":
                            types[1] = "PhlF"
                        rbs = pdict[types[0]]
                        cds = pdict[types[1]]
                        
                        r = FeatureLocation(start if s==1 else end-len(rbs), start+len(rbs) if s==1 else end ,strand=s)
                        c = FeatureLocation(start+len(rbs) if s==1 else start, end if s==1 else end-len(rbs),strand=s)
                        print(pdict[types[0]] == str((r.extract(record).seq)))
                        print(pdict[types[1]] == str((c.extract(record).seq)))
                        SequenceDictionary  += [[types[0], str((r.extract(record).seq)), "rbs", (r.strand == -1), r.start + 1, r.end]]
                        SequenceDictionary  += [[types[1], str((c.extract(record).seq)), "cds", (c.strand == -1), c.start + 1, c.end]]
                elif "scar" in q['label'][0].lower() and (end-start) <= 4 :
                        scars = True
                        label = q['label'][0]
                        label = re.sub(r'\W+', '', label)
                        SequenceDictionary  += [[scardict[str((f.extract(record).seq))], str((f.extract(record).seq)), str(f.type).lower(), (loc.strand == -1), start+1, end+0]]
                elif str(f.type).lower() in ["misc_feature"] and 'promoter' in re.sub(r'\W+', '', q['label'][0]):
                        label = re.sub(r'\W+', '', q['label'][0])
                        SequenceDictionary  += [[label, str((f.extract(record).seq)), "promoter", (loc.strand == -1), start+1, end+0]]
                elif str(f.type).lower() not in ["source", "primer", "primer_bind", "rep_origin", "misc_feature", "repressor"]:
                        
                        label = q['label'][0]
                        label = re.sub(r'\W+', '', label)
                        SequenceDictionary  += [[label, str((f.extract(record).seq)), str(f.type).lower(), (loc.strand == -1), start+1, end+0]]
                SequenceDictionary = sorted(SequenceDictionary,key=lambda x: (x[-1]))
        if not scars:
            SequenceDictionary = []
        for i in range(len(SequenceDictionary)):
                if "2 scar" in SequenceDictionary[i][2]:
                       SequenceDictionary = SequenceDictionary[i:]
                       break
##                elif "0 scar" in SequenceDictionary[i][2]:
##                       SequenceDictionary = SequenceDictionary[i:]
##                       break
        for j in range(len(SequenceDictionary)-1, 0, -1):
                if "2 scar" in SequenceDictionary[j][2]:
                    SequenceDictionary = SequenceDictionary[:j+1]
                    break
##                elif "0 scar" in SequenceDictionary[j][2]:
##                    SequenceDictionary = SequenceDictionary[:j+1]
##                    break
        SequenceDictionary += [["plasmid_vector" if "Vector" in strdirectory else "source",str(record.seq),"plasmid_vector" if "Vector" in strdirectory else "source",(f.location.strand == -1),1,len(str(record.seq))]]
        print(SequenceDictionary[:-1])
        df = pd.DataFrame(SequenceDictionary)
        
        df.to_csv(strdirectory + filename[:-3].replace(" ", "_") + ".csv", header=False, index=False)
    else:
        continue

    
    


