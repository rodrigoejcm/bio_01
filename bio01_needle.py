import pandas as pd
import math
import bio01_api as napi
import textwrap

import warnings
warnings.filterwarnings('ignore')


def pairwise_align_needle(seq1,seq2,smatrix,api,dna):

        score_matrix = define_score_matrix(smatrix)

        GAP_COST = -10

        SEQ_1 = list(seq1.upper())
        SEQ_2 = list(seq2.upper())

        SEQ_1.insert(0,'')
        SEQ_2.insert(0,'')

        GA = pd.DataFrame(columns=SEQ_1, index=SEQ_2)

        for i in range(0, len(SEQ_1)):
                GA.iloc[0, i]=[i*GAP_COST,'LR']

        for j in range(0, len(SEQ_2)):
                GA.iloc[j, 0]=[j*GAP_COST, 'TD']


        GC = GA.copy()
        row = GC.shape[0]
        col = GC.shape[1]


        for j in range(1,col):
                for i in range(1,row):
                        pw_col = GC.columns[j]
                        pw_row = GC.index[i]
                        cost = find_in_substitution_matrix(score_matrix,pw_col,pw_row)

                        #TOP DOWN
                        TD = GC.iloc[i -1,j][0] + GAP_COST
                        #LEFT RIGHT
                        LR = GC.iloc[i ,j - 1][0] + GAP_COST
                        # DIAGONAL
                        DD = GC.iloc[i -1,j - 1][0] + cost

                        ops = {'DD': DD, 'LR': LR, 'TD': TD  }

                        GC.iloc[i,j] = [ max(ops.values()) , max(ops, key=ops.get)]

        ##GC.to_csv('example.csv')

        i = row-1
        j = col-1

        align_x = ''
        align_y = ''
        score = 0

        while True:
                if(i == 0 and j == 0):
                        break

                inicio = GC.iloc[i,j]

                if inicio[1] == 'DD':
                        align_x += str(GC.columns[j])
                        align_y += str(GC.index[i])
                        score += find_in_substitution_matrix(score_matrix,GC.columns[j],GC.index[i])
                        i-=1
                        j-=1

                elif inicio[1] == 'LR':
                        align_x += str(GC.columns[j])
                        align_y += '-'
                        j-=1
                        #score += GAP_COST
                elif inicio[1] == 'TD':
                        align_x += '-'
                        align_y += str(GC.index[i])
                        i-=1
                        #score += GAP_COST

        #Seq invertida
        local_seq_1 = align_x[::-1]
        local_seq_2 = align_y[::-1]

        print("\n\n--- NEEDLE ---")
        print("Type {}".format("DNA" if dna  else "Protein"))
        print("Gap Cost {}".format(GAP_COST))
        print("Score Matrix {}".format(smatrix))
        print("Score {}".format(score))

        textwrap.wrap(align_x[::-1], width=100)
        print(textwrap.fill(align_x[::-1], width=100))
        textwrap.wrap(align_y[::-1], width=100)
        print(textwrap.fill(align_y[::-1], width=100))

        ## if api is set....
        if api:
                napi.pairwise_align_needle_api(seq1.upper(),seq2.upper(),smatrix, dna)


def define_score_matrix(sco_matrix):
        # https://github.com/JGI-Bioinformatics/DASW/blob/master/data/substitution_matrices/EBLOSUM50
        if sco_matrix == "blosum50":
                arquivo = pd.read_csv('matrix/EBLOSUM50.txt', sep=" ", index_col=0)
        else:
                arquivo = pd.read_csv('matrix/'+sco_matrix+'.txt', sep="  ")
        return  arquivo

## Function to find values on the table
def find_in_substitution_matrix(smatrix,row,column):
        return smatrix.get_value(row, column)


#####
# HBA_HUMAN Hemoglobin subunit alpha OS=Homo sapiens OX=9606 GN=HBA1 PE=1 SV=2
### MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR

# HBB_HUMAN Hemoglobin subunit beta OS=Homo sapiens OX=9606 GN=HBB PE=1 SV=2
### MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH

### MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH
#python3 bio01.py --seq1 "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR" --seq2 "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH" -a





