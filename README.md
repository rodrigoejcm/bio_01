# Bio01 - Minitrabalho

Implementação do algoritimo de alinhamento global de needleman-wunsch para a disciplina de bioinformática

### Files:

- bio01.py
- bio01_needle.py
- bio01_needle.py

### Instructions
```
pip install "the required packages listed in requirements.txt"
```

### Parameters details: 

```
python bio01.py -h
```
```
  --seq1=SEQ1        First input sequence
  --seq2=SEQ2        Second input sequence
  -a                 Includes Pairwise alignment from EMBOSS Api (Default: False )
  -d                 DNA pairwise align (Default: Protein)
  --smatrix=SMATRIX  Scoring matrix: pam250/pam30/blosum50/blosum62/dndfull (Default= blosum50)
 ```
 
 ### Features:
 
 - In order to compare two sequences, define in eache paramenter ( seq1 & seq2 ) a string with the protein/dna sequence
 - It is possible to also run the Emboss Needle API alongside the local algorithm with the flag '-a' for API.
 - If DNA comparison is needed, the flag '-d' should be present
 - Besides the default scoring matrix (blosum50), there are others to choose from the above list.
 
 ### Examples

Protein:
```
# Default Scoring matrix and Emboss API comparison
python bio01.py --seq1 "HEAGAWGHEE" --seq2 "PAWHEAE" -a

# Blosum62 Scoring matrix
python bio01.py --seq1 "HEAGAWGHEE" --seq2 "PAWHEAE" --smatris="blosum62"

```

DNA:
```
#DNA with dnafull scoring matrix
python bio01.py --seq1 "HEAGAWGHEE" --seq2 "PAWHEAE" -d --smatris="dnafull"
```


 
 
