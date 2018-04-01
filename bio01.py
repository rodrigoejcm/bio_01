from optparse import OptionParser
import bio01_needle as needle

usage = "Usage: %prog [options...] [seqFile]"
parser = OptionParser(usage=usage)

parser.add_option('--seq1', type="string", dest="seq1", help='First input sequence')
parser.add_option('--seq2', type="string", dest="seq2", help='Second input sequence')
parser.add_option("-a", action="store_true", dest="api" , help='Includes Pairwise alignment from EMBOSS Api (Default: False )' , default=False)
parser.add_option("-d", action="store_true", dest="dna" , help='DNA pairwise align (Default: Protein)', default=False)
parser.add_option('--smatrix', help='Scoring matrix: pam250/pam30/blosum50/blosum62/dndfull', type='choice', action='store', dest='smatrix',  choices=['pam250', 'pam30', 'blosum50', 'blosum62', 'dnafull'],default="blosum50")

(options, args) = parser.parse_args()


if not options.seq1 or not options.seq2 :
        parser.print_help()
elif options.dna and options.smatrix != "dnafull":
        print("If DNA pairwise align, then scoring matrix should be dnafull")
        parser.print_help()
else:
        needle.pairwise_align_needle(options.seq1,options.seq2,options.smatrix,options.api, options.dna)