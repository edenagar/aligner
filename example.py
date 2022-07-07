import json
from Bio import SeqIO
from Bio.Seq import Seq

from aligner import ab1_anelyzer

references_path = 'samples/references.json'
referenceData = {}
with open(references_path) as fp:
    referenceData = json.load(fp)

referenceData_nhr_100 = Seq(referenceData["nhr_100"])
N2_em_nhr_100 = 'samples/N2_em_nhr_100.ab1'
C52_em_117 = 'samples/C52_em_117.ab1'


def align_example():
    example = ab1_anelyzer(N2_em_nhr_100)
    example.align_with(referenceData_nhr_100)
    example.print_best_alignment()
    example.txt_best_alignment()
    example.plot_heat_map(0, 613)


align_example()

def missmatch_example():
    example = ab1_anelyzer(N2_em_nhr_100)
    example.align_with(referenceData_nhr_100)
    example.plot_matter_sequence(example.changed[4], 550)

missmatch_example()