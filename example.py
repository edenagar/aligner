import json

with open(
        '/Users/edenn/bio-info/Sanger/Sanger Image Processor - application/source_code/sanger_alignment/examples/references.json',
        'r') as fp:
    referenceData = json.load(fp)

def aligner_example():
    dataSanger_hr110 = SeqIO.read(
        '/Users/edenn/bio-info/Sanger/Sanger Image Processor - application/source_code/sanger_alignment/examples/N2_em_nhr_100.ab1',
        "abi")
    dataSanger_C52_em_117 = SeqIO.read(
        "/Users/edenn/bio-info/Sanger/Sanger Image Processor - application/source_code/sanger_alignment/examples/C52_em_117.ab1",
        "abi")
    dataSanger_ALM505_F48_88 = SeqIO.read(
        "/Users/edenn/bio-info/Sanger/Sanger Image Processor - application/source_code/sanger_alignment/examples/ALM505_F48_88.ab1",
        "abi")
    dataSanger_F28_1_52 = SeqIO.read(
        "/Users/edenn/bio-info/Sanger/Sanger Image Processor - application/source_code/sanger_alignment/examples/F28-1_52.ab1",
        "abi")

    seq1 = matter_to_seq(dataSanger_hr110)[0]
    seq2 = matter_to_seq(dataSanger_C52_em_117)[0]
    seq3 = matter_to_seq(dataSanger_ALM505_F48_88)[0]
    seq4 = matter_to_seq(dataSanger_F28_1_52)[0]

    seq_ref1 = Seq(referenceData['nhr_110']).reverse_complement()
    seq_ref2 = Seq(referenceData['C52_em_117'])

    temp1 = sanger_alignment_v2(seq_ref1, seq1)
    print("dataSanger_hr110:")
    print(temp1[0])
    temp2 = sanger_alignment_v2(seq_ref2, seq2)
    print("dataSanger_C52_em_117:")
    print(temp2[0])



def heatmap_example():
    sanger_heatmap(SeqIO.read(
        "/Users/edenn/bio-info/Sanger/Sanger Image Processor - application/source_code/sanger_alignment/examples/N2_em_nhr_100.ab1",
        "abi"), 0, 613)




def missmatch_example():
    dataSanger_hr110 = SeqIO.read(
        '/Users/edenn/bio-info/Sanger/Sanger Image Processor - application/source_code/sanger_alignment/examples/N2_em_nhr_100.ab1',
        "abi")
    [seq1, index1] = matter_to_seq(dataSanger_hr110)
    seq_ref1 = Seq(referenceData['nhr_110']).reverse_complement()
    temp1 = sanger_alignment_v2(seq_ref1, seq1)
    print(temp1[0])
    changed = get_index_of_changed_nuc(temp1)
    print("arrey of missmatch nucleotides: ")
    print(index1)
    print("we will find the epsilon around missmatch number 7:")
    print(index1[changed[4]])
    area_by_index(index1[changed[4]], dataSanger_hr110, 50)  # SHOWS THAT THERE RESIDUE BETWEEN TWO PEAKS



test = ab1_anelyzer('/Users/edenn/bio-info/Sanger/Sanger Image Processor - application/source_code/sanger_alignment/examples/N2_em_nhr_100.ab1')

test.alignWith(Seq(referenceData['nhr_110']))


# aligner_example()
#heatmap_example()
# missmatch_example()
#if theres a miss match find all amplitude in amplitude or the changed amplitude
#waight by location in sanger