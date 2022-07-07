import pandas as pd
from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.Align import PairwiseAlignment
from scipy.signal import find_peaks
import numpy as np
import matplotlib.pyplot as plt



def max_peaks(dataSanger, threshold=10):
	data9 = list(dataSanger.annotations["abif_raw"]["DATA9"])
	data10 = list(dataSanger.annotations["abif_raw"]["DATA10"])
	data11 = list(dataSanger.annotations["abif_raw"]["DATA11"])
	data12 = list(dataSanger.annotations["abif_raw"]["DATA12"])
	array_of_peaks = [data9, data10, data11, data12]

	data9_peak = (find_peaks(data9, threshold)[0])
	data10_peak = (find_peaks(data10, threshold)[0])
	data11_peak = (find_peaks(data11, threshold)[0])
	data12_peak = (find_peaks(data12, threshold)[0])
	data9_peak = list(map(lambda x: [x, 9], data9_peak))
	data10_peak = list(map(lambda x: [x, 10], data10_peak))
	data11_peak = list(map(lambda x: [x, 11], data11_peak))
	data12_peak = list(map(lambda x: [x, 12], data12_peak))
	peaks = data9_peak + data10_peak + data12_peak + data11_peak
	peaks.sort(key=lambda x: x[0])

	max_peaks = filter(
			lambda peak:
			max(data9[peak[0]], data10[peak[0]], data11[peak[0]], data12[peak[0]]) == array_of_peaks[peak[1] - 9][peak[0]],
			peaks)
	max_peaks = [[a[0],number_to_nucleotide((a[1]))] for a in list(max_peaks)]
	return { "data9_peak":data9_peak, 
						"data10_peak" : data10_peak,
						"data11_peak": data11_peak, 
						"data12_peak" : data12_peak,
						"max_peaks" : max_peaks}


def peaks_to_seq(arrey_of_peaks):
	array_of_nucleotides = [a[1] for a in arrey_of_peaks]
	pred_seq = ''.join(array_of_nucleotides)
	return pred_seq

class ab1_anelyzer:
	def __init__(self, path):
			self.seqIO = SeqIO.read(path, "abi")
			self.peaks = max_peaks(self.seqIO)
			self.seq = peaks_to_seq(self.peaks[max_peaks])


	def alignWith(self, reference : str):
			alignment = sanger_alignment_v2(Seq(reference).reverse_complement(),self.seq)
			alignment_reverse = sanger_alignment_v2(Seq(reference),self.seq)
			return alignment.score >= alignment_reverse.score if alignment else alignment_reverse

	def sanger_heatmap(ab1_file, start=0, finish=0):
		G = ab1_file.annotations["abif_raw"]["DATA9"]
		A = ab1_file.annotations["abif_raw"]["DATA10"]
		T = ab1_file.annotations["abif_raw"]["DATA11"]
		C = ab1_file.annotations["abif_raw"]["DATA12"]
		if finish == 0:
				finish = len(G)
		diff_from_max = []
		for i in range(len(G)):
				list = [G[i], A[i], T[i], C[i]]
				max_parameter = max(list)
				list.remove(max(list))
				second_max_parameter = max(list)
				diff_from_max.append(max_parameter - second_max_parameter)
		plt.rcParams["figure.figsize"] = 5, 2

		x = np.array(range(start, finish))
		y = np.array(diff_from_max[start:finish])

		fig, (ax, ax2) = plt.subplots(nrows=2, sharex=True)

		extent = [x[0] + 1 - (x[1] - x[0]) / 2., x[-1] + (x[1] - x[0]) / 2., 0, 1]
		ax.imshow(y[np.newaxis, :], cmap="plasma", aspect="auto", extent=extent)
		ax.set_yticks([])
		ax.set_xlim(extent[0], extent[1])

		ax2.plot(x, y)

		plt.tight_layout()
		plt.show()


referenceData = {}






def get_index_of_changed_nuc(sequence: PairwiseAlignment):
	seq_aligned = sequence[0].aligned[1]
	reverted_list = [a[1] for a in seq_aligned]
	return reverted_list



def sanger_alignment_v2(reference: Seq, seq_sanger: Seq, match_score: int = 20, mismatch_score: int = -20,
											extend_gap_score: int = -10, open_gap_score: int = -10) -> PairwiseAlignment:

	aligner = Align.PairwiseAligner()
	aligner.mode = 'global'
	aligner.match_score = match_score
	aligner.mismatch_score = mismatch_score
	aligner.query_gap_score = extend_gap_score
	aligner.open_gap_score = open_gap_score

	aligner.query_right_open_gap_score = 0
	aligner.query_right_extend_gap_score = 0
	aligner.query_left_open_gap_score = 0
	aligner.query_left_extend_gap_score = 0

	# aligner.algorithm = 'Smith-Waterman' TO KNOW: if not set

	return  aligner.align(reference, seq_sanger)


def number_to_nucleotide(number):
	if number == 9: return 'G'
	if number == 10: return 'A'
	if number == 11: return 'T'
	if number == 12: return 'C'


def matter_to_seq(dataSanger, threshold=10):
	data9 = list(dataSanger.annotations["abif_raw"]["DATA9"])
	data10 = list(dataSanger.annotations["abif_raw"]["DATA10"])
	data11 = list(dataSanger.annotations["abif_raw"]["DATA11"])
	data12 = list(dataSanger.annotations["abif_raw"]["DATA12"])
	array_of_peaks = [data9, data10, data11, data12]

	data9_peak = (find_peaks(data9, threshold)[0])
	data10_peak = (find_peaks(data10, threshold)[0])
	data11_peak = (find_peaks(data11, threshold)[0])
	data12_peak = (find_peaks(data12, threshold)[0])
	data9_peak = list(map(lambda x: [x, 9], data9_peak))
	data10_peak = list(map(lambda x: [x, 10], data10_peak))
	data11_peak = list(map(lambda x: [x, 11], data11_peak))
	data12_peak = list(map(lambda x: [x, 12], data12_peak))
	peaks = data9_peak + data10_peak + data12_peak + data11_peak
	peaks.sort(key=lambda x: x[0])

	max_peaks = filter(
			lambda peak:
			max(data9[peak[0]], data10[peak[0]], data11[peak[0]], data12[peak[0]]) == array_of_peaks[peak[1] - 9][peak[0]],
			peaks)
	max_peaks = list(max_peaks)

	array_of_nucleotides = [number_to_nucleotide(a[1]) for a in max_peaks]
	pred_seq = ''.join(array_of_nucleotides)
	return [pred_seq, [a[0] for a in max_peaks]]


def area_by_index(index: int, ab1_file, epsilon: int = 10):
	array_of_peaks = matter_to_seq(ab1_file)[1]
	min_index = max(index - epsilon, 0)
	max_index = min(len(array_of_peaks) - 1, index + epsilon)

	G = ab1_file.annotations["abif_raw"]["DATA9"][min_index:max_index]
	A = ab1_file.annotations["abif_raw"]["DATA10"][min_index:max_index]
	T = ab1_file.annotations["abif_raw"]["DATA11"][min_index:max_index]
	C = ab1_file.annotations["abif_raw"]["DATA12"][min_index:max_index]
	table=[G,A,T,C]
	df = pd.DataFrame(table,columns = range(index-epsilon,index+ epsilon),index=['G: ','A: ','T: ','C: '])
	print(df.to_string())
	X = range(min_index, max_index)
	plt.plot(X, G, 'black')
	plt.plot(X, A, "g")
	plt.plot(X, T, 'r')
	plt.plot(X, C, 'b')
	plt.show()




