import argparse
import numpy as np

from itertools import groupby

BASE_TO_INDEX = {"A":0, "C":1, "G":2, "T":3}
INDEX_TO_BASE = {0:"A", 1:"C", 2:"G", 3:"T"}


def base_to_index(seq):
    index_seq = []
    for base in seq:
        index_seq.append(BASE_TO_INDEX[base])
    return index_seq


def index_to_base(seq):
    base_seq = ""
    for index in seq:
        base_seq += INDEX_TO_BASE[index]
    return base_seq


def init_alignment_matrix(seq1, seq2, score):
    n = len(seq1)
    m = len(seq2)
    alignment_matrix = np.zeros((n+1,m+1))
    # fill first row and column
    alignment_matrix[0, 0] = 0
    for j in range(1,m+1):
        alignment_matrix[0,j] = alignment_matrix[0,j-1] + score[4, seq2[j-1]]
    for i in range(1,n+1):
        alignment_matrix[i,0] = alignment_matrix[i-1,0] + score[4, seq1[i-1]]
    return alignment_matrix


def fill_alignment_matrix(alignment_matrix, seq1, seq2, score, can_be_neg):
    n = len(seq1)
    m = len(seq2)
    path = np.zeros((n+1,m+1))
    min_value = - np.PINF if can_be_neg else 0
    vals = np.full((4, ), min_value)
    # fill rest of array
    for i in range(1,n+1):
        for j in range(1,m+1):
            vals[0] = alignment_matrix[i-1,j] + score[seq1[i-1], 4]
            vals[1] = alignment_matrix[i,j-1] + score[4, seq2[j-1]]
            vals[2] = alignment_matrix[i-1,j-1] + score[seq1[i-1], seq2[j-1]]
            alignment_matrix[i,j] = np.amax(vals)
            path[i,j] = np.argmax(vals)
    return alignment_matrix, path


def traceback(path, i, j, seq1, seq2, condition):
    # traceback path to reconstruct the alignment
    trace1, trace2 = "", ""
    while condition(i, j):
        if path[i,j] == 0:
            trace1 += INDEX_TO_BASE[seq1[i-1]]
            trace2 += '-'
            i -= 1
        elif path[i,j] == 1:
            trace1 += '-'
            trace2 += INDEX_TO_BASE[seq2[j-1]]
            j -= 1
        elif path[i,j] == 2:
            trace1 += INDEX_TO_BASE[seq1[i-1]]
            trace2 += INDEX_TO_BASE[seq2[j-1]]
            i -= 1
            j -= 1
        else:
            break
    # reverse the aligned sequences
    return trace1[::-1], trace2[::-1], i, j


def print_result(trace1, trace2, alignment_type, max_score):
    # print the aligned sequences (50 chars width) and score
    i,j = 0,0
    while i < len(trace1)-1:
        j = min(i+50, len(trace1))
        print(trace1[i:j])
        print(trace2[i:j],'\n')
        i = j
    print("%s:%d" % (alignment_type, max_score))


def global_alignment(seq1, seq2, score):
    n = len(seq1)
    m = len(seq2)

    alignment_matrix = init_alignment_matrix(seq1, seq2, score)
    
    alignment_matrix, path = fill_alignment_matrix(alignment_matrix, seq1, seq2, score, True)
    
    trace1, trace2, _, _ = traceback(path, n, m, seq1, seq2, lambda x, y: x+y > 0)
    
    print_result(trace1, trace2, 'global', alignment_matrix[n,m])


def local_alignment(seq1, seq2, score):
    n = len(seq1)
    m = len(seq2)

    alignment_matrix, path = fill_alignment_matrix(np.zeros((n+1,m+1)), seq1, seq2, score, False)
    
    i, j = np.unravel_index(np.argmax(alignment_matrix), alignment_matrix.shape)
    
    trace1, trace2, _, _ = traceback(path, i, j, seq1, seq2, lambda x, y: x > 0 and y > 0)
    
    print_result(trace1, trace2, 'local', alignment_matrix[i,j])


def overlap_alignment(seq1, seq2, score):
    n = len(seq1)
    m = len(seq2)

    alignment_matrix, path = fill_alignment_matrix(np.zeros((n+1,m+1)), seq1, seq2, score, True)
    
    i, j = n, np.argmax(alignment_matrix[n,:])
    
    trace1, trace2, i1, _ = traceback(path, i, j, seq1, seq2, lambda x, y: x > 0 and y > 0)
    
    trace1 = index_to_base(seq1[:i1]) + trace1 + '-' * (m - j)
    trace2 = '-' * i1 + trace2 + index_to_base(seq2[j:m])
    
    print_result(trace1, trace2, 'overlap', alignment_matrix[i,j])


def fastaread(fasta_name):
    f = open(fasta_name)
    faiter = (x[1] for x in groupby(f, lambda line: line.startswith(">")))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('seq_a', help='Path to first FASTA file (e.g. fastas/HomoSapiens-SHH.fasta)')
    parser.add_argument('seq_b', help='Path to second FASTA file')
    parser.add_argument('--align_type', help='Alignment type (e.g. local)', required=True)
    parser.add_argument('--score', help='Score matrix in.tsv format (default is score_matrix.tsv) ', default='score_matrix.tsv')
    command_args = parser.parse_args()

    _, seq1 = next(fastaread(command_args.seq_a))
    _, seq2 = next(fastaread(command_args.seq_b))
    score = np.genfromtxt(command_args.score, delimiter='\t')[1:,1:]
    
    if command_args.align_type == 'global':
        alignment_fn = global_alignment
    elif command_args.align_type == 'local':
        alignment_fn = local_alignment
    elif command_args.align_type == 'overlap':
        alignment_fn = overlap_alignment
    else:
        raise ValueError('no such alignment')
    alignment_fn(base_to_index(seq1), base_to_index(seq2), score)


if __name__ == '__main__':
    main()
