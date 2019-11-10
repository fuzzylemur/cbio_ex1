import argparse
import numpy as np

from itertools import groupby


def global_alignment(seq1, seq2, score):
    n = len(seq1)
    m = len(seq2)
    
    d = {"A":0, "C":1, "G":2, "T":3}

    arr = np.zeros((n+1,m+1))
    path = np.zeros((n+1,m+1))

    # fill first row and column
    arr[0, 0] = 0
    for j in range(1,m+1):
        arr[0,j] = arr[0,j-1] + score[4, d[seq2[j-1]]]
    for i in range(1,n+1):
        arr[i,0] = arr[i-1,0] + score[4, d[seq1[i-1]]]

    # fill rest of array
    for i in range(1,n+1):
        for j in range(1,m+1):
            val1 = arr[i-1,j] + score[d[seq1[i-1]], 4]
            val2 = arr[i,j-1] + score[4, d[seq2[j-1]]]
            val3 = arr[i-1,j-1] + score[d[seq1[i-1]], d[seq2[j-1]]]
            vals = [val1, val2, val3]
            arr[i,j] = max(vals)
            path[i,j] = np.argmax(vals)

    #print(arr)
    #print(path)

    # traceback path to reconstruct the alignment
    trace1, trace2, i, j = "", "", n, m
    while i+j > 0:
        if path[i,j] == 0:
            trace1 += seq1[i-1]
            trace2 += '-'
            i -= 1
        elif path[i,j] == 1:
            trace1 += '-'
            trace2 += seq2[j-1]
            j -= 1
        else:
            trace1 += seq1[i-1]
            trace2 += seq2[j-1]
            i -= 1
            j -= 1
    # reverse the aligned sequences
    trace1 = trace1[::-1]
    trace2 = trace2[::-1]

    # print the aligned sequences (50 chars width) and score
    i,j = 0,0
    while i < len(trace1)-1:
        j = min(i+50, len(trace1))
        print(trace1[i:j])
        print(trace2[i:j],'\n')
        i = j
    print("global:%d" % arr[n,m])


def local_alignment(seq1, seq2, score):
    n = len(seq1)
    m = len(seq2)


def overlap_alignment(seq1, seq2, score):
    n = len(seq1)
    m = len(seq2)


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
    alignment_fn(seq1, seq2, score)


if __name__ == '__main__':
    main()
