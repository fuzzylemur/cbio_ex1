import argparse
import numpy as np

from itertools import groupby

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

    h1, seq1 = next(fastaread(command_args.seq_a))
    h2, seq2 = next(fastaread(command_args.seq_b))
    seq1 = "GGCC"
    seq2 = "GGTTCC"
    n = len(seq1)
    m = len(seq2)
    score = np.genfromtxt(command_args.score, delimiter='\t')[1:,1:]

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
            #print(i, j, vals)
            arr[i,j] = max(vals)
            path[i,j] = np.argmax(vals)

    max_score = arr[i,j]

    trace1, trace2 = "",""
    print(arr)
    print(path)


    if command_args.align_type == 'global':
        raise NotImplementedError
    elif command_args.align_type == 'local':
        raise NotImplementedError
    elif command_args.align_type == 'overlap':
        raise NotImplementedError
    # print the best alignment and score


if __name__ == '__main__':
    main()
