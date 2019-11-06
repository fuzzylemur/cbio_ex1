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

    h1, s = next(fastaread(command_args.seq_a))
    h2, t = next(fastaread(command_args.seq_b))
    s = "TTT"
    t = "AA"
    n = len(s)
    m = len(t)
    score = np.genfromtxt(command_args.score, delimiter='\t')[1:,1:]

    d = {"A":0, "C":1, "G":2, "T":3}

    arr = np.zeros((n,m))
    path = np.zeros((n,m))

    arr[0, 0] = score[d[s[0]], d[t[0]]]
    for i in range(1,n):
        arr[i,0] = arr[i-1,0] + score[4, d[s[i]]]

    for i in range(n):
        for j in range(1,m):
            val1 = arr[i-1,j] + score[d[s[i]], 4]
            val2 = arr[i,j-1] + score[4, d[t[j]]]
            val3 = arr[i-1,j-1] + score[d[s[i]], d[t[j]]]
            vals = [val1,val2,val3]
            print(i, j, vals)
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
