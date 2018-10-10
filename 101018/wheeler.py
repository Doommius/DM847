#! /usr/bin/env python
"""
A simple Burrows-Wheeler transform function in python.
Algorithm presented in:
Burrows M, Wheeler DJ: A Block Sorting Lossless Data Compression Algorithm.
    Technical Report 124. Palo Alto, CA: Digital Equipment Corporation; 1994.
USAGE: bwt.py [-h] [-i INDEX] STRING
"""

import argparse

def bw_transform(s):
    n = len(s)
    m = sorted([s[i:n]+s[0:i] for i in range(n)])
    I = m.index(s)
    L = ''.join([q[-1] for q in m])
    return (I, L)

from operator import itemgetter

def bw_restore(I, L):
    n = len(L)
    X = sorted([(i, x) for i, x in enumerate(L)], key=itemgetter(1))

    T = [None for i in range(n)]
    for i, y in enumerate(X):
        j, _ = y
        T[j] = i

    Tx = [I]
    for i in range(1, n):
        Tx.append(T[Tx[i-1]])

    S = [L[i] for i in Tx]
    S.reverse()
    return ''.join(S)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Transforms strings using Burrows-Wheeler transform')
    parser.add_argument('-i', '--index', metavar="INDEX", help='Performs ' +
        'inverse BWT. A transform index must be specified as an INDEX value.')
    parser.add_argument('STRING', type=str, help='Transform this string.')
    args = parser.parse_args()

    if args.index:
        print(bw_restore(args.index, args.STRING))
    else:
        print("%d %s" % bw_transform(args.STRING))