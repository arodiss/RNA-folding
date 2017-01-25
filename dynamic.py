import sys
import re
import numpy as np


def get_pairings_num(sequence):
    sequence = str(sequence).strip()
    length = len(sequence)

    if length == 0:
        print ("Please provide sequence to analyze")
        sys.exit(1)

    if len(re.findall("^([UACG]+)$", sequence)) == 0:
        print ("RNA sequence can only contain literals A, U, C, G")
        sys.exit(1)

    results = np.full(fill_value=None, shape=(length, length))
    for i in range(0, length):
        for j in range(0, length):
            if j < i + 5:
                results[i, j] = 0

    for k in range(5, length):
        for i in range(0, length-k):
            j = i + k
            if pairs(sequence[i], sequence[j]):
                possible_ts = []
                for t in range(i + 1, j - 1):
                    possible_ts.append(
                        1 + results[i, t - 1] + results[t + 1, j]
                    )
                results[i, j] = max(
                    results[i, j-1],
                    max(possible_ts)
                )
            else:
                results[i, j] = results[i, j - 1]

    return results[0, length - 1]


def pairs(literal, other_literal):
    if literal == "A":
        return other_literal == "U"
    if literal == "U":
        return other_literal == "A"
    if literal == "C":
        return other_literal == "G"
    if literal == "G":
        return other_literal == "C"
    raise ValueError('Invalid literal: ' + literal)


if __name__ == "__main__":
    if len(sys.argv) <= 1:
        print ("Please provide sequence to analyze")
        sys.exit(1)
    print (get_pairings_num(sys.argv[1]))