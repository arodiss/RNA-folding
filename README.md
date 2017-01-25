# RNA-folding
Simple RNA folding model.
Created for test purposes, not usable for any real-life application.

## Supported rules
Upon finding solution following rules are considered (not exctly what real RNA does):
- **complementarity**: U pairs A, C pairs G
- **no sharp corners**: any pair should be separated by at least 4 literals (nukleotids)
- **no crossings**: if position A pairs with position B (B > A) and position Y pairs with position Z (Z > Y), it's either `A > Z`, `Y > B`, `B > Z && Y > A` or `Z > B && A > Y`
- one literal (nukleotide) can be paired no more than one time
From all possible pairings, only one with maximum number of connections is considered

## Running
```
python3 pairings.py --method=dynamic --sequence=AAAUUCCGG
python3 pairings.py --method=recursive --sequence=AAAUUCCGG
```
