import random
from math import log2, ceil

h = list(range(17,34))
N = [2048]
REPS = 100000

def build_B_hist(N, h):
  B_hist = {}
  for i in range(REPS):
    key = sorted([random.randint(0,N) for _ in range(h)])
    diff_set = [(key[i] - key[i-1]) for i in range(1, h)] + [N - key[h-1]]
    B = ceil(log2(max(diff_set)))
    B_hist[B] = B_hist.get(B, 0) + 1
  for b in B_hist:
    B_hist[b] /= REPS
  print("N: %d h: %d" % (N, h), end=" ")
  sum = 0
  for i in sorted(B_hist.items()):
    sum += i[1]
    print("(B = 2^%d, δ = %.2f)" % (i[0], -log2(sum)), end=", ")
  print()

for Ni in N:
  for hi in h:
    build_B_hist(Ni, hi)