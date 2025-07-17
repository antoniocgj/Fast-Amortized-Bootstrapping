import random
from math import log2, ceil

h = list(range(22,43))
N = [2048, 4096, 8192]
B_values = list(range(7, 12))
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
  print(N, h, end="\t", sep="\t")
  sum = 0
  for i in B_values:
    if(i not in B_hist):
      print("-inf" if sum == 0 else "0", end="\t")
      continue
    sum += B_hist[i]
    print("%.2f" % (log2(sum)), end="\t")
  print()

header = ["$B = 2^%d$" % i for i in B_values]
print("$N$", "$h$", *header , sep="\t")
for Ni in N:
  for hi in h:
    build_B_hist(Ni, hi)
    