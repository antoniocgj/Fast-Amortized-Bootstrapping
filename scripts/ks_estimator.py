import random
from math import log2, sqrt

N = 2048
h = 64


def pass1(v:list, val):
  res = v.copy()
  for i in range(0, N, val):
    for j in range(i, i + val):
      if(j in res):
        res.remove(j)
        break
  return res

def pass_all(B, v:list):
  val = B
  k = 0
  while(len(v) > 0):
    v = pass1(v, val)
    k+=1
    # print(len(v), val, v)
    # val *= 2
  return k
  
attempts = 100
def get_k_50(B):
  k_vec = [0]*h*2
  for _ in range(attempts):
    s_idx = [random.randint(0,N-1) for _ in range(h)]
    k_vec[pass_all(B, s_idx)] += 1

  for i in range(2*h):
    if(sum(k_vec[:i]) > attempts//2):
      return i

print("B","k","cost", "noise", "size", sep="\t\t")
for i in range(1, int(log2(N)+1)):
  B = 1 << i
  k = get_k_50(B)
  cost = 2*k*(N/B)*B
  noise = N * sqrt(2* (k/B) * log2(B))
  size = 4 * (N**2)*(k/B)* log2(B)
  print(B, k, "%.2f" % log2(cost), "%.2f" % log2(noise), "%.2f" % log2(size), sep="\t\t")