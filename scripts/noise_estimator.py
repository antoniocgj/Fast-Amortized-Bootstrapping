from mpmath import mpf, mp, nstr
import math
mp.dps = 500

class LeveledRGSW:
  def __init__(self, N_in, h, h_packing, h_out, key_base, sigma_in, sigma_packing, N, sigma, k, ell, bg_bit, ell_packing, b_packing, ell_ks, b_ks, msg_prec, r_prec, verbose=False, key_type="binary_hw"):
    self.msg_prec, self.sigma_packing, self.sigma_in, self.sigma, self.N_in, self.N, self.k, self.ell, self.bg_bit, self.verbose, self.h, self.h_packing, self.h_out, self.key_base, self.ell_packing, self.b_packing, self.ell_ks, self.b_ks, self.r_prec, self.key_type = mpf(msg_prec), mpf(sigma_packing), mpf(sigma_in), mpf(sigma), (N_in), (N), (k), (ell), (bg_bit), verbose, h, h_packing, h_out, key_base, ell_packing, b_packing, ell_ks, b_ks, r_prec, key_type
    self.bg = mpf(1 << self.bg_bit)
    self.beta = self.bg/2
    self.epsilon = 1 /(2*(self.bg**self.ell))
    self.sigma_var = sigma**2

  def vertical_packing(self, in_prec, LUT_var, suppress=False):
    term1 = (self.k + 1)*self.ell*self.N*(self.beta**2)*self.sigma_var
    term2 = (1 + self.k*self.h_out/3)*self.epsilon**2
    if(self.verbose and not suppress):
      print("\n## Vertical Packing ##")
      print("Variance <= LUT_var + precison * (Gaussian + Mod Switching)")
      print("Variance <= %d + %d * (2^%f + 2^%f) = 2^%f" %(LUT_var, in_prec, mp.log(term1,2), mp.log(term2,2), mp.log(in_prec*(term1 + term2), 2)))
      print("Failure Probability: 2^%s" % (self.failure_rate(in_prec*(term1 + term2))))
    return in_prec * (term1 + term2) + LUT_var

  def GSW_product(self, A_var, B_var, suppress=False):
    term1 = (self.k + 1)*self.ell*self.N*(self.beta**2)*A_var
    term2 = (1 + self.k*self.h_out/3)*self.epsilon**2
    if(self.verbose and not suppress):
      print("\n## GSW product ##")
      print("Variance <= Gaussian + Mod Switching + Variance(B)")
      print("Variance <= 2^%f + 2^%f + 2^%f = 2^%f" %(mp.log(term1, 2), mp.log(term2, 2), mp.log(B_var, 2), mp.log(term1 + term2 + B_var, 2)))
      print("Failure Probability: 2^%s" % (self.failure_rate(term1 + term2 + B_var)))
    return term1 + term2 + B_var
  
  def CMUX(self, C_var, A_var, B_var, suppress=False):
    term1 = (self.k + 1)*self.ell*self.N*(self.beta**2)*C_var
    term2 = (1 + self.k*self.h_out/3)*self.epsilon**2
    if(self.verbose and not suppress):
      print("\n## CMUX ##")
      print("Variance <= Gaussian + Mod Switching + Variance(max(A,B))")
      print("Variance <= 2^%f + 2^%f + 2^%f = 2^%f" %(mp.log(term1, 2), mp.log(term2, 2), mp.log(max(A_var, B_var), 2), mp.log(term1 + term2 + max(A_var, B_var), 2)))
      print("Failure Probability: 2^%s" % (self.failure_rate(term1 + term2 + max(A_var, B_var))))
    return term1 + term2 + max(A_var, B_var)
  
  def auto(self, A_var, suppress=False):
    term1 = self.k*self.ell*self.N*(self.beta**2)*self.sigma_var
    term2 = (self.k*self.h_out/3)*self.epsilon**2
    return term1 + term2 + A_var

  def NCMUX(self, C_var, A_var, B_var, suppress=False):
    return self.CMUX(C_var, A_var, self.auto(B_var, suppress), suppress)

  def RGSW_monomial_mul(self, in_var, r_prec, suppress=False):
    res_var = in_var
    res_var = self.NCMUX(self.sigma_var, res_var, res_var, False)
    one_NCMUX = res_var
    one_CMUX = self.CMUX(self.sigma_var, res_var, res_var, True)
    for _ in range(r_prec - 1):
      res_var = self.CMUX(self.sigma_var, res_var, res_var, True)
    if(self.verbose and not suppress):
      print("\n## RGSW monomial mul ##")
      print("Variance <= NCMUX + (r_prec - 1) * CMUX(in_var)")
      print("Variance <= 2^%f + %d * 2^%f = 2^%f" % (mp.log(one_NCMUX,2), r_prec - 1, mp.log(one_CMUX,2), mp.log(res_var,2)))
      print("Failure Probability: 2^%s" % (self.failure_rate(res_var)))
    return res_var

  def sub_a(self, in_var, key_type, suppress=False):
    out_var = in_var
    if(key_type == "binary"): 
      out_var = self.GSW_product(self.sigma_var, in_var, suppress)
    if(key_type == "ternary"): 
      out_var = self.CMUX(self.sigma_var, in_var, in_var, suppress)
    if(key_type == "gaussian"):
      out_var =  self.auto(self.GSW_product(self.sigma_var, self.auto(in_var, suppress), suppress), suppress)
    return out_var

  def sparse_mul(self, h, r_prec, suppress=False):
    in_var = 0
    for _ in range(h):
      in_var = self.RGSW_monomial_mul(in_var, r_prec, suppress)
      in_var = self.sub_a(in_var, self.key_type, suppress)
    in_var = self.RGSW_monomial_mul(in_var, r_prec, suppress)
    if(self.verbose and not suppress):
      print("\n## Sparse Mul ##")
      print("Variance <= 2^%f = 2^%f" %(mp.log(in_var, 2), mp.log(in_var, 2)))
      print("Failure Probability: 2^%s" % (self.failure_rate(in_var)))
    return in_var
  
  def full_bootstrap(self, suppress=False, force_print=False):
    var_bs = self.sparse_mul(self.h, self.r_prec)
    var_bs_repack = self.packing_ks(var_bs)
    var_res = self.hw_reduction_ks(var_bs_repack)
    var_res_fbs = (((2**(-mp.log(2*self.N, 2))))**2)*(self.h/12)*self.key_base**2
    if(force_print or (self.verbose and not suppress)):
      print("\n## Full RLWE Bootstrapping ##")
      print("Failure Probability: 2^%s" % (self.failure_rate(var_res)))
      print("Variance <= Input var + Mod Switching")
      print("Variance <= 2^%f + 2^%f = 2^%f" % (mp.log(var_res,2), mp.log(var_res_fbs,2), mp.log(var_res + var_res_fbs,2)))
      print("Failure Probability FBS: 2^%s" % (self.failure_rate(var_res + var_res_fbs)))
    return var_res + var_res_fbs

  def print_bootstrap_fr(self, var):
    print("\n## Failure Probability FBS ##")
    print("Ptxt space\t\tFailure Rate")
    for msg in range(16):
      prob = self.failure_rate(var, msg_size=msg)
      print("\t%d\t\t2^%s" % (msg, prob))
      if(float(prob) > -10): 
        break
  
  def failure_rate(self, var, msg_size=None):
    msg_size = self.msg_prec if (msg_size == None) else msg_size
    interval = 2**(- msg_size - 1)
    sigma = mp.sqrt(var)
    return nstr((mp.log(1 - mp.erf(interval/(sigma*mp.sqrt(2))), 2)))
  
  def packing_ks(self, in_var, suppress=False):
    beta_packing = (1 << self.b_packing)/2
    epsilon_packing =  1 /(((1 << self.b_packing)**self.ell_packing))
    mul_noise = self.N * (self.N_in * (beta_packing**2) * self.ell_packing) * (self.sigma_packing**2)
    ms_noise = (self.h_out/12)*self.N_in*epsilon_packing**2
    if(self.verbose and not suppress):
      print("\n## Packing KS ##")
      print("Variance <= Gaussian + Mod Switching + Input var")
      print("Variance <= 2^%f + 2^%f + 2^%f = 2^%f" %(mp.log(mul_noise, 2), mp.log(ms_noise, 2), mp.log(in_var, 2), mp.log(mul_noise + ms_noise + in_var, 2)))
      print("Failure Probability: 2^%s" % (self.failure_rate(mul_noise + ms_noise + in_var)))
    return mul_noise + ms_noise + in_var

  def hw_reduction_ks(self, in_var, suppress=False):
    mul_noise = self.ell_ks*(((1<<self.b_ks)/2)**2)*self.N_in*(self.sigma_in)**2
    ms_noise = (self.h_packing/12)*((1 /(((1<<self.b_ks)**self.ell_ks)))**2)
    if(self.verbose and not suppress):
      print("\n## HW reduction KS ##")
      print("Variance <= Gaussian + Mod Switching + Input var")
      print("Variance <= 2^%f + 2^%f + 2^%f = 2^%f" %(mp.log(mul_noise, 2), mp.log(ms_noise, 2), mp.log(in_var, 2), mp.log(mul_noise + ms_noise + in_var, 2)))
      print("Failure Probability: 2^%s" % (self.failure_rate(mul_noise + ms_noise + in_var)))
    return mul_noise + ms_noise + in_var


sigma = 2**-51 #2.2148688116005568e-16
prec = 4


# main ones binary
# 2/3
scheme_0_1 = LeveledRGSW(2048, 12, 256, 512, 1, 2**-13, 2**-44, 2048, 2**-50, 1, ell=1, bg_bit=23, ell_packing=2, b_packing=14, ell_ks=9, b_ks=1, msg_prec=3, r_prec=9, verbose=True)

scheme_2_3 = LeveledRGSW(2048, 18, 256, 512, 1, 2**-15, 2**-44, 2048, 2**-50, 1, ell=1, bg_bit=23, ell_packing=2, b_packing=14, ell_ks=12, b_ks=1, msg_prec=3, r_prec=9, verbose=True)
# 4/5
scheme_4_5 = LeveledRGSW(2048, 21, 256, 512, 1, 2**-17, 2**-44, 2048, 2**-50, k=1, ell=1, bg_bit=23, ell_packing=2, b_packing=14, ell_ks=14, b_ks=1, msg_prec=5, r_prec=9, verbose=True) 
# 6/7
scheme_6_7 = LeveledRGSW(2048, 28, 256, 512, 1, 2**-20, 2**-44, 2048, 2**-50, k=1, ell=1, bg_bit=23, ell_packing=2, b_packing=14, ell_ks=16, b_ks=1, msg_prec=7, r_prec=8, verbose=True)
# 8/9
scheme_8_9_HFR = LeveledRGSW(2048, 29, 256, 512, 1, 2**-21, 2**-44, 4096, 2**-50, k=1, ell=1, bg_bit=23, ell_packing=2, b_packing=14, ell_ks=18, b_ks=1,  msg_prec=9, r_prec=8, verbose=True)

scheme_8_9_a = LeveledRGSW(2048, 32, 256, 512, 1, 2**-23, 2**-44, 8192, 2**-50, k=1, ell=1, bg_bit=22, ell_packing=2, b_packing=14, ell_ks=19, b_ks=1,  msg_prec=9, r_prec=9, verbose=True)
scheme_8_9_b = LeveledRGSW(2048, 33, 256, 512, 1, 2**-23, 2**-44, 8192, 2**-50, k=1, ell=1, bg_bit=22, ell_packing=2, b_packing=14, ell_ks=19, b_ks=1,  msg_prec=9, r_prec=9, verbose=True)

# main ones ternary
# 2/3
scheme_t_2_3_a = LeveledRGSW(2048, 16, 256, 512, 1, 2**-15, 2**-44, 2048, 2**-50, 1, ell=1, bg_bit=23, ell_packing=2, b_packing=14, ell_ks=10, b_ks=1, msg_prec=3, r_prec=8, verbose=True, key_type="ternary")
scheme_t_2_3_b = LeveledRGSW(2048, 17, 256, 512, 1, 2**-15, 2**-44, 2048, 2**-50, 1, ell=1, bg_bit=23, ell_packing=2, b_packing=14, ell_ks=12, b_ks=1, msg_prec=3, r_prec=8, verbose=True, key_type="ternary")
# 4/5
scheme_t_4_5 = LeveledRGSW(2048, 21, 256, 512, 1, 2**-17, 2**-44, 2048, 2**-50, k=1, ell=1, bg_bit=23, ell_packing=2, b_packing=14, ell_ks=14, b_ks=1, msg_prec=5, r_prec=9, verbose=True,key_type="ternary") 
# 6/7
scheme_t_6_7 = LeveledRGSW(2048, 27, 256, 512, 1, 2**-20, 2**-44, 2048, 2**-50, k=1, ell=1, bg_bit=23, ell_packing=2, b_packing=14, ell_ks=16, b_ks=1, msg_prec=7, r_prec=8, verbose=True, key_type="ternary")
# 8/9
scheme_t_8_9 = LeveledRGSW(2048, 30, 256, 512, 1, 2**-23, 2**-44, 8192, 2**-50, k=1, ell=1, bg_bit=22, ell_packing=2, b_packing=14, ell_ks=19, b_ks=1,  msg_prec=9, r_prec=9, verbose=True, key_type="ternary")

## increased dimension section

# 4096 -> 2048
scheme_L1 = LeveledRGSW(4096, 32, 256, 512, 1, 2**-23, 2**-44, 2048, 2**-50, k=1, ell=1, bg_bit=23, ell_packing=2, b_packing=14, ell_ks=16, b_ks=1, msg_prec=7, r_prec=8, verbose=True)
# 4096 -> 8192
scheme_L2 = LeveledRGSW(4096, 33, 256, 512, 1, 2**-24, 2**-44, 8192, 2**-50, k=1, ell=1, bg_bit=22, ell_packing=2, b_packing=14, ell_ks=19, b_ks=1,  msg_prec=9, r_prec=9, verbose=True)



# 8192 -> 2048
scheme_L3 = LeveledRGSW(8192, 32, 256, 512, 1, 2**-23, 2**-44, 2048, 2**-50, k=1, ell=1, bg_bit=23, ell_packing=2, b_packing=14, ell_ks=16, b_ks=1, msg_prec=7, r_prec=8, verbose=True)
# 8192 -> 8192
scheme_L4 = LeveledRGSW(8192, 32, 256, 512, 1, 2**-24, 2**-44, 8192, 2**-50, k=1, ell=1, bg_bit=22, ell_packing=2, b_packing=14, ell_ks=19, b_ks=1,  msg_prec=9, r_prec=9, verbose=True)


# 4096 -> 8192 ternary
scheme_t_L2 = LeveledRGSW(4096, 32, 256, 512, 1, 2**-24, 2**-44, 8192, 2**-50, k=1, ell=1, bg_bit=22, ell_packing=2, b_packing=14, ell_ks=19, b_ks=1,  msg_prec=9, r_prec=9, verbose=True, key_type="ternary")


## Sparse-generic section

scheme_A2 = LeveledRGSW(2048, 17, 256, 512, 4, 2**-15, 2**-44, 2048, 2**-50, k=1, ell=1, bg_bit=23, ell_packing=2, b_packing=14, ell_ks=12, b_ks=1, msg_prec=3, r_prec=9, verbose=True, key_type="gaussian")
scheme_A4 = LeveledRGSW(2048, 24, 256, 512, 4, 2**-18, 2**-44, 2048, 2**-50, k=1, ell=1, bg_bit=23, ell_packing=2, b_packing=14, ell_ks=15, b_ks=1, msg_prec=5, r_prec=8, verbose=True, key_type="gaussian")
scheme_A6 = LeveledRGSW(2048, 30, 256, 512, 4, 2**-23, 2**-44, 8192, 2**-50, k=1, ell=1, bg_bit=22, ell_packing=2, b_packing=14, ell_ks=19, b_ks=1, msg_prec=7, r_prec=8, verbose=True, key_type="gaussian")




# 4096 -> 2048
scheme_G2 = LeveledRGSW(4096, 32, 256, 512, 4, 2**-23, 2**-44, 2048, 2**-50, k=1, ell=1, bg_bit=23, ell_packing=2, b_packing=14, ell_ks=16, b_ks=1, msg_prec=7, r_prec=8, verbose=True, key_type="gaussian")
# 4096 -> 8192
scheme_G3 = LeveledRGSW(4096, 32, 256, 512, 4, 2**-23, 2**-44, 8192, 2**-50, k=1, ell=1, bg_bit=22, ell_packing=2, b_packing=14, ell_ks=18, b_ks=1,  msg_prec=9, r_prec=9, verbose=True, key_type="gaussian")



# 8192 -> 2048
scheme_G4 = LeveledRGSW(8192, 32, 256, 512, 4, 2**-23, 2**-44, 2048, 2**-50, k=1, ell=1, bg_bit=23, ell_packing=2, b_packing=14, ell_ks=16, b_ks=1, msg_prec=7, r_prec=8, verbose=True, key_type="gaussian")
# 8192 -> 8192
scheme_G5 = LeveledRGSW(8192, 32, 256, 512, 4, 2**-24, 2**-44, 8192, 2**-50, k=1, ell=1, bg_bit=22, ell_packing=2, b_packing=14, ell_ks=19, b_ks=1,  msg_prec=9, r_prec=9, verbose=True, key_type="gaussian")

scheme = scheme_t_8_9
bs_var = scheme.full_bootstrap(force_print=True)
scheme.print_bootstrap_fr(bs_var)
