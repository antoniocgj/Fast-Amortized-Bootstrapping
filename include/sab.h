#include <mosfhet.h>

typedef struct _tmp_pool{
  TRLWE_DFT rlwe_dft;
  TRLWE rlwe, rlwe_in, * rlwe_poly1, * rlwe_poly2;
  TLWE * extracted_poly;
  TRGSW rgsw;
} * tmp_pool;

typedef struct _SAB_Key{
  uint64_t in_N, in_k, out_N, out_k, h, r_prec, b_prec;
  bool include_zeros, gaussian_secret, ternary_secret;
  TRLWE_KS_Key aut_minus1;
  TRLWE_KS_Key * aut_ksk;
  TRGSW_DFT *** s;
  TRGSW_DFT ** s_coff;
  TRGSW_DFT ** s_sign;
  TRLWE_KS_Key packing_key;
  TRLWE_KS_Key hw_reducing_key;
  tmp_pool tmp;
} * SAB_Key;

uint64_t RS_sparse_binary_key(TRLWE_Key * output_key, int N, int k, int h, double sigma, uint64_t target_r_prec);
uint64_t RS_sparse_ternary_key(TRLWE_Key * output_key, int N, int k, int h, double sigma, uint64_t target_r_prec);
uint64_t RS_sparse_arbitrary_key(TRLWE_Key * output_key, int N, int k, int h, uint64_t key_bound, double sigma, uint64_t target_r_prec);
  

bool check_key(TRLWE_Key key, uint64_t h, uint64_t r_prec, SAB_Key sab);
void RGSW_encrypt_bits(TRGSW_DFT * out, TRGSW tmp, TRGSW_Key key, uint64_t in, uint64_t prec);
void RGSW_encrypt(TRGSW_DFT out, TRGSW tmp, TRGSW_Key key, uint64_t c, uint64_t e);
SAB_Key new_sparse_amortized_bootstrapping(TRLWE_Key input_key, TRLWE_Key repacking_key, TRGSW_Key output_key, uint64_t b_prec, uint64_t b_packing, uint64_t ell_packing, uint64_t t_ks, uint64_t b_ks, uint64_t h, uint64_t r_prec, bool include_zeros, bool ternary, bool gaussian);
void NCMUX(TRLWE out, TRLWE in1, TRLWE in2, TRGSW_DFT selector, SAB_Key sab);
void RGSW_monomial_mul(TRLWE * p0, TRGSW_DFT * e, SAB_Key sab);
void sub_a_ga(TRLWE * p, uint64_t * a, uint64_t key_idx, SAB_Key sab);
void sub_a(TRLWE * p, uint64_t * a, uint64_t key_idx, SAB_Key sab);
void sparse_mul(TRLWE * p, uint64_t * a, uint64_t a_idx, SAB_Key sab);  
void sab_blind_rotate(TRLWE * out, TRLWE in, SAB_Key sab);
void sab_rlwe_bootstrap_wo_extract(TRLWE * out, TRLWE in, TRLWE tv, SAB_Key sab);
uint64_t get_min_prec(TRLWE_Key key);
TRLWE * setup_single_tv(uint64_t * b, TRLWE tv, SAB_Key sab);
void sab_rlwe_to_lwe_bootstrap(TLWE * out, TRLWE in, TRLWE tv, SAB_Key sab);
void sab_rlwe_bootstrap(TRLWE out, TRLWE in, TRLWE tv, SAB_Key sab);
void sab_LUT_packing(TRLWE out, uint64_t * in, SAB_Key sab);