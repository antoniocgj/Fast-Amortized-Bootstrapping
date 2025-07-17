#include <sab.h>

uint64_t get_min_prec(TRLWE_Key key){
  const uint64_t N = key->s[0]->N;
  uint64_t r_max = 0;
  for (size_t i = 0; i < key->k; i++){
    uint64_t previous = N;
    for (size_t j = 0; j < N; j++){
      const uint64_t coeff = key->s[i]->coeffs[N - j - 1];
      if(coeff != 0){
        const uint64_t r_diff = previous - (N - j - 1);
        if(r_max < r_diff) r_max = r_diff;
        previous = (N - j - 1);
      }
    }  
  }
  return (uint64_t)(log2(r_max) + 1);
}

bool check_key(TRLWE_Key key, uint64_t h, uint64_t r_prec, SAB_Key sab){
  const uint64_t N = key->s[0]->N;
  const uint64_t r_max = (1ULL<<r_prec);
  for (size_t i = 0; i < key->k; i++){
    uint64_t cnt_h = 0;
    uint64_t previous = N;
    for (size_t j = 0; j < N; j++){
      const uint64_t coeff = key->s[i]->coeffs[N - j - 1];
      if(coeff != 0){
        cnt_h++;
        const uint64_t r_diff = previous - (N - j - 1);
        if(r_diff >= r_max){
          if(sab->include_zeros) cnt_h += r_diff/r_max;
          else return false;
        }
        if(coeff != 1 && !sab->gaussian_secret && !sab->ternary_secret) return false;
        if(coeff != 1 && coeff != -1 && !sab->gaussian_secret) return false;
        previous = (N - j - 1);
      }
    }  
    if(cnt_h > h) return false;
    if(cnt_h < h && !sab->include_zeros) return false;
  }
  return true;
}

// rejection sampling for new key
uint64_t RS_sparse_binary_key(TRLWE_Key * output_key, int N, int k, int h, double sigma, uint64_t target_r_prec){
  uint64_t attempts = 0;
  for (size_t i = 0; i < 1<<15; i++){
    attempts++;
    *output_key = trlwe_new_sparse_binary_key(N, k, h, sigma);
    uint64_t r_prec = get_min_prec(*output_key);
    if(r_prec <= target_r_prec) break;
    free_trlwe_key(*output_key);
  }
  return attempts;
}

uint64_t RS_sparse_ternary_key(TRLWE_Key * output_key, int N, int k, int h, double sigma, uint64_t target_r_prec){
  uint64_t attempts = 0;
  for (size_t i = 0; i < 1<<15; i++){
    attempts++;
    *output_key = trlwe_new_ternary_key(N, k, h, sigma);
    uint64_t r_prec = get_min_prec(*output_key);
    if(r_prec <= target_r_prec) break;
    free_trlwe_key(*output_key);
  }
  return attempts;
}

uint64_t RS_sparse_arbitrary_key(TRLWE_Key * output_key, int N, int k, int h, uint64_t key_bound, double sigma, uint64_t target_r_prec){
  uint64_t attempts = 0;
  for (size_t i = 0; i < 128; i++){
    attempts++;
    *output_key = trlwe_new_sparse_generic_key(N, k, h, key_bound, sigma);
    uint64_t r_prec = get_min_prec(*output_key);
    if(r_prec <= target_r_prec) break;
    free_trlwe_key(*output_key);
  }
  return attempts;
}

void RGSW_encrypt_bits(TRGSW_DFT * out, TRGSW tmp, TRGSW_Key key, uint64_t in, uint64_t prec){
  for (size_t b = 0; b < prec; b++){
    const uint64_t val = (in >> b)&1;
    trgsw_monomial_sample(tmp, val, 0, key);
    trgsw_to_DFT(out[b], tmp);
  }
}

void RGSW_encrypt(TRGSW_DFT out, TRGSW tmp, TRGSW_Key key, uint64_t c, uint64_t e){
  trgsw_monomial_sample(tmp, c, e, key);
  trgsw_to_DFT(out, tmp);
}

TLWE_Key __gbl_extracted_key;
TRLWE_Key __gbl_rlwe_key;
TRLWE_Key __gbl_rlwe_key_in;
TRLWE_Key __gbl_rlwe_key_packing;

// #define MEASURE_NOISE
#ifdef MEASURE_NOISE
#define DECRYPTION_FUNCTION(X) tlwe_phase(X, __gbl_extracted_key)
#include <noise_util.h>
#endif

SAB_Key new_sparse_amortized_bootstrapping(TRLWE_Key input_key, TRLWE_Key repacking_key, TRGSW_Key output_key, uint64_t b_prec, uint64_t b_packing, uint64_t ell_packing, uint64_t t_ks, uint64_t b_ks, uint64_t h, uint64_t r_prec, bool include_zeros, bool ternary, bool gaussian){
  SAB_Key res = (SAB_Key) safe_malloc(sizeof(*res));
  const uint64_t in_N = input_key->s[0]->N;
  const uint64_t in_k = input_key->k;
  const uint64_t out_N = output_key->trlwe_key->s[0]->N;
  const uint64_t out_k = output_key->trlwe_key->k;
  const uint64_t r_max = (1ULL<<r_prec);
  uint64_t m1[1] = {2*out_N - 1};
  TRGSW tmp = trgsw_alloc_new_sample(output_key->l, output_key->Bg_bit, out_k, out_N);
  // auto key for -1
  TRLWE_KS_Key * aut_ks = trlwe_new_automorphism_KS_keyset_2(output_key->trlwe_key, m1, 1, output_key->l, output_key->Bg_bit);
  res->aut_minus1 = aut_ks[0];
  free(aut_ks);

  // auto keys for gaussian secrets
  if(gaussian) res->aut_ksk = trlwe_new_automorphism_KS_keyset(output_key->trlwe_key, true, output_key->l, output_key->Bg_bit);
  // include zeros or ternary
  if(include_zeros || gaussian) res->s_coff = (TRGSW_DFT **) safe_malloc(sizeof(TRGSW_DFT *)*in_k);
  if(ternary) res->s_sign = (TRGSW_DFT **) safe_malloc(sizeof(TRGSW_DFT *)*in_k);
  // repacking
  TLWE_Key extracted_key = tlwe_alloc_key(out_N*out_k, output_key->trlwe_key->sigma);
  trlwe_extract_tlwe_key(extracted_key, output_key->trlwe_key);
  __gbl_rlwe_key = output_key->trlwe_key;
  __gbl_extracted_key = extracted_key;
  __gbl_rlwe_key_in = input_key;
  __gbl_rlwe_key_packing = repacking_key;
  res->packing_key = trlwe_new_full_packing_KS_key(repacking_key, extracted_key, ell_packing, b_packing);
  // free_tlwe_key(extracted_key);
  // HW reducing KS
  res->hw_reducing_key = trlwe_new_KS_key(input_key, repacking_key, t_ks, b_ks);

  // encrypt sparse representation
  uint64_t previous = in_N;
  res->s = (TRGSW_DFT ***) safe_malloc(sizeof(TRGSW_DFT ***)*in_k);
  for (size_t i = 0; i < in_k; i++){    
    res->s[i] = (TRGSW_DFT **) safe_malloc(sizeof(TRGSW_DFT **)*(h+1));
    uint64_t cnt_h = 0;
    if(ternary) res->s_sign[i] = trgsw_alloc_new_DFT_sample_array(h, output_key->l, output_key->Bg_bit, out_k, out_N);
    if(include_zeros || gaussian) res->s_coff[i] = trgsw_alloc_new_DFT_sample_array(h, output_key->l, output_key->Bg_bit, out_k, out_N);
    for (size_t j = 0; j < in_N; j++){
      const uint64_t coeff = input_key->s[i]->coeffs[in_N - j - 1];
      if(coeff != 0){
        const uint64_t current = in_N - j - 1;
        uint64_t r_diff = previous - current;
        // encrypt zero coeffs (if needed)
        while (r_diff >= r_max){
          assert(false);
          res->s[i][cnt_h] = trgsw_alloc_new_DFT_sample_array(r_prec, output_key->l, output_key->Bg_bit, out_k, out_N);
          RGSW_encrypt_bits(res->s[i][cnt_h], tmp, output_key, (r_max - 1), r_prec);
          if(gaussian) RGSW_encrypt(res->s_coff[i][cnt_h], tmp, output_key, 1, 0);  // 1*x^0
          else if(include_zeros) RGSW_encrypt(res->s_coff[i][cnt_h], tmp, output_key, 0, 0); // 0
          else assert(false); // 0 valued coefficient without saving coefficients 
          if(ternary) RGSW_encrypt(res->s_sign[i][cnt_h], tmp, output_key, 0, 0); // 0
          r_diff -= (r_max - 1);
          cnt_h++;
        }
        // encrypt non zero coeff
        res->s[i][cnt_h] = trgsw_alloc_new_DFT_sample_array(r_prec, output_key->l, output_key->Bg_bit, out_k, out_N);
        RGSW_encrypt_bits(res->s[i][cnt_h], tmp, output_key, r_diff, r_prec);
        if(gaussian){
          RGSW_encrypt(res->s_coff[i][cnt_h], tmp, output_key, 1, coeff); // 1*x^coeff
        }else{
          if(include_zeros) RGSW_encrypt(res->s_coff[i][cnt_h], tmp, output_key, 1, 0); // 1
          if(ternary) RGSW_encrypt(res->s_sign[i][cnt_h], tmp, output_key, coeff == -1, 0); // 1 if negative, 0 otherwise
        }
        previous = current;
        cnt_h++;
      }
    }
    assert(cnt_h == h);
    res->s[i][cnt_h] = trgsw_alloc_new_DFT_sample_array(r_prec, output_key->l, output_key->Bg_bit, out_k, out_N);
    RGSW_encrypt_bits(res->s[i][cnt_h], tmp, output_key, previous, r_prec);
  }
  // save values
  res->gaussian_secret = gaussian;
  res->include_zeros = include_zeros;
  res->ternary_secret = ternary;
  res->in_N = in_N;
  res->in_k = in_k;
  res->out_N = out_N;
  res->out_k = out_k;
  res->h = h;
  res->r_prec = r_prec;
  res->b_prec = b_prec;
  // alloc tmps
  res->tmp = (tmp_pool) safe_malloc(sizeof(*res->tmp));
  res->tmp->rlwe_dft = trlwe_alloc_new_DFT_sample(out_k, out_N);
  res->tmp->rlwe = trlwe_alloc_new_sample(out_k, out_N);
  res->tmp->rlwe_in = trlwe_alloc_new_sample(in_k, in_N);
  res->tmp->rlwe_poly1 = trlwe_alloc_new_sample_array(in_N, out_k, out_N);
  res->tmp->rlwe_poly2 = trlwe_alloc_new_sample_array(in_N, out_k, out_N);
  res->tmp->extracted_poly = tlwe_alloc_sample_array(in_N, out_N*out_k);
  res->tmp->rgsw = tmp;
  return res;
}

SAB_Key copy_SAB_key(SAB_Key sab_key){
  assert(false); // not implemented
  return sab_key;
}

void CMUX(TRLWE out, TRLWE in1, TRLWE in2, TRGSW_DFT selector, SAB_Key sab){
  trlwe_sub(sab->tmp->rlwe, in2, in1); // B - A
  trgsw_mul_trlwe_DFT(sab->tmp->rlwe_dft, sab->tmp->rlwe, selector);  // S(B - A)
  trlwe_from_DFT(sab->tmp->rlwe, sab->tmp->rlwe_dft); 
  trlwe_add(out, sab->tmp->rlwe, in1); // S(B - A) + A 
}

void NCMUX(TRLWE out, TRLWE in1, TRLWE in2, TRGSW_DFT selector, SAB_Key sab){
  trlwe_eval_automorphism(sab->tmp->rlwe, in2, 2*in2->b->N - 1, sab->aut_minus1);
  CMUX(out, in1, sab->tmp->rlwe, selector, sab);
}

// p0 <- p0 * X^e
void RGSW_monomial_mul(TRLWE * p0, TRGSW_DFT * e, SAB_Key sab){
  const uint32_t r_prec = sab->r_prec, in_N = sab->in_N;  
  TRLWE * p[2] = {p0, sab->tmp->rlwe_poly2};
  for (size_t i = 0; i < r_prec; i++){
    const uint64_t power = 1ULL << i;
    const uint64_t out = (i+1)&1, in = out^1;
    for (size_t j = 0; j < power; j++){
      NCMUX(p[out][j], p[in][j], p[in][in_N - power + j], e[i], sab);
    }
    for (size_t j = 0; j < in_N - power; j++){
      CMUX(p[out][j + power], p[in][j + power], p[in][j], e[i], sab);
    }
  }
  if(p[r_prec&1]!=p0){
    for (size_t i = 0; i < in_N; i++){
      trlwe_copy(p0[i], p[r_prec&1][i]);
    }
  }
}

void sub_a_ga(TRLWE * p, uint64_t * a, uint64_t key_idx, SAB_Key sab){
  for (size_t i = 0; i < sab->in_N; i++){
    assert(a[i] != 0);
    const uint64_t w_inv = inverse_mod_2N(a[i], sab->out_N);
    assert(w_inv != 0);
    trlwe_eval_automorphism(sab->tmp->rlwe, p[i], w_inv, sab->aut_ksk[(w_inv - 1)>>1]);
    trgsw_mul_trlwe_DFT(sab->tmp->rlwe_dft, sab->tmp->rlwe, sab->s_coff[0][key_idx]);
    trlwe_from_DFT(sab->tmp->rlwe, sab->tmp->rlwe_dft);
    trlwe_eval_automorphism(p[i], sab->tmp->rlwe, a[i], sab->aut_ksk[(a[i] - 1)>>1]);
  }
}

void sub_a(TRLWE * p, uint64_t * a, uint64_t key_idx, SAB_Key sab){
  if(sab->gaussian_secret) return sub_a_ga(p, a, key_idx, sab);
  for (size_t i = 0; i < sab->in_N; i++){
    if(sab->include_zeros){
      trlwe_mul_by_xai_minus_1(sab->tmp->rlwe, p[i], a[i]);
      trgsw_mul_trlwe_DFT(sab->tmp->rlwe_dft, sab->tmp->rlwe, sab->s_coff[0][key_idx]);
      trlwe_from_DFT(sab->tmp->rlwe, sab->tmp->rlwe_dft);
      trlwe_addto(p[i], sab->tmp->rlwe);
    }else if (sab->ternary_secret){
      trlwe_mul_by_xai(sab->tmp->rlwe, p[i], a[i]);
      trlwe_copy(p[i], sab->tmp->rlwe);
      trlwe_mul_by_xai_minus_1(sab->tmp->rlwe, p[i], -2*a[i]);
      trgsw_mul_trlwe_DFT(sab->tmp->rlwe_dft, sab->tmp->rlwe, sab->s_sign[0][key_idx]);
      trlwe_from_DFT(sab->tmp->rlwe, sab->tmp->rlwe_dft);
      trlwe_addto(p[i], sab->tmp->rlwe);
    }else{
      trlwe_mul_by_xai(sab->tmp->rlwe, p[i], a[i]);
      trlwe_copy(p[i], sab->tmp->rlwe);
    }
  }
}

// p = p * x^{-as}
void sparse_mul(TRLWE * p, uint64_t * a, uint64_t a_idx, SAB_Key sab){ 
  #ifdef MEASURE_NOISE
  TorusPolynomial __debug_poly = polynomial_new_torus_polynomial(p[0]->b->N);
  #endif
  for (size_t i = 0; i < sab->h; i++){
    RGSW_monomial_mul(p, sab->s[a_idx][i], sab);
    #ifdef MEASURE_NOISE
    for (size_t j = 0; j < sab->in_N; j++){
      trlwe_phase(__debug_poly, p[j], __gbl_rlwe_key);
      ARRAY_LOG_BY_PREC(__debug_poly->coeffs, sab->b_prec, __debug_poly->N);
    }
    printf("Sparse mul -- h: %ld - ", i);
    print_reset_noise();
    #endif
    sub_a(p, a, i, sab);
  }
  RGSW_monomial_mul(p, sab->s[a_idx][sab->h], sab);
}

uint64_t odd_mod_switch(uint64_t in, uint64_t prec){
  const uint64_t bit_size = sizeof(uint64_t) * 8;
  const uint64_t round_offset = 1UL << (bit_size - prec - 1);
  const uint64_t val = (in + round_offset)>>(bit_size - prec);
  if(val&1){
    return val;
  }else{
    return (val - 1)&((1ULL<<prec) - 1);
  }
}

void mod_switch_a(uint64_t * out, uint64_t * in, uint64_t prec, uint64_t size, bool round_to_odd){
  for (size_t j = 0; j < size; j++){
    if(round_to_odd) out[j] = odd_mod_switch(in[j], prec);
    else out[j] = torus2int(in[j], prec);
  } 
}


void sab_blind_rotate(TRLWE * out, TRLWE in, SAB_Key sab){
  uint64_t * a = (uint64_t *) safe_malloc(sizeof(uint64_t)*sab->in_N); 
  const uint64_t log_N2 = (uint64_t) log2(2*sab->out_N);
  assert(sab->in_k == 1); // TODO: add mul X^N for k > 1
  for (size_t i = 0; i < sab->in_k; i++){
    // mod switch a
    mod_switch_a(a, in->a[i]->coeffs, log_N2, sab->in_N, sab->gaussian_secret);
    // compute -s[i]*a[i]
    sparse_mul(out, a, i, sab);
  }
  free(a);
}

TRLWE * setup_single_tv(uint64_t * b, TRLWE tv, SAB_Key sab){
  const int N = sab->out_N, log_N2 = (int) log2(N*2);
  const Torus prec_offset = double2torus(1./(2*sab->b_prec));
  TRLWE * res = trlwe_alloc_new_sample_array(sab->in_N, sab->out_k, sab->out_N);
  for (size_t i = 0; i < sab->in_N; i++){
    trlwe_mul_by_xai(res[i], tv, torus2int(b[i] + prec_offset, log_N2));
  }
  return res;
}

// compute acc = tv * X^b * X^N 
void setup_tv_xb(TRLWE * acc, uint64_t * b, TRLWE tv, SAB_Key sab){
  const int N = sab->out_N, log_N2 = (int) log2(N*2);
  const uint64_t prec_offset = 1ULL << (64 - sab->b_prec - 1);
  for (size_t i = 0; i < sab->in_N; i++){
    trlwe_mul_by_xai(acc[i], tv, torus2int(b[i] + prec_offset, log_N2));
  }
}

void sab_rlwe_bootstrap_wo_extract(TRLWE * out, TRLWE in, TRLWE tv, SAB_Key sab){
  setup_tv_xb(out, in->b->coeffs, tv, sab);
  sab_blind_rotate(out, in, sab);
}

void sab_rlwe_to_lwe_bootstrap(TLWE * out, TRLWE in, TRLWE tv, SAB_Key sab){
  sab_rlwe_bootstrap_wo_extract(sab->tmp->rlwe_poly1, in, tv, sab);
  for (size_t i = 0; i < sab->in_N; i++){
    trlwe_extract_tlwe(out[i], sab->tmp->rlwe_poly1[i], 0);
  }
}

void sab_rlwe_bootstrap(TRLWE out, TRLWE in, TRLWE tv, SAB_Key sab){
  #ifdef MEASURE_NOISE
  TorusPolynomial __debug_poly = polynomial_new_torus_polynomial(in->b->N);
  trlwe_phase(__debug_poly, in, __gbl_rlwe_key_in);
  ARRAY_LOG_BY_PREC(__debug_poly->coeffs, sab->b_prec, __debug_poly->N);
  printf("Bootstrapping input -- ");
  print_reset_noise();
  #endif
  sab_rlwe_bootstrap_wo_extract(sab->tmp->rlwe_poly1, in, tv, sab);
  #ifdef MEASURE_NOISE
  reset_noise();
  #endif
  for (size_t i = 0; i < sab->in_N; i++){
    trlwe_extract_tlwe(sab->tmp->extracted_poly[i], sab->tmp->rlwe_poly1[i], 0);
    #ifdef MEASURE_NOISE
    const uint64_t mod_mask = (1ULL<<(sab->b_prec - 1)) - 1;
    // if(torus2int(DECRYPTION_FUNCTION(sab->tmp->extracted_poly[i]), sab->b_prec) != (i&mod_mask)){
    //   printf("%lu: %lu != %lu\n", i, torus2int(DECRYPTION_FUNCTION(sab->tmp->extracted_poly[i]), sab->b_prec), i&mod_mask);
    // }
    DECRYPT_LOG_BY_PREC(sab->tmp->extracted_poly[i], sab->b_prec);
    #endif
  }
  #ifdef MEASURE_NOISE
  printf("After bootstrapping (LWE noise) -- ");
  print_reset_noise();
  #endif
  trlwe_full_packing_keyswitch(sab->tmp->rlwe_in, sab->tmp->extracted_poly, sab->in_N, sab->packing_key);
  #ifdef MEASURE_NOISE
  trlwe_phase(__debug_poly, sab->tmp->rlwe_in, __gbl_rlwe_key_packing);
  ARRAY_LOG_BY_PREC(__debug_poly->coeffs, sab->b_prec, __debug_poly->N);
  printf("After repacking -- ");
  print_reset_noise();
  #endif
  trlwe_keyswitch(out, sab->tmp->rlwe_in, sab->hw_reducing_key);
  #ifdef MEASURE_NOISE
  trlwe_phase(__debug_poly, out, __gbl_rlwe_key_in);
  ARRAY_LOG_BY_PREC(__debug_poly->coeffs, sab->b_prec, __debug_poly->N);
  printf("After HW KS -- ");
  print_reset_noise();
  #endif
}

void sab_LUT_packing(TRLWE out, uint64_t * in, SAB_Key sab){
  const uint64_t size = 1ULL << (sab->b_prec - 1);
  trlwe_noiseless_trivial_sample(out, 0);
  out->b->coeffs[0] = int2torus(in[0], sab->b_prec);
  for (size_t i = 1; i < out->b->N; i++){
    out->b->coeffs[out->b->N - i] = -int2torus(in[i/(out->b->N/size)], sab->b_prec);
  }
}