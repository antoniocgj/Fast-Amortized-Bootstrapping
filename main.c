#include <sab.h>
#include <benchmark_util.h>

// #define PRINT_POLY

void tlwe_print(TLWE c, TLWE_Key key, uint64_t prec){
  printf("%lu", torus2int(tlwe_phase(c, key), prec));
}

void trlwe_print(TRLWE c, TRLWE_Key key, uint64_t prec){
  TorusPolynomial tmp = polynomial_new_torus_polynomial(c->b->N);
  trlwe_phase(tmp, c, key);
#ifdef PRINT_POLY
  for (int64_t i = c->b->N - 1; i >= 0; i--){
    const uint64_t val = torus2int(tmp->coeffs[i], prec);
    if(val){
      printf("+ %lux^%ld ", val, i);
    }
  }
#else
  for (int64_t i = 0; i < c->b->N; i++){
    const uint64_t val = torus2int(tmp->coeffs[i], prec);
    printf("%lu, ", val);
  }
#endif
  printf("\n");
  free_polynomial(tmp);
}

void print_exp_poly(TRLWE * c, TRLWE_Key key, uint64_t n, uint64_t prec){
  TorusPolynomial tmp = polynomial_new_torus_polynomial(c[0]->b->N);
  for (size_t i = 0; i < n; i++){
    int64_t idx = -1;
    int64_t sign = 1;
    trlwe_phase(tmp, c[i], key);
    for (size_t j = 0; j < c[i]->b->N; j++){
      const uint64_t val = torus2int(tmp->coeffs[j], prec);
      if(val){
        if(idx != -1){
          printf("Decryption failed\n");
        }
        if(val > 1ULL << (prec - 1)) sign = -1;
        else sign = 1;
        idx = j;
      }
    }
    if(idx == -1){
      printf("Decryption failed\n");
    }
    #ifdef PRINT_POLY
    printf("+ %ldx^%ld ", idx*sign, i);
    #else
    printf("%ld, ", idx*sign);
    #endif
  }
  printf("\n");
  free_polynomial(tmp);
}

void test_NCMUX(){
  const uint64_t in_N = 1024, in_k = 1, out_N = 2048, out_k = 1, l = 1, bg_bit = 23, b_packing = 12, ell_packing = 2, t_aut = l, b_aut = bg_bit, h_in = 64, msg_prec = 4;
  TRLWE_Key input_key = trlwe_new_sparse_binary_key(in_N, in_k, h_in, pow(2, -20));
  TRLWE_Key out_key = trlwe_new_sparse_binary_key(out_N, out_k, h_in, pow(2, -53));
  TRGSW_Key output_key = trgsw_new_key(out_key, l, bg_bit);
  const uint64_t r_prec = get_min_prec(input_key);
  printf("Min precision: %lu\n", r_prec);
  SAB_Key sab = new_sparse_amortized_bootstrapping(input_key, out_key, output_key, msg_prec, b_packing, ell_packing, t_aut, b_aut, h_in, r_prec, false, false, false);
  TRLWE * rlwe = trlwe_alloc_new_sample_array(4, out_k, out_N);
  TRGSW_DFT * sel = trgsw_alloc_new_DFT_sample_array(2, l, bg_bit, out_k, out_N);
  trgsw_monomial_DFT_sample(sel[0], 0, 0, output_key);
  trgsw_monomial_DFT_sample(sel[1], 1, 0, output_key);
  trlwe_noiseless_trivial_sample(rlwe[1], NULL);
  trlwe_noiseless_trivial_sample(rlwe[2], NULL);
  rlwe[1]->b->coeffs[1] = int2torus(7, 4);
  rlwe[2]->b->coeffs[1] = int2torus(7, 4);
  NCMUX(rlwe[0], rlwe[1], rlwe[2], sel[0], sab);
  NCMUX(rlwe[3], rlwe[1], rlwe[2], sel[1], sab);
  trlwe_print(rlwe[0], out_key, 4);
  trlwe_print(rlwe[3], out_key, 4);
}

void test_monomial_mul(){
  const uint64_t in_N = 1024, in_k = 1, out_N = 2048, out_k = 1, l = 1, bg_bit = 23, b_packing = 12, ell_packing = 2, t_aut = l, b_aut = bg_bit, h_in = 64, msg_prec = 4;
  TRLWE_Key input_key = trlwe_new_sparse_binary_key(in_N, in_k, h_in, pow(2, -20));
  TRLWE_Key out_key = trlwe_new_sparse_binary_key(out_N, out_k, h_in, pow(2, -53));
  TRGSW_Key output_key = trgsw_new_key(out_key, l, bg_bit);
  const uint64_t r_prec = get_min_prec(input_key);
  printf("Min precision: %lu\n", r_prec);
  SAB_Key sab = new_sparse_amortized_bootstrapping(input_key, out_key, output_key, msg_prec, b_packing, ell_packing, t_aut, b_aut, h_in, r_prec, false, false, false);

  TRGSW_DFT * sel = trgsw_alloc_new_DFT_sample_array(r_prec, l, bg_bit, out_k, out_N);
  
  TRLWE rlwe_in = trlwe_new_sample(NULL, input_key);
  TRLWE rlwe_tv = trlwe_new_noiseless_trivial_sample(NULL, out_k, out_N);
  rlwe_tv->b->coeffs[1] += int2torus(1, 4);


  TRLWE * rlwe_acc = setup_single_tv(rlwe_in->b->coeffs, rlwe_tv, sab);
  
  print_exp_poly(rlwe_acc, out_key, sab->in_N, 4);
  RGSW_encrypt_bits(sel, sab->tmp->rgsw, output_key, 10, r_prec);
  RGSW_monomial_mul(rlwe_acc, sel, sab);
  print_exp_poly(rlwe_acc, out_key, sab->in_N, 4);
}

void test_sab_br(){
  const uint64_t in_N = 1024, in_k = 1, out_N = 2048, out_k = 1, l = 1, bg_bit = 23, b_packing = 12, ell_packing = 2, t_aut = l, b_aut = bg_bit, h_in = 64, msg_prec = 4;
  TRLWE_Key input_key = trlwe_new_sparse_binary_key(in_N, in_k, h_in, pow(2, -20));
  TRLWE_Key out_key = trlwe_new_sparse_binary_key(out_N, out_k, h_in, pow(2, -53));
  TRGSW_Key output_key = trgsw_new_key(out_key, l, bg_bit);
  const uint64_t r_prec = get_min_prec(input_key);
  printf("Min precision: %lu\n", r_prec);
  SAB_Key sab = new_sparse_amortized_bootstrapping(input_key, out_key, output_key, msg_prec, b_packing, ell_packing, t_aut, b_aut, h_in, r_prec, false, false, false);
  
  TorusPolynomial poly_in = polynomial_new_torus_polynomial(in_N);

  const uint64_t mod_mask = (1ULL<<(msg_prec - 1)) - 1;
  for (size_t i = 0; i < in_N; i++) poly_in->coeffs[i] = int2torus(i&mod_mask, msg_prec);

  TRLWE rlwe_in = trlwe_new_sample(poly_in, input_key);
  TRLWE rlwe_tv = trlwe_new_noiseless_trivial_sample(NULL, out_k, out_N);
  rlwe_tv->b->coeffs[1] += int2torus(1, 4);


  TRLWE * rlwe_acc = trlwe_alloc_new_sample_array(sab->in_N, sab->out_k, sab->out_N);
  MEASURE_TIME("", 10, "SAB", 
    sab_rlwe_bootstrap_wo_extract(rlwe_acc, rlwe_in, rlwe_tv, sab);
  );

  print_exp_poly(rlwe_acc, out_key, sab->in_N, 4);
}

void test_sab(){
  const uint64_t reps = 3;
#if defined(SET_2_3_2048)
  const uint64_t in_N = 2048, in_k = 1, out_N = 2048, out_k = 1, l = 1, bg_bit = 23, b_packing = 14, ell_packing = 2, t_ks = 12, b_ks = 1, h_in = 39, h_out = 512, msg_prec = 3;
  const double sigma_in = pow(2, -15);
  const double sigma_out = pow(2, -50);
  const uint64_t target_r_prec = 7;
#elif defined(SET_2_3_4096)
  const uint64_t in_N = 4096, in_k = 1, out_N = 2048, out_k = 1, l = 1, bg_bit = 23, b_packing = 14, ell_packing = 2, t_ks = 12, b_ks = 1, h_in = 32, h_out = 512, msg_prec = 3;
  const double sigma_in = pow(2, -15);
  const double sigma_out = pow(2, -50);
  const uint64_t target_r_prec = 8; 
#elif defined(SET_2_3_8192)
  const uint64_t in_N = 8192, in_k = 1, out_N = 2048, out_k = 1, l = 1, bg_bit = 23, b_packing = 14, ell_packing = 2, t_ks = 12, b_ks = 1, h_in = 25, h_out = 512, msg_prec = 3;
  const double sigma_in = pow(2, -15);
  const double sigma_out = pow(2, -50);
  const uint64_t target_r_prec = 10; 
#elif defined(SET_4_5_2048)
  const uint64_t in_N = 2048, in_k = 1, out_N = 2048, out_k = 1, l = 1, bg_bit = 23, b_packing = 14, ell_packing = 2, t_ks = 14, b_ks = 1, h_in = 42, h_out = 512, msg_prec = 5;
  const double sigma_in = pow(2, -17);
  const double sigma_out = pow(2, -50);
  const uint64_t target_r_prec = 7;
#elif defined(SET_4_5_4096)
  const uint64_t in_N = 4096, in_k = 1, out_N = 2048, out_k = 1, l = 1, bg_bit = 23, b_packing = 14, ell_packing = 2, t_ks = 14, b_ks = 1, h_in = 34, h_out = 512, msg_prec = 5;
  const double sigma_in = pow(2, -18);
  const double sigma_out = pow(2, -50);
  const uint64_t target_r_prec = 8;
#elif defined(SET_4_5_8192)
  const uint64_t in_N = 8192, in_k = 1, out_N = 2048, out_k = 1, l = 1, bg_bit = 23, b_packing = 14, ell_packing = 2, t_ks = 14, b_ks = 1, h_in = 26, h_out = 512, msg_prec = 5;
  const double sigma_in = pow(2, -18);
  const double sigma_out = pow(2, -50);
  const uint64_t target_r_prec = 10;
#elif defined(SET_6_7_4096)
  const uint64_t in_N = 4096, in_k = 1, out_N = 2048, out_k = 1, l = 1, bg_bit = 23, b_packing = 14, ell_packing = 2, t_ks = 17, b_ks = 1, h_in = 33, h_out = 512, msg_prec = 7;
  const double sigma_in = pow(2, -21);
  const double sigma_out = pow(2, -50);
  const uint64_t target_r_prec = 9;
#elif defined(SET_6_7_8192)
  const uint64_t in_N = 8192, in_k = 1, out_N = 2048, out_k = 1, l = 1, bg_bit = 23, b_packing = 14, ell_packing = 2, t_ks = 17, b_ks = 1, h_in = 27, h_out = 512, msg_prec = 7;
  const double sigma_in = pow(2, -21);
  const double sigma_out = pow(2, -50);
  const uint64_t target_r_prec = 10;
#elif defined(SET_8_9_4096)
  const uint64_t in_N = 4096, in_k = 1, out_N = 8192, out_k = 1, l = 1, bg_bit = 22, b_packing = 14, ell_packing = 2, t_ks = 20, b_ks = 1, h_in = 34, h_out = 512, msg_prec = 9;
  const double sigma_in = pow(2, -24);
  const double sigma_out = pow(2, -51);
  const uint64_t target_r_prec = 9;
#elif defined(SET_8_9_8192)
  const uint64_t in_N = 8192, in_k = 1, out_N = 8192, out_k = 1, l = 1, bg_bit = 22, b_packing = 14, ell_packing = 2, t_ks = 20, b_ks = 1, h_in = 28, h_out = 512, msg_prec = 9;
  const double sigma_in = pow(2, -24);
  const double sigma_out = pow(2, -51);
  const uint64_t target_r_prec = 10;
#elif defined(SET_8_9_HIGH_FR)
  const uint64_t in_N = 8192, in_k = 1, out_N = 4096, out_k = 1, l = 1, bg_bit = 23, b_packing = 14, ell_packing = 2, t_ks = 17, b_ks = 1, h_in = 28, h_out = 512, msg_prec = 9;
  const double sigma_in = pow(2, -22);
  const double sigma_out = pow(2, -50);
  const uint64_t target_r_prec = 10;
#else
  const uint64_t in_N = 2048, in_k = 1, out_N = 2048, out_k = 1, l = 1, bg_bit = 23, b_packing = 14, ell_packing = 2, t_ks = 12, b_ks = 1, h_in = 39, h_out = 512, msg_prec = 3;
  const double sigma_in = pow(2, -15);
  const double sigma_out = pow(2, -50);
  const uint64_t target_r_prec = 7;
#endif
  printf("Sparse bootstrapping with binary keys\n");
  printf("Input: (N=%ld, h=%ld, binary, σ=2^%ld)\n", in_N, h_in, (int64_t) round(log2(sigma_in)));
  printf("Packing: (N=%ld, h=%ld, ternary, σ=2^%ld)\n", in_N, (uint64_t) 256, (int64_t) -44);
  printf("Output: (N=%ld, h=%ld, ternary, σ=2^%ld)\n", out_N, h_out, (int64_t) round(log2(sigma_out)));
  printf("Decomposition: BS(ℓ=%ld, β=2^%ld) PCK(ℓ=%ld, β=2^%ld) KS(ℓ=%ld, β=2^%ld)\n", 
        l, bg_bit, ell_packing, b_packing, t_ks, b_ks);
  printf("Message precision: %ld - Repetitions: %ld\n", msg_prec, reps);

  // end of parameters
  TRLWE_Key input_key;
  uint64_t rs_attempts = RS_sparse_binary_key(&input_key, in_N, in_k, h_in, sigma_in, target_r_prec);
  TRLWE_Key out_key = trlwe_new_ternary_key(out_N, out_k, h_out, sigma_out);
  TRLWE_Key packing_key = trlwe_new_ternary_key(in_N, in_k, 256, pow(2, -44));
  // TRLWE_Key packing_key = out_key;
  TRGSW_Key output_key = trgsw_new_key(out_key, l, bg_bit);
  const uint64_t r_prec = get_min_prec(input_key);
  printf("\nMax monomial distance (log B): %lu\n", r_prec);
  printf("Rejection Sampling Attempts: %lu\n", rs_attempts);
  SAB_Key sab = new_sparse_amortized_bootstrapping(input_key, packing_key, output_key, msg_prec, b_packing, ell_packing, t_ks, b_ks, h_in, r_prec, false, false, false);
  
  TorusPolynomial poly_in = polynomial_new_torus_polynomial(in_N);

  const uint64_t mod_mask = (1ULL<<(msg_prec - 1)) - 1;
  for (size_t i = 0; i < in_N; i++) poly_in->coeffs[i] = int2torus(i&mod_mask, msg_prec);

  TRLWE rlwe_in = trlwe_new_sample(poly_in, input_key);
  TRLWE rlwe_tv = trlwe_new_noiseless_trivial_sample(NULL, out_k, out_N);
  uint64_t LUT[1ULL << msg_prec];
  generate_random_bytes(sizeof(uint64_t)*(1ULL << msg_prec), (uint8_t *) LUT);
  for (size_t i = 0; i < (1ULL << msg_prec); i++) LUT[i] &= mod_mask;
  sab_LUT_packing(rlwe_tv, LUT, sab);
  // rlwe_tv->b->coeffs[0] += int2torus(1, 3);

  // TRLWE rlwe_out = trlwe_alloc_new_sample(in_k, in_N);
  // printf("in: "); trlwe_print(rlwe_in, input_key, msg_prec);
  // printf("tv: "); trlwe_print(rlwe_tv, out_key, msg_prec);
  MEASURE_TIME("", reps, "Bootstrapping time", 
    sab_rlwe_bootstrap(rlwe_in, rlwe_in, rlwe_tv, sab);
  );
  TorusPolynomial res_poly = polynomial_new_torus_polynomial(in_N);
  trlwe_phase(res_poly, rlwe_in, input_key);
  bool pass = true;
  for (size_t i = 0; i < in_N; i++){
    const uint64_t res = torus2int(res_poly->coeffs[i], msg_prec);
    uint64_t expected = LUT[torus2int(poly_in->coeffs[i], msg_prec)];
    for (size_t j = 1; j < reps; j++) expected = LUT[expected];
    if(res != expected){
      printf("\nFail %lu: %lu != %lu\n", i, res, expected);
      pass = false;
    }
  }
  if(pass) printf("Pass");
  printf("\n");
  // printf("out: "); trlwe_print(rlwe_in, input_key, msg_prec);
}

void test_sab_tern(){
  const uint64_t reps = 3;
#if defined(SET_2_3_2048)
  const uint64_t in_N = 2048, in_k = 1, out_N = 2048, out_k = 1, l = 1, bg_bit = 23, b_packing = 14, ell_packing = 2, t_ks = 12, b_ks = 1, h_in = 35, h_out = 512, msg_prec = 3;
  const double sigma_in = pow(2, -15);
  const double sigma_out = pow(2, -50);
  const uint64_t target_r_prec = 7;
#elif defined(SET_2_3_4096)
  const uint64_t in_N = 4096, in_k = 1, out_N = 2048, out_k = 1, l = 1, bg_bit = 23, b_packing = 14, ell_packing = 2, t_ks = 12, b_ks = 1, h_in = 26, h_out = 512, msg_prec = 3;
  const double sigma_in = pow(2, -15);
  const double sigma_out = pow(2, -50);
  const uint64_t target_r_prec = 9;
#elif defined(SET_2_3_8192)
  const uint64_t in_N = 8192, in_k = 1, out_N = 2048, out_k = 1, l = 1, bg_bit = 23, b_packing = 14, ell_packing = 2, t_ks = 12, b_ks = 1, h_in = 23, h_out = 512, msg_prec = 3;
  const double sigma_in = pow(2, -16);
  const double sigma_out = pow(2, -50);
  const uint64_t target_r_prec = 10;
#elif defined(SET_4_5_2048)
  const uint64_t in_N = 2048, in_k = 1, out_N = 2048, out_k = 1, l = 1, bg_bit = 23, b_packing = 14, ell_packing = 2, t_ks = 14, b_ks = 1, h_in = 38, h_out = 512, msg_prec = 5;
  const double sigma_in = pow(2, -17);
  const double sigma_out = pow(2, -50);
  const uint64_t target_r_prec = 7;
#elif defined(SET_4_5_4096)
  const uint64_t in_N = 4096, in_k = 1, out_N = 2048, out_k = 1, l = 1, bg_bit = 23, b_packing = 14, ell_packing = 2, t_ks = 15, b_ks = 1, h_in = 28, h_out = 512, msg_prec = 5;
  const double sigma_in = pow(2, -18);
  const double sigma_out = pow(2, -50);
  const uint64_t target_r_prec = 9;
#elif defined(SET_4_5_8192)
  const uint64_t in_N = 8192, in_k = 1, out_N = 2048, out_k = 1, l = 1, bg_bit = 23, b_packing = 14, ell_packing = 2, t_ks = 15, b_ks = 1, h_in = 24, h_out = 512, msg_prec = 5;
  const double sigma_in = pow(2, -18);
  const double sigma_out = pow(2, -50);
  const uint64_t target_r_prec = 10;
#elif defined(SET_6_7_4096)
  const uint64_t in_N = 4096, in_k = 1, out_N = 2048, out_k = 1, l = 1, bg_bit = 23, b_packing = 14, ell_packing = 2, t_ks = 17, b_ks = 1, h_in = 30, h_out = 512, msg_prec = 7;
  const double sigma_in = pow(2, -21);
  const double sigma_out = pow(2, -50);
  const uint64_t target_r_prec = 9;
#elif defined(SET_6_7_8192)
  const uint64_t in_N = 8192, in_k = 1, out_N = 2048, out_k = 1, l = 1, bg_bit = 23, b_packing = 14, ell_packing = 2, t_ks = 17, b_ks = 1, h_in = 25, h_out = 512, msg_prec = 7;
  const double sigma_in = pow(2, -21);
  const double sigma_out = pow(2, -50);
  const uint64_t target_r_prec = 10;
#elif defined(SET_8_9_4096)
  const uint64_t in_N = 4096, in_k = 1, out_N = 8192, out_k = 1, l = 1, bg_bit = 22, b_packing = 14, ell_packing = 2, t_ks = 19, b_ks = 1, h_in = 32, h_out = 512, msg_prec = 9;
  const double sigma_in = pow(2, -24);
  const double sigma_out = pow(2, -50);
  const uint64_t target_r_prec = 9;
#elif defined(SET_8_9_8192)
  const uint64_t in_N = 8192, in_k = 1, out_N = 8192, out_k = 1, l = 1, bg_bit = 22, b_packing = 14, ell_packing = 2, t_ks = 19, b_ks = 1, h_in = 26, h_out = 512, msg_prec = 9;
  const double sigma_in = pow(2, -23);
  const double sigma_out = pow(2, -50);
  const uint64_t target_r_prec = 10;
#else 
  const uint64_t in_N = 2048, in_k = 1, out_N = 2048, out_k = 1, l = 1, bg_bit = 23, b_packing = 14, ell_packing = 2, t_ks = 10, b_ks = 1, h_in = 17, h_out = 512, msg_prec = 3;
  const double sigma_in = pow(2, -15);
  const double sigma_out = pow(2, -50);
  const uint64_t target_r_prec = 8;
  #endif
  printf("Sparse bootstrapping with ternary keys\n");
  printf("Input: (N=%ld, h=%ld, ternary, σ=2^%ld)\n", in_N, h_in, (int64_t) round(log2(sigma_in)));
  printf("Packing: (N=%ld, h=%ld, ternary, σ=2^%ld)\n", in_N, (uint64_t) 256, (int64_t) -44);
  printf("Output: (N=%ld, h=%ld, ternary, σ=2^%ld)\n", out_N, h_out, (int64_t) round(log2(sigma_out)));
  printf("Decomposition: BS(ℓ=%ld, β=2^%ld) PCK(ℓ=%ld, β=2^%ld) KS(ℓ=%ld, β=2^%ld)\n", 
        l, bg_bit, ell_packing, b_packing, t_ks, b_ks);
  printf("Message precision: %ld - Repetitions: %ld\n", msg_prec, reps);

  // end of parameters
  TRLWE_Key input_key;
  const uint64_t rs_attempts = RS_sparse_ternary_key(&input_key, in_N, in_k, h_in, sigma_in, target_r_prec);
  TRLWE_Key out_key = trlwe_new_ternary_key(out_N, out_k, h_out, sigma_out);
  TRLWE_Key packing_key = trlwe_new_ternary_key(in_N, in_k, 256, pow(2, -44));
  // TRLWE_Key packing_key = out_key;
  TRGSW_Key output_key = trgsw_new_key(out_key, l, bg_bit);
  const uint64_t r_prec = get_min_prec(input_key);
  printf("\nMonomial distance (log B): %lu\n", r_prec);
  printf("Rejection Sampling Attempts: %lu\n", rs_attempts);
  SAB_Key sab = new_sparse_amortized_bootstrapping(input_key, packing_key, output_key, msg_prec, b_packing, ell_packing, t_ks, b_ks, h_in, r_prec, false, true, false);
  
  TorusPolynomial poly_in = polynomial_new_torus_polynomial(in_N);

  const uint64_t mod_mask = (1ULL<<(msg_prec - 1)) - 1;
  for (size_t i = 0; i < in_N; i++) poly_in->coeffs[i] = int2torus(i&mod_mask, msg_prec);

  TRLWE rlwe_in = trlwe_new_sample(poly_in, input_key);
  TRLWE rlwe_tv = trlwe_new_noiseless_trivial_sample(NULL, out_k, out_N);
  uint64_t LUT[1ULL << msg_prec];
  generate_random_bytes(sizeof(uint64_t)*(1ULL << msg_prec), (uint8_t *) LUT);
  for (size_t i = 0; i < (1ULL << msg_prec); i++) LUT[i] &= mod_mask;
  sab_LUT_packing(rlwe_tv, LUT, sab);
  // rlwe_tv->b->coeffs[0] += int2torus(1, 3);

  // TRLWE rlwe_out = trlwe_alloc_new_sample(in_k, in_N);
  // printf("in: "); trlwe_print(rlwe_in, input_key, msg_prec);
  // printf("tv: "); trlwe_print(rlwe_tv, out_key, msg_prec);
  MEASURE_TIME("", reps, "Bootstrapping time", 
    sab_rlwe_bootstrap(rlwe_in, rlwe_in, rlwe_tv, sab);
  );
  TorusPolynomial res_poly = polynomial_new_torus_polynomial(in_N);
  trlwe_phase(res_poly, rlwe_in, input_key);
  bool pass = true;
  for (size_t i = 0; i < in_N; i++){
    const uint64_t res = torus2int(res_poly->coeffs[i], msg_prec);
    uint64_t expected = LUT[torus2int(poly_in->coeffs[i], msg_prec)];
    for (size_t j = 1; j < reps; j++) expected = LUT[expected];
    if(res != expected){
      printf("\nFail %lu: %lu != %lu\n", i, res, expected);
      pass = false;
    }
  }
  if(pass) printf("Pass");
  printf("\n");
  // printf("out: "); trlwe_print(rlwe_in, input_key, msg_prec);
}


void test_sab_arbitrary(){
  const uint64_t reps = 3;
  const uint64_t key_bound = 8; // key in Z_8 = [-3,+4]
#if defined(SET_A2)
  const uint64_t in_N = 4096, in_k = 1, out_N = 2048, out_k = 1, l = 1, bg_bit = 23, b_packing = 14, ell_packing = 2, t_ks = 19, b_ks = 1, h_in = 32, h_out = 512, msg_prec = 5;
  const double sigma_in = pow(2, -24);
  const double sigma_out = pow(2, -50);
  const uint64_t target_r_prec = 9;
#elif defined(SET_A3)
  const uint64_t in_N = 8192, in_k = 1, out_N = 2048, out_k = 1, l = 1, bg_bit = 23, b_packing = 14, ell_packing = 2, t_ks = 19, b_ks = 1, h_in = 32, h_out = 512, msg_prec = 5;
  const double sigma_in = pow(2, -24);
  const double sigma_out = pow(2, -50);
  const uint64_t target_r_prec = 10;
#elif defined(SET_A4)
  const uint64_t in_N = 4096, in_k = 1, out_N = 8192, out_k = 1, l = 1, bg_bit = 22, b_packing = 14, ell_packing = 2, t_ks = 19, b_ks = 1, h_in = 32, h_out = 512, msg_prec = 7;
  const double sigma_in = pow(2, -24);
  const double sigma_out = pow(2, -50);
  const uint64_t target_r_prec = 9;
#elif defined(SET_A5)
  const uint64_t in_N = 8192, in_k = 1, out_N = 8192, out_k = 1, l = 1, bg_bit = 22, b_packing = 14, ell_packing = 2, t_ks = 19, b_ks = 1, h_in = 32, h_out = 512, msg_prec = 7;
  const double sigma_in = pow(2, -24);
  const double sigma_out = pow(2, -50);
  const uint64_t target_r_prec = 10;
#else 
  const uint64_t in_N = 4096, in_k = 1, out_N = 2048, out_k = 1, l = 1, bg_bit = 23, b_packing = 14, ell_packing = 2, t_ks = 19, b_ks = 1, h_in = 32, h_out = 512, msg_prec = 5;
  const double sigma_in = pow(2, -24);
  const double sigma_out = pow(2, -50);
  const uint64_t target_r_prec = 9;
#endif
  printf("Sparse bootstrapping with arbitrary keys\n");
  printf("Input: (N=%ld, h=%ld, θ=%ld, σ=2^%ld)\n", in_N, h_in, key_bound, (int64_t) round(log2(sigma_in)));
  printf("Packing: (N=%ld, h=%ld, ternary, σ=2^%ld)\n", in_N, (uint64_t) 256, (int64_t) -44);
  printf("Output: (N=%ld, h=%ld, ternary, σ=2^%ld)\n", out_N, h_out, (int64_t) round(log2(sigma_out)));
  printf("Decomposition: BS(ℓ=%ld, β=2^%ld) PCK(ℓ=%ld, β=2^%ld) KS(ℓ=%ld, β=2^%ld)\n", 
         l, bg_bit, ell_packing, b_packing, t_ks, b_ks);
  printf("Message precision: %ld - Repetitions: %ld\n", msg_prec, reps);

  // end of parameters
  TRLWE_Key input_key;
  const uint64_t rs_attempts = RS_sparse_arbitrary_key(&input_key, in_N, in_k, h_in, key_bound, sigma_in, target_r_prec);
  TRLWE_Key out_key = trlwe_new_ternary_key(out_N, out_k, h_out, sigma_out);
  TRLWE_Key packing_key = trlwe_new_ternary_key(in_N, in_k, 256, pow(2, -44));
  // TRLWE_Key packing_key = out_key;
  TRGSW_Key output_key = trgsw_new_key(out_key, l, bg_bit);
  const uint64_t r_prec = get_min_prec(input_key);
  printf("\nMonomial distance (log B): %lu\n", r_prec);
  printf("Rejection Sampling Attempts: %lu\n", rs_attempts);
  SAB_Key sab = new_sparse_amortized_bootstrapping(input_key, packing_key, output_key, msg_prec, b_packing, ell_packing, t_ks, b_ks, h_in, r_prec, false, false, true);
  
  TorusPolynomial poly_in = polynomial_new_torus_polynomial(in_N);

  const uint64_t mod_mask = (1ULL<<(msg_prec - 1)) - 1;
  for (size_t i = 0; i < in_N; i++) poly_in->coeffs[i] = int2torus(i&mod_mask, msg_prec);

  TRLWE rlwe_in = trlwe_new_sample(poly_in, input_key);
  TRLWE rlwe_tv = trlwe_new_noiseless_trivial_sample(NULL, out_k, out_N);
  uint64_t LUT[1ULL << msg_prec];
  for (size_t i = 0; i < (1ULL << msg_prec); i++) LUT[i] = ((i*i) &mod_mask);
  sab_LUT_packing(rlwe_tv, LUT, sab);
  // rlwe_tv->b->coeffs[0] += int2torus(1, 3);

  // TRLWE rlwe_out = trlwe_alloc_new_sample(in_k, in_N);
  // printf("in: "); trlwe_print(rlwe_in, input_key, msg_prec);
  // printf("tv: "); trlwe_print(rlwe_tv, out_key, msg_prec);
  MEASURE_TIME("", reps, "Bootstrapping time", 
    sab_rlwe_bootstrap(rlwe_in, rlwe_in, rlwe_tv, sab);
  );
  TorusPolynomial res_poly = polynomial_new_torus_polynomial(in_N);
  trlwe_phase(res_poly, rlwe_in, input_key);
  bool pass = true;
  for (size_t i = 0; i < in_N; i++){
    const uint64_t res = torus2int(res_poly->coeffs[i], msg_prec);
    uint64_t expected = LUT[torus2int(poly_in->coeffs[i], msg_prec)];
    for (size_t j = 1; j < reps; j++) expected = LUT[expected];
    if(res != expected){
      printf("\nFail %lu: %lu != %lu\n", i, res, expected);
      pass = false;
    }
  }
  if(pass) printf("Pass");
  printf("\n");
  // printf("out: "); trlwe_print(rlwe_in, input_key, msg_prec);
}

void test_sab_lwe(){
  const uint64_t in_N = 1024, in_k = 1, out_N = 2048, out_k = 1, l = 1, bg_bit = 23, b_packing = 12, ell_packing = 2, t_aut = 20, b_aut = 1, h_in = 64, msg_prec = 4;
  TRLWE_Key input_key = trlwe_new_sparse_binary_key(in_N, in_k, h_in, pow(2, -20));
  TRLWE_Key out_key = trlwe_new_sparse_binary_key(out_N, out_k, h_in, pow(2, -53));
  TRGSW_Key output_key = trgsw_new_key(out_key, l, bg_bit);
  const uint64_t r_prec = get_min_prec(input_key);
  printf("Min precision: %lu\n", r_prec);
  SAB_Key sab = new_sparse_amortized_bootstrapping(input_key, out_key, output_key, msg_prec, b_packing, ell_packing, t_aut, b_aut, h_in, r_prec, false, false, false);
  
  TorusPolynomial poly_in = polynomial_new_torus_polynomial(in_N);

  const uint64_t mod_mask = (1ULL<<(msg_prec - 1)) - 1;
  for (size_t i = 0; i < in_N; i++) poly_in->coeffs[i] = int2torus(i&mod_mask, msg_prec);

  TRLWE rlwe_in = trlwe_new_sample(poly_in, input_key);
  TRLWE rlwe_tv = trlwe_new_noiseless_trivial_sample(NULL, out_k, out_N);
  uint64_t LUT[1ULL << msg_prec];
  for (size_t i = 0; i < (1ULL << msg_prec); i++) LUT[i] = i;
  sab_LUT_packing(rlwe_tv, LUT, sab);

  printf("in: "); trlwe_print(rlwe_in, input_key, msg_prec);
  printf("tv: "); trlwe_print(rlwe_tv, out_key, msg_prec);
  MEASURE_TIME("", 1, "SAB LWE", 
    sab_rlwe_to_lwe_bootstrap(sab->tmp->extracted_poly, rlwe_in, rlwe_tv, sab);
  );
  TLWE_Key extracted_key = tlwe_alloc_key(out_N*out_k, output_key->trlwe_key->sigma);

  trlwe_extract_tlwe_key(extracted_key, output_key->trlwe_key);
  printf("out: ");
  for (size_t i = 0; i < in_N; i++){
    const uint64_t res = torus2int(tlwe_phase(sab->tmp->extracted_poly[i], extracted_key), msg_prec);
    const uint64_t expected = LUT[torus2int(poly_in->coeffs[i], msg_prec)];
    if(res != expected){
      printf("\nFail %lu: %lu != %lu\n", i, res, expected);
    }
    printf("%lu, ", res);
  }
  printf("\n");

  // trlwe_print(rlwe_out, input_key, msg_prec);
}

int main(int argc, char const *argv[])
{
#if defined(TERNARY)
  test_sab_tern();
#elif defined(ARBITRARY)
  test_sab_arbitrary();
#else
  test_sab();
#endif
  return 0;
}
