#include <stdint.h>
#include <math.h>

static uint64_t mod_diff(uint64_t a, uint64_t b, uint64_t q){
  const uint64_t half_q = q==0? 1ULL<<63 : q/2;
  const uint64_t diff_abs = a < b ? b-a : a-b;
  return diff_abs > half_q ? q - diff_abs : diff_abs;
}

#define MAX_POINTS (1<<16)
static double __variance_array[MAX_POINTS];
static double __variance_sum = 0;
static double __variance_cnt = 0;
static uint64_t __var_idx = 0;

#define PRINT_VAR 
#define PRINT_TORUS

#ifndef __Q
#define __Q 0
#endif
#ifndef DECRYPTION_FUNCTION
  #error "Define a decryption function using '#define DECRYPTION_FUNCTION' before importing this code"
#endif

static uint64_t mod_switch(uint64_t v, uint64_t p, uint64_t q) {
  const double double_q = q == 0? pow(2,64) : ((double)q); 
  const double double_p = p == 0? pow(2,64) : ((double)p); 
  uint64_t val = (uint64_t) round((((double) v) * double_q) / double_p);
  return val < q? val : val - q;
}

static void __log_noise(uint64_t value, uint64_t expected){
  const uint64_t diff = mod_diff(value, expected, __Q);
  __variance_sum += pow(diff, 2);
  __variance_cnt++;
}

static void __log_noise_prec(uint64_t value, uint64_t prec){
  const uint64_t mod_mask = (1ULL << (64 - prec)) - 1;
  const uint64_t noise_mod = (1ULL << (64 - prec));
  const uint64_t diff = value & mod_mask;
  const uint64_t mod_diff = diff > noise_mod/2 ? noise_mod - diff : diff; 
  __variance_sum += pow(mod_diff, 2);
  __variance_cnt++;
}

static void __log_noise_array(uint64_t * value, uint64_t * expected, uint64_t size){
  for (size_t i = 0; i < size; i++){
    __log_noise(value[i], expected[i]);
  }
}

static void __log_noise_array_prec(uint64_t * value, uint64_t prec, uint64_t size){
  for (size_t i = 0; i < size; i++){
    __log_noise_prec(value[i], prec);
  }
}

static void print_noise(){
  double log_sigma = log2(sqrt(__variance_sum/((double)__variance_cnt)));

  #ifdef PRINT_TORUS
  log_sigma -= 64;
  #endif

  #ifdef PRINT_VAR
  printf("Sigma^2: 2^%lf\n", 2*log_sigma);
  #else
  printf("Sigma: 2^%lf\n", log2(sigma));
  #endif
}

static void reset_noise(){
  __variance_sum = 0;
  __variance_cnt = 0;
}

static void print_reset_noise(){
  print_noise();
  reset_noise();
}

#define LOG_NOISE(value, expected_value) __log_noise(value, expected_value)
#define DECRYPT_LOG_NOISE(ciphertext, expected_value) __log_noise(DECRYPTION_FUNCTION(ciphertext), expected_value)
#define DECRYPT_LOG_BY_PREC(ciphertext, prec) __log_noise_prec(DECRYPTION_FUNCTION(ciphertext), prec)
#define ARRAY_LOG_BY_PREC(value, prec, size) __log_noise_array_prec(value, prec, size)


