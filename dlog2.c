#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

#define MAX_PRIMES 20
#define MAX_EXPONENT 10

typedef struct {
    uint64_t module;
    uint64_t generator;
    uint64_t order;
} AlgebraicGroup;

typedef struct {
    uint64_t primes[MAX_PRIMES];
    uint64_t exponents[MAX_PRIMES];
    int count;
} Factorization;

typedef struct {
    uint64_t j_values[1000];
    uint64_t r_values[1000];
    int size;
} LookupTable;

typedef struct {
    uint64_t coeffs[MAX_EXPONENT];
    int size;
} Coefficients;

void init_group(AlgebraicGroup* group, uint64_t module, uint64_t generator);
uint64_t ring_pow(uint64_t base, uint64_t exp, uint64_t mod);
uint64_t gcd(uint64_t a, uint64_t b);
void EEA(uint64_t a, uint64_t b, int64_t* gcd, int64_t* x, int64_t* y);
uint64_t ring_inv(uint64_t a, uint64_t mod);
void factorize(uint64_t n, Factorization* fact);
void generate_lookup(uint64_t prime, uint64_t alpha, uint64_t n, LookupTable* table);
void solve(uint64_t alpha, uint64_t beta, uint64_t prime, uint64_t exponent, uint64_t n, LookupTable* table, Coefficients* coefs);
uint64_t solve_crt(uint64_t* remainders, uint64_t* moduli, int count);
uint64_t sph_log(uint64_t alpha, uint64_t beta, uint64_t n);
uint64_t brute_force_dl(uint64_t a, uint64_t b, uint64_t m);

int main() {
    // Test cases
    printf("Testing Silver-Pohlig-Hellman Algorithm\n");
    printf("========================================\n\n");
    
    // Test 1: Small example
    uint64_t alpha1 = 2, beta1 = 8, n1 = 10;
    clock_t start = clock();
    uint64_t result1 = sph_log(alpha1, beta1, n1);
    clock_t end = clock();
    double time1 = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Test 1: alpha=%llu, beta=%llu, n=%llu\n", alpha1, beta1, n1);
    printf("Result: x = %llu\n", result1);
    printf("Verification: %llu^%llu mod %llu = %llu\n", alpha1, result1, n1+1, ring_pow(alpha1, result1, n1+1));
    printf("Time: %.6f seconds\n\n", time1);
    
    // Test 2: From Python code - type 1, p=2
    uint64_t alpha2 = 15, beta2 = 38, n2 = 46;
    start = clock();
    uint64_t result2 = sph_log(alpha2, beta2, n2);
    end = clock();
    double time2 = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Test 2: alpha=%llu, beta=%llu, n=%llu\n", alpha2, beta2, n2);
    printf("Result: x = %llu\n", result2);
    printf("Verification: %llu^%llu mod %llu = %llu\n", alpha2, result2, n2+1, ring_pow(alpha2, result2, n2+1));
    printf("Time: %.6f seconds\n\n", time2);
    
    // Test 3: From Python code - type 1, p=3
    uint64_t alpha3 = 798, beta3 = 165, n3 = 910;
    start = clock();
    uint64_t result3 = sph_log(alpha3, beta3, n3);
    end = clock();
    double time3 = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Test 3: alpha=%llu, beta=%llu, n=%llu\n", alpha3, beta3, n3);
    printf("Result: x = %llu\n", result3);
    printf("Verification: %llu^%llu mod %llu = %llu\n", alpha3, result3, n3+1, ring_pow(alpha3, result3, n3+1));
    printf("Time: %.6f seconds\n\n", time3);
    
    return 0;
}

//Always init before usage
void init_group(AlgebraicGroup* group, uint64_t module, uint64_t generator) {
    group->generator = generator;
    group->module = module;
    group->order = module - 1;
}

//Pow mod
uint64_t ring_pow(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t result = 1;
    base = base % mod;
    
    while (exp > 0) {
        if (exp % 2 == 1) {
            result = (__uint128_t)result * base % mod;
        }
        exp = exp >> 1;
        base = (__uint128_t)base * base % mod;
    }
    
    return result;
}

//fact gcd
uint64_t gcd(uint64_t a, uint64_t b) {
    while (b != 0) {
        uint64_t temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

//Extended Euclidean alogrithm
void EEA(uint64_t a, uint64_t b, int64_t* gcd_val, int64_t* x, int64_t* y) {
    int64_t u0 = 1, u1 = 0;
    int64_t v0 = 0, v1 = 1;
    int64_t a_signed = (int64_t)a;
    int64_t b_signed = (int64_t)b;
    
    while (b_signed != 0) {
        int64_t q = a_signed / b_signed;
        int64_t r = a_signed % b_signed;
        
        int64_t u_temp = u0 - q * u1;
        int64_t v_temp = v0 - q * v1;
        
        u0 = u1;
        u1 = u_temp;
        v0 = v1;
        v1 = v_temp;
        
        a_signed = b_signed;
        b_signed = r;
    }
    
    *gcd_val = a_signed;
    *x = u0;
    *y = v0;
}

//Default a^-1
uint64_t ring_inv(uint64_t a, uint64_t mod) {
    int64_t gcd_val, x, y;
    EEA(a, mod, &gcd_val, &x, &y);
    
    if (gcd_val != 1) {
        return 0; // No inverse exists
    }
    
    // Make x positive
    x = x % (int64_t)mod;
    if (x < 0) {
        x += (int64_t)mod;
    }
    
    return (uint64_t)x;
}

//Simple factorizer
void factorize(uint64_t n, Factorization* fact) {
    fact->count = 0;
    
    //Check for 2
    if (n % 2 == 0) {
        fact->primes[fact->count] = 2;
        fact->exponents[fact->count] = 0;
        while (n % 2 == 0) {
            fact->exponents[fact->count]++;
            n /= 2;
        }
        fact->count++;
    }
    //Check odd
    for (uint64_t i = 3; i * i <= n; i += 2) {
        if (n % i == 0) {
            fact->primes[fact->count] = i;
            fact->exponents[fact->count] = 0;
            while (n % i == 0) {
                fact->exponents[fact->count]++;
                n /= i;
            }
            fact->count++;
        }
    }
    
    // If n > 1, then it's a prime factor
    if (n > 1) {
        fact->primes[fact->count] = n;
        fact->exponents[fact->count] = 1;
        fact->count++;
    }
}

//Lookup to use
void generate_lookup(uint64_t prime, uint64_t alpha, uint64_t n, LookupTable* table) {
    table->size = 0;
    
    for (uint64_t j = 0; j < prime; j++) {
        uint64_t exponent = (n * j) / prime;
        uint64_t r_j = ring_pow(alpha, exponent, n + 1);
        
        table->j_values[table->size] = j;
        table->r_values[table->size] = r_j;
        table->size++;
    }
}

void solve(uint64_t alpha, uint64_t beta, uint64_t prime, uint64_t exponent,uint64_t n, LookupTable* table, Coefficients* coefs) {
    coefs->size = 0;
    uint64_t inverse_alpha = ring_inv(alpha, n + 1);
    
    for (uint64_t l = 0; l < exponent; l++) {
        if (l == 0) {
            uint64_t exp = n / prime;
            uint64_t res = ring_pow(beta, exp, n + 1);
            
            // Find res in table
            int found = -1;
            for (int i = 0; i < table->size; i++) {
                if (table->r_values[i] == res) {
                    found = i;
                    break;
                }
            }
            
            if (found >= 0) {
                coefs->coeffs[coefs->size++] = table->j_values[found];
            } else {
                coefs->coeffs[coefs->size++] = 0;
            }
        } else {
            // Calculate alpha_pow = sum(coefs[m] * prime^m)
            uint64_t alpha_pow = 0;
            uint64_t prime_power = 1;
            for (int m = 0; m < coefs->size; m++) {
                alpha_pow += coefs->coeffs[m] * prime_power;
                prime_power *= prime;
            }
            
            uint64_t prime_pow_l_plus_1 = 1;
            for (uint64_t i = 0; i <= l; i++) {
                prime_pow_l_plus_1 *= prime;
            }
            
            uint64_t main_power = n / prime_pow_l_plus_1;
            uint64_t inv_alpha_pow = ring_pow(inverse_alpha, alpha_pow, n + 1);
            uint64_t temp = (__uint128_t)beta * inv_alpha_pow % (n + 1);
            uint64_t res = ring_pow(temp, main_power, n + 1);
            
            // Find res in table
            int found = -1;
            for (int i = 0; i < table->size; i++) {
                if (table->r_values[i] == res) {
                    found = i;
                    break;
                }
            }
            
            if (found >= 0) {
                coefs->coeffs[coefs->size++] = table->j_values[found];
            } else {
                coefs->coeffs[coefs->size++] = 1;
            }
        }
    }
}

uint64_t solve_crt(uint64_t* remainders, uint64_t* moduli, int count) {
    uint64_t product = 1;
    for (int i = 0; i < count; i++) {
        product *= moduli[i];
    }
    
    uint64_t x = 0;
    for (int j = 0; j < count; j++) {
        uint64_t sup = product / moduli[j];
        uint64_t inv = ring_inv(sup, moduli[j]);
        uint64_t add = (__uint128_t)remainders[j] * sup % product;
        add = (__uint128_t)add * inv % product;
        x = (x + add) % product;
    }
    
    return x;
}

uint64_t sph_log(uint64_t alpha, uint64_t beta, uint64_t n) {
    // Factorize n
    Factorization fact;
    factorize(n, &fact);
    
    uint64_t remainders[MAX_PRIMES];
    uint64_t moduli[MAX_PRIMES];
    
    // For each prime factor
    for (int i = 0; i < fact.count; i++) {
        uint64_t prime = fact.primes[i];
        uint64_t exp = fact.exponents[i];
        
        // Create lookup table
        LookupTable table;
        generate_lookup(prime, alpha, n, &table);
        
        // Find coefficients
        Coefficients coefs;
        solve(alpha, beta, prime, exp, n, &table, &coefs);
        
        // Calculate y_i = sum(coefs[m] * prime^m)
        uint64_t y_i = 0;
        uint64_t prime_power = 1;
        for (int m = 0; m < coefs.size; m++) {
            y_i += coefs.coeffs[m] * prime_power;
            prime_power *= prime;
        }
        
        // Calculate modulus
        uint64_t mod_i = 1;
        for (uint64_t j = 0; j < exp; j++) {
            mod_i *= prime;
        }
        
        remainders[i] = y_i % mod_i;
        moduli[i] = mod_i;
    }
    
    // Solve using Chinese Remainder Theorem
    return solve_crt(remainders, moduli, fact.count);
}

uint64_t brute_force_dl(uint64_t a, uint64_t b, uint64_t m) {
    for (uint64_t x = 2; x < m; x++) {
        if (ring_pow(a, x, m) == b) {
            return x;
        }
    }
    return 0;
}
