#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <omp.h>
//#include <cstring>
#include <algorithm>
#include <chrono>
#include <gmpxx.h> // large integers
#include <string>  // for std::stoull
#include <cstring>
#include <vector>
using namespace std;
//#define ITERATIONS 1000000000

// gcd(a,b)
inline unsigned long long gcd(unsigned long long a, unsigned long long b) {
    while (b != 0) {
        unsigned long long rem = a % b;
        a = b;
        b = rem;
    }
    return a;
}

// computes a*b modulo n
inline unsigned long long mod_product(unsigned long long a, unsigned long long b, unsigned long long n) {
    return ((((__uint128_t) a) * ((__uint128_t) b)) % n);
}

// computes a^b modulo n
inline unsigned long long powmod(unsigned long long a, unsigned long long b, unsigned long long n) {
    unsigned long long res = 1;
    while (b > 0) {
        if (b & 1) {
            res = mod_product(a, res, n);
        }
        a = mod_product(a, a, n);
        b >>= 1;
    }
    return res;
}

void mpz_cofactor(mpz_t mpz_product, vector<unsigned long long> factors) {
    mpz_set_ui(mpz_product, 1);
    for (unsigned long long f : factors) {
        mpz_t fac;
        mpz_init(fac);
        mpz_import(fac,          // GMP variable
           1,                // one "word"
           -1,               // most significant word first
           sizeof(f), // size of a word
           0,                // native endianness
           0,                // no nails
           &f);     // pointer to data
        mpz_lcm(mpz_product, mpz_product, fac);
        mpz_clear(fac);
    }
}

inline int factors(unsigned long long base, unsigned long long n, unsigned long long lower, unsigned long long iters, int threads, bool quiet) {
    int success = 1; // false
    if (!quiet) cout << "Starting Trial Division..." << endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    omp_set_num_threads(threads);
    vector<unsigned long long> factor_threads[threads];
    #pragma omp parallel
    {
        unsigned long long q;
        int tid = omp_get_thread_num();
        #pragma omp for
        for (unsigned long long i = lower; i < iters; i++) {
            q = 6*n*i+1;
            if( q%5!=0 && q%7!=0 && q%11!=0 && q%13!=0 &&  
                q%17!=0 && q%19!=0 && q%23!=0 && q%29!=0 && q%31!=0 && q%37!=0 && 
                q%41!=0 && q%43!=0 && q%47!=0 && q%53!=0 && q%59!=0 && 
                powmod(base, n, q) == (q-1) ) {
                factor_threads[tid].push_back(q);
                    //#pragma omp critical 
                    //{
                    //    //cout << "Factor of " << base << "^" << n << "-1 found : " << q << endl;
                    //}
            }
        }
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    if (!quiet) std::cout << "\nTime: " << elapsed.count() << " seconds\n";
    vector<unsigned long long> factors;
    factors.push_back(base+1);
    for (int i = 0; i < threads; i++) {
        int szz = factor_threads[i].size();
        for (int j = 0; j < szz; j++) {
            factors.push_back(factor_threads[i][j]);
        }
    }
    sort(factors.begin(), factors.end());
    if (factors.size() > 1) success = 0;
    mpz_t mpz_product;
    mpz_init(mpz_product); // Initialize mpz_t
    mpz_cofactor(mpz_product, factors);
    char *cf = mpz_get_str(NULL, 10, mpz_product);
    string cofactor = "(" + std::to_string(base) + "^" + std::to_string(n) + "+1)/" + cf;
    if (!quiet) cout << "Cofactor Result : ";
    cout << cofactor << endl;
    free(cf);
    return success;
}

// compile with  g++ -O3 -fopenmp trial-const.cpp -lgmp -lgmpxx -o trial-const
// i.e.
// sudo apt-get update && sudo apt-get install g++ && sudo apt-get install libgmp3-dev && g++ -O3 -fopenmp trial-const-p3.cpp -lgmp -lgmpxx -o trial-const-p3
int main(int argc, char* argv[]) {
    if (argc < 4 || argc > 9) {
        return 0;
    }
    int counter = 1;
    bool quiet = false;
    if (!strcmp(argv[counter],"-q")) {
        counter++;
        quiet = true;
    }
    unsigned long long base = stoull(argv[counter]);
    counter++;
    unsigned long long n = stoull(argv[counter]);
    counter++;
    unsigned long long lower = 1;
    unsigned long long iters = 0;
    int t = 1;
    if (!strcmp(argv[counter],"-lo")) {
        counter++;
        lower = stoul(argv[counter]);
        counter++;
    }
    iters = stoul(argv[counter]);
    counter++;
    if (!strcmp(argv[counter],"-t")) {
        counter++;
        t = stoi(argv[counter]);
        counter++;
    }
    base = 3;
    return factors(base, n, lower, iters, t, quiet);
}