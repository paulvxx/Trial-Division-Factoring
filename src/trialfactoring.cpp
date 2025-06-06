#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <omp.h>
#include <cstring>
#include <chrono>
#include <gmpxx.h> // large integers
#include <string>  // for std::stoull
#include <vector>
using namespace std;
#define size_t WHEEL_SIZE = 5760;
#define UPDATE 10000

// COMPILE WITH g++ -O3 trialdiv.cpp -o trialdiv
// g++ -O3 trialdiv.cpp -o trialdiv -lgmp -lgmpxx

// gcd function
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

// wheel sieve for n*k+1 factors
// mult = multiplier (n) in n*k+1 factors
// use 2,3,5,7,11,13 in the wheel sieve
// offset = first i such that n*i+1 is not divisible by 2,3,5,7,11,13
std::vector<unsigned long long> wheel(unsigned long long mult, unsigned long long offset) { 
    unsigned long long wh = 2*3*5*7*11*13;
    unsigned long long wheel_range = wh / gcd(mult, wh);
    //unsigned long long wheel_range = 5760 / gcd(mult, 5760);
    std::vector<unsigned long long> wheelMultipliers;
    unsigned long long prev = offset;
    for (int i = offset+1; i <= wheel_range; i++) {
        unsigned long long q = mult * i + 1;
        if ((q%12==1 || q%12==11) && q%2 && q%3 && q%5!=0 && q%7!=0 && q%11!=0 && q%13!=0) { 
                wheelMultipliers.push_back((i - prev) * mult);
            prev = i;
        }
    }
    return wheelMultipliers;
}

// output product of factors as single factor
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
        mpz_mul(mpz_product, mpz_product, fac);
        mpz_clear(fac);
    }
}

// perform repunit (b^n +- 1)/(b +- 1) sieve : 
// b : base
// n : exponent (should be odd prime)
// limit : sieve limit for factorization
vector<unsigned long long> sieve_repunit(
    unsigned int b, unsigned int n, unsigned long long lower, unsigned long long limit, bool plus, int num_theads, string filename="") {
    vector<unsigned long long> factors;
    factors.push_back(b-1);
    // threading
    //omp_set_num_threads(num_theads);
    // is prime check here
    // if not prime, return nothing
    // use factors of the form 2*k*n+1 (default)
    // do gmp check to see size
    // if size is small enough (i.e. less than 64 bits) do factoring
    // or if size is < sieve_limit ^ 2,   output factored number?
    // check if n == 2;
    // else perform sieve
    //unsigned long long increment = 2*n;
    unsigned long long q; // = 1 + 2*increment;
    std::vector<unsigned long long> wheel_vec = wheel(2*n, (lower) / (2*n));
    int wrap = wheel_vec.size();
    //q = 1 + 2*n*wheel_vec[0];
    q = 1 + wheel_vec[0];
    // for progress status
    unsigned long long update = limit / UPDATE;
    unsigned long long milestone = update;
    int index = 1;
    auto start_time = std::chrono::high_resolution_clock::now();
    cout << "Starting Trial Division..." << endl;
    while (q < limit) {
        // small prime filter
        if (q%17!=0 && q%19!=0 && q%23!=0 && q%29!=0 && q%31!=0 && q%37!=0 && q%41!=0 && q%43!=0 && q%47!=0 
            && q%53!=0 && q%59!=0 && q%61!=0 && q%67!=0 && q%71!=0 && q%73!=0 && q%79!=0 && q%83!=0
            //&& q%89!=0 && q%97!=0 && q%101!=0 && q%103!=0 && q%107!=0 && q%109!=0 && q%113!=0
            //&& q%127!=0 && q%131!=0 && q%137!=0 && q%139!=0 && q%149!=0 && q%151!=0 && q%157!=0
            //&& q%163!=0 && q%167!=0 && q%173!=0 && q%179!=0 && q%181!=0 && q%191!=0 && q%193!=0
            ) {
            //unsigned int q28 = q % 12;
            //if ((q28==1 || q28==3 || q28==9 || q28==27 || q28==25 || q28==19)) {
                //unsigned int q28 = q % 12;
                //(q28==1 || q28==3 || q28==9 || q28==27 || q28==25 || q28==19)
                if ( powmod((unsigned long long) b, (unsigned long long) n, q)==1) {
                    int sz = factors.size();
                    bool inlist = false;
                    for (int i = 0; i < sz; i++) {
                        if (q%factors[i]==0) {
                            inlist=true; 
                            break;
                        }
                    }
                    if (!inlist) {
                        factors.push_back(q);
                        cout << "\n" << q << " divides (" << b << "^" << n; 
                        if (plus) cout << "+1)/";
                        else cout << "-1)/";
                        cout << b-1 << "!" << endl;
                    }
                }
            //}
        }

        //q += 2*n*wheel_vec[index];
        q += wheel_vec[index];
        index = index + 1;
        if (index == wrap) index = 0;

        if (q >= milestone) {
            double percentage = (static_cast<double>(q) / limit) * 100.0;
            std::cout << "\rProgress : "<< std::fixed << std::setprecision(2) << std::setw(6) << percentage << "%" << std::flush;
            milestone += update;
        }
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << "\nTime taken: " << elapsed.count() << " seconds\n";
    std::cout << "Factors Found: [";
    //string cofactor = "(" + std::to_string(b) + "^" + std::to_string(n) + "-1)/" + "(" + std::to_string(b-1);
    int sz = factors.size();
    for (int i = 0; i < sz-1; i++) {
        cout << factors[i];
        cout << ", ";
        //cofactor = cofactor + "*" + std::to_string(factors[i]);
    }
    if (sz > 0) {
        cout << factors[sz-1] << "]" << endl;
        //cofactor = cofactor + "*" + std::to_string(factors[sz-1]);
    }
    //cofactor += ")";
    mpz_t mpz_product;
    mpz_init(mpz_product); // Initialize mpz_t
    mpz_cofactor(mpz_product, factors);
    char *cf = mpz_get_str(NULL, 10, mpz_product);
    string cofactor = "(" + std::to_string(b) + "^" + std::to_string(n) + "-1)/" + cf;
    cout << "Remaining Cofactor:  " << cofactor << endl;
    if (filename!="") {
        std::ofstream output_file;
        output_file.open(filename, std::ios::app);
        if (!output_file) {
            std::cerr << "Error opening file for writing.\n";
            return factors;
        }
        output_file << cofactor << "\n";
        output_file.close();
    }
    free(cf);
    return factors;
}

void usage() {
    cout << "Sieving utility to trial divide numbers k*b^n+c [(k,b,n,c) < 64 bits]" << endl;
    cout << "Usage:" << endl;
    cout << "Program MUST follow the following order:" << endl;
    cout << "./sieve-program <type REQUIRED> <factoring limit REQUIRED> <options REQUIRED> <extended options>" << endl;
    cout << "OR run ./sieve-program --help to get this usage" << endl;
}

int main(int argc, char* argv[]) {
    //unsigned long long result = powmod(3, 10001, 67);
    //cout << result << endl;    
    if (argc < 7) {
        if (argc == 1 || (argc == 2 && strcmp(argv[1], "--help") == 0)) {
            usage();
            return 0;
        }
        cout << "Error : Command line arguments <type>, <factoring limit>, and <options> are required! For usage and options run ./sieve --help" << endl;
        return -1;
    }
    int i = 1;
    bool plus = false;
    if (!strcmp(argv[i], "-r")) {
        i++;
    } 
    /// ... elif for -c and -g
    else {
        cout << "Error." << endl;
        return -1;
    }

    if (argc > 7 && !strcmp(argv[i], "-p")) {
        plus = true;
        i++;
    }

    unsigned long long lower = 0;
    if (!strcmp(argv[i], "-low")) {
        lower = std::stoull(argv[i]);
        i++;
    }
    // ...
    if (strcmp(argv[i], "-l")) {
        cout << "Error." << endl;
    }
    i++;
    unsigned long long limit = std::stoull(argv[i]);
    //cout << limit << endl;
    i++;
    // if (strcmp("-e", argv[i]) == 0) Not yet implemented
    // else scan arguments
    unsigned long long base, n;
    if (strcmp("-s", argv[i]) == 0) {
        i++;
        // base
        base = std::stoull(argv[i]); 
        i++;
        // n
        n = std::stoull(argv[i]);
        i++;
        // filename
        string filename = "";
        if (i < (argc+1) && strcmp("-f", argv[i]) == 0) {
            i++;
            filename = argv[i];
        }
        cout << "Factoring (" << base << "^" << n;
        if (plus) cout << "+1)/" << base+1 << endl; 
        else cout << "-1)/" << base-1 << endl;
        sieve_repunit(base,n,lower,limit,plus,2,filename);
        //sieve_repunit
    }
    // not yet implemented
    //else if (strcmp("-m", argv[i]) == 0) {
    //    i++;
    //}
    //else if (strcmp("-f", argv[i]) == 0) {
    //    i++;
    //
    //}
    else {
        cout << "Error." << endl;
        return -1;
    }

    //unsigned long long t = 0;
    //std::vector<unsigned long long> v = wheel(17);
    //for (i = 0; i < v.size(); i++) {
    //    cout << v[i] << ",";
    //}
    //cout << endl;
    return 0;
}
