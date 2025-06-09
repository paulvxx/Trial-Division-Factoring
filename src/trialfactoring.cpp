#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <omp.h>
#include <cstring>
#include <algorithm>
#include <chrono>
#include <gmpxx.h> // large integers
#include <string>  // for std::stoull
#include <vector>
using namespace std;
#define size_t WHEEL_SIZE = 5760;
#define UPDATE 10000

inline bool is_quadratic_residue_switch(unsigned int base, unsigned long long q, bool plus) {
    // current hard-coded implementation
    if (base > 12) return true;
    if (!plus) {
        switch (base) {
            case 2: {
                int r = q & 7;
                return r == 1 || r == 7;
            }
            case 3: {
                int r = q % 12;
                return r == 1 || r == 11;
            }
            case 5: {
                int r = q % 10;
                return r == 1 || r == 9;
            }
            case 6: {
                int r = q % 24;
                return r == 1 || r == 5 || r == 19 || r == 23;
            }
            case 7: {
                int r = q % 28;
                return r == 1 || r == 3 || r == 9 || r == 19 || r == 25 || r == 27;
            }
            case 10: {
                int r8 = q & 7;
                int r = q % 5;
                return ((r8==1 || r8==7) && (r==1 || r==4)) || ((r8==3 || r8==5) && (r==2 || r==3));
            }
            case 11: {
                int r4 = q & 3;
                int r = q % 11;
                return (r4==1 && (r==1 || r==3 || r==4 || r==5 || r==9)) || (r4==3 && (r==2 || r==6 || r==7 || r==8 || r==10));
            }
            case 12: {
                int r = q % 12;
                return r == 1 || r == 11;
            }
            default:
                // If the base isn't in our hardcoded list, we can't optimize.
                return true;
        }
    }
    else {
        switch (base) {
            case 2: {
                int r = q & 7;
                return r == 1 || r == 3;
            }
            case 3: {
                int r = q % 3;
                return r == 1;
            }
            case 5: {
                int r = q % 20;
                return r == 1 || r == 3 || r == 7 || r == 9;
            }
            case 6: {
                int r = q % 24;
                return r == 1 || r == 5 || r == 7 || r == 11;
            }
            case 7: {
                int r = q % 7;
                return r == 1 || r == 2 || r == 4;
            }
            case 9: {
                int r = q & 3;
                return r == 1;
            }            
            case 10: {
                int r8 = q & 7;
                int r = q % 5;
                return ((r8==1 || r8==3) && (r==1 || r==4)) || ((r8==5 || r8==7) && (r==2 || r==3)); 
            }
            case 11: {
                int r = q % 11;
                return (r==1 || r==3 || r==4 || r==5 || r==9);
            }
            case 12: {
                int r = q % 3;
                return r == 1;
            }
            default:
                // If the base isn't in our hardcoded list, we can't optimize.
                return true;
        }        
    }
}

struct SieveState {
    unsigned long long start; // The first valid candidate >= the thread's lower bound
    unsigned int index;            // The index of the *next* increment to use
    unsigned long long iters;                  // Number of iterations
};

// COMPILE WITH
// g++ -O3 -fopenmp trialdiv-demo.cpp -o trialdiv-demo -lgmp -lgmpxx
// gcd function

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

// wheel sieve for n*k+1 factors
// mult = multiplier (n) in n*k+1 factors
// use 2,3,5,7,11,13 in the wheel sieve
// offset = first i such that n*i+1 is not divisible by 2,3,5,7,11,13
std::vector<unsigned int> wheel(unsigned long long mult) { 
    // quick wheel to skip through numbers of the form mult*k+1 which are divisible by small primes
    unsigned long long wh = 2*3*5*7*11*13; //*7*11*13;
    // how many numbers the cycle iterates through before it repeats (i.e. 1,5,7,11,13,17,19,23,25,... skips through 6 numbers before repeating)
    unsigned long long wheel_range = wh / gcd(mult, wh);
    //unsigned long long wheel_range = 5760 / gcd(mult, 5760);
    std::vector<unsigned int> wheelMultipliers;
    unsigned long long prev = 0; // offset is supposed to work for lower bound limits
    for (int i = 1; i <= wheel_range; i++) {
        unsigned long long q = mult * i + 1;
        // TODO LATER : Implement Quadratic Residue (QR) Sieve 
        // i.e. b is a quadratic residue all divisors of (b^n-1)
        if (q%2!=0 && q%3!=0 && q%5!=0 && q%7!=0 && q%11!=0 && q%13!=0) { 
            //if (q < 100000000) cout << q << endl;
            wheelMultipliers.push_back((i - prev) * mult);
            prev = i;
        }
    }
    return wheelMultipliers;
}

unsigned long long get_starting_value(unsigned long long offset, 
    unsigned long long modulus, int tid=0) {
        const unsigned int wheel_mod = 30030; // default
        unsigned int n = modulus;
        unsigned int mod_round = offset % n;
        if (mod_round == 0) {
            offset += 1;
        }
        else if (mod_round != 1) {
            offset = offset + (n - mod_round) + 1;
        }
        while (gcd(offset,wheel_mod)!=1) offset += n;
        //cout << tid << "," << offset << endl;
        return offset;
}

SieveState get_sieve_state(std::vector<unsigned int> wheel, 
    unsigned long long lower, 
    unsigned long long upper,  
    unsigned long long modulus, int tid=0) {
        //cout << "Called! Begin " << tid << endl;
        const unsigned int wheel_mod = 30030; // default
        SieveState s;
        unsigned long long n = modulus;
        s.start = get_starting_value(lower, n, tid);
        s.index = 0;
        unsigned long long upto = s.start % (30030*n);
        unsigned long long counter = 1;
        //cout << upto << endl;
        //cout << counter << endl;
        while (counter != upto) {
            counter += wheel[s.index];
            s.index++;
            if (s.index == 5760) s.index = 0;
        }
        unsigned long long div_iters = (upper - lower) / (30030*n);
        unsigned long long remainder = upper - lower - div_iters*n;
        unsigned long long num_iters = div_iters * 5760;
        counter = lower + div_iters*30030*n + (s.start - lower);
        unsigned int idx = s.index;
        while (counter <= upper) {
            num_iters++;
            counter += wheel[idx];
            idx++;
            if (idx == 5760) idx = 0;
        }
        s.iters = num_iters;
        //cout << "Called End! " << tid << endl;
        return s;
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
    unsigned int b, unsigned int n, unsigned long long lower, unsigned long long limit, bool plus, int threads, string filename = "") {
    vector<unsigned long long> factors;
    if (lower < 2) lower = 2;
    if (limit == 0 || lower > limit) return factors; // Return early if limit is invalid
    int pm1 = plus ? 1 : -1;
    factors.push_back(b+pm1);
    // TODOs
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
    //unsigned long long q; // = 1 + 2*increment;
    //cout << (lower / n) << endl;
    std::vector<unsigned int> wheel_vec = wheel(n);
    // size of the wheel
    int wrap = wheel_vec.size();
    // increment if needed
    // for progress status
    unsigned long long update = limit / UPDATE;
    unsigned long long milestone = update;
    //
    //if (q == 1) {
    //    q = q + wheel_vec[index];
    //    index++;
    //    if (index == wrap) index = 0;
    //    num_iterations--;
    //}
    //
    cout << "Starting Trial Division..." << endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    // Begin trial division
    SieveState state = get_sieve_state(wheel_vec, lower, limit, n);
    unsigned long long q = state.start;
    unsigned int index = state.index;
    unsigned long long num_iterations = state.iters;
    //cout << lower << endl;
    //cout << q << endl;
    //cout << index << endl;
    //cout << num_iterations << endl;
    // no multithreading due to tradeoffs
    if (threads == 1 || num_iterations <= 10000000) {
        for (unsigned long long i = 0; i < num_iterations; i++) {
            // small prime filter sieve
            if (q%17!=0 && q%19!=0 && q%23!=0 && q%29!=0 && q%31!=0 && q%37!=0 && q%41!=0 && q%43!=0 && q%47!=0 
                && q%53!=0 && q%59!=0 && q%61!=0 && q%67!=0 && q%71!=0 && q%73!=0 && q%79!=0 && q%83!=0
                ) {
                    // proper residue to check for (+-1)
                    unsigned long long res = 1;
                    if (plus) res = q-1;
                    // compute modular exponentiation COSTLY
                    if ( is_quadratic_residue_switch(b,q,plus) && powmod((unsigned long long) b, (unsigned long long) n, q) == res ) {
                        int sz = factors.size();
                        bool inlist = false;
                        for (int i = 0; i < sz; i++) {
                            if (q%factors[i]==0) {
                                inlist=true; 
                                break;
                            }
                        }
                        // only store a factor if it is a prime factor
                        if (!inlist) {
                            factors.push_back(q);
                            cout << "\n" << q << " divides (" << b << "^" << n; 
                            if (plus) cout << "+1)/";
                            else cout << "-1)/";
                            cout << b+pm1 << "!" << endl;
                        }
                    }
            }
            //q += 2*n*wheel_vec[index];
            q += wheel_vec[index];
            index++;
            if (index == wrap) index = 0;

            if (q >= milestone) {
                double percentage = (static_cast<double>(q) / limit) * 100.0;
                std::cout << "\rProgress : "<< std::fixed << std::setprecision(2) << std::setw(6) << percentage << "%" << std::flush;
                milestone += update;
            }
        }
    }
    // threading here
    else {
        omp_set_num_threads(threads);
        long long total_range = limit - lower + 1;
        vector<unsigned long long> factor_threads[threads];
        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            long long block = total_range / threads;
            long long start = tid * block + lower;
            long long end;
            // Last thread picks up the remainder
            if (tid == threads - 1) {
                end = limit;
            } else {
                end = (tid + 1) * block + lower - 1;
            }
            SieveState state_thread = get_sieve_state(wheel_vec, start, end, n, tid);
            unsigned long long q = state_thread.start;
            unsigned long long index = state_thread.index;
            unsigned long long num_iterations_t = state_thread.iters;
            #pragma omp critical 
            {
                cout << "Thread #" << tid << " starting with : q=" << q << ", index=" << index << ", iterations=" << num_iterations_t << endl;
            }
            for (unsigned long long i = 0; i < num_iterations_t; i++) {
                //cout << i << endl;
                // small prime filter sieve
                if (q%17!=0 && q%19!=0 && q%23!=0 && q%29!=0 && q%31!=0 && q%37!=0 && q%41!=0 && q%43!=0 && q%47!=0 
                    && q%53!=0 && q%59!=0 && q%61!=0 && q%67!=0 && q%71!=0 && q%73!=0 && q%79!=0 && q%83!=0
                    ) {
                        // proper residue to check for (+-1)
                        unsigned long long res = 1;
                        if (plus) res = q-1;
                        // compute modular exponentiation COSTLY
                        if ( is_quadratic_residue_switch(b,q,plus) && powmod((unsigned long long) b, (unsigned long long) n, q) == res ) {
                                int sz = factor_threads[tid].size();
                                bool inlist = false;
                                for (int i = 0; i < sz; i++) {
                                    if (q%factor_threads[tid][i]==0) {
                                        inlist=true; 
                                        break;
                                    }
                                }
                                // only store a factor if it is a prime factor
                                if (!inlist) {
                                    factor_threads[tid].push_back(q);
                                    #pragma omp critical 
                                    {
                                        cout << "\n" << q << " divides (" << b << "^" << n; 
                                        if (plus) cout << "+1)/";
                                        else cout << "-1)/";
                                        cout << b+pm1 << "!" << endl;
                                    }
                                }
                        }
                }
                //q += 2*n*wheel_vec[index];
                q += wheel_vec[index];
                index++;
                if (index == wrap) index = 0;

                if (tid == 0 && q >= milestone) {
                    double percentage = (static_cast<double>(i) / num_iterations_t) * 100.0;
                    std::cout << "\rProgress : "<< std::fixed << std::setprecision(2) << std::setw(6) << percentage << "%" << std::flush;
                    milestone += update;
                }
            }
        }
        for (int i = 0; i < threads; i++) {
            int szz = factor_threads[i].size();
            for (int j = 0; j < szz; j++) {
                unsigned long long curr_q = factor_threads[i][j];
                bool add = true;
                for (int k = 0; k < factors.size(); k++) {
                    if (curr_q % factors[k] == 0) {
                        add = false;
                        break;
                    }
                }
                if (add) factors.push_back(curr_q);
            }
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
    }
    // Print the cofactor to the screen, and possibly to output file if specified
    mpz_t mpz_product;
    mpz_init(mpz_product); // Initialize mpz_t
    mpz_cofactor(mpz_product, factors);
    char *cf = mpz_get_str(NULL, 10, mpz_product);
    string pm1_str = "-1)/";
    if (plus) pm1_str = "+1)/";
    string cofactor = "(" + std::to_string(b) + "^" + std::to_string(n) + pm1_str + cf;
    cout << "Remaining Cofactor:  " << cofactor << endl;
    if (filename!="") {
        std::ofstream output_file;
        output_file.open(filename, std::ios::app);
        if (!output_file) {
            std::cerr << "Error opening file for writing.\n";
        }
        else {
            output_file << cofactor << "\n";
            output_file.close();
        }
    }
    free(cf);
    return factors;
}

// Command line usage
void usage() {
    cout << "Lightweight Utility program to trial divide numbers of the form (b^n +- 1)/(b +- 1)\n"
         << "Usage: ./sieve -b <base> -n <exponent> -l <limit> [options]\n"
         << "\nRequired arguments:\n"
         << "  -b <base>          The base of the expression.\n"
         << "  -n <exponent>      The exponent of the expression.\n"
         << " OR alternatively, specify multiple base/exponent pairs in a file with -i\n"
         << "  -l <limit>         The upper limit for trial division.\n"
         << "\nOptions:\n"
         << "  -p                 Factor (b^n+1)/(b+1) instead of (b^n-1)/(b-1).\n"
         << "  -low <lower_bound> Start sieving from this number (default: 1).\n"
         << "  -t <threads>       Number of threads to use (default: 1).\n"
         << "  -i <filename>      Read from file instead of command line.\n"
         << "  One candidate per line : <base> <n>\n"
         << "  -f <filename>      Append the remaining cofactor to a file.\n"
         << "  --help             Show this help message.\n"
         << "Compile with : g++ -O3 -fopenmp trialdiv-demo.cpp -o trialdiv-demo -lgmp -lgmpxx\n"
         << "OR : g++ -std=c++17 -O3 -fopenmp trialdiv-demo.cpp -lgmpxx -lgmp -o trialdiv-demo\n"
         << "Example: ./trialdiv-demo -b 3 -n 200003 -l 100000000000000 -t 4 -f \"factor.txt\" \n";
}

// driver function
int main(int argc, char* argv[]) {
    unsigned long long base = 0, n = 0, limit = 0, lower = 1;
    bool plus = false;
    int num_threads = 1;
    string filename = "";
    string input_filename = "";
    bool inputf = false;

    //vector<unsigned int> mults = wheel(200003);
    //int x = calculate_iterations(0, 1000, 0, 510, mults);
    //for (unsigned int v : mults) {
    //    cout << v << ","; // << endl;
    //}
    //cout << "\n";
    //SieveState state = get_sieve_state(mults, 100000, 10000000000, 200003);
    //cout << "Start trial dividing at : " << state.start << endl;
    //cout << "Next Index to Increment : " << state.index << endl;
    //cout << "Number of Iterations : " << state.iters << endl;
    if (argc == 1) {
        usage();
        return 0;
    }

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--help") {
            usage();
            return 0;
        } else if (arg == "-b") {
            base = stoull(argv[++i]);
            if (inputf) {
                cerr << "Error: Cannot Specify both input file and command-line base/exponent" << endl;
                return 1;
            }
        } else if (arg == "-n") {
            if (inputf) {
                cerr << "Error: Cannot Specify both input file and command-line base/exponent" << endl;
                return 1;
            }
            n = stoull(argv[++i]);
        } else if (arg == "-l") {
            limit = stoull(argv[++i]);
        } else if (arg == "-p") {
            plus = true;
        } else if (arg == "-i") {
            if (base !=0 || n != 0) {
                cerr << "Error: Cannot Specify both input file and command-line base/exponent" << endl;
                return 1;                
            }
            input_filename = argv[++i];
            inputf = true;
        } else if (arg == "-low") {
            lower = stoull(argv[++i]);
        } else if (arg == "-t") {
            num_threads = stoi(argv[++i]);
        } else if (arg == "-f") {
            filename = argv[++i];
        } else {
            cerr << "Error: Unknown option '" << arg << "'. Use --help for usage." << endl;
            return 1;
        }
    }
    if (((base == 0 || n == 0) && !inputf) || limit == 0) {
        cerr << "Error: Missing required arguments. Use --help for usage." << endl;
        return 1;
    }

    if (inputf) {
        // open input file name for processing
        std::ifstream file(input_filename);
        if (!file.is_open()) {
            std::cerr << "Error opening input file." << std::endl;
            return 1; // Or handle the error appropriately
        }
        string line;
        while (std::getline(file, line)) {
            std::istringstream stream(line);
            std::string token;
            std::vector<std::string> args;
            // The extraction operator '>>' reads word by word,
            // automatically using whitespace as a delimiter.
            while (stream >> token) {
                args.push_back(token);
            }
            if (args.size() != 2) {
                cerr << "Error: Line " << line << " does not match the required format <base> <exponent>" << endl;
                continue;
            }
            base = stoull(args[0]);
            n = stoull(args[1]);
            cout << "\nFactoring (" << base << "^" << n;
            if (plus) cout << "+1)/(" << base + 1 << ")\n";
            else cout << "-1)/(" << base - 1 << ")\n";
            cout << "Sieve limit: " << limit << "\n";
            vector<unsigned long long> factors = sieve_repunit(base, n, lower, limit, plus, num_threads, filename);
        }
    }
    else {
        cout << "Factoring (" << base << "^" << n;
        if (plus) cout << "+1)/(" << base + 1 << ")\n";
        else cout << "-1)/(" << base - 1 << ")\n";
        cout << "Sieve limit: " << limit << "\n";
        // call the sieving function
        ////
        vector<unsigned long long> factors = sieve_repunit(base, n, lower, limit, plus, num_threads, filename);
        //vector<unsigned long long> factors = sieve_repunit(base, n, lower, limit, plus, filename);
    }
    return 0;
}
