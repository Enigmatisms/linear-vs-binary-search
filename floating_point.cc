#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <stdint.h>
#include <algorithm>
#include <random>

#ifdef _MSC_VER
    #define FORCEINLINE __forceinline
    #define NOINLINE __declspec(noinline)
    #define ALIGN(n) __declspec(align(n))
    FORCEINLINE uint32_t bsr(uint32_t x) {
        unsigned long res;
        _BitScanReverse(&res, x);
        return res;
    }
    FORCEINLINE uint32_t bsf(uint32_t x) {
        unsigned long res;
        _BitScanForward(&res, x);
        return res;
    }
#else
    #define FORCEINLINE __attribute__((always_inline)) inline
    #define NOINLINE __attribute__((noinline))
    #define ALIGN(n) __attribute__((aligned(n)))
    FORCEINLINE uint32_t bsr(uint32_t x) {
        return 31 - __builtin_clz(x);
    }
    FORCEINLINE uint32_t bsf(uint32_t x) {
        return __builtin_ctz(x);
    }
#endif

//if true, then average latency of one search is measured
//if false, then average throughput performance of one search is measured
//implementation-wise, setting it makes the index of the next array/key
//to be searched dependent on the answer of the current search
#define MEASURE_LATENCY true

//controls inlining of all search functions being tested
//NOINLINE means that inlining is forbidden:
//in this case benchmarking code is less likely to influence search performance
#define TESTINLINE NOINLINE

//======================= search implementations =======================

static TESTINLINE int binary_search_std (const float *arr, int n, float key) {
    return std::lower_bound(arr, arr + n, key) - arr;
}

static TESTINLINE int binary_search_simple (const float *arr, int n, float key) {
    intptr_t left = -1;
    intptr_t right = n;
    while (right - left > 1) {
        intptr_t middle = (left + right) >> 1;
        if (arr[middle] < key)
            left = middle;
        else
            right = middle;
    }
    return right;
}

intptr_t MINUS_ONE = -1;
static TESTINLINE int binary_search_branchless (const float *arr, int n, float key) {
    assert((n & (n+1)) == 0); //n = 2^k - 1
    //intptr_t pos = -1;            //generates "or r9, -1" on MSVC -- false dependency harms throughput
    intptr_t pos = MINUS_ONE;       //workaround for MSVC: generates mov without dependency
    intptr_t logstep = bsr(n);
    intptr_t step = intptr_t(1) << logstep;
    while (step > 0) {
        pos = (arr[pos + step] < key ? pos + step : pos);
        step >>= 1;
    }
    return pos + 1;
}

template<intptr_t MAXN> static TESTINLINE int binary_search_branchless_UR (const float *arr, int n, float key) {
    assert((n & (n+1)) == 0); //n = 2^k - 1
    assert(n+1 == MAXN);

    //intptr_t pos = -1;
    intptr_t pos = MINUS_ONE;
    #define STEP(logstep) \
        pos = (arr[pos + (1<<logstep)] < key ? pos + (1<<logstep) : pos);
    STEP(7)
    STEP(6)
    STEP(5)
    STEP(4)
    STEP(3)
    STEP(2)
    STEP(1)
    STEP(0)
    #undef STEP
    return pos + 1;
}

static TESTINLINE int binary_search_branchlessM (const float *arr, int n, float key) {
    assert((n & (n+1)) == 0); //n = 2^k - 1
    //intptr_t pos = -1;            //generates "or r9, -1" on MSVC -- false dependency harms throughput
    intptr_t pos = MINUS_ONE;       //workaround for MSVC: generates mov without dependency
    intptr_t logstep = bsr(n);
    intptr_t step = intptr_t(1) << logstep;
    while (step > 0) {
        pos += (arr[pos + step] < key) * step;
        step >>= 1;
    }
    return pos + 1;
}

static TESTINLINE int binary_search_branchlessA (const float *arr, int n, float key) {
    assert((n & (n+1)) == 0); //n = 2^k - 1
    //intptr_t pos = -1;            //generates "or r9, -1" on MSVC -- false dependency harms throughput
    intptr_t pos = MINUS_ONE;       //workaround for MSVC: generates mov without dependency
    intptr_t logstep = bsr(n);
    intptr_t step = intptr_t(1) << logstep;
    while (step > 0) {
        pos += (-(arr[pos + step] < key)) & step;
        step >>= 1;
    }
    return pos + 1;
}

static TESTINLINE int binary_search_branchless_pre (const float *arr, int n, float key) {
    assert((n & (n+1)) == 0); //n = 2^k - 1
    //intptr_t pos = -1;            //generates "or r9, -1" on MSVC -- false dependency harms throughput
    intptr_t pos = MINUS_ONE;       //workaround for MSVC: generates mov without dependency
    intptr_t logstep = bsr(n);
    intptr_t step = intptr_t(1) << logstep;
    float pivot = arr[pos + step];
    while (step > 1) {
        intptr_t nextstep = step >> 1;
        float pivotL = arr[pos + nextstep];
        float pivotR = arr[pos + step + nextstep];
        pos = (pivot < key ? pos + step : pos);
        pivot = (pivot < key ? pivotR : pivotL);
        step = nextstep;
    }
    pos = (pivot < key ? pos + step : pos);
    return pos + 1;
}

static TESTINLINE int quaternary_search_branchless (const float *arr, int n, float key) {
    assert((n & (n+1)) == 0); //n = 2^k - 1
    //intptr_t pos = -1;            //generates "or r9, -1" on MSVC -- false dependency harms throughput
    intptr_t pos = MINUS_ONE;       //workaround for MSVC: generates mov without dependency
    intptr_t logstep = bsr(n) - 1;
    intptr_t step = intptr_t(1) << logstep;
    while (step > 0) {
        float pivotL = arr[pos + step * 1];
        float pivotM = arr[pos + step * 2];
        float pivotR = arr[pos + step * 3];
        pos = (pivotL < key ? pos + step : pos);
        pos = (pivotM < key ? pos + step : pos);
        pos = (pivotR < key ? pos + step : pos);
        step >>= 2;
    }
    pos = (arr[pos + 1] < key ? pos + 1 : pos);
    return pos + 1;
}

//======================= testing code =======================

//length of each input array (including one sentinel element)
//must be power of two
static const int SIZE = 256;
//how many searches are done in benchmark
static const int TRIES = (1<<30) / 128;
//number of pre-generated input arrays (rotated cyclically)
static const int ARR_SAMPLES = (4<<20) / SIZE;
//number of pre-generated keys for search (rotated cyclically)
static const int KEY_SAMPLES = (4<<16);

//number of elements in every array to be searched (excluding sentinel element)
int n = SIZE - 1;

//input arrays for search (aligned for AVX)
ALIGN(32) float input[ARR_SAMPLES][SIZE];
//keys to be searched
float keys[KEY_SAMPLES];


int main() {
    //used RNG
    std::mt19937 rnd;
    std::uniform_real_distribution<float> distr(0, 1);

    //generate all input arrays
    for (int s = 0; s < ARR_SAMPLES; s++) {
        for (int i = 0; i < n; i++)
            input[s][i] = distr(rnd);
        std::sort(input[s], input[s] + n);
        //set sentinel element to INT_MAX
        for (int i = n; i < (n+15)/16*16; i++)
            input[s][i] = 10.f;
    }

    //test correctness of searches AND generate all keys to be searched
    const int ITERS = std::max(TRIES / 10, std::max(KEY_SAMPLES, ARR_SAMPLES)) + 10;
    for (int t = 0; t < ITERS; t++) {
        const float *arr = input[t % ARR_SAMPLES];
        float key = keys[t % KEY_SAMPLES] = distr(rnd);

        int res[16];
        int sk = 0;
        res[sk++] = binary_search_std(arr, n, key);
        res[sk++] = binary_search_simple(arr, n, key);
        res[sk++] = binary_search_branchless(arr, n, key);
        res[sk++] = binary_search_branchless_UR<SIZE>(arr, n, key);
        //some experimental implementations:
        res[sk++] = binary_search_branchlessM(arr, n, key);
        res[sk++] = binary_search_branchlessA(arr, n, key);
        res[sk++] = binary_search_branchless_pre(arr, n, key);
        res[sk++] = quaternary_search_branchless(arr, n, key);

        //program terminates if any search gives different answer
        for (int i = 1; i < sk; i++)
            if (res[i-1] != res[i]) {
                printf("ERROR: ");
                for (int j = 0; j < sk; j++)
                    printf(" %d", res[j]);
                printf("\n");
                exit(0);
            }
    }

    //print some info about current benchmark parameters
    static const int DARR = 1779033703 & (ARR_SAMPLES - 1);
    static const int DKEY = 2654435769 & (KEY_SAMPLES - 1);
    printf("Arrays: %d x %d   (+%d)\n", ARR_SAMPLES, SIZE, DARR);
    printf("Keys: %d        (+%d)\n", KEY_SAMPLES, DKEY);
    printf("Memory: %d B\n", int(sizeof(input) + sizeof(keys)));

    //benchmark environment
    //note: it could had been a function instead of macro
    //but originally I wanted to benchmark inlined code of every search
    #define TEST_SEARCH(func) \
    { \
        int start = clock(); \
        int check = 0; \
        for (int t = 0; t < TRIES; t++) { \
            int i = (t * DARR + (MEASURE_LATENCY ? check&1 : 0)) & (ARR_SAMPLES - 1); \
            int j = (t * DKEY + (MEASURE_LATENCY ? check&1 : 0)) & (KEY_SAMPLES - 1); \
            const float *arr = input[i]; \
            float key = keys[j]; \
            int res = func(arr, n, key); \
            check += res; \
        } \
        double elapsed = double(clock() - start) / CLOCKS_PER_SEC; \
        printf("%8.1lf ns : %40s   (%d)\n", 1e+9 * elapsed / TRIES, #func, check); \
    }

    //run performance benchmark and print formatted results
    //TEST_SEARCH(linearX_search_scalar);
    //TEST_SEARCH(linear_search_scalar);

    TEST_SEARCH(binary_search_std);
    TEST_SEARCH(binary_search_simple);
    TEST_SEARCH(binary_search_branchless);
    TEST_SEARCH(binary_search_branchless_UR<SIZE>);

    //some experimental implementations:

    TEST_SEARCH(binary_search_branchlessM);
    TEST_SEARCH(binary_search_branchlessA);

    TEST_SEARCH(binary_search_branchless_pre);
    TEST_SEARCH(quaternary_search_branchless);

    return 0;
}
