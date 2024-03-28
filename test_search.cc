#include <chrono>
#include <random>
#include <algorithm>

using namespace std;

/**
 * runtime:
 * memory:
 */

#define ALIGN(n) __attribute__((aligned(n)))

int binary_search_branchless_UR (const float *arr, float key) {
    // best latency performance, twice faster than std::lower_bound
    intptr_t pos = -1;
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

class TicTocMicro {
public:
    void tic() {
        start = std::chrono::system_clock::now();
    }

    int toc() const noexcept {
        auto end = std::chrono::system_clock::now();
        return std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    }
private:
    std::chrono::system_clock::time_point start;
};

int main(int argc, char** argv) {

    ALIGN(32) float arr[256];

    std::mt19937 rnd;
    std::uniform_real_distribution<float> distr(0, 1);

    for (int i = 0; i < 256; i++)
        arr[i] = distr(rnd);
    std::sort(arr, arr + 256);
    arr[255] = 1;


    TicTocMicro timer;
    float query = 0;
    for (int i = 0; i < 20; i++) {
        timer.tic();
        int result = binary_search_branchless_UR(arr, query);
        int micro_s = timer.toc();
        printf("(%d) Floating point: %f, searched: [%f, %f], index: %d\n", micro_s, query, arr[result], arr[result + 1], result);
        timer.tic();
        int std_result = std::lower_bound(arr, arr + 255, query) - arr;
        micro_s = timer.toc();
        printf("(%d) std: [%f, %f], %d \n", micro_s, arr[result], arr[result + 1], result);
        query = distr(rnd);
    }
    return 0;
}