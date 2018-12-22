#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "fibonacci.h"
#define ERR_0 "Error: not enough input arguments\n"
#define ERR_1 "Error: the index must be positive\n"
#define ERR_2 "Error: not enough bits to store result\n"
#define PRINT_ERR(str, err) { printf(str);  return err; }

void test_rec(int n);
void test_dp(int n);
void test_fd(int n);

int main(int argc, char **argv) {
    int n;

    if (argc > 1)  n = atoi(argv[1]);
    else  PRINT_ERR(ERR_0, -1)
    if (n < 1)  PRINT_ERR(ERR_1, 1)
    else if (n > 47)  PRINT_ERR(ERR_2, 2)

    test_rec(n);
    test_dp(n);
    test_fd(n);
    return 0;
}


// Tests
void test_rec(int n) {
    clock_t start, end;
    double cpu_time;
    unsigned long f;

    start = clock();
    f = fib_rec(n);
    end = clock();
    cpu_time = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Recursive: %lu\nElapsed time: %gs\n\n", f, cpu_time);
}

void test_dp(int n) {
    clock_t start, end;
    double cpu_time;
    unsigned long f;

    start = clock();
    f = fib_dp(n);
    end = clock();
    cpu_time = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Dynamic prog.: %lu\nElapsed time: %gs\n\n", f, cpu_time);
}

void test_fd(int n) {
    clock_t start, end;
    double cpu_time;
    unsigned long f;

    start = clock();
    f = fib_fd(n);
    end = clock();
    cpu_time = (end - start) / (double) CLOCKS_PER_SEC;
    printf("Fast doubling: %lu\nElapsed time: %gs\n\n", f, cpu_time);
}