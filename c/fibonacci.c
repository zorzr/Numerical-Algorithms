#include "fibonacci.h"
// NOTE: n must be a positive integer

/*
    Recursive calculation: really bad solution.
    O(n) stack occupation
    O(a^n) complexity (with a = 1.618)
*/
unsigned long fib_rec(int n) {
    if (n < 2)  return n;
    return fib_rec(n - 1) + fib_rec(n - 2);
}

/*
    Much nicer approach using dynamic programming.
    O(n) computational complexity
    O(1) memory usage (constant)
*/
unsigned long fib_dp(int n) {
    unsigned long fib[] = {1, 1};

    for (int i = 1; i < n; i++)
        fib[i % 2] = fib[0] + fib[1];
    return fib[n % 2];
}

/*
    Fast doubling: awesome for n big, slow for n small.
    Noticeable improvements from DP for n > 150
    O(log(n)) computational complexity
    O(1) memory usage (constant)
*/
unsigned long fib_fd(int n) {
    unsigned long a = 0, b = 1, t;

    for (int i = 31; i >= 0; i--) {
		t = a * a + b * b;
		a = a * (2 * b - a);
		b = t;

		if (((n >> i) & 1) != 0) {
			t = a + b;
			a = b;   b = t;
		}
	}
    return a;
}