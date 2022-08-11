from math import fabs, ceil, sqrt, exp, log, inf
import random
from itertools import chain
from sympy import *
import time
import numpy as np
import matplotlib.pyplot as plt


def prime_generation(n):

    if n < 2:
        return []

    nums = []
    isPrime = []

    for i in range(0, n + 1):
        nums.append(i)
        isPrime.append(True)

    isPrime[0] = False
    isPrime[1] = False

    for j in range(2, int(n / 2)):
        if isPrime[j] == True:
            for i in range(2 * j, n + 1, j):
                isPrime[i] = False

    primes = []
    for i in range(0, n + 1):
        if isPrime[i] == True:
            primes.append(nums[i])

    return primes


def quadratic_residue(a, n):

    l = 1
    q = (n - 1) // 2
    x = q ** l
    if x == 0:
        return 1

    a %= n
    z = 1
    while x != 0:
        if x % 2 == 0:
            a = (a ** 2) % n
            x //= 2
        else:
            x -= 1
            z = (z * a) % n

    return z


def shanks_tonelli(n, p):

    assert quadratic_residue(n, p) == 1, "not a square (mod p)"
    q = p - 1
    s = 0

    while q % 2 == 0:
        q //= 2
        s += 1

    if s == 1:
        r = pow(n, (p + 1) // 4, p)
        return r, p - r

    for z in range(2, p):
        if p - 1 == quadratic_residue(z, p):
            break

    c = pow(z, q, p)
    r = pow(n, (q + 1) // 2, p)
    t = pow(n, q, p)
    m = s
    t2 = 0

    while (t - 1) % p != 0:
        t2 = (t * t) % p
        for i in range(1, m):
            if (t2 - 1) % p == 0:
                break
            t2 = (t2 * t2) % p
        b = pow(c, 1 << (m - i - 1), p)
        r = (r * b) % p
        c = (b * b) % p
        t = (t * c) % p
        m = i

    return r, p - r


def get_factor_base(N, B):

    factor_base = []
    primes = prime_generation(B)

    for p in primes:
        if quadratic_residue(N, p) == 1:
            factor_base.append(p)

    return factor_base


def get_smooth_numbers(factor_base, N, I):

    sieve_seq = [x ** 2 - N for x in range(root, root + I)]

    sieve_list = sieve_seq.copy()

    if factor_base[0] == 2:
        i = 0
        while sieve_list[i] % 2 != 0:
            i += 1
        for j in range(i, len(sieve_list), 2):
            while sieve_list[j] % 2 == 0:
                sieve_list[j] //= 2

    for p in factor_base[1:]:
        residues = shanks_tonelli(N, p)

        for r in residues:
            for i in range((r - root) % p, len(sieve_list), p):
                while sieve_list[i] % p == 0:
                    sieve_list[i] //= p

    xlist = []
    smooth_numbers = []
    indices = []

    for i in range(len(sieve_list)):
        if len(smooth_numbers) >= len(factor_base) + T:
            break
        if sieve_list[i] == 1:
            smooth_numbers.append(sieve_seq[i])
            xlist.append(i + root)
            indices.append(i)

    return smooth_numbers, xlist, indices


def factor(n, factor_base):

    factors = []
    if n < 0:
        factors.append(-1)

    for p in factor_base:
        if p == -1:
            pass
        else:
            while n % p == 0:
                factors.append(p)
                n //= p

    return factors


def generate_matrix(smooth_numbers, factor_base):

    M = []
    factor_base.insert(0, -1)

    for n in smooth_numbers:
        exponent_vector = [0] * (len(factor_base))
        n_factors = factor(n, factor_base)

        for i in range(len(factor_base)):
            if factor_base[i] in n_factors:
                exponent_vector[i] = (
                    exp_vector[i] + n_factors.count(factor_base[i])
                ) % 2

        if 1 not in exponent_vector:
            return True, n
        else:
            pass

        M.append(exponent_vector)

    return False, transpose(M)


def transpose(matrix):

    transposed_matrix = []
    for i in range(len(matrix[0])):
        new_row = []
        for row in matrix:
            new_row.append(row[i])
        transposed_matrix.append(new_row)

    return transposed_matrix


def gauss_elimination(M):

    marks = [False] * len(M[0])

    for i in range(len(M)):
        row = M[i]

        for num in row:
            if num == 1:
                j = row.index(num)
                marks[j] = True

                for k in chain(range(0, i), range(i + 1, len(M))):
                    if M[k][j] == 1:
                        for i in range(len(M[k])):
                            M[k][i] = (M[k][i] + row[i]) % 2
                break

    M = transpose(M)

    row_solution = []
    for i in range(len(marks)):
        if marks[i] == False:
            free_row = [M[i], i]
            row_solution.append(free_row)

    return row_solution, marks, M


def solve_row(row_solution, M, marks, K=0):

    vector_solution, indices = [], []
    free_row = row_solution[K][0]

    for i in range(len(free_row)):
        if free_row[i] == 1:
            indices.append(i)

    for r in range(len(M)):
        for i in indices:
            if M[r][i] == 1 and marks[r]:
                vector_solution.append(r)
                break

    vector_solution.append(row_solution[K][1])

    return vector_solution


def solve(vector_solution, smooth_numbers, xlist, N):

    solution_nums = [smooth_numbers[i] for i in vector_solution]
    x_nums = [xlist[i] for i in vector_solution]

    square = 1
    for n in solution_nums:
        square *= n

    b = 1
    for n in x_nums:
        b *= n

    a = sqrt(square)

    factor = gcd(b - a, N)

    return factor


# Pollard p-1 algorithm
def pollard(N):

    a = 2
    B = floor(sqrt(N)) + 1

    for i in range(2, B):
        a = (a ** i) % N
        d = gcd(a - 1, N)
        if d > 1 and d < N:
            return d, N // d


# Quadratic sieve
def quadratic_sieve(n, B, I):

    global N
    global root
    global T
    N, root, K, T = n, int(sqrt(n)), 0, 1

    factor_base = get_factor_base(N, B)

    global F
    F = len(factor_base)

    smooth_numbers, xlist, indices = get_smooth_numbers(factor_base, N, I)

    is_square, t_matrix = generate_matrix(smooth_numbers, factor_base)

    if is_square == True:
        x = smooth_numbers.index(t_matrix)
        factor = gcd(xlist[x] + sqrt(t_matrix), N)
        return factor, N / factor

    row_solution, marks, M = gauss_elimination(t_matrix)
    vector_solution = solve_row(row_solution, M, marks, 0)

    factor = solve(vector_solution, smooth_numbers, xlist, N)

    for K in range(1, len(row_solution)):
        if factor == 1 or factor == N:
            vector_solution = solve_row(row_solution, M, marks, K)
            factor = solve(vector_solution, smooth_numbers, xlist, N)
        else:
            return factor, int(N / factor)


def main_calc(N):

    B = floor(sqrt(N)) + 1
    I = B * 3
    print("Factoring {}".format(N))
    print("...")

    start1 = time.time()
    x1 = pollard(N)

    end1 = time.time()
    print("Pollard result {}".format(x1))
    print("Pollard time: {0:.4f}".format(end1 - start1))

    if N > 10 ** 14:
        return end1 - start1, inf

    start2 = time.time()
    x2 = quadratic_sieve(N, B, I)
    end2 = time.time()
    print("QS result: {}".format(x2))
    print("QS time: {0:.4f}".format(end2 - start2))

    print("====================")
    return end1 - start1, end2 - start2


if __name__ == "__main__":

    N1 = 83177 * 28909
    N2 = 970061 * 931097
    N3 = 51199 * 25541
    N4 = 61051 * 95261
    N5 = 422621 * 18149
    N6 = 44101 * 200153
    N7 = 144271 * 72869
    N8 = 84047 * 216851
    N9 = 148639 * 58211
    N10 = 144611 * 74297
    N11 = 766916207 * 85560361
    N12 = 7280783129 * 78263737

    a1, b1 = main_calc(N1)
    a2, b2 = main_calc(N2)
    a3, b3 = main_calc(N3)
    a4, b4 = main_calc(N4)
    a5, b5 = main_calc(N5)
    a6, b6 = main_calc(N6)
    a7, b7 = main_calc(N7)
    a8, b8 = main_calc(N8)
    a9, b9 = main_calc(N9)
    a10, b10 = main_calc(N10)
    a11, b11 = main_calc(N11)
    a12, b12 = main_calc(N12)

    pollard = np.mean([a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12])
    quadratic_sieve = np.mean([b1, b2, b3, b4, b5, b6, b7, b8, b9, b10])

    print(pollard)
    print(quadratic_sieve)

    plt.plot(
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
        np.log([a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12]),
        linestyle="-",
        label="Pollard's p-1",
        color="c",
    )
    plt.plot(
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
        np.log([b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12]),
        linestyle="--",
        label="Quadratic sieve",
        color="m",
    )
    plt.legend()
    plt.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)
    plt.ylabel("Logarithmic time (s)")
    plt.xlabel("Integers factored")
    plt.show()
