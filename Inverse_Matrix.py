import math


def main():
    n = int(input())
    A = [[0.0 for i in range(n)] for j in range(n)]
    U = [[0.0 for i in range(n)] for j in range(n)]
    L = [[0.0 for i in range(n)] for j in range(n)]
    for i in range(n):
        A[i] = [float(x) for x in input().split()]
    print("A:")
    print_matrix(A)

    U[0] = A[0].copy()
    recursive_LU_factorization(A, 0, 0, L, U)
    reverse_matrix(A, L, U, n)

    print("L:")
    print_matrix(L)
    print("U:")
    print_matrix(U)
    print("A Inverse:")
    print_matrix(A)


def recursive_LU_factorization(matrix, p, q, L, U):
    L[p][q] = 1.0
    if p + 1 == len(matrix):
        return
    for i in range(p + 1, len(matrix)):
        scale = matrix[i][q] / matrix[p][q]
        L[i][q] = scale
        for j in range(q, len(matrix[0])):
            U[i][j] = matrix[i][j] - (matrix[p][j] * scale)
    recursive_LU_factorization(matrix, p + 1, q + 1, L, U)


def recursive_forward_substitution(L, p, q, b):
    if p + 1 == len(L):
        return
    for i in range(p + 1, len(L)):
        b[i] -= L[i][q] * b[p]
    recursive_forward_substitution(L, p + 1, q + 1, b)


def recursive_backward_substitution(U, p, q, y):
    if p < 0:
        return

    y[p] /= U[p][q]
    for i in reversed(range(0, p)):
        y[i] -= y[p] * U[i][q]
    recursive_backward_substitution(U, p - 1, q - 1, y)


def reverse_matrix(matrix, L, U, n):
    e = [[0.0 for i in range(n)] for j in range(n)]
    for k in range(n):
        e[k][k] = 1.0
    for i in range(n):
        recursive_forward_substitution(L, 0, 0, e[i])
        recursive_backward_substitution(U, n - 1, n - 1, e[i])
    for i in range(n):
        for j in range(n):
            matrix[i][j] = e[j][i]


def print_matrix(matrix):
    print("=================================")
    for i in range(len(matrix)):
        print("[ ", end="")
        for j in range(len(matrix[i])):
            if j != 0:
                if int(matrix[i][j]) == matrix[i][j] or abs(math.floor(matrix[i][j]) - matrix[i][j]) <= 0.00000000001:
                    print("%8.0f" % (matrix[i][j]), end="")
                else:
                    print("%8.2f" % (matrix[i][j]), end="")
            else:
                if int(matrix[i][j]) == matrix[i][j] or abs(math.floor(matrix[i][j]) - matrix[i][j]) <= 0.00000000001:
                    print("%3.0f" % (matrix[i][j]), end="")
                else:
                    print("%3.2f" % (matrix[i][j]), end="")
        print("  ]")
    print("=================================\n")


if __name__ == '__main__':
    main()
