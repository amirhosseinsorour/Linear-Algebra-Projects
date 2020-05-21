import math


def main():
    n = int(input())
    A = [[0.0 for i in range(n + 1)] for j in range(n)]
    for i in range(n):
        A[i] = [float(x) for x in input().split()] + [0.0]
    b = [float(x) for x in input().split()]
    for i in range(n):
        A[i][n] = b[i]
    print_matrix(A)
    recursive_echelon(A, 0, leftmost_nonzero_column(A))
    recursive_reduced_echelon(A, n - 1, n - 1)
    print_results(A)
    # print(b)


step = 0


def print_matrix(matrix):
    print("Step %d :" % step)
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


def print_results(matrix):
    print("Results:")
    print("=================================")
    for n in range(len(matrix)):
        print("x%d = %g" % (n + 1, matrix[n][len(matrix[0]) - 1]))


def interchange_row(matrix, i, j):
    tmp_list = matrix[i]
    matrix[i] = matrix[j]
    matrix[j] = tmp_list


def leftmost_nonzero_column(matrix):
    isZero = True
    for j in range(len(matrix[0])):
        for i in range(len(matrix)):
            if matrix[i][j] != 0:
                isZero = False
                break
        if not isZero:
            return j
    return -1


def recursive_echelon(matrix, p, q):
    if p + 1 == len(matrix):
        return
    row = p
    while matrix[row][q] == 0:
        row += 1
        if row == len(matrix):
            print("No Solution or Many Solutions!")
            exit(0)
    if row != p:
        interchange_row(matrix, p, row)
    for i in range(p + 1, len(matrix)):
        scale = matrix[i][q] / matrix[p][q]
        for j in range(q, len(matrix[0])):
            matrix[i][j] -= matrix[p][j] * scale
        global step
        step += 1
        print_matrix(matrix)
    recursive_echelon(matrix, p + 1, q + 1)


def recursive_reduced_echelon(matrix, p, q):
    if p < 0:
        return
    if matrix[p][q] == 0:
        print("No Solution or Many Solutions!")
        exit(0)

    m = len(matrix[0])
    matrix[p][m - 1] /= matrix[p][q]
    matrix[p][q] = 1
    for i in reversed(range(0, p)):
        matrix[i][m - 1] -= matrix[p][m - 1] * matrix[i][q]
        matrix[i][q] = 0
    global step
    step += 1
    print_matrix(matrix)
    recursive_reduced_echelon(matrix, p - 1, q - 1)


if __name__ == '__main__':
    main()
