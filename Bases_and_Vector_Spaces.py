import copy
import math


def main():
    n, m = input().split()
    n = int(n)
    m = int(m)
    A = [[0.0 for i in range(m + 1)] for j in range(n)]
    for i in range(n):
        A[i] = [float(x) for x in input().split()] + [0.0]
    print("Matrix A:")
    print_matrix(A, n, m)
    reduced_echelon_A = copy.deepcopy(A)
    recursive_reduced_echelon(reduced_echelon_A, 0, leftmost_nonzero_column(A))
    global first_echelon
    first_echelon = False
    print("Reduced Echelon Matrix [ A | 0 ]:")
    print_matrix(reduced_echelon_A, n, m + 1)
    null_bases(reduced_echelon_A, m)
    row_bases(reduced_echelon_A, n, m)
    col_bases(A, n)
    linear_comb(A)


pivot_col = []
non_pivot_col = []


def print_matrix(matrix, n, m):
    print("=================================")
    for i in range(n):
        print("[ ", end="")
        for j in range(m):
            if j != 0:
                if int(matrix[i][j]) == matrix[i][j] or abs(
                        math.floor(abs(matrix[i][j])) - abs(matrix[i][j])) <= 0.00000000001:
                    print("%8.0d" % (int(matrix[i][j])), end="")
                else:
                    print("%8.2f" % (matrix[i][j]), end="")
            else:
                if int(matrix[i][j]) == matrix[i][j] or abs(
                        math.floor(abs(matrix[i][j])) - abs(matrix[i][j])) <= 0.00000000001:
                    print("%3.0d" % (int(matrix[i][j])), end="")
                else:
                    print("%3.2f" % (matrix[i][j]), end="")
        print("  ]")
    print("=================================\n")


def print_bases(matrix):
    print("=================================")
    for j in range(len(matrix[0])):
        print("[ ", end="")
        for i in range(len(matrix)):
            if i != 0:
                if int(matrix[i][j]) == matrix[i][j] or abs(
                        math.floor(abs(matrix[i][j])) - abs(matrix[i][j])) <= 0.00000000001:
                    print("]    [%8.0d" % (int(matrix[i][j])), end="")
                else:
                    print("]    [%8.2f" % (matrix[i][j]), end="")
            else:
                if int(matrix[i][j]) == matrix[i][j] or abs(
                        math.floor(abs(matrix[i][j])) - abs(matrix[i][j])) <= 0.00000000001:
                    print("%3.0d" % (int(matrix[i][j])), end="")
                else:
                    print("%3.2f" % (matrix[i][j]), end="")
        print("  ]")
    print("=================================\n")


def print_results(matrix):
    print("Results:")
    print("=================================")
    for n in range(len(matrix)):
        print("x%d = %g" % (n + 1, matrix[n][len(matrix[0]) - 1]))


def print_linear_comb(matrix, scales):
    global pivot_col
    print_vector(matrix[len(matrix) - 1])
    print(" = ", end="")
    for j in range(len(pivot_col)):
        print("%g" % (scales[j]), end="")
        print_vector(matrix[j])
        if j == len(pivot_col) - 1:
            print("")
        else:
            print(" + ", end="")


def print_vector(vector):
    print("(", end="")
    for i in range(len(vector)):
        if i != len(vector) - 1:
            print("%g, " % vector[i], end="")
        else:
            print("%g)" % vector[i], end="")


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


first_echelon = True


def recursive_reduced_echelon(matrix, p, q):
    if int(math.modf(matrix[p][q])[1]) == 0:
        matrix[p][q] = 0.0
    if p + 1 == len(matrix) and q + 1 == len(matrix[0]):
        return
    global first_echelon
    row = p
    while matrix[row][q] == 0:
        row += 1
        if row == len(matrix):
            if first_echelon:
                global non_pivot_col
                non_pivot_col.append(q)
            recursive_reduced_echelon(matrix, p, q + 1)
            return
    if row != p:
        interchange_row(matrix, p, row)
    scale = matrix[p][q]
    for j in range(q, len(matrix[0])):
        matrix[p][j] /= scale
    for i in range(0, len(matrix)):

        if i == p:
            continue
        scale = matrix[i][q]
        for j in range(q, len(matrix[0])):
            matrix[i][j] -= matrix[p][j] * scale
    if first_echelon:
        global pivot_col
        pivot_col.append(q)
    recursive_reduced_echelon(matrix, p + 1, q + 1)


def null_bases(matrix, m):
    global non_pivot_col, pivot_col
    null_base = [[0.0 for i in range(m)] for j in range(len(non_pivot_col))]
    for j in range(len(non_pivot_col)):
        null_base[j][non_pivot_col[j]] = 1.0
    for i in range(len(non_pivot_col)):
        for (j) in range(len(pivot_col)):
            null_base[i][pivot_col[j]] = matrix[j][non_pivot_col[i]] * -1.0
    print("Bases for null space:")
    print_bases(null_base)


def row_bases(matrix, n, m):
    row_base = []
    for i in range(n):
        isZero = True
        for j in range(m):
            if int(matrix[i][j]) != 0:
                isZero = False
        if not isZero:
            row_base.append(matrix[i])
            row_base[i].pop(m)
    print("Bases for row space:")
    print_bases(row_base)


def col_bases(matrix, n):
    global pivot_col
    col_base = [[0.0 for i in range(n)] for j in range(len(pivot_col))]
    for i in range(len(pivot_col)):
        for j in range(n):
            col_base[i][j] = matrix[j][pivot_col[i]]
    print("Bases for column space:")
    print_bases(col_base)


def linear_comb(matrix):
    print("Non pivot columns as a linear combination of pivot columns :")
    print("=================================")
    global pivot_col, non_pivot_col
    pivot_col_matrix = []
    for j in pivot_col:
        pivot_col_matrix.append(transpose(matrix)[j])
    for k in non_pivot_col:
        tmp_matrix = copy.deepcopy(pivot_col_matrix)
        tmp_matrix.append(transpose(matrix)[k])
        tmp_matrix_t = transpose(tmp_matrix)
        recursive_reduced_echelon(tmp_matrix_t, 0, 0)
        print_linear_comb(tmp_matrix, transpose(tmp_matrix_t)[len(tmp_matrix_t[0]) - 1])


def transpose(matrix):
    transposed_matrix = [[0.0 for i in range(len(matrix))] for j in range(len(matrix[0]))]
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            transposed_matrix[j][i] = matrix[i][j]
            if int(math.modf(transposed_matrix[j][i])[1]) == 0:
                transposed_matrix[j][i] = 0.0
    return transposed_matrix


if __name__ == '__main__':
    main()
