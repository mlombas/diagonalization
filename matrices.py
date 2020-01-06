from polynomials import *

def minor(mat, i, j):
    result = []
    for index,row in enumerate(mat):
        if index != i:
            result.append(row[:j] + row[j + 1:])

    return result


def determinant(mat):
    if len(mat) != len(mat[0]): raise RuntimeError("Cannot take the determinant of a non-square matrix")

    if len(mat) == 1: return mat[0][0]

    result = 0
    row = mat[0]
    for j,num in enumerate(row):
        result += (-1)**j * num * determinant(minor(mat, 0, j))
    
    return result


def characteristic_polynomial(mat):
    if len(mat) != len(mat[0]): raise RuntimeError("Cannot calculate characteristic polynomial of a non-square matrix")

    result = [[el for el in row] for row in mat]

    for i in range(len(result)):
        result[i][i] = Poly([Term(-1,1), Term.numerical(mat[i][i])])

    return determinant(result)


def gauss(mat):
    copy = [[el for el in row] for row in mat]

    n = len(copy)
    m = len(copy[0])

    row = 0
    while row < n:
        if copy[row][row] != 0: 
            for other in range(row + 1, n): 
                multiplier = (copy[other][row] / copy[row][row]) 
                for i in range(m): 
                    copy[other][i] -= (multiplier * copy[row][i]) 

            row += 1
        else: 
            found = False
            for i in range(row + 1, m): 
                if copy[i][row] != 0: 
                    temp = copy[row]
                    copy[row] = copy[i]
                    copy[i] = temp

                    found = True
                    break

            if not found: 
                return copy
    
    return copy


def rank(mat):
    gaussed = gauss(mat)
    return len([row for row in gaussed if not all(el == 0 for el in row)])


def diagonalize(mat):
    if len(mat) != len(mat[0]): raise RuntimeError("Cannot diagonalize non-square matrix")

    positions = len(mat[0])
    result = []

    roots = calc_roots(characteristic_polynomial(mat))
    
    #If its not diagonalizable
    if len(roots) != len(mat):
        return None

    for pos,eigenvalue in enumerate(roots):
        row = [0] * positions
        row[pos] = eigenvalue

        result.append(row)

    return result


