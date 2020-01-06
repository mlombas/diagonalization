from polynomials import *
from matrices import *

def print_matrix(mat):
    max_length = max(
            max(len(str(element)) for element in row) for row in mat
    )

    for row in mat:
        print("".join(f"{str(el):<{max_length + 1}}" for el in row))


def read_matrix():
    matrix = []
    print("Enter the elements of the matrix, separated by spaces. Enter an empty line to end")

    while (line := input()):
        matrix.append([float(num) for num in line.split(" ")])

    return matrix

def main():
    mat = read_matrix()

    if any(len(mat) != len(row) for row in mat):
        print("Cannot diagonalize a non-square matrix")
        return

    poly = characteristic_polynomial(mat)
    roots = calc_roots(poly)
    multiplicities = {root: roots.count(root) for root in set(roots)}

    print(f"Let us call A to the matrix you entered")
    print(f"Calculate the characteristic polynomial, pA(x) = |A - xI|, which is {poly}, has degree {poly.get_degree()}")
    
    roots_text = ", ".join(f"{root} (multiplicity {mult})" for root, mult in multiplicities.items())
    print(f"Calculate its roots, which are {roots_text}")
    print(f"Thus, it is diagonalizable over {'R' if all(not isinstance(root, complex) for root in roots) else 'C'}")

    print(f"Let us check the dimension of the eigenspaces")
    for root,mult in multiplicities.items():
        if mult == 1:
            print(f"dim(S({root})) = 1 as its multiplicity as a root of pA(x) is 1")
        else:
            rank_matrix = [[el for el in row] for row in mat]
            for i in range(len(mat)):
                rank_matrix[i][i] -= root

            r = rank(rank_matrix)
            dim = len(mat) - r

            if dim == mult:
                print(f"dim(S({root})) = {len(mat)} - rank(A - {root}I) = {len(mat)} - {r} = {dim}, which is its multiplicity as a root of pA(x)")
            else:
                print(f"dim(S({root})) = {len(mat)} - rank(A - {root}I) = {len(mat)} - {r} = {dim}, but its multiplicity as a root of pA(x) is {mult}, so A is not diagonalizable")
                return

    print(f"Then, A is diagonalizable and D is:")
    print_matrix(diagonalize(mat))

main()

