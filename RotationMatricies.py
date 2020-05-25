import math


def rotation_matrix(xyz,phi):
    matrix = []
    if xyz == "x":
        matrix.append([1,0,0])
        matrix.append([0,math.cos(phi),-math.sin(phi)])
        matrix.append([0,math.sin(phi),math.cos(phi)])
    elif xyz == "y":
        matrix.append([math.cos(phi),0,math.sin(phi)])
        matrix.append([0,1,0])
        matrix.append([-math.sin(phi),0,math.cos(phi)])
    elif xyz == "z":
        matrix.append([math.cos(phi),-math.sin(phi),0])
        matrix.append([math.sin(phi),math.cos(phi),0])
        matrix.append([0,0,1])  
    return round_to_zero(matrix)

def round_to_zero(matrix):
    row_num = 0
    for row in matrix:
        col_num = 0
        for num in row:
            if abs(num) < 0.0001:
                matrix[row_num][col_num] = 0
            col_num +=1
        row_num +=1
    return matrix

def multiply_matricies(matrix_one,matrix_two):
    product = []
    for row_one in matrix_one:
        product_row = []
        for col_two_num in range(0,len(matrix_two[1])):
            value = 0
            for  col_one_num in range(0,len(matrix_one[1])):
                value += row_one[col_one_num]*matrix_two[col_one_num][col_two_num]
            product_row.append(value)
        product.append(product_row)
    return product  