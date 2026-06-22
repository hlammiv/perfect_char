import sympy as sym
import numpy as np
import itertools


c_i_j_dict = read_dictionary = np.load(
    './data/c_for_i_j_result.npy', allow_pickle='TRUE').item()

base_dim = 3
nb_character = 17
array_dim = (nb_character,) * base_dim
base_tensor = np.zeros(array_dim)  # store c_i_j

for i in range(17):
    for j in range(17):
        base_tensor[i][j] = c_i_j_dict[str(i+1) + '_' + str(j+1)]


def beta(i):
    string = 'beta_'+str(i)
    if i == 1:
        return sym.Symbol(string)-1
    else:
        return sym.Symbol(string)

# nth order c_i_j_...
def nth_order_c(n):
    if n == 2:
        return base_tensor
    elif n > 2:
        target_dim = n+1
        previous_tesnor = base_tensor
        for dim in range(3+1, target_dim+1):
            next_array_dim = (nb_character,) * dim
            next_tensor = np.zeros(next_array_dim)
            for idx in itertools.product(*[range(s) for s in (nb_character,) * (dim-1)]):
                index_except_last = idx[:dim-2]
                index_last = idx[-1]
                temp_array = np.zeros(nb_character)
                for counter in range(nb_character):
                    temp_array += previous_tesnor[index_except_last][counter] * \
                        base_tensor[counter][index_last]
                next_tensor[idx] = temp_array
            previous_tesnor = next_tensor
        return next_tensor

# nth order coefficients for all chi_i
def nth_order_coefficients(n):
    target_order = n
    result = []

    # base result
    for i in range(17):
        result.append(beta(i+1))
    result = np.array(result)

    for i in range(1, target_order):
        order = i+1
        c_tensor = nth_order_c(order)
        for idx in itertools.product(*[range(s) for s in (nb_character, ) * order]):
            b = 1/order * (-1)**(order+1)
            for i in idx:
                b *= beta(i+1)
            result += b * c_tensor[idx]

    final_result=[]
    for unsimplified in result:
        final_result.append(sym.simplify(sym.expand(unsimplified)))

    return final_result
