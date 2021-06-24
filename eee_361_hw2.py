import cv2 as cv
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import linalg
import math
import copy
import time

size_A_array = [1000, 500, 100]
tau_array = [0.01, 0.1]
a_array = []

for size in range(3):
    size_A = size_A_array[size]

    a_diag = np.zeros((size_A, size_A), int)
    np.fill_diagonal(a_diag, 1)
    a_temp = np.zeros((size_A, size_A), float)
    for x in range(size_A):
        for y in range(size_A):
            if x != y:
                a_temp[x, y] = np.random.uniform(-1, 1, 1)
    a = a_diag + a_temp

    for tau_location in range(2):
        tau = tau_array[tau_location]

        for x in range(size_A):
            for y in range(size_A):
                if x != y:
                    if np.abs(a[x, y]) > tau:
                        a[x, y] = 0
        a_array.append(a)

def find_errors(a, b, e_ij, x):
    a = a
    b = b
    e_ij = e_ij
    x = x
    Es_i = np.zeros(3)
    Eo_i = np.zeros(3)

    for i in range(3):
        for j in range(10):
            Es_i[i] += np.square(np.linalg.norm(e_ij[i, j]))
        Es_i[i] = np.sqrt(1/10 * Es_i[i])

    for i in range(3):
        for j in range(10):
            Eo_i[i] += np.square(np.linalg.norm(b[i, j] - a @ x[i, j]))
        Eo_i[i] = np.sqrt(1/10 * Eo_i[i])

    return Es_i, Eo_i

def plot_error_estimated_sol(Es_i):
    Es_i = Es_i
    plt.plot([1, 2, 3],Es_i[:3], marker='o')
    plt.plot([1, 2, 3],Es_i[3:6], marker='o')
    plt.plot([1, 2, 3],Es_i[6:9], marker='o')
    plt.plot([1, 2, 3],Es_i[9:12], marker='o')
    plt.plot([1, 2, 3],Es_i[12:15], marker='o')
    plt.plot([1, 2, 3],Es_i[15:18], marker='o')
    plt.legend(["A=10000, tau=0.01", "A=10000, tau=0.1", "A=500, tau=0.01", "A=500, tau=0.1", "A=100, tau=0.01",
                "A=100, tau=0.1"])
    plt.xlabel("Noise Level Index")
    plt.ylabel("Es_i")
    plt.yscale('log')
    plt.title("Es_i Plot")
    plt.show()

def plot_error_fit_observ(Eo_i):
    Eo_i = Eo_i
    plt.plot([1, 2, 3], Eo_i[0:3], marker='o')
    plt.plot([1, 2, 3], Eo_i[3:6], marker='o')
    plt.plot([1, 2, 3], Eo_i[6:9], marker='o')
    plt.plot([1, 2, 3], Eo_i[9:12], marker='o')
    plt.plot([1, 2, 3], Eo_i[12:15], marker='o')
    plt.plot([1, 2, 3], Eo_i[15:18], marker='o')
    plt.legend(["A=10000, tau=0.01", "A=10000, tau=0.1", "A=500, tau=0.01", "A=500, tau=0.1", "A=100, tau=0.01", "A=100, tau=0.1"])
    plt.xlabel("Noise Level Index")
    plt.ylabel("Eo_i")
    plt.title("Eo_i Plot")
    # plt.yscale('log')
    plt.show()

def a_pseudo_method(a, b):
    a_pseudo = np.linalg.pinv(a)
    x = a_pseudo @ b
    return x


def conjugate_gradient(S, b, x0):
    x = copy.deepcopy(x0)
    r = b
    d = r
    r_k_norm = np.dot(r.T, r)
    for i in range(S[0].size):
        Sd = np.dot(S, d)
        alpha = r_k_norm / np.dot(d.T, Sd)
        x += alpha * d
        r += -1 * alpha * Sd
        r_kplus1_norm = np.dot(r.T, r)
        beta = r_kplus1_norm / r_k_norm
        r_k_norm = r_kplus1_norm
        if r_kplus1_norm < 1e-15:
            # print('Itr:', i)
            break
        d = r + beta * d
    return x

def gmres_method(S, b, x0):
    q = b / np.linalg.norm(b)
    x = []
    h = np.zeros((S[0].size + 1, S[0].size))
    v = np.zeros([S[0].size, 1])
    # q = np.array(q).reshape([-1, 1])
    x_first = np.random.randn(S[0].size)
    e1 = np.zeros(S[0].size+1)
    e1[0] = 1
    for i in range(S[0].size):
        if i == 0:
            v = np.outer(S[:, i], q[i])
        else:
            v = S @ q

        for j in range(i):
            h[j, i] = np.dot(q.T, v)
            v = v - h[j, i] * q

        h[i + 1, i] = np.linalg.norm(v)

        if h[i + 1, i] != 0:
            q = v / h[i + 1, i]

        Q, R = QR_decomposition(h)
        res = Q.T @ e1
        y_k = np.linalg.norm(b)*np.linalg.inv(R) @ res
        if i == 30:
            break
    return x

# taken from https://github.com/danbar/qr_decomposition/blob/master/qr_decomposition/qr_decomposition.py
def QR_decomposition(A):
    """Gram-schmidt orthogonalization"""
    Q=np.zeros_like(A)
    for i in range(A.shape[1]):
        u = np.copy(A[:,i])
        for j in range(i):
            u -= np.dot(np.dot(Q[:, j].T, A[:,i]), Q[:, j])
        e = u / np.linalg.norm(u)
        Q[:, i] = e
    R = np.dot(Q.T, A)
    return (Q,R)


Es_i = np.zeros(3*6)
Eo_i = np.zeros(3*6)
r_arr = np.zeros(180)
for a_size in range(6):
    a = np.array(a_array[a_size])
    size_A = a[0].size
    sigma = [1, 0.01, 0.0001]
    x_0 = np.ndarray([10, size_A])
    x = np.ndarray([3, 10, size_A])
    x_intihal = np.ndarray([3, 10, size_A])
    e_ij = np.ndarray([3, 10, size_A])
    e_ij_intihal = np.ndarray([3, 10, size_A])
    b = np.ndarray([3, 10, size_A])
    for j in range(10):
        x_0[j, :] = np.random.randn(size_A)
        b_0 = a @ x_0[j]
        for i in range(3):
            w = sigma[i] * np.random.randn(size_A)
            b[i, j, :] = b_0 + w
            # x[i, j, :] = a_pseudo_method(a, b[i, j, :])
            # check = linalg.cg(A=a, b=b[i, j, :], x0=np.zeros([size_A])) // used for check purposes
            # check = linalg.gmres(A=a, b=b[i, j, :], x0=np.zeros([size_A])) // used for plotting and check purposes
            # x[i, j, :] = conjugate_gradient(S=a, b=b[i, j, :], x0=np.zeros([size_A]))
            x[i, j, :] = gmres_method(S=a, b=b[i, j, :], x0=x_0[j])
            # x[i, j, :] = gmres_method(S=a, b=b[i, j, :], x0=np.zeros([size_A])) // does not work properly so i used scipy library gmres function to plot
            # check = np.asarray(check[0])
            # x_intihal[i, j, :] = check
            e_ij[i, j, :] = x_0[j] - x[i, j]
            # e_ij_intihal[i, j, :] = x_0[j] - x_intihal[i, j]

    Es_i[3*a_size:3*a_size+3], Eo_i[3*a_size: 3*a_size+3] = find_errors(a, b, e_ij_intihal, x)
# plot_error_estimated_sol(Es_i)
plot_error_fit_observ(Eo_i)