import matplotlib.pyplot as plt
import numpy as np
import math


def save_current_sample(sample, file_path):
    file = open(file_path, 'w')
    for i in range(sample.shape[0]):
        for j in range(sample.shape[1]):
            file.write(str(sample[i][j]))
            if j != (sample.shape[1] - 1): file.write(',')
        if i != (sample.shape[0] - 1): file.write('\n')
    file.close()


N = 100

# 1
sample_file = open('./stat_range.csv', 'r')
lines = sample_file.readlines()
stat_range = np.ndarray((len(lines), 11))
for line in range(len(lines)):
    elements = lines[line].split(',')
    for i in range(stat_range.shape[1]): stat_range[line][i] = float(elements[i])
sample_file.close()

data_ind_seq = np.lexsort((stat_range[:, 10], stat_range[:, 9], stat_range[:, 8], stat_range[:, 7],
                           stat_range[:, 6], stat_range[:, 5], stat_range[:, 4], stat_range[:, 3],
                           stat_range[:, 2], stat_range[:, 0], stat_range[:, 1]))
ranked_sample = stat_range.copy()
for i in range(stat_range.shape[0]):
    for j in range(stat_range.shape[1]):
        ranked_sample[i][j] = stat_range[int(data_ind_seq[i])][j]
save_current_sample(ranked_sample, './ranked_range_2.csv')

var_sample = np.ndarray((0, ranked_sample.shape[1]))
var_sample_freq = np.ndarray((0, 2))
for i in range(ranked_sample.shape[0]):
    if i == 0 or not np.array_equal(ranked_sample[i], ranked_sample[i - 1]):
        var_sample = np.vstack((var_sample, ranked_sample[i]))
        var_sample_freq = np.vstack((var_sample_freq, [1, 1 / N]))
    else:
        var_sample_freq[var_sample_freq.shape[0] - 1][0] += 1
        var_sample_freq[var_sample_freq.shape[0] - 1][1] = var_sample_freq[var_sample_freq.shape[0] - 1][0] / N
save_current_sample(var_sample, './var_range_2.csv')
save_current_sample(var_sample_freq, './var_range_freq_2.csv')

interval_range_y = np.ndarray((0, 5))
k = 1 + 3.31 * np.log10(N)
if int(k) % 2 == 1:
    k = int(k)
else:
    k = int(k) + 1
h = (max(var_sample[:, 1]) - min(var_sample[:, 1])) / k
print('k =', k, 'h =', h)
for i in range(1, k + 1):
    int_sample = np.vstack((interval_range_y, [min(var_sample[:, 1]) + (h * (i - 1)), min(var_sample[:, 1]) + (h * i),
                                         (min(var_sample[:, 1]) + (h * (i - 1)) + min(var_sample[:, 1]) + (h * i)) / 2,
                                         0, 0]))
current_int_ind = 0
for i in range(len(var_sample)):
    if var_sample[i][1] >= interval_range_y[current_int_ind][1]: current_int_ind += 1
    if current_int_ind >= k: current_int_ind = k - 1
    interval_range_y[current_int_ind][3] += var_sample_freq[i][0]
    interval_range_y[current_int_ind][4] += var_sample_freq[i][1]
save_current_sample(interval_range_y, './interval_range_2.csv')

plt.clf()
plt.title('Abs. frequency polygon')
plt.plot(interval_range_y[:, 2], interval_range_y[:, 3], 'bo-')
plt.xticks(interval_range_y[:, 2], np.around(interval_range_y[:, 2], 2))
plt.xlabel('Middle of interval'); plt.ylabel('Abs. frequency')
plt.grid()
plt.show()

plt.clf()
plt.title('Rel. frequency polygon')
plt.plot(interval_range_y[:, 2], interval_range_y[:, 4], 'bo-')
plt.xticks(interval_range_y[:, 2], np.around(interval_range_y[:, 2], 2))
plt.xlabel('Middle of interval'); plt.ylabel('Rel. frequency')
plt.grid()
plt.show()

plt.clf()
plt.title('Abs. frequency histogram')
plt.bar(interval_range_y[:, 2], height=interval_range_y[:, 3] / h, width=h, color='b', tick_label=np.around(interval_range_y[:, 2], 2))
plt.xlabel('Middle of interval'); plt.ylabel('n / h')
plt.grid()
plt.show()

plt.clf()
plt.title('Rel. frequency histogram')
plt.bar(interval_range_y[:, 2], height=interval_range_y[:, 4] / h, width=h, color='b', tick_label=np.around(interval_range_y[:, 2], 2))
plt.xlabel('Middle of interval'); plt.ylabel('n_rel. / h')
plt.grid()
plt.show()

plt.clf()
plt.title('Empirical distribution function')
points = np.ndarray((0, 2))
for i in range(k):
    accumulated_freq = 0
    for j in range(i): accumulated_freq += (interval_range_y[j][3] / N)
    points = np.vstack((points, [[interval_range_y[i][2], accumulated_freq],
                                 [interval_range_y[i][2], accumulated_freq + (interval_range_y[i][3] / N)]]))
plt.plot(points[:, 0], points[:, 1], 'bo-')
plt.xticks(points[:, 0], np.around(points[:, 0], 2))
plt.xlabel('Middle of interval'); plt.ylabel('Accumulated frequency')
plt.grid()
plt.show()

plt.clf()
plt.title('Empirical distribution function')
points = np.ndarray((0, 2))
for i in range(k):
    accumulated_freq = 0
    for j in range(i): accumulated_freq += (interval_range_y[j][4] / N)
    points = np.vstack((points, [[interval_range_y[i][2], accumulated_freq],
                                 [interval_range_y[i][2], accumulated_freq + (interval_range_y[i][4] / N)]]))
plt.plot(points[:, 0], points[:, 1], 'bo-')
plt.xticks(points[:, 0], np.around(points[:, 0], 2))
plt.xlabel('Middle of interval'); plt.ylabel('Accumulated frequency')
plt.grid()
plt.show()

# 2
conditional_options = np.arange(k)
C = interval_range_y[k // 2][2]
for interval in range(k):
    u = int(((interval_range_y[interval][2] - C) / h + 0.1).round())
    conditional_options[interval] = u
print('Условные варианты:', conditional_options)

M = np.zeros(4)
for i in range(4):
    for j in range(k):
        M[i] += interval_range_y[j][3] * pow(conditional_options[j], i + 1) / N
print('Условные эмпирические моменты: M1 =', M[0], 'M2 =', M[1], 'M3 =', M[2], 'M4 =', M[3])

m1 = M[0] * h + C
m2 = (M[1] - pow(M[0], 2)) * pow(h, 2)
m3 = (M[2] - 3 * M[1] * M[0] + 2 * pow(M[0], 3)) * pow(h, 3)
m4 = (M[3] - 4 * M[2] * M[0] + 6 * M[1] * pow(M[0], 2) - 3 * pow(M[0], 4)) * pow(h, 4)
print('Центральные эмпирические моменты: m1 =', m1, 'm2 =', m2, 'm3 =', m3, 'm4 =', m4)

y_sample_mean_u = m1; D_u = m2
print('Выборочное среднее, вычисленное с помощью усл. вариант:', y_sample_mean_u)
print('Дисперсия, вычисленная с помощью усл. вариант:', D_u)

y_sample_mean = 0.0
for i in range(k): y_sample_mean += interval_range_y[i][3] * interval_range_y[i][2] / N
D = 0.0
for i in range(k): D += interval_range_y[i][3] * pow(interval_range_y[i][2] - y_sample_mean, 2) / N
print('Выборочное среднее, вычисленное с помощью стандартной формулы:', y_sample_mean)
print('Дисперсия, вычисленная с помощью стандартной формулы:', D)

s_2 = N * D / (N - 1); sigma = math.sqrt(D); s_y = math.sqrt(s_2)
print('Исправленная оценка дисперсии:', s_2, '\nСтат. оценки СКО: sigma =', sigma, 's =', s_y)

As = m3 / pow(s_y, 3); E = m4 / pow(s_y, 4) - 3
print('Стат. оценка коэфф. ассиметрии:', As, '\nСтат. оценка эксцесса:', E)

M0 = interval_range_y[np.argmax(interval_range_y[:, 3])][2]; m_e = C
print('Мода:', M0, 'Медиана:', m_e)

# 3.a
interval_range_file = open('./interval_range.csv', 'r')
lines = interval_range_file.readlines()
interval_range_x = np.ndarray((len(lines), 5))
for line in range(len(lines)):
    elements = lines[line].split(',')
    for i in range(interval_range_x.shape[1]): interval_range_x[line][i] = float(elements[i])
interval_range_file.close()
x_sample_mean = 17.979884799999997
s_x = 8.079444486616898

corr_table = np.zeros((interval_range_y.shape[0], interval_range_x.shape[0]))
for i in range(interval_range_y.shape[0]):
    for j in range(interval_range_x.shape[0]):
        for sample_element in range(N):
            if (interval_range_x[j][0] <= stat_range[sample_element][0] <= interval_range_x[j][1]) \
                    and (interval_range_y[i][0] <= stat_range[sample_element][1] <= interval_range_y[i][1]):
                corr_table[i][j] += 1
save_current_sample(corr_table, './corr_table.csv')

# 3.b
r_xy = 0.0
for i in range(interval_range_y.shape[0]):
    for j in range(interval_range_x.shape[0]):
        r_xy += (corr_table[i][j] * interval_range_y[i][2] * interval_range_x[j][2])
r_xy -= (N * x_sample_mean * y_sample_mean)
r_xy /= (N * s_x * s_y)
print('r_ =', r_xy, '\n', r_xy - 3 * ((1 - pow(r_xy, 2)) / (math.sqrt(N))), '<= r <=', r_xy + 3 * ((1 + pow(r_xy, 2)) / (math.sqrt(N))))

# 3.c
z_ = 1.1513 * np.log10((1 + r_xy) / (1 - r_xy))
sigma_z = 1 / math.sqrt(N - 3)
print('z_ =', z_, 'sigma_z =', sigma_z)
lambda_gamma = 1.96
z_left = z_ - lambda_gamma * sigma_z; z_right = z_ + lambda_gamma * sigma_z
print('Дов. инт-л для ген. значения: (', z_left, ';', z_right, ')')
r_left = (np.exp(2 * z_left) - 1) / (np.exp(2 * z_left) + 1); r_right = (np.exp(2 * z_right) - 1) / (np.exp(2 * z_right) + 1)
print('Дов. инт-л для ген. значения коэфф. корреляции: (', r_left, ';', r_right, ')')

# 3.d
T_watch = (r_xy * math.sqrt(N - 2)) / math.sqrt(1 - pow(r_xy, 2))
print('Т_набл =', T_watch)

# 3.e
sample_file = open('./var_range.csv', 'r')
lines = sample_file.readlines()
sample = np.ndarray((len(lines), 11))
for line in range(len(lines)):
    elements = lines[line].split(',')
    for i in range(sample.shape[1]): sample[line][i] = float(elements[i])
sample_file.close()

plt.clf()
plt.title('Sample graph')
plt.plot(sample[:, 0], sample[:, 1], 'bo', label='Sample')
plt.xticks(interval_range_x[:, 2], np.around(interval_range_x[:, 2], 2))
plt.yticks(interval_range_y[:, 2], np.around(interval_range_y[:, 2], 2))
plt.xlabel('Middle of x-interval'); plt.ylabel('Middle of y-interval')
plt.legend()
plt.grid()
plt.show()

plt.clf()
plt.title('Sample graph')
plt.plot(sample[:, 0], sample[:, 1], 'bo', label='Sample')
plt.plot((sample[:, 1] - y_sample_mean) * r_xy * s_x / s_y + x_sample_mean, sample[:, 1],
         'g-', label='Root mean square regression line from X to Y')
plt.plot(sample[:, 0], (sample[:, 0] - x_sample_mean) * r_xy * s_y / s_x + y_sample_mean,
         'r-', label='Root mean square regression line from Y to X')
plt.xticks(interval_range_x[:, 2], np.around(interval_range_x[:, 2], 2))
plt.yticks(interval_range_y[:, 2], np.around(interval_range_y[:, 2], 2))
plt.xlabel('Middle of x-interval'); plt.ylabel('Middle of y-interval')
plt.legend()
plt.grid()
plt.show()

# 3.f
k_x = 7; k_y = 7
x_y = np.zeros(k_y); y_x = np.zeros(k_x)
n_y = np.array([2, 1, 23, 28, 26, 15, 5]); n_x = np.array([8, 9, 17, 16, 25, 15, 10])
for i in range(k_y):
    for j in range(k_x):
        x_y[i] += corr_table[i][j] * interval_range_x[j][2] / n_y[i]
for i in range(k_x):
    for j in range(k_y):
        y_x[i] += corr_table[j][i] * interval_range_y[j][2] / n_x[i]
print('x_y =', x_y, '\ny_x =', y_x)

Dxy = np.zeros(k_y); Dyx = np.zeros(k_x)
for i in range(k_y):
    for j in range(k_x):
        Dxy[i] += corr_table[i][j] * pow(interval_range_x[j][2] - x_y[i], 2) / n_y[i]
for i in range(k_x):
    for j in range(k_y):
        Dyx[i] += corr_table[j][i] * pow(interval_range_y[j][2] - y_x[i], 2) / n_x[i]
print('D_xy =', Dxy, '\nD_yx =', Dyx)

D_in_gr_xy = 0.0; D_in_gr_yx = 0.0
for i in range(k_y): D_in_gr_xy += Dxy[i] * n_y[i] / N
for i in range(k_x): D_in_gr_yx += Dyx[i] * n_x[i] / N
print('D_in_gr_xy =', D_in_gr_xy, '\nD_in_gr_yx =', D_in_gr_yx)

D_between_gr_xy = 0.0; D_between_gr_yx = 0.0
for i in range(k_y): D_between_gr_xy += n_y[i] * pow(x_y[i] - x_sample_mean, 2) / N
for i in range(k_x): D_between_gr_yx += n_x[i] * pow(y_x[i] - y_sample_mean, 2) / N
print('D_between_gr_xy =', D_between_gr_xy, '\nD_between_gr_yx =', D_between_gr_yx)

D_general_xy = D_in_gr_xy + D_between_gr_xy; D_general_yx = D_in_gr_yx + D_between_gr_yx
print('D_general_xy =', D_general_xy, '\nD_general_yx =', D_general_yx)

sigma_xy = math.sqrt(D_between_gr_xy); sigma_x = math.sqrt(D_general_xy); eta_xy = sigma_xy / sigma_x
sigma_yx = math.sqrt(D_between_gr_yx); sigma_y = math.sqrt(D_general_yx); eta_yx = sigma_yx / sigma_y
print('eta_xy =', eta_xy, '\neta_yx =', eta_yx)

# 3.g
matrix = np.zeros((3, 4))
for i in range(k_x):
    matrix[0][0] += n_x[i] * pow(interval_range_x[i][2], 4)
    matrix[0][1] += n_x[i] * pow(interval_range_x[i][2], 3)
    matrix[0][2] += n_x[i] * pow(interval_range_x[i][2], 2)
    matrix[1][2] += n_x[i] * interval_range_x[i][2]
    matrix[0][3] += n_x[i] * pow(interval_range_x[i][2], 2) * y_x[i]
    matrix[1][3] += n_x[i] * interval_range_x[i][2] * y_x[i]
    matrix[2][3] += n_x[i] * y_x[i]
matrix[1][0] = matrix[0][1]; matrix[1][1] = matrix[0][2]; matrix[2][0] = matrix[0][2]; matrix[2][1] = matrix[1][2]; matrix[2][2] = N
print(matrix)

determinant = np.linalg.det(matrix[:, 0:3])
a = np.linalg.det(matrix[:, [3, 1, 2]]) / determinant
b = np.linalg.det(matrix[:, [0, 3, 2]]) / determinant
c = np.linalg.det(matrix[:, [0, 1, 3]]) / determinant
print(a, b, c)
paraboloid_curve = a * sample[:, 0] * sample[:, 0] + b * sample[:, 0] + c

matrix = np.zeros((2, 3))
for i in range(k_x):
    matrix[0][0] += n_x[i] * pow(interval_range_x[i][2], 2)
    matrix[0][1] += n_x[i] * interval_range_x[i][2]
    matrix[0][2] += n_x[i] * interval_range_x[i][2] / y_x[i]
    matrix[1][2] += n_x[i] / y_x[i]
matrix[1][0] = matrix[0][1]; matrix[1][1] = N
print(matrix)

b = (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) / (matrix[0][1] * matrix[1][0] - matrix[1][1] * matrix[0][0])
a = (matrix[1][2] - b * matrix[1][1]) / matrix[1][0]
print(a, b)

plt.clf()
plt.title('Sample graph')
plt.plot(sample[:, 0], sample[:, 1], 'bo', label='Sample')
plt.plot(sample[:, 0], paraboloid_curve, 'r-', label='Paraboloid correlation curve')
plt.plot(sample[:, 0], 1 / (a * sample[:, 0] + b), 'g-', label='Custom correlation curve')
plt.xticks(interval_range_x[:, 2], np.around(interval_range_x[:, 2], 2))
plt.yticks(interval_range_y[:, 2], np.around(interval_range_y[:, 2], 2))
plt.xlabel('Middle of x-interval'); plt.ylabel('Middle of y-interval')
plt.legend()
plt.grid()
plt.show()
