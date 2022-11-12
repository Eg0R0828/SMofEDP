import matplotlib.pyplot as plt
import numpy as np
import random, math


def save_current_sample(sample, file_path):
    file = open(file_path, 'w')
    for i in range(sample.shape[0]):
        for j in range(sample.shape[1]):
            file.write(str(sample[i][j]))
            if j != (sample.shape[1] - 1): file.write(',')
        if i != (sample.shape[0] - 1): file.write('\n')
    file.close()


N = 100
data_matrix = np.ndarray((N, 11))

# 1
gen_data_file_path = input('Enter gen. data file path: ')
data_sampling_mode = input('Enter data sampling mode (random / random_unique / mechanical): ')
gen_data_file = open(gen_data_file_path, 'r')
lines = gen_data_file.readlines()
gen_data_matrix = np.ndarray((len(lines) - 1, data_matrix.shape[1]))
for line in range(len(lines)):
    if line == 0: continue
    elements = lines[line].split(',')
    for i in range(data_matrix.shape[1]): gen_data_matrix[line - 1][i] = float(elements[i])
gen_data_file.close()

data_ind_seq = np.ndarray([N])
if data_sampling_mode == 'mechanical':
    for i in range(N): data_ind_seq[i] = int(i * gen_data_matrix.shape[0] / N)
elif data_sampling_mode == 'random':
    for i in range(N): data_ind_seq[i] = random.randint(0, gen_data_matrix.shape[0])
elif data_sampling_mode == 'random_unique':
    for i in range(N):
        data_ind_seq[i] = random.randint(0, gen_data_matrix.shape[0])
        if i > 0:
            while data_ind_seq[i] in data_ind_seq[0:i]:
                data_ind_seq[i] = random.randint(0, gen_data_matrix.shape[0])
for i in range(data_matrix.shape[0]):
    for j in range(data_matrix.shape[1]):
        data_matrix[i][j] = gen_data_matrix[int(data_ind_seq[i])][j]
save_current_sample(data_matrix, './sample.csv')

# 2.a
data_ind_seq.sort()
for i in range(data_matrix.shape[0]):
    for j in range(data_matrix.shape[1]):
        data_matrix[i][j] = gen_data_matrix[int(data_ind_seq[i])][j]
save_current_sample(data_matrix, './stat_range.csv')

data_ind_seq = np.lexsort((data_matrix[:, 10], data_matrix[:, 9], data_matrix[:, 8], data_matrix[:, 7],
                           data_matrix[:, 6], data_matrix[:, 5], data_matrix[:, 4], data_matrix[:, 3],
                           data_matrix[:, 2], data_matrix[:, 1], data_matrix[:, 0]))
ranked_sample = data_matrix.copy()
for i in range(data_matrix.shape[0]):
    for j in range(data_matrix.shape[1]):
        ranked_sample[i][j] = data_matrix[int(data_ind_seq[i])][j]
save_current_sample(ranked_sample, './ranked_range.csv')

# 2.b
var_sample = np.ndarray((0, ranked_sample.shape[1]))
var_sample_freq = np.ndarray((0, 2))
for i in range(ranked_sample.shape[0]):
    if i == 0 or not np.array_equal(ranked_sample[i], ranked_sample[i - 1]):
        var_sample = np.vstack((var_sample, ranked_sample[i]))
        var_sample_freq = np.vstack((var_sample_freq, [1, 1 / N]))
    else:
        var_sample_freq[var_sample_freq.shape[0] - 1][0] += 1
        var_sample_freq[var_sample_freq.shape[0] - 1][1] = var_sample_freq[var_sample_freq.shape[0] - 1][0] / N
save_current_sample(var_sample, './var_range.csv')
save_current_sample(var_sample_freq, './var_range_freq.csv')

k = 1 + 3.31 * np.log10(N)
if int(k) % 2 == 1: k = int(k)
else: k = int(k) + 1
h = (max(var_sample[:, 0]) - min(var_sample[:, 0])) / k
print('k =', k, 'h =', h)

# 2.c
interval_range = np.ndarray((0, 5))
for i in range(1, k + 1):
    interval_range = np.vstack((interval_range, [min(var_sample[:, 0])+(h*(i-1)), min(var_sample[:, 0])+(h*i),
                                        (min(var_sample[:, 0])+(h*(i-1)) + min(var_sample[:, 0])+(h*i)) / 2, 0, 0]))
current_int_ind = 0
for i in range(len(var_sample)):
    if var_sample[i][0] >= interval_range[current_int_ind][1]: current_int_ind += 1
    if current_int_ind >= k: current_int_ind = k - 1
    interval_range[current_int_ind][3] += var_sample_freq[i][0]
    interval_range[current_int_ind][4] += var_sample_freq[i][1]
save_current_sample(interval_range, './interval_range.csv')

# 2.d
plt.clf()
plt.title('Abs. frequency polygon')
plt.plot(interval_range[:, 2], interval_range[:, 3], 'bo-')
plt.xticks(interval_range[:, 2], np.around(interval_range[:, 2], 2))
plt.xlabel('Middle of interval'); plt.ylabel('Abs. frequency')
plt.grid()
plt.show()

plt.clf()
plt.title('Rel. frequency polygon')
plt.plot(interval_range[:, 2], interval_range[:, 4], 'bo-')
plt.xticks(interval_range[:, 2], np.around(interval_range[:, 2], 2))
plt.xlabel('Middle of interval'); plt.ylabel('Rel. frequency')
plt.grid()
plt.show()

# 2.e
plt.clf()
plt.title('Abs. frequency histogram')
plt.bar(interval_range[:, 2], height=interval_range[:, 3] / h, width=h, color='b', tick_label=np.around(interval_range[:, 2], 2))
plt.xlabel('Middle of interval'); plt.ylabel('n / h')
plt.grid()
plt.show()

plt.clf()
plt.title('Rel. frequency histogram')
plt.bar(interval_range[:, 2], height=interval_range[:, 4] / h, width=h, color='b', tick_label=np.around(interval_range[:, 2], 2))
plt.xlabel('Middle of interval'); plt.ylabel('n_rel. / h')
plt.grid()
plt.show()

# 2.f
plt.clf()
plt.title('Empirical distribution function')
points = np.ndarray((0, 2))
for i in range(k):
    accumulated_freq = 0
    for j in range(i): accumulated_freq += (interval_range[j][3] / N)
    points = np.vstack((points, [[interval_range[i][2], accumulated_freq],
                                 [interval_range[i][2], accumulated_freq + (interval_range[i][3] / N)]]))
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
    for j in range(i): accumulated_freq += (interval_range[j][4] / N)
    points = np.vstack((points, [[interval_range[i][2], accumulated_freq],
                                 [interval_range[i][2], accumulated_freq + (interval_range[i][4] / N)]]))
plt.plot(points[:, 0], points[:, 1], 'bo-')
plt.xticks(points[:, 0], np.around(points[:, 0], 2))
plt.xlabel('Middle of interval'); plt.ylabel('Accumulated frequency')
plt.grid()
plt.show()

# 2.g
conditional_options = np.arange(interval_range.shape[0])
C = interval_range[interval_range.shape[0] // 2][2]
N = 100.0
h = 0.0
for interval in range(interval_range.shape[0]):
    h = interval_range[interval][1] - interval_range[interval][0]
    u = int(((interval_range[interval][2] - C) / h + 0.1).round())
    conditional_options[interval] = u
print('Условные варианты:', conditional_options)

M = np.zeros(4)
for k in range(4):
    for j in range(interval_range.shape[0]):
        M[k] += interval_range[j][3] * pow(conditional_options[j], k + 1) / N
print('Условные эмпирические моменты: M1 =', M[0], 'M2 =', M[1], 'M3 =', M[2], 'M4 =', M[3])
m1 = M[0] * h + C
m2 = (M[1] - pow(M[0], 2)) * pow(h, 2)
m3 = (M[2] - 3 * M[1] * M[0] + 2 * pow(M[0], 3)) * pow(h, 3)
m4 = (M[3] - 4 * M[2] * M[0] + 6 * M[1] * pow(M[0], 2) - 3 * pow(M[0], 4)) * pow(h, 4)
print('Центральные эмпирические моменты: m1 =', m1, 'm2 =', m2, 'm3 =', m3, 'm4 =', m4)

x_sample_mean_u = m1; D_u = m2
print('Выборочное среднее, вычисленное с помощью усл. вариант:', x_sample_mean_u)
print('Дисперсия, вычисленная с помощью усл. вариант:', D_u)
x_sample_mean = 0.0
for i in range(interval_range.shape[0]): x_sample_mean += interval_range[i][3] * interval_range[i][2] / N
D = 0.0
for i in range(interval_range.shape[0]): D += interval_range[i][3] * pow(interval_range[i][2] - x_sample_mean, 2) / N
print('Выборочное среднее, вычисленное с помощью стандартной формулы:', x_sample_mean)
print('Дисперсия, вычисленная с помощью стандартной формулы:', D)
s_2 = N * D / (N - 1); sigma = math.sqrt(D); S = math.sqrt(s_2)
print('Исправленная оценка дисперсии:', s_2, '\nСтат. оценки СКО: sigma =', sigma, 's =', S)

As = m3 / pow(S, 3); E = m4 / pow(S, 4) - 3
print('Стат. оценка коэфф. ассиметрии:', As, '\nСтат. оценка эксцесса:', E)

# 2.h
M0 = interval_range[np.argmax(interval_range[:, 3])][2]; m_e = C
print('Мода:', M0, 'Медиана:', m_e)

# 2.i
t_Student = 1.984
print('Дов. инт-л для мат. ожидания: (', x_sample_mean - t_Student * S / math.sqrt(N), ';',
      x_sample_mean + t_Student * S / math.sqrt(N), ')')
q = 0.143
print('Дов. инт-л для оценки СКВО: (', S * (1 - q), ';', S * (1 + q), ')')

# 2.j
Phi_z = [-0.5, -0.4452, -0.3461, -0.17, 0.0557, 0.2642, 0.4032, 0.5]
Hi_2_watch = 0.0
print('x_i | x_i+1 | n_i | z_i | z_i+1 | Ф(z_i) | Ф(z_i+1) | p_i | n`_i')
for i in range(interval_range.shape[0]):
    x_i = interval_range[i][0]; x_next = interval_range[i][1]; n_i = interval_range[i][3]
    if i == 0: z_i = '-inf'
    else: z_i = (x_i - x_sample_mean) / S
    if i == interval_range.shape[0] - 1: z_next = '+inf'
    else: z_next = (x_next - x_sample_mean) / S
    p_i = Phi_z[i + 1] - Phi_z[i]; n_i_ = N * p_i
    print(x_i, ' | ', x_next, ' | ', n_i, ' | ', z_i, ' | ', z_next, ' | ', Phi_z[i], ' | ', Phi_z[i + 1], ' | ', p_i, ' | ', n_i_)
    Hi_2_watch += pow(n_i - n_i_, 2) / n_i_
print('X2набл. =', Hi_2_watch)
k = interval_range.shape[0] - 3
print('X2крит.( alpha = 0.05; k =', k, ') = 9.5')
