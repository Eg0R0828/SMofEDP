import matplotlib.pyplot as plt
import numpy as np
import random, math


def distance(a, b, m=2):
    s = 0
    for k_m in range(m): s += pow(a[k_m] - b[k_m], 2)
    return math.sqrt(s)


def change_centroids(groups, centroids, k):
    for j in range(k):
        s = [0, 0]; count = 0
        for n in range(N):
            if groups[j][n] == -1: break
            s[0] += sample[groups[j][n]][0]; s[1] += sample[groups[j][n]][1]; count += 1
        if count != 0: centroids[j] = [s[0] / count, s[1] / count]


def change_centroid(group_array):
    s = [0, 0]; count = 0
    for n in group_array:
        if n == -1: break
        s[0] += sample[n][0]; s[1] += sample[n][1]
        count += 1
    return s[0] / count, s[1] / count


def groups_equal(arr1, arr2):
    for i in range(arr1.shape[0]):
        for j in range(arr1.shape[1]):
            if arr1[i][j] != arr2[i][j]: return False
    return True


def k_means(k, freq_centroids_changing):
    print('k =', k, 'freq_centroids_changing =', freq_centroids_changing)
    centroids = np.ndarray((k, 2))
    groups = np.full((k, N), -1)
    old_groups = np.full((k, N), -1)
    first_stage = True
    for i in range(k):
        index = random.randint(0, N - 1)
        while index in groups: index = random.randint(0, N - 1)
        groups[i][0] = index
        centroids[i] = sample[index]

    while first_stage or not groups_equal(old_groups, groups):
        if first_stage: first_stage = False
        else:
            old_groups = groups.copy()
            groups = np.full((k, N), -1)
        for i in range(N):
            if i not in groups:
                min_distance = None; group_number = -1
                for j in range(k):
                    if group_number == -1 or distance(sample[i], centroids[j]) < min_distance:
                        min_distance = distance(sample[i], centroids[j])
                        group_number = j
                for j in range(N):
                    if groups[group_number][j] == -1:
                        groups[group_number][j] = i
                        break
                if freq_centroids_changing: change_centroids(groups, centroids, k)
        if not freq_centroids_changing: change_centroids(groups, centroids, k)

        F1 = 0; F2 = 0; F3 = 0
        for i in range(k):
            cluster_n = 0
            for j in range(N):
                if groups[i][j] == -1: break
                cluster_n += sample_freq[groups[i][j]][0]
                F1 += pow(distance(sample[groups[i][j]], centroids[i]), 2)
                F3 += pow(np.sum((sample[groups[i][j]] - centroids[i]) * (sample[groups[i][j]] - centroids[i]) *
                                 sample_freq[groups[i][j]][0]) / cluster_n, 2)
                for jj in range(N):
                    if j == jj: continue
                    if groups[i][jj] == -1: break
                    F2 += pow(distance(sample[groups[i][j]], sample[groups[i][jj]]), 2)
        print('F1 =', F1, 'F2 =', F2, 'F3 =', F3)

    styles = ['b', 'g', 'y', 'r', 'm', 'k', 'c']
    plt.clf()
    plt.title('Grouped sample graphic')
    for i in range(k):
        plt.plot([sample[j][0] for j in groups[i][:] if j != -1], [sample[j][1] for j in groups[i][:] if j != -1],
                 styles[i] + '+', label=str(i + 1) + ' group')
        plt.plot(centroids[i][0], centroids[i][1], styles[i] + 'o', label=str(i + 1) + ' group centroid')
    plt.xticks((interval_range_x[:, 2] - x_sample_mean) / s_x, np.around((interval_range_x[:, 2] - x_sample_mean) / s_x, 2))
    plt.yticks((interval_range_y[:, 2] - y_sample_mean) / s_y, np.around((interval_range_y[:, 2] - y_sample_mean) / s_y, 2))
    plt.xlabel('Middle of x-interval (normalized)'); plt.ylabel('Middle of y-interval (normalized)')
    plt.legend()
    plt.grid()
    plt.show()


def concentration_searching(radius):
    groups = np.full((N, N), -1)
    centroids = np.ndarray((N, 2))
    k = 0

    while True:
        start_center = random.randint(0, N - 1); local_radius = radius
        while start_center in groups: start_center = random.randint(0, N - 1)

        while True:
            centroids[k][0] = sample[start_center][0]; centroids[k][1] = sample[start_center][1]
            first_stage = True
            old_group = groups[k].copy()
            while first_stage or not groups_equal(old_group, groups[k]):
                if not first_stage:
                    old_group = groups[k].copy()
                    groups[k] = np.full(N, -1)
                else: first_stage = False
                for i in range(N):
                    if distance(sample[i], centroids[k]) <= local_radius:
                        for j in range(N):
                            if groups[k][j] == -1:
                                groups[k][j] = i
                                break
                centroids[k][0], centroids[k][1] = change_centroid(groups[k])

            crossing = False
            for i in range(N):
                if groups[k][i] != -1 and ((k > 1 and groups[k][i] in groups[0:k - 1, :]) or (k == 1 and groups[k][i] in groups[0])):
                    if local_radius > delta: local_radius -= delta
                    groups[k] = np.full(N, -1)
                    crossing = True
                    break
            if not crossing: break

        k += 1
        complete = True
        for i in range(N):
            if i not in groups: complete = False
        if complete: break

    F1 = 0; F2 = 0; F3 = 0
    for i in range(k):
        cluster_n = 0
        for j in range(N):
            if groups[i][j] == -1: break
            cluster_n += sample_freq[groups[i][j]][0]
            F1 += pow(distance(sample[groups[i][j]], centroids[i]), 2)
            F3 += pow(np.sum((sample[groups[i][j]] - centroids[i]) * (sample[groups[i][j]] - centroids[i]) *
                             sample_freq[groups[i][j]][0]) / cluster_n, 2)
            for jj in range(N):
                if j == jj: continue
                if groups[i][jj] == -1: break
                F2 += pow(distance(sample[groups[i][j]], sample[groups[i][jj]]), 2)
    print('F1 =', F1, 'F2 =', F2, 'F3 =', F3)

    styles = ['b', 'g', 'y', 'r', 'm', 'k', 'c']
    plt.clf()
    plt.title('Grouped sample graphic')
    for i in range(k):
        plt.plot([sample[j][0] for j in groups[i][:] if j != -1], [sample[j][1] for j in groups[i][:] if j != -1],
                 styles[i % len(styles)] + '+', label=str(i + 1) + ' group')
        plt.plot(centroids[i][0], centroids[i][1], styles[i % len(styles)] + 'o', label=str(i + 1) + ' group centroid')
    plt.xticks((interval_range_x[:, 2] - x_sample_mean) / s_x,
               np.around((interval_range_x[:, 2] - x_sample_mean) / s_x, 2))
    plt.yticks((interval_range_y[:, 2] - y_sample_mean) / s_y,
               np.around((interval_range_y[:, 2] - y_sample_mean) / s_y, 2))
    plt.xlabel('Middle of x-interval (normalized)'); plt.ylabel('Middle of y-interval (normalized)')
    plt.legend()
    plt.grid()
    plt.show()


# 1
N = 100
k_ = int(math.sqrt(N / 2))
x_sample_mean = 17.979884799999997
s_x = 8.079444486616898
y_sample_mean = 1013.9642857142858
s_y = 7.355282831323163

sample_file = open('./var_range.csv', 'r')
lines = sample_file.readlines()
sample = np.ndarray((len(lines), 2))
for line in range(len(lines)):
    elements = lines[line].split(',')
    for i in range(2): sample[line][i] = float(elements[i])
sample_file.close()

sample_file = open('./var_range_freq.csv', 'r')
lines = sample_file.readlines()
sample_freq = np.ndarray((len(lines), 2))
for line in range(len(lines)):
    elements = lines[line].split(',')
    for i in range(2): sample_freq[line][i] = float(elements[i])
sample_file.close()

sample_file = open('./interval_range.csv', 'r')
lines = sample_file.readlines()
interval_range_x = np.ndarray((len(lines), 5))
for line in range(len(lines)):
    elements = lines[line].split(',')
    for i in range(interval_range_x.shape[1]): interval_range_x[line][i] = float(elements[i])
sample_file.close()

sample_file = open('./interval_range_2.csv', 'r')
lines = sample_file.readlines()
interval_range_y = np.ndarray((len(lines), 5))
for line in range(len(lines)):
    elements = lines[line].split(',')
    for i in range(interval_range_y.shape[1]): interval_range_y[line][i] = float(elements[i])
sample_file.close()

plt.clf()
plt.title('Sample graphic')
plt.plot(sample[:, 0], sample[:, 1], 'bo', label='Sample')
plt.xticks(interval_range_x[:, 2], np.around(interval_range_x[:, 2], 2))
plt.yticks(interval_range_y[:, 2], np.around(interval_range_y[:, 2], 2))
plt.xlabel('Middle of x-interval'); plt.ylabel('Middle of y-interval')
plt.legend()
plt.grid()
plt.show()

# 2
sample[:, 0] = (sample[:, 0] - x_sample_mean) / s_x; sample[:, 1] = (sample[:, 1] - y_sample_mean) / s_y

plt.clf()
plt.title('Normalized sample graphic')
plt.plot(sample[:, 0], sample[:, 1], 'bo', label='Normalized sample')
plt.xticks((interval_range_x[:, 2] - x_sample_mean) / s_x, np.around((interval_range_x[:, 2] - x_sample_mean) / s_x, 2))
plt.yticks((interval_range_y[:, 2] - y_sample_mean) / s_y, np.around((interval_range_y[:, 2] - y_sample_mean) / s_y, 2))
plt.xlabel('Middle of x-interval (normalized)'); plt.ylabel('Middle of y-interval (normalized)')
plt.legend()
plt.grid()
plt.show()

# 3 - 5
for i in range(2, k_ + 1):
    k_means(i, True); k_means(i, False)

# 6 - 8
D = np.zeros((N, N))
for i in range(N):
    for j in range(N):
        if i != j: D[i][j] = distance(sample[i], sample[j])
delta = 0.01
concentration_searching(float(input('Enter a 1st cluster sphere radius (' + str(np.min(D)) + '; ' + str(np.max(D)) + '): ')))
