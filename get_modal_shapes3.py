import sys

import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

from pprint import pprint

import os

def get_labels(X, n_clusters_min, n_clusters_max):
    scores = []
    labels = []
    preds = []

    n_clusters_range = range(n_clusters_min, n_clusters_max+1)
    # find first minimum in the score
    # scores are from -1 to 1
    previous_score = None
    for n in n_clusters_range:
        print("Test {} clusters".format(n))

        model = KMeans(n_clusters=n).fit(X)
        label = model.labels_
        labels.append(label)

        # TODO: test other metric as well
        score = silhouette_score(X, label, metric='euclidean')
        # score = silhouette_score(X, label, metric='cityblock')
        # if previous_score is not None and previous_score > score:
        #     break

        print("Score: {}".format(score))
        scores.append(score)
        previous_score = score

    idx = scores.index(max(scores))
    # TODO: uncomment this
    #idx = 6
    n_clusters = n_clusters_range[idx]
    print("Max score gave {} clusters with {} score".format(n_clusters, scores[idx]))

    labels = labels[idx]

    return labels

def sort_by_labels(amps, freqs, damps, channels, labels):
    sorted_by_label = {}

    for i in range(len(labels)):
        # print("{}:\t{}".format(freqs[i], labels[i]))
        l = labels[i]
        if l not in sorted_by_label.keys():
            sorted_by_label[l] = [(amps[i], freqs[i], damps[i], channels[i])]
        else:
            sorted_by_label[l].append((amps[i], freqs[i], damps[i], channels[i]))

    return sorted_by_label

# calculate weighted averages
def calc_weighted_avg(grouped_by_label):
    avg_freqs = []
    avg_damps = []
    for k in grouped_by_label.keys():
        f_avg = 0.0
        d_avg = 0.0
        a_sum = 0.0
        for x in grouped_by_label[k]:
            # weigthed by amplitude or by energy of peak
            #f_avg += x[0]*x[0]*x[1]
            #d_avg += x[0]*x[0]*x[2]
            #a_sum += x[0]*x[0]
            f_avg += x[0]*x[1]
            d_avg += x[0]*x[2]
            a_sum += x[0]
        f_avg = f_avg/a_sum if a_sum > 0.0 else f_avg
        d_avg = d_avg/a_sum if a_sum > 0.0 else d_avg
        avg_freqs.append(f_avg)
        avg_damps.append(d_avg)
    return (avg_freqs, avg_damps)

def get_sorted_indices(l):
    return [i[0] for i in sorted(enumerate(l), key=lambda x:x[1])]

def main():

    dir_name = sys.argv[1]
#    files = (file for file in os.listdir(dir_name) if (not ("_0_" in file and "_1_" in file) and file.endswith(sys.argv[2])))
    files = (file for file in os.listdir(dir_name) if file.endswith(sys.argv[2]))

    # print(files)

    modes = []

    for file_name in files:
        file_name = os.path.join(dir_name, file_name)
        print(file_name)
        # file_name = file_name.replace('//', "\\")
        # file_name = file_name.replace('/', "\\")
        # print(file_name.split("\\")[-1])
        channel = int(file_name.split("\\")[-1].split("_")[-2])
        # print("Channel: %d" % channel)
        num_modes = 0
        num_dims = 0
        with open(file_name, 'r') as file:
            num_modes, num_dims = next(file).split()
            temp_modes = []
            for i in range(int(num_modes)):
                amp, freq, damp, phase, l_idx, r_idx = next(file).split()
                # TODO: mind that samples only from 100-600 range are taken into account
#                if float(freq) < 600.0:
                temp_modes.append((amp, freq, damp, channel))
            # if all(float(m[1]) < 850.0 for m in temp_modes) and all(float(m[0]) < 98.0 for m in temp_modes) and all(float(m[2]) < 0.2 for m in temp_modes):
            if (
                all(float(m[1]) < 850.0 for m in temp_modes) and
                all(float(m[0]) < 198.0 for m in temp_modes) and
                all(float(m[2]) < 0.3 for m in temp_modes)
            ):
                if temp_modes:
                    for mode in temp_modes:
                        modes.append(mode)

    amps = [float(mode[0]) for mode in modes]
    max_a = max(amps)
    freqs = [float(mode[1]) for mode in modes]
    damps = [float(mode[2]) for mode in modes]
    channels = [int(mode[3]) for mode in modes]

    n_samples = len(amps)
    print("Number of samples: %s" % n_samples)

#    plt.figure(figsize=(max_a, max_a))

    X = np.array(freqs)
    X = X.reshape(-1, 1)

#    labels = get_labels(X, 2, int(n_samples/2))
    # labels = get_labels(X, 2, 20)
    labels = get_labels(X, 2, 10)

    sorted_by_label = sort_by_labels(amps, freqs, damps, channels, labels)

    # pprint(sorted_by_label)
    for k in sorted_by_label.keys():
        print("For mode no {}, there is {} samples.".format(k, len(sorted_by_label[k])))

#    plt.scatter(X, amps, c=labels)

#    plt.title("")

    avg_freqs, avg_damps = calc_weighted_avg(sorted_by_label)

    sorted_idx = get_sorted_indices(avg_freqs)


    print("mode\tavg freq\t\tavg damping\t\tamp per channel")
    
    # for idx in sorted_idx:
    #     print("{}\t{}\t{}\t{}".format(
    #         idx,
    #         avg_freqs[idx],
    #         avg_damps[idx],
    #         [{m[3]: m[0]} for m in sorted_by_label[idx]]
    #     ))


    pprint(
        [
            {
                'mode': idx,
                # 'amp per channel': [{m[3]: m[0]} for m in sorted_by_label[idx]],
                'freq': avg_freqs[idx],
                'damp': avg_damps[idx]
            } for idx in sorted_idx
        ]
    )
#    plt.show()
    # return

if __name__ == "__main__":
    main()

