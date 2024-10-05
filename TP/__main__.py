from TP.loading import load_directory
from TP.kmers import stream_kmers, make_dict
from time import time

def jaccard(fileA, fileB, k):   # as seen in lecture
    """ Computes the Jaccard similarity between two list of sequences """
    dictA = make_dict(fileA, k)
    intersection = 0
    union = sum(dictA.values())
    for kmer in stream_kmers(fileB, k):
        if kmer in dictA and dictA[kmer] > 0:
            intersection += 1
            dictA[kmer] -= 1
        else:
            union += 1
    return intersection / union

def dict_intersection(dictA:dict, dictB:dict):
    """ Computes the intersection of two dictionaries """
    intersection = 0
    for kmer, countB in dictB.items():
        countA = dictA.get(kmer, 0)
        intersection += min(countA, countB)
    return intersection


if __name__ == "__main__":
    # Load all the files in a dictionary
    files = load_directory("data")
    filenames = list(files.keys())
    k = 21
    
    print("Computing Jaccard similarity for all pairs of samples")
    # start = time()
    # print("    naive version")
    # for i in range(len(files)):
    #     for j in range(i+1, len(files)):
    #         dist_j = jaccard(files[filenames[i]], files[filenames[j]], k)
    #         print(filenames[i], filenames[j], dist_j)
    # print("Time:", time()-start)
    # print()

    start = time()
    print("    precomputing the dictionnaries")
    dicts = [make_dict(files[filename], k) for filename in filenames]
    dicts_size = [sum(d.values()) for d in dicts]
    print("time:", time()-start)
    print("    computing the pairwise similarities")
    for i in range(len(files)):
        for j in range(i+1, len(files)):
            intersection = dict_intersection(dicts[i], dicts[j])
            dist_j = intersection / (dicts_size[i] + dicts_size[j] - intersection)
            print(filenames[i], filenames[j], dist_j)
    print("total time:", time()-start)
