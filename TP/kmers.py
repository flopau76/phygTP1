def encode_nucl_old(letter):
    """ Encodes a nucleotide on two bits """
    encode = {'A': 0b00, 'C': 0b01, 'T': 0b10, 'G': 0b11}
    return encode[letter]

def encode_nucl_rev_old(letter):
    """ Encodes the complementary of a nucleotide on two bits """
    encode = {'T': 0b00, 'G': 0b01, 'A': 0b10, 'C': 0b11}
    return encode[letter]

def encode_nucl(letter):
    """ Encodes a nucleotide on two bits
    using the ascii code"""
    return (ord(letter) >> 1) & 0b11

def encode_nucl_rev(letter):
    """ Encodes the complementary of a nucleotide on two bits
    using the ascii code"""
    return encode_nucl(letter) ^ 0b10

def encode_kmer(seq, k):
    """ Encodes the first k-mer, and its reverse complementary """
    kmer = 0
    rev_kmer = 0
    for letter in seq[0:k]:
        kmer <<= 2
        rev_kmer >>= 2
        kmer += encode_nucl(letter)
        rev_kmer += encode_nucl_rev(letter) << (2*(k-1))
    return kmer, rev_kmer

def stream_kmers(file, k):
    """ Enumerates all k-mers present in file (list of sequences)
    Note: we keep only the lowest value of each kmer and its reverse complementary """
    for seq in file:
        mask = (1 << (2*k)) - 1 # to keep only the k rightmost nucleotides
        kmer, rev_kmer  = encode_kmer(seq, k)
        yield  min(kmer, rev_kmer)
        for i in range(len(seq)-k):
            kmer <<= 2
            kmer &= mask    # bitwise AND
            kmer += encode_nucl(seq[i+k])
            rev_kmer >>= 2
            rev_kmer += encode_nucl_rev(seq[i+k]) << (2*(k-1))
            yield min(kmer, rev_kmer)

def make_dict(file, k):
    """ Creates a dictionary representing the multi-set of k-mers present in file
    Note: file is a list of sequences"""
    dict = {}
    for kmer in stream_kmers(file, k):
        if kmer in dict:
            dict[kmer] += 1
        else:
            dict[kmer] = 1
    return dict


def kmer2str(val, k):
    """ Transform a kmer integer into a its string representation
    :param int val: An integer representation of a kmer
    :param int k: The number of nucleotides involved into the kmer.
    :return str: The kmer string formatted
    """
    letters = ['A', 'C', 'T', 'G']
    str_val = []
    for _ in range(k):
        str_val.append(letters[val & 0b11])
        val >>= 2

    str_val.reverse()
    return "".join(str_val)