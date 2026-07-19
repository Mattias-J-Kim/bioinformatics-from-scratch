def count_kmers(sequence, k):
    a = {}
    for i in range(0, len(sequence)-k+1):
        b = sequence[i:i+k]
        if b not in a:
            a[b] = 1
        else:
            a[b] += 1
    return a


def most_frequent_kmers(sequence, k):
    a = []
    b = count_kmers(sequence, k)
    c = max(b.values())
    for i in b:
        if b[i] == c:
            a.append(i)
    return sorted(a)


def find_clumps(sequence, k, L, t):              # Before Optimization
    climps_list = []
    for i in range(0, len(sequence)-L+1):
        a = sequence[i:i+L]
        b = count_kmers(a, k)
        for j in b:
            if b[j] >= t:
                if j not in climps_list:
                    climps_list.append(j)
    return sorted(climps_list)


def find_clumps_optimized(sequence, k, L, t):
    climps_list = []
    counts = count_kmers(sequence[0:L], k)

    for kmer, c in counts.items():
        if c >= t:
            climps_list.append(kmer)

    for i in range(1, len(sequence) - L + 1):
        outgoing = sequence[i-1 : i-1+k]
        incoming = sequence[i+L-k : i+L]

        counts[outgoing] -= 1
        counts[incoming] = counts.get(incoming, 0) + 1

        if counts[incoming] >= t and incoming not in climps_list:
            climps_list.append(incoming)

    return sorted(climps_list)


seq = "CGGACTCGACAGATGTGAAGAACGAAAATGTGACGACGGACGACGGCATGAACGAAAACTGACAGAAGTGACAAAGTGA"
print(find_clumps(seq, 5, 50, 4))
print(find_clumps(seq, 5, 50, 2))
