"""
De Bruijn Graph Genome Assembler
 
Reconstructs a genome sequence from a set of k-mers by building a De Bruijn
graph and finding an Eulerian path through it.
 
References:
    Compeau, P., Pevzner, P. & Tesler, G. (2011).
        "How to apply de Bruijn graphs to genome assembly."
        Nature Biotechnology 29:987-991.
    Hierholzer, C. (1873). Mathematische Annalen 6:30-32.
 
The general Eulerian path pattern used here is documented separately at:
    github.com/Mattias-J-Kim/algorithms-practice -> eulerian-path/
"""
 
import copy
 
 
def build_debruijn_graph(kmers):
    """
    Build a De Bruijn graph as an adjacency list.
 
    Args:
        kmers: list of str, all of length k
    Returns:
        dict mapping (k-1)-mer -> list of (k-1)-mers
 
    Learning note:
        Each k-mer is an EDGE, not a node. "ATG" becomes an edge from its
        prefix "AT" to its suffix "TG". The nodes are therefore (k-1)-mers,
        and the graph has one edge per input k-mer.
 
        Repeated k-mers must produce repeated entries in the destination
        list. Writing {"TG": ["GC"]} where {"TG": ["GC", "GC"]} is required
        silently deletes an edge, and the assembly will come out short.
    """
    a = {}
    for i in range(0, len(kmers)):
        if kmers[i][:-1] not in a:
            a[kmers[i][:-1]] = [kmers[i][1:]]
        else:
            a[kmers[i][:-1]].append(kmers[i][1:])
    return a
 
 
def find_start_node(graph):
    """
    Find the node an Eulerian path must start from.
 
    Returns:
        the node where out_degree - in_degree == 1
    Raises:
        ValueError if no such node exists or more than one does
 
    Learning note:
        Any node interior to a path must depart once per arrival, so
        out == in. Only the start (+1) and the end (-1) break this symmetry.
        Picking an arbitrary key instead fails on graphs with repeats,
        because traversal from the wrong node strands unused edges.
 
        In-degree cannot be read off the adjacency list directly, since the
        list only stores outgoing edges. It must be counted by scanning
        every destination list for occurrences of the node. Counting
        occurrences rather than testing membership is what makes duplicate
        edges come out right.
    """
    edge_lists = graph.values()
    nodes = graph.keys()
    diffs = {}
    for i in nodes:
        out_degree = len(graph[i])
        in_degree = 0
        for j in edge_lists:
            for l in range(0, len(j)):
                if i == j[l]:
                    in_degree += 1
        diff = out_degree - in_degree
        diffs[i] = diff
 
    diff_values = diffs.values()
    start_candidates_count = 0
    for k in diff_values:
        if k == 1:
            start_candidates_count += 1
 
    if start_candidates_count == 1:
        for t in diffs:
            if diffs[t] == 1:
                return t
    else:
        raise ValueError("No Eulerian path exists")
 
 
def eulerian_path(graph):
    """
    Hierholzer traversal. Return the node visit order.
 
    Learning note:
        The graph is deep-copied because traversal consumes edges by popping
        them. A shallow copy would still share the inner destination lists,
        so the caller's graph would be emptied and verification afterwards
        would be impossible.
 
        stack[-1] is inspected, not stack[0]: the structure must be LIFO so
        that the most recently entered branch is explored to exhaustion
        first. Popping then returns control to exactly the junction it
        descended from, which makes backtracking implicit -- the stack
        itself stores the junction history.
 
        A node is appended to `circuit` only when it has no outgoing edges
        left, so the first node appended lies near the END of the path.
        `circuit` records the order in which nodes became stuck, not the
        order they are visited, hence the final reverse().
 
        graph.get(current, []) is needed because the terminal node appears
        as an edge target but never as a dictionary key.
    """
    local_graph = copy.deepcopy(graph)
    start = find_start_node(local_graph)
    stack = [start]
    circuit = []
 
    while len(stack) > 0:
        current = stack[-1]
        if len(local_graph.get(current, [])) > 0:
            next_node = local_graph[current].pop()
            stack.append(next_node)
        else:
            circuit.append(stack.pop())
    circuit.reverse()
    return circuit
 
 
def reconstruct_genome(path):
    """
    Collapse a node path into the assembled sequence.
 
    Learning note:
        Consecutive nodes overlap by k-2 characters, so only the final
        character of each subsequent node is new. The first node is taken
        whole; every node after it contributes exactly one character.
        Output length is therefore len(path) + k - 2, which equals
        len(kmers) + k - 1.
    """
    route = ""
    route += path[0]
    for i in range(1, len(path)):
        route += path[i][-1]
    return route
 
 
def verify_assembly(genome, k, original_kmers):
    """
    Return True if the assembled genome decomposes back into the input k-mers.
 
    Learning note:
        String equality against the original genome is the WRONG test when
        repeats are present -- see the Test 2 result below. The correct test
        is whether the output is consistent with the evidence that was
        actually supplied, i.e. whether re-decomposing it reproduces the
        input k-mer multiset. sorted() comparison performs the multiset
        check without importing Counter.
    """
    rebuilt_kmers = []
    for i in range(0, len(genome) - k + 1):
        rebuilt_kmers.append(genome[i:i + k])
    return sorted(rebuilt_kmers) == sorted(original_kmers)
 
 
# Test 1: no repeated (k-1)-mers
kmers = ["ATG", "TGC", "GCG"]
graph = build_debruijn_graph(kmers)
path = eulerian_path(graph)
genome = reconstruct_genome(path)
k = len(kmers[0])
print(genome)
print(verify_assembly(genome, k, kmers))
 
# Test 2: repeated (k-1)-mers present
# Assembled output differs from the source string, but is a valid assembly.
# See README for why this is expected rather than a defect.
source = "TAATGCCATGGGATGTT"
kmers_repeat = []
for i in range(0, len(source) - k + 1):
    kmers_repeat.append(source[i:i + k])
 
graph_repeat = build_debruijn_graph(kmers_repeat)
path_repeat = eulerian_path(graph_repeat)
genome_repeat = reconstruct_genome(path_repeat)
print(source)
print(genome_repeat)
print(genome_repeat == source)
print(len(genome_repeat) == len(kmers_repeat) + k - 1)
print(verify_assembly(genome_repeat, k, kmers_repeat))
