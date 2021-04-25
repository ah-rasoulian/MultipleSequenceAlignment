import math


class Node:
    def __init__(self, node_id, sequence_number, child_1=None, child_2=None, child_1_distance=None,
                 child_2_distance=None):
        self.id = node_id
        self.sequence_number = sequence_number
        self.child_1 = child_1
        self.child_2 = child_2
        self.child_1_distance = child_1_distance
        self.child_2_distance = child_2_distance


def global_align(x, y, s_match=1, s_mismatch=-1, s_gap=-2):
    A = []
    for i in range(len(y) + 1):
        A.append([0] * (len(x) + 1))
    for i in range(len(y) + 1):
        A[i][0] = s_gap * i
    for i in range(len(x) + 1):
        A[0][i] = s_gap * i
    for i in range(1, len(y) + 1):
        for j in range(1, len(x) + 1):
            A[i][j] = max(
                A[i][j - 1] + s_gap,
                A[i - 1][j] + s_gap,
                A[i - 1][j - 1] + (s_match if (y[i - 1] == x[j - 1] and y[i - 1] != '-') else 0) + (
                    s_mismatch if (y[i - 1] != x[j - 1] and y[i - 1] != '-' and x[j - 1] != '-') else 0) + (
                    s_gap if (y[i - 1] == '-' or x[j - 1] == '-') else 0)
            )
    align_X = ""
    align_Y = ""
    i = len(x)
    j = len(y)

    while i > 0 or j > 0:
        current_score = A[j][i]

        if i > 0 and j > 0 and (
                ((x[i - 1] == y[j - 1] and y[j - 1] != '-') and current_score == A[j - 1][i - 1] + s_match) or
                ((y[j - 1] != x[i - 1] and y[j - 1] != '-' and x[i - 1] != '-') and current_score == A[j - 1][
                    i - 1] + s_mismatch) or
                ((y[j - 1] == '-' or x[i - 1] == '-') and current_score == A[j - 1][i - 1] + s_gap)
        ):
            align_X = x[i - 1] + align_X
            align_Y = y[j - 1] + align_Y
            i = i - 1
            j = j - 1
        elif i > 0 and (current_score == A[j][i - 1] + s_gap):
            align_X = x[i - 1] + align_X
            align_Y = "-" + align_Y
            i = i - 1
        else:
            align_X = "-" + align_X
            align_Y = y[j - 1] + align_Y
            j = j - 1
    return align_X, align_Y, A[len(y)][len(x)]


def main():
    global tree
    n = int(input())  # getting counts of sequences
    sequences = []
    for i in range(n):
        new_sequence = input().upper()  # getting sequences that we are going to align
        sequences.append(new_sequence)
        tree[i] = Node(len(tree), i)

    distance_matrix = {}
    alignment_matrix = {}

    for seq1 in tree:
        distance_column = {}
        align_column = {}
        for seq2 in tree:
            if seq1.id != seq2.id:
                align_x, align_y, distance = global_align(sequences[seq1.sequence_number], sequences[seq2.
                                                          sequence_number])
                distance_column[seq2.id] = distance
                align_column[seq2.id] = [align_x, align_y]
        distance_matrix[seq1.id] = distance_column
        alignment_matrix[seq1.id] = align_column

    # creating guide tree
    N = n
    while N > 2:
        net_divergence_r = {}
        # calculating r values
        for seq1_id in distance_matrix.keys():
            _sum = 0
            for seq2_id in distance_matrix.get(seq1_id).keys():
                _sum += distance_matrix[seq1_id][seq2_id]
            net_divergence_r[seq1_id] = _sum

        # calculating new distance matrix
        distance_prime_matrix = {}
        for seq1_id in distance_matrix.keys():
            distance_column = {}
            for seq2_id in distance_matrix.get(seq1_id).keys():
                distance_column[seq2_id] = distance_matrix[seq1_id][seq2_id] - (
                            net_divergence_r.get(seq1_id) + net_divergence_r.get(seq2_id)) / (N - 2)
            distance_prime_matrix[seq1_id] = distance_column

        # find minimum distance index
        minimum_distance_index1 = -1
        minimum_distance_index2 = -1
        minimum_distance_value = math.inf
        for seq1_id in distance_prime_matrix.keys():
            for seq2_id in distance_prime_matrix.get(seq1_id).keys():
                if distance_prime_matrix[seq1_id][seq2_id] < minimum_distance_value:
                    minimum_distance_index1, minimum_distance_index2 = seq1_id, seq2_id
                    minimum_distance_value = distance_prime_matrix[seq1_id][seq2_id]

        # create parent node of 2 minimum distanced nodes
        child_1_distance = distance_matrix[minimum_distance_index1][minimum_distance_index2] / 2 + (
                    net_divergence_r[minimum_distance_index1] - net_divergence_r[minimum_distance_index2]) / (
                                       2 * (N - 2))
        child_2_distance = distance_matrix[minimum_distance_index1][minimum_distance_index2] - child_1_distance
        new_node = Node(len(tree), None, minimum_distance_index1, minimum_distance_index2, child_1_distance,
                        child_2_distance)
        tree[new_node.id] = new_node

        # create new distance matrix
        new_distance_matrix = {}
        for seq1_id in distance_matrix.keys():
            distance_column = {}
            if seq1_id not in [minimum_distance_index1, minimum_distance_index2]:
                for seq2_id in distance_matrix[seq1_id].keys():
                    if seq2_id not in [minimum_distance_index1, minimum_distance_index2]:
                        distance_column[seq2_id] = distance_matrix[seq1_id][seq2_id]
                new_node_distance = (distance_matrix[minimum_distance_index1][seq1_id] + distance_matrix[minimum_distance_index2][seq1_id] - distance_matrix[minimum_distance_index1][minimum_distance_index2]) / 2
                distance_column[new_node.id] = new_node_distance
                new_distance_matrix[seq1_id] = distance_column

        distance_matrix = new_distance_matrix
        N -= 1


tree = {}

if __name__ == '__main__':
    main()
