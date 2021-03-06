import math
from collections import Counter


class GuideTreeNode:
    def __init__(self, node_id, sequence_number, priority, child_1=None, child_2=None, child_1_distance=None,
                 child_2_distance=None):
        self.id = node_id
        self.sequence_number = sequence_number
        self.priority = priority
        self.child_1 = child_1
        self.child_2 = child_2
        self.child_1_distance = child_1_distance
        self.child_2_distance = child_2_distance


def print_tree(input_tree):
    for key in input_tree.keys():
        print(vars(input_tree.get(key)))


def delete_children(parent_id):
    global node_list
    if parent_id is None:
        return
    else:
        node_list.remove(parent_id)
        delete_children(tree.get(parent_id).child_1)
        delete_children(tree.get(parent_id).child_2)


def find_multiple_sequence_alignment(parent_node):
    parent_node: GuideTreeNode
    if parent_node.child_1 is None:
        return [parent_node.sequence_number]
    else:
        left_sequences = find_multiple_sequence_alignment(tree.get(parent_node.child_1))
        right_sequences = find_multiple_sequence_alignment(tree.get(parent_node.child_2))

        left_consensus_sequence = find_consensus_sequence(left_sequences)
        right_consensus_sequence = find_consensus_sequence(right_sequences)

        align_left, align_right, distance = global_align(left_consensus_sequence, right_consensus_sequence)

        align_based_on_consensus_alignment(left_sequences, left_consensus_sequence, align_left)
        align_based_on_consensus_alignment(right_sequences, right_consensus_sequence, align_right)

        left_sequences.extend(right_sequences)
        return left_sequences


def find_consensus_sequence(sequences_index):
    sequences_index: []
    consensus_sequence = ""
    for i in range(len(sequences[sequences_index[0]])):
        ith_chars = []
        for seq_id in sequences_index:
            ith_chars.append(sequences[seq_id][i])
        consensus_sequence += Counter(ith_chars).most_common(1)[0][0]
    return consensus_sequence


def align_based_on_consensus_alignment(sequences_index, consensus_sequence, aligned_consensus):
    global sequences
    for i in range(len(aligned_consensus)):
        if aligned_consensus[i] == '-':
            if len(consensus_sequence) < i + 1:
                consensus_sequence += '-'
                for seq_id in sequences_index:
                    sequences[seq_id] += '-'
            elif consensus_sequence[i] != '-':
                consensus_sequence = consensus_sequence[0:i] + '-' + consensus_sequence[i:len(consensus_sequence)]
                for seq_id in sequences_index:
                    sequences[seq_id] = sequences[seq_id][0:i] + '-' + sequences[seq_id][i:len(sequences[seq_id])]
    return sequences_index


def calculate_score(final_sequences):
    double_gap_score = -2
    gap_score = -2
    match_score = 1
    mismatch_score = -1
    score = 0
    for i in range(len(final_sequences[0])):
        ith_chars = []
        for seq in final_sequences:
            ith_chars.append(seq[i])
        for first in range(len(ith_chars)):
            for second in range(first + 1, len(ith_chars)):
                if ith_chars[first] == '-' and ith_chars[second] == '-':
                    score += double_gap_score
                elif ith_chars[first] == '-' or ith_chars[second] == '-':
                    score += gap_score
                elif ith_chars[first] == ith_chars[second]:
                    score += match_score
                else:
                    score += mismatch_score
    return score


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


def create_guide_tree(distance_matrix, N):
    global tree
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
        minimum_distance_index1 = 0
        minimum_distance_index2 = 0
        minimum_distance_value = math.inf
        for seq1_id in distance_prime_matrix.keys():
            for seq2_id in distance_prime_matrix.get(seq1_id).keys():
                if distance_prime_matrix[seq1_id][seq2_id] < minimum_distance_value or (distance_prime_matrix[seq1_id][seq2_id] == minimum_distance_value and (min(tree.get(seq1_id).priority, tree.get(seq2_id).priority) < min(tree.get(minimum_distance_index1).priority, tree.get(minimum_distance_index2).priority))):
                    minimum_distance_index1, minimum_distance_index2 = seq1_id, seq2_id
                    minimum_distance_value = distance_prime_matrix[seq1_id][seq2_id]

        # create parent node of 2 minimum distanced nodes
        child_1_distance = distance_matrix[minimum_distance_index1][minimum_distance_index2] / 2 + (
                net_divergence_r[minimum_distance_index1] - net_divergence_r[minimum_distance_index2]) / (
                                   2 * (N - 2))
        child_2_distance = distance_matrix[minimum_distance_index1][minimum_distance_index2] - child_1_distance

        new_node = GuideTreeNode(len(tree), None, min(tree.get(minimum_distance_index1).priority,
                                                      tree.get(minimum_distance_index2).priority),
                                 minimum_distance_index1, minimum_distance_index2, child_1_distance,
                                 child_2_distance)
        tree[new_node.id] = new_node

        # create new distance matrix
        new_distance_matrix = {}
        new_node_column = {}
        for seq1_id in distance_matrix.keys():
            distance_column = {}
            if seq1_id not in [minimum_distance_index1, minimum_distance_index2]:
                for seq2_id in distance_matrix[seq1_id].keys():
                    if seq2_id not in [minimum_distance_index1, minimum_distance_index2]:
                        distance_column[seq2_id] = distance_matrix[seq1_id][seq2_id]
                new_node_distance = (distance_matrix[minimum_distance_index1][seq1_id] +
                                     distance_matrix[minimum_distance_index2][seq1_id] -
                                     distance_matrix[minimum_distance_index1][minimum_distance_index2]) / 2
                distance_column[new_node.id] = new_node_distance
                new_node_column[seq1_id] = new_node_distance
                new_distance_matrix[seq1_id] = distance_column
        new_distance_matrix[new_node.id] = new_node_column
        distance_matrix = new_distance_matrix

        N -= 1


def main():
    global tree, node_list, sequences
    n = int(input())  # getting counts of sequences
    for i in range(n):
        new_sequence = input().upper()  # getting sequences that we are going to align
        sequences.append(new_sequence)
        tree[i] = GuideTreeNode(len(tree), i, i)

    distance_matrix = {}

    for seq1_id in tree.keys():
        distance_column = {}
        for seq2_id in tree.keys():
            if seq1_id != seq2_id:
                align_x, align_y, distance = global_align(sequences[tree.get(seq1_id).sequence_number],
                                                          sequences[tree.get(seq2_id).sequence_number])
                distance_column[seq2_id] = distance
        distance_matrix[seq1_id] = distance_column

    # creating guide tree
    create_guide_tree(distance_matrix, n)

    # delete first child of root and all its descendants to find second child
    node_list = list(tree.keys())
    rootChild1ID = tree.get(len(tree) - 1).id
    delete_children(rootChild1ID)
    # getting second child of root node
    rootChild2ID = node_list[len(node_list) - 1]

    # create root node and add it to the tree
    root = GuideTreeNode(len(tree), None, min(tree.get(rootChild1ID).priority, tree.get(rootChild2ID).priority),
                         rootChild1ID, rootChild2ID)
    tree[root.id] = root

    find_multiple_sequence_alignment(root)

    # printing MSA
    for sequence in sequences:
        print(sequence)
    # printing MSA score
    print(calculate_score(sequences))


tree = {}
node_list = []
sequences = []

if __name__ == '__main__':
    main()
