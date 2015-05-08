import itertools


class Node(object):
    def __init__(self, node_id):
        self.node_id = node_id
        self.children = []
    def create_labeled_tree_string(self):
        tree_string = str(self.node_id) + ':['
        for child in self.children:
            tree_string += child.create_labeled_tree_string()
        tree_string += ']'
        return tree_string
    def create_unlabled_tree_string(self):
        tree_string = '['
        for child in self.children:
            tree_string += child.create_unlabled_tree_string()
        tree_string += ']'
        return tree_string
    def get_tree_node_ids(self):
        node_ids = [self.node_id]
        for child in self.children:
            node_ids.extend(child.get_tree_node_ids())
        return node_ids
    def fill_gamma_matrix(self, gamma_matrix):
        for subtree_node_id in self.get_tree_node_ids():
            gamma_matrix[self.node_id][subtree_node_id] = 1
        for child in self.children:
            child.fill_gamma_matrix(gamma_matrix)

def _create_subtree(node_string, current_id):
    node = Node(next(current_id))
    cnt = 0
    start_idx = 0
    for idx in xrange(len(node_string)):
        if node_string[idx] == '(':
            cnt += 1
        else:
            cnt -= 1
        if cnt == 0:
            node.children.append(_create_subtree(node_string[start_idx+1:idx], current_id))
            start_idx = idx + 1
    return node

def create_subtree(tree_string):
    return _create_subtree(tree_string, itertools.count(-1)).children[0]

def _create_from_parent_array(parent_array, current_id):
    node = Node(next(current_id) - 1)
    while len(parent_array) > 0:
        if parent_array[0] == node.node_id + 1:
            child, parent_array = _create_from_parent_array(parent_array[1:], current_id)
            node.children.append(child)
        else:
            break
    return node, parent_array
    
def create_from_parent_array(parent_array):
    return _create_from_parent_array(parent_array, itertools.count(0))[0].children[0]

