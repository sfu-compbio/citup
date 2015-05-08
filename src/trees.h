#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <assert.h>


struct Node
{
public:
    Node(int node_id) : node_id(node_id)
    {
    }

    int count_nodes() const
    {
        int cnt = 0;
        
        for (std::vector<Node*>::const_iterator child = children.begin(); child != children.end(); child++)
        {
            cnt += (*child)->count_nodes();
        }

        return cnt + 1;
    }

    std::string create_tree_string() const
    {
        std::string tree_string = "(";
        
        for (std::vector<Node*>::const_iterator child = children.begin(); child != children.end(); child++)
        {
            tree_string += (*child)->create_tree_string();
        }
        
        tree_string += "):";

        std::stringstream ss;
        ss << node_id;
        tree_string += ss.str();

        return tree_string;
    }

    std::vector<int> get_tree_node_ids() const
    {
        std::vector<int> node_ids;

        node_ids.push_back(node_id);

        for (std::vector<Node*>::const_iterator child = children.begin(); child != children.end(); child++)
        {
            std::vector<int> child_node_ids = (*child)->get_tree_node_ids();

            node_ids.insert(node_ids.end(), child_node_ids.begin(), child_node_ids.end());
        }

        return node_ids;
    }

    void fill_gamma_matrix(std::vector<std::vector<int> >& gamma_matrix) const
    {
        std::vector<int> tree_node_ids = get_tree_node_ids();

        for (std::vector<int>::const_iterator node_id_iter = tree_node_ids.begin(); node_id_iter != tree_node_ids.end(); node_id_iter++)
        {
            assert(node_id < gamma_matrix.size());
            assert(*node_id_iter < gamma_matrix.size());

            gamma_matrix[node_id][*node_id_iter] = 1;
        }

        for (std::vector<Node*>::const_iterator child = children.begin(); child != children.end(); child++)
        {
            (*child)->fill_gamma_matrix(gamma_matrix);
        }
    }

    std::vector<std::vector<int> > get_gamma_matrix() const
    {
        std::vector<std::vector<int> > gamma_matrix = get_zero_matrix();

        fill_gamma_matrix(gamma_matrix);

        return gamma_matrix;
    }

    void fill_adjacency_matrix(std::vector<std::vector<int> >& adjacency_matrix) const
    {
        for (std::vector<Node*>::const_iterator child = children.begin(); child != children.end(); child++)
        {
            assert(node_id < adjacency_matrix.size());
            assert((*child)->node_id < adjacency_matrix.size());

            adjacency_matrix[node_id][(*child)->node_id] = 1;
            adjacency_matrix[(*child)->node_id][node_id] = 1;

            (*child)->fill_adjacency_matrix(adjacency_matrix);
        }
    }

    std::vector<std::vector<int> > get_adjacency_matrix() const
    {
        std::vector<std::vector<int> > adjacency_matrix = get_zero_matrix();

        fill_adjacency_matrix(adjacency_matrix);

        return adjacency_matrix;
    }

    void fill_adjacency_list(std::vector<std::vector<int> >& adjacency_list) const
    {
        for (std::vector<Node*>::const_iterator child = children.begin(); child != children.end(); child++)
        {
            adjacency_list.push_back(std::vector<int>());
            adjacency_list.back().push_back(node_id);
            adjacency_list.back().push_back((*child)->node_id);

            (*child)->fill_adjacency_list(adjacency_list);
        }
    }

    std::vector<std::vector<int> > get_adjacency_list() const
    {
        std::vector<std::vector<int> > adjacency_list;

        fill_adjacency_list(adjacency_list);

        return adjacency_list;
    }

    std::vector<std::vector<int> > get_zero_matrix() const
    {
        std::vector<std::vector<int> > zero_matrix;

        int num_nodes = count_nodes();

        zero_matrix.resize(num_nodes);
        for(int i = 0; i < num_nodes; i++)
        {
            zero_matrix[i].resize(num_nodes, 0);
        }

        return zero_matrix;
    }

    int node_id;
    std::vector<Node*> children;
};
    
Node* create_subtree(const std::string& node_string, char open_node_char, int& node_id)
{
    Node* node = new Node(node_id++);
    int cnt = 0;
    int start_idx = 0;
    for (size_t idx = 0; idx < node_string.size(); idx++)
    {
        if (node_string[idx] == open_node_char)
        {
            cnt += 1;
        }
        else
        {
            cnt -= 1;
        }
        if (cnt == 0)
        {
            node->children.push_back(create_subtree(node_string.substr(start_idx+1, idx-start_idx-1), open_node_char, node_id));
            start_idx = idx + 1;
        }
    }
    return node;
}

Node* interpret_tree_string(const std::string& node_string)
{
    char open_node_char = node_string[0];

    // The create_subtree function creates an extra root node for the whole string
    int node_id = -1;
    Node* tree = create_subtree(node_string, open_node_char, node_id);
    assert(tree->children.size() == 1);

    return tree->children[0];
}

