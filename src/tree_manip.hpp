#pragma once

#include <cassert>
#include <memory>
#include <stack>
#include <boost/format.hpp>
#include <queue>
#include <set>
#include <regex>
#include <boost/range/adaptor/reversed.hpp>
#include <utility>
#include "tree.hpp"
#include "xop.hpp"
#include "lot.hpp"

using namespace std;
using namespace boost;

namespace op {
    class TreeManip {
    public:
        TreeManip();
        explicit TreeManip(const Tree::SharedPtr& t);
        ~TreeManip() = default;

        void                        createLeafNodeMap(map<string, Node *> & leafmap) const;
        Node *                      getNodeWithSplit(const Split & s) const;

        void                        setTree(const Tree::SharedPtr & t);
        Tree::SharedPtr             getTree();
        void                        setIsRooted(bool is_rooted);
        bool                        getIsRooted() const;
        double                      calcTreeLength() const;
        void                        scaleAllEdgeLengths(double scaler) const;
        void                        createTestTree();
        unsigned                    randomStep(Lot & lot, double mu, double sigma);

        string                      debugMakeNewick(unsigned precision, bool use_names) const;
        string                      makeNewick(unsigned precision, bool use_names) const;
        void                        buildFromNewick(const string &newick, bool rooted, bool allow_polytomies);
        void                        setLeafNames(const vector<string> & leafnames) const;
        void                        storeSplits(set<Split> & internal_splits, set<Split> & leaf_splits) const;
        void                        setEdgeLength(const Split & s, double edge_length) const;
        void                        rerootAt(int node_number) const;
        void                        nniNodeSwap(Node * a, Node * b);
        void                        dropSplit(const Split & s) const;
        void                        addSplit(const Split & s) const;
        void                        debugCheckSplits() const;
        void                        debugCheckForSinNombreLeaves(const string & msg) const;
        unsigned                    debugCountUnusedNodes() const;

        static string               nexusTranslateCommand();    // creates translate command from _taxon_names
        static string               renumberNewick(const string & newick, map<unsigned, unsigned> taxon_index_map);
        static void                 stripComments(string & newick);

        static vector<string>           _taxon_names; // all trees stored by any TreeManip use this taxon ordering
        static map<string, unsigned>    _taxon_map;   // maps taxon name found in the treefile to the index into _taxon_names

    private:

        void                        refreshPreorder() const;
        void                        refreshLevelorder() const;

        static void                 rerootHelper(Node * m, Node * t);

        static void                 extractNodeNumberFromName(Node * nd, set<unsigned> & used, unsigned nlvs);
        static void                 extractEdgeLen(Node * nd, string edge_length_string);
        static unsigned             countNewickLeaves(const string &newick);
        static void                 stripOutNexusComments(string & newick);
        static bool                 canHaveSibling(const Node * nd, bool rooted, bool allow_polytomies);

        Tree::SharedPtr             _tree;
        bool                        _is_rooted;

    public:

        typedef std::shared_ptr< TreeManip > SharedPtr;
    };

    inline TreeManip::TreeManip() {
        //cerr << "Constructing a TreeManip" << endl;
        _is_rooted = true;
        _tree.reset();
    }

    inline TreeManip::TreeManip(const Tree::SharedPtr & t) {
        //cerr << "Constructing a TreeManip with a supplied tree" << endl;
        _is_rooted = true;
        _tree.reset();
        setTree(t);
    }

    inline Node * TreeManip::getNodeWithSplit(const Split & s) const {
        Node * thend = nullptr;
        for (auto nd : _tree->_preorder) {
            if (nd->getSplit() == s) {
                thend = nd;
                break;
            }
        }
        return thend;
    }

    inline void TreeManip::createLeafNodeMap(map<string, Node *> & leafmap) const {
        leafmap.clear();
        for (auto nd : _tree->_preorder) {
            if (!nd->getLeftChild()) {
                // leaf node
                leafmap[nd->getName()] = nd;
            }
        }
    }

    inline void TreeManip::setTree(const Tree::SharedPtr & t) {
        assert(t);
        _tree = t;
    }

    inline Tree::SharedPtr TreeManip::getTree() {
        return _tree;
    }

    inline void TreeManip::setIsRooted(bool is_rooted) {
        _is_rooted = is_rooted;
    }

    inline bool TreeManip::getIsRooted() const {
        return _is_rooted;
    }

    inline double TreeManip::calcTreeLength() const {
        double TL = 0.0;
        for (auto nd : _tree->_preorder) {
            TL += nd->_edge_length;
        }
        return TL;
    }

    inline void TreeManip::scaleAllEdgeLengths(double scaler) const {
        for (auto nd : _tree->_preorder) {
            nd->_edge_length *= scaler;
        }
    }

    inline unsigned TreeManip::randomStep(Lot & lot, double mu, double sigma) {
        // Take a random step by adding a Normal(mu, sigma) variate to each edge
        // Internal edges that drop below zero result in a random change in topology.
        // Leaf edges that drop below zero are reflected back onto the positive real line.
        // Note that sigma is the standard deviation (not the variance) of the normal variate.

        // Choose random normal deviates used to modify edge lengths
        // Note: decide new edge lengths first, adding them to vector v,
        // then use v to make changes to the tree. This is because changes
        // to the tree would disrupt the preorder sequence, making it very
        // difficult to figure out which edge lengths still need to be
        // modified.
        typedef pair<double, Node *> edge_node_pair_t;
        vector<edge_node_pair_t> v;
        for (auto nd : _tree->_preorder) {
            // skip internal node whose parent is NULL (that is the root node)
            if (nd->_parent->_parent == nullptr)
                continue;

            if (nd->_left_child) {
                // nd is an internal node
                double edge_length = nd->_edge_length;
                double normal_variate = mu + sigma*lot.normal();
                double new_edge_length = edge_length + normal_variate;
                v.emplace_back(new_edge_length, nd);
            }
            else {
                // nd is a leaf node
                double edge_length = nd->_edge_length;
                double normal_variate = mu + sigma*lot.normal();
                double new_edge_length = edge_length + normal_variate;
                if (new_edge_length < Node::_smallest_edge_length) {
                    // Case 1: v = new_edge_length > 0
                    // ----|----------|--------------|--------------|------
                    //     0          v           smallest    new_edge_length
                    //                |<-----x------>|<-----x------>|
                    // Case 2: v = new_edge_length < 0
                    // ----|----------|--------------|-------------------------|
                    //     v          0           smallest              new_edge_length
                    //     |<-----------x----------->|<-----------x----------->|
                    double v = new_edge_length;
                    double smallest = Node::_smallest_edge_length;
                    double x = smallest - v;
                    new_edge_length = smallest + x;
                    nd->_edge_length = new_edge_length;
                }
            }
        }

        // Update tree based on new edge lengths in v
        unsigned cumulative_nni = 0;
        for (auto & p : v) {
            // Let x equal p.second. Because only non-root internal nodes were included in v, we know that
            // x has two children (a is the left child of x and b is the right child of x), and because
            // the sole descendant of the root is excluded, we also know that x has parent y. That is,
            // one of the three configurations below pertains if the tree is binary:
            //
            //  a     b               a     b
            //   \   /                 \   /
            //    \ /                   \ /
            //     x     d         d     x
            //      \   /           \   /
            //       \ /             \ /
            //        y               y
            //        |               |
            //        |               |
            //        c               c
            //
            bool same_orthant = (p.first >= 0.0);
            if (same_orthant) {
                p.second->_edge_length = p.first;
                if (p.second->_edge_length < Node::_smallest_edge_length)
                    p.second->_edge_length = Node::_smallest_edge_length;
            }
            else {
                cumulative_nni++;

                // Flip coin to decide whether a or b will be involved in the move
                double u = lot.uniform();
                bool a_travels_with_x = (u < 0.5);
                bool x_is_left_child_of_y = (p.second->_parent->_left_child == p.second);

                if (a_travels_with_x && x_is_left_child_of_y) {
                    // |-starting-|     |--ending--|
                    //
                    //  a     b          a     d
                    //  \\   /           \\   /
                    //   \\ /             \\ /
                    //     x     d          x     b
                    //      \   /            \   /
                    //       \ /    -->       \ /
                    //        y                y
                    //        |                |
                    //        |                |
                    //        c                c
                    Node * x = p.second;
                    Node * b = x->_left_child->_right_sib;
                    Node * d = x->_right_sib;
                    nniNodeSwap(b, d);
                    x->_edge_length = -p.first;
                    if (x->_edge_length < Node::_smallest_edge_length)
                        x->_edge_length = Node::_smallest_edge_length;
                }
                else if (a_travels_with_x && !x_is_left_child_of_y) {
                    // |-starting-|   |--ending--|
                    //
                    //     a     b       a     d
                    //     \\   /        \\   /
                    //      \\ /          \\ /
                    //  d     x       b     x
                    //   \   /         \   /
                    //    \ /    -->    \ /
                    //     y             y
                    //     |             |
                    //     |             |
                    //     c             c
                    Node * x = p.second;
                    Node * y = x->_parent;
                    Node * b = x->_left_child->_right_sib;
                    Node * d = y->_left_child;

                    // //temporary!
                    // cerr << debugMakeNewick(5, false) << endl;

                    nniNodeSwap(b, d);
                    x->_edge_length = -p.first;
                    if (x->_edge_length < Node::_smallest_edge_length)
                        x->_edge_length = Node::_smallest_edge_length;
                }
                else if (!a_travels_with_x && x_is_left_child_of_y) {
                    // |-starting-|     |--ending--|
                    //
                    //  a     b          d     b
                    //   \   //           \   //
                    //    \ //             \ //
                    //     x     d          x     a
                    //      \   /            \   /
                    //       \ /    -->       \ /
                    //        y                y
                    //        |                |
                    //        |                |
                    //        c                c
                    Node * x = p.second;
                    Node * a = x->_left_child;
                    Node * d = x->_right_sib;
                    nniNodeSwap(a, d);
                    x->_edge_length = -p.first;
                    if (x->_edge_length < Node::_smallest_edge_length)
                        x->_edge_length = Node::_smallest_edge_length;
                }
                else {
                    // |-starting-|   |--ending--|
                    //
                    //     a     b       d     b
                    //      \   //        \   //
                    //       \ //          \ //
                    //  d     x       a     x
                    //   \   /         \   /
                    //    \ /    -->    \ /
                    //     y             y
                    //     |             |
                    //     |             |
                    //     c             c
                    Node * x = p.second;
                    Node * y = x->_parent;
                    Node * a = x->_left_child;
                    Node * d = y->_left_child;

                    // //temporary!
                    // Node * b = a->_right_sib;
                    // Node * c = y->_parent;
                    // x->describeNode("x");
                    // y->describeNode("y");
                    // a->describeNode("a");
                    // b->describeNode("b");
                    // c->describeNode("c");
                    // d->describeNode("d");
                    // cerr << makeNewick(5,false) << endl;

                    nniNodeSwap(a, d);
                    x->_edge_length = -p.first;
                    if (x->_edge_length < Node::_smallest_edge_length)
                        x->_edge_length = Node::_smallest_edge_length;
                }
            }
        }
        return cumulative_nni;
    }

    inline void TreeManip::createTestTree() {
        _tree.reset();
        _tree = std::make_shared<Tree>();
        _tree->_nodes.resize(6);

        Node * root_node       = &_tree->_nodes[0];
        Node * first_internal  = &_tree->_nodes[1];
        Node * second_internal = &_tree->_nodes[2];
        Node * first_leaf      = &_tree->_nodes[3];
        Node * second_leaf     = &_tree->_nodes[4];
        Node * third_leaf      = &_tree->_nodes[5];

        // Here is the structure of the tree (numbers in
        // parentheses are node numbers, other numbers
        // are edge lengths):
        //
        // first_leaf (0)   second_leaf (1)   third_leaf (2)
        //      \              /                  /
        //       \ 0.1        / 0.1              /
        //        \          /                  /
        //     second_internal (3)             / 0.2
        //             \                      /
        //              \ 0.1                /
        //               \                  /
        //                first_internal (4)
        //                        |
        //                        | 0.1
        //                        |
        //                    root_node (5)
        //
        root_node->_parent = nullptr;
        root_node->_left_child = first_internal;
        root_node->_right_sib = nullptr;
        root_node->_number = 5;
        root_node->_name = "root node";
        root_node->_edge_length = 0.0;

        first_internal->_parent = root_node;
        first_internal->_left_child = second_internal;
        first_internal->_right_sib = nullptr;
        first_internal->_number = 4;
        first_internal->_name = "first internal node";
        first_internal->_edge_length = 0.1;

        second_internal->_parent = first_internal;
        second_internal->_left_child = first_leaf;
        second_internal->_right_sib = third_leaf;
        second_internal->_number = 3;
        second_internal->_name = "second internal node";
        second_internal->_edge_length = 0.1;

        first_leaf->_parent = second_internal;
        first_leaf->_left_child = nullptr;
        first_leaf->_right_sib = second_leaf;
        first_leaf->_number = 0;
        first_leaf->_name = "first leaf";
        first_leaf->_edge_length = 0.1;

        second_leaf->_parent = second_internal;
        second_leaf->_left_child = nullptr;
        second_leaf->_right_sib = nullptr;
        second_leaf->_number = 1;
        second_leaf->_name = "second leaf";
        second_leaf->_edge_length = 0.1;

        third_leaf->_parent = first_internal;
        third_leaf->_left_child = nullptr;
        third_leaf->_right_sib = nullptr;
        third_leaf->_number = 2;
        third_leaf->_name = "third leaf";
        third_leaf->_edge_length = 0.1;

        _tree->_is_rooted = true;
        _tree->_root = root_node;
        _tree->_nleaves = 3;

        // Note that root node is not included in _preorder
        _tree->_preorder.push_back(first_internal);
        _tree->_preorder.push_back(second_internal);
        _tree->_preorder.push_back(first_leaf);
        _tree->_preorder.push_back(second_leaf);
        _tree->_preorder.push_back(third_leaf);

        _tree->_levelorder.push_back(first_internal);
        _tree->_levelorder.push_back(second_internal);
        _tree->_levelorder.push_back(third_leaf);
        _tree->_levelorder.push_back(first_leaf);
        _tree->_levelorder.push_back(second_leaf);
    }

    inline string TreeManip::makeNewick(unsigned precision, bool use_names) const {
        if (use_names && _taxon_names.empty()) {
            throw Xop("Cannot use taxon names in makeNewick when no taxon names have been saved");
        }
        string newick;
        const format tip_node_format( str(format("%%d:%%.%df") % precision) );
        const format tip_node_format_using_names( str(format("%%s:%%.%df") % precision) );
        const format internal_node_format( str(format("):%%.%df") % precision) );
        stack<Node *> node_stack;

        Node * root_tip = (_tree->_is_rooted ? nullptr : _tree->_root);
        for (auto nd : _tree->_preorder) {
            if (nd->_left_child) {
                newick += "(";
                node_stack.push(nd);
                if (root_tip) {
                    if (use_names) {
                        string with_underscores = std::regex_replace(_taxon_names[root_tip->_number], std::regex(" "), "_");
                        newick += str(format(tip_node_format_using_names) % with_underscores % nd->_edge_length);
                    }
                    else {
                        newick += str(format(tip_node_format) % (root_tip->_number + 1) % nd->_edge_length);
                    }
                    newick += ",";
                    root_tip = nullptr;
                }
            }
            else {
                if (use_names) {
                    string with_underscores = std::regex_replace(_taxon_names[nd->_number], std::regex(" "), "_");
                    newick += str(format(tip_node_format_using_names) % with_underscores % nd->_edge_length);
                }
                else
                    newick += str(format(tip_node_format) % (nd->_number + 1) % nd->_edge_length);
                if (nd->_right_sib)
                    newick += ",";
                else {
                    Node * popped = (node_stack.empty() ? nullptr : node_stack.top());
                    while (popped && !popped->_right_sib) {
                        node_stack.pop();
                        if (node_stack.empty()) {
                            newick += ")";
                            popped = nullptr;
                        }
                        else {
                            newick += str(format(internal_node_format) % popped->_edge_length);
                            popped = node_stack.top();
                        }
                    }
                    if (popped && popped->_right_sib) {
                        node_stack.pop();
                        newick += str(format(internal_node_format) % popped->_edge_length);
                        newick += ",";
                    }
                }
            }
        }

        return newick;
    }

    inline string TreeManip::debugMakeNewick(unsigned precision, bool use_names) const {
        if (use_names && _taxon_names.empty()) {
            throw Xop("Cannot use taxon names in makeNewick when no taxon names have been saved");
        }
        string newick;
        const format tip_node_format( str(format("%%d[&index=%%d]:%%.%df") % precision) );
        const format tip_node_format_using_names( str(format("%%s[&index=%%d]:%%.%df") % precision) );
        const format internal_node_format( str(format(")[&index=%%d]:%%.%df") % precision) );
        stack<Node *> node_stack;

        Node * root_tip = (_tree->_is_rooted ? nullptr : _tree->_root);
        for (auto nd : _tree->_preorder) {
            if (nd->_left_child) {
                newick += "(";
                node_stack.push(nd);
                if (root_tip) {
                    if (use_names) {
                        string with_underscores = std::regex_replace(_taxon_names[root_tip->_number], std::regex(" "), "_");
                        newick += str(format(tip_node_format_using_names) % with_underscores % nd->_number % nd->_edge_length);
                    }
                    else {
                        newick += str(format(tip_node_format) % (root_tip->_number + 1) % nd->_number % nd->_edge_length);
                    }
                    newick += ",";
                    root_tip = nullptr;
                }
            }
            else {
                if (use_names) {
                    string with_underscores = std::regex_replace(_taxon_names[nd->_number], std::regex(" "), "_");
                    newick += str(format(tip_node_format_using_names) % with_underscores % nd->_number % nd->_edge_length);
                }
                else
                    newick += str(format(tip_node_format) % (nd->_number + 1) % nd->_number % nd->_edge_length);
                if (nd->_right_sib)
                    newick += ",";
                else {
                    Node * popped = (node_stack.empty() ? nullptr : node_stack.top());
                    while (popped && !popped->_right_sib) {
                        node_stack.pop();
                        if (node_stack.empty()) {
                            newick += ")";
                            popped = nullptr;
                        }
                        else {
                            newick += str(format(internal_node_format) % popped->_number % popped->_edge_length);
                            popped = node_stack.top();
                        }
                    }
                    if (popped && popped->_right_sib) {
                        node_stack.pop();
                        newick += str(format(internal_node_format) % popped->_number % popped->_edge_length);
                        newick += ",";
                    }
                }
            }
        }

        return newick;
    }

    inline void TreeManip::extractNodeNumberFromName(Node * nd, set<unsigned> & used, unsigned nlvs) {
        assert(nd);
        bool success = true;
        unsigned x = 0;
        try {
            x = stoi(nd->_name);
            assert(x > 0);
            --x;
        } catch(invalid_argument &) {
            success = false;
        }

        if (!success) {
            // Node name could not be converted to an integer value
            // Assume node name is an actual taxon name
            string no_underscores = std::regex_replace(nd->_name, std::regex("_"), " ");
            try {
                x = _taxon_map.at(no_underscores);
            } catch(out_of_range &) {
                if (_taxon_map.size() == nlvs) {
                    throw Xop(str(format("taxon name \"%s\" is not a valid taxon name") % nd->_name));
                }
                // Add this taxon name to _taxon_names and _taxon_map
                x = static_cast<unsigned>(_taxon_names.size());
                _taxon_names.emplace_back(no_underscores);
                _taxon_map[no_underscores] = x;
            }
        }

        // attempt to insert x into the set of node numbers already used
        pair<set<unsigned>::iterator, bool> insert_result = used.insert(x);
        if (insert_result.second) {
            // insertion was made, so x has NOT already been used
            nd->_number = static_cast<int>(x);
        } else {
            // insertion was not made, so the set already contained x
            throw Xop(str(format("leaf number %d used more than once") % x));
        }
    }

    inline void TreeManip::extractEdgeLen(Node * nd, string edge_length_string) {
        assert(nd);
        bool success = true;
        double d = 0.0;
        try {
            d = stof(edge_length_string);
        }
        catch(invalid_argument &) {
            // edge_length_string could not be converted to a double value
            success = false;
        }

        if (success) {
            // conversion succeeded
            nd->_edge_length = (d < 0.0 ? 0.0 : d);
        }
        else
            throw Xop(str(format("%s is not interpretable as an edge length") % edge_length_string));
    }

    inline unsigned TreeManip::countNewickLeaves(const string & newick) {
        regex taxonexpr(R"([(,]\s*(\d+|\S+?|['].+?['])\s*(?=[,):]))");
        sregex_iterator m1(newick.begin(), newick.end(), taxonexpr);
        sregex_iterator m2;
        return static_cast<unsigned>(std::distance(m1, m2));
    }

    inline void TreeManip::stripOutNexusComments(string & newick) {
        regex commentexpr("\\[.*?\\]");
        newick = regex_replace(newick, commentexpr, string(""));
    }

    inline void TreeManip::debugCheckSplits() const {
        cerr << "debugCheckSplits: postorder traversal of tree..." << endl;
        for (auto nd : adaptors::reverse(_tree->_preorder)) {
            cerr << "  split = " << nd->_split.createPatternRepresentation() << "  number = " << nd->_number << "  name = \"" << nd->_name << "\"" << endl;
        }
        cerr << endl;
    }

    inline void TreeManip::debugCheckForSinNombreLeaves(const string & msg) const {
        for (auto nd : _tree->_preorder) {
            if (!nd->_left_child && nd->_name.size() == 0) {
                cerr << "***** sin nombre leaf *****" << endl;
                cerr << "***** " << msg << ": split = " << nd->_split.createPatternRepresentation() << "  number = " << nd->_number << "  name = \"" << nd->_name << "\"" << endl;
                cerr << "***** sin nombre leaf *****" << endl;
                throw Xop("sin nombre leaf detected");
            }
        }
    }

    inline unsigned TreeManip::debugCountUnusedNodes() const {
        // Find the number of nodes not currently in use
        unsigned nunused = 0;
        for (auto & node : _tree->_nodes) {
            Node * nd = &node;
            if (nd->_left_child == nullptr && nd->_right_sib == nullptr && nd->_parent == nullptr) {
                nunused++;
            }
        }
        return nunused;
    }

    inline void TreeManip::refreshPreorder() const {
        // Create the vector of node pointers in the preorder sequence
        _tree->_preorder.clear();
        _tree->_preorder.reserve(_tree->_nodes.size() - 1); // _preorder does not include root node

        if (!_tree->_root)
            return;

        Node * first_preorder = _tree->_root->_left_child;

        // sanity check: first preorder node should be the only child of the root node
        assert(first_preorder->_right_sib == nullptr);

        Node * nd = first_preorder;
        _tree->_preorder.push_back(nd);

        while (true)
        {
            if (!nd->_left_child && !nd->_right_sib)
            {
                // nd has no children and no siblings, so the next preorder is the right sibling of
                // the first ancestral node that has a right sibling.
                Node * anc = nd->_parent;
                while (anc && !anc->_right_sib)
                    anc = anc->_parent;
                if (anc)
                {
                    // We found an ancestor with a right sibling
                    _tree->_preorder.push_back(anc->_right_sib);
                    nd = anc->_right_sib;
                }
                else
                {
                    // nd is the last preorder node in the tree
                    break;
                }
            }
            else if (nd->_right_sib && !nd->_left_child)
            {
                // nd has no children (it is a tip), but does have a sibling on its right
                _tree->_preorder.push_back(nd->_right_sib);
                nd = nd->_right_sib;
            }
            else if (nd->_left_child && !nd->_right_sib)
            {
                // nd has children (it is an internal node) but no siblings on its right
                _tree->_preorder.push_back(nd->_left_child);
                nd = nd->_left_child;
            }
            else
            {
                // nd has both children and siblings on its right
                _tree->_preorder.push_back(nd->_left_child);
                nd = nd->_left_child;
            }

        }   // end while loop

        // renumber internal nodes in postorder sequence
        int curr_internal = static_cast<int>(_tree->_nleaves);
        for (auto nd : adaptors::reverse(_tree->_preorder))
        {
            if (nd->_left_child)
            {
                // nd is an internal node

                // node numbers for internal nodes start at _nleaves and go up
                nd->_number = curr_internal;
                curr_internal++;
            }
        }

        if (_tree->_is_rooted)
            _tree->_root->_number = curr_internal;
    }

    //                            1. start by adding only descendant of root node to buffer queue
    //                               queue = [1], stack = []
    //                            2. move node at front of buffer queue to back of stack vector
    //                               queue = [], stack = [1]
    //                 8    9     3. add this node's immediate children to back of buffer queue
    //                  \  /         queue = [2,3], stack = [1]
    //                   \/       4. repeat 2 and 3 until all nodes are processed
    //      4  5    6    7           (2) queue = [3], stack = [1,2]
    //       \ |     \  /            (3) queue = [3,4,5], stack = [1,2]
    //        \|      \/             (2) queue = [4,5], stack = [1,2,3]
    //         2      3              (3) queue = [4,5,6,7], stack = [1,2,3]
    //          \    /               (2) queue = [5,6,7], stack = [1,2,3,4]
    //           \  /                (3) no-op: 4 has no children
    //            \/                 (2) queue = [6,7], stack = [1,2,3,4,5]
    //            1                  (3) no-op: 5 has no children
    //            |                  (2) queue = [7], stack = [1,2,3,4,5,6]
    //            0                  (3) no-op: 6 has no children
    //                               (2) queue = [], stack = [1,2,3,4,5,6,7]
    //                               (3) queue = [8,9], stack = [1,2,3,4,5,6,7]
    //                               (2) queue = [9], stack = [1,2,3,4,5,6,7,8]
    //                               (3) no-op: 8 has no children
    //                               (2) queue = [], stack = [1,2,3,4,5,6,7,8,9]
    //                               (3) no-op: 9 has no children
    //                            5. stack vector is now in level order
    inline void TreeManip::refreshLevelorder() const {
        if (!_tree->_root)
            return;

        // q is the buffer queue
        queue<Node *> q;

        // _tree->_levelorder is the stack vector
        _tree->_levelorder.clear();
        _tree->_levelorder.reserve(_tree->_nodes.size() - 1);

        Node * nd = _tree->_root->_left_child;

        // sanity check: first node should be the only child of the root node
        assert(nd->_right_sib == nullptr);

        // Push nd onto the back of the queue
        q.push(nd);

        while (!q.empty())
        {
            // pop nd off front of queue
            nd = q.front(); q.pop();

            // and push it onto the stack
            _tree->_levelorder.push_back(nd);

            // add all children of nd to the back of the queue
            Node * child = nd->_left_child;
            if (child)
            {
                q.push(child);
                child = child->_right_sib;
                while (child)
                {
                    q.push(child);
                    child = child->_right_sib;
                }
            }
        }   // end while loop
    }

    inline bool TreeManip::canHaveSibling(const Node * nd, bool rooted, bool allow_polytomies) {
        assert(nd);
        if (!nd->_parent) {
            // trying to give the root node a sibling
            return false;
        }

        if (allow_polytomies)
            return true;

        bool nd_can_have_sibling = true;
        if (nd != nd->_parent->_left_child) {
            if (nd->_parent->_parent) {
                // trying to give a sibling to a sibling of nd, and nd's parent is not the root
                nd_can_have_sibling = false;
            }
            else {
                if (rooted) {
                    // the root node has exactly 2 children in rooted trees
                    nd_can_have_sibling = false;
                }
                else if (nd != nd->_parent->_left_child->_right_sib) {
                    // trying to give the root node more than 3 children
                    nd_can_have_sibling = false;
                }
            }
        }

        return nd_can_have_sibling;
    }

    inline void TreeManip::rerootAt(int node_number) const {
        // Locate the node having a _number equal to node_number
        Node * nd = nullptr;
        for (auto & curr : _tree->_nodes) {
            if (curr._number == node_number) {
                nd = &curr;
                break;
            }
        }
        if (!nd)
            throw Xop(str(format("no node found with number equal to %d") % node_number));

        if (nd->_left_child)
            throw Xop(str(format("cannot currently root trees at internal nodes (e.g. node %d)") % nd->_number));

        Node * t = nd;
        Node * m = nd->_parent;
        while (nd->_parent) {
            // Begin by swapping the mover's edge length with nd's edge length
            double tmp_edgelen = m->_edge_length;
            m->_edge_length = nd->_edge_length;
            nd->_edge_length = tmp_edgelen;

            // Make the "mover" node m (along with all of its children except nd) the rightmost child of the "target" node t
            rerootHelper(m, t);

            // The next target node is always the previous mover, and the next mover node is always nd's parent
            t = m;
            m = nd->_parent;
        }
        _tree->_root = nd;
    }

    inline void TreeManip::rerootHelper(Node * m, Node * t) {
        assert(m);
        assert(t);

        // Save nodes to which m attaches
        Node * m_left_child    = m->_left_child;
        Node * m_right_sib     = m->_right_sib;
        Node * m_parent        = m->_parent;

        // Starting with t, walk down the tree to identify x, the child of m that is on the path from m to t
        Node * x = t;
        while (x->_parent != m) {
            x = x->_parent;
            assert(x);
        }
        Node * x_right_sib = x->_right_sib;

        // Identify x_left_sib, the immediate left sibling of x (will be NULL if x is _left_child of m)
        Node * x_left_sib = nullptr;
        if (x != m_left_child) {
            x_left_sib = m_left_child;
            while (x_left_sib->_right_sib != x) {
                x_left_sib = x_left_sib->_right_sib;
                assert(x_left_sib);
            }
        }

        // identify m_left_sib, the immediate left sibling of m (will be NULL if m is the root node or is _left_child of its parent)
        Node * m_left_sib = nullptr;
        if (m_parent && m != m_parent->_left_child) {
            m_left_sib = m_parent->_left_child;
            while (m_left_sib->_right_sib != m) {
                m_left_sib = m_left_sib->_right_sib;
                assert(m_left_sib);
            }
        }

        // Put x where m is now
        if (!m_parent) {
            // m is the root node
            assert(!m_right_sib);
            assert(!m_left_sib);
            x->_right_sib = nullptr;
            x->_parent = nullptr;
            if (x == m_left_child)
                m->_left_child = x_right_sib;
            else
                x_left_sib->_right_sib = x_right_sib;
        }
        else if (m == m_parent->_left_child) {
            // m is the leftmost child of its parent
            x->_right_sib = m_right_sib;
            x->_parent = m_parent;
            m->_right_sib = nullptr;
            m->_parent = nullptr;
            m_parent->_left_child = x;
            if (x == m_left_child)
                m->_left_child = x_right_sib;
            else
                x_left_sib->_right_sib = x_right_sib;
        }
        else {
            // m is not leftmost child of its parent
            m_left_sib->_right_sib = x;
            x->_right_sib = m_right_sib;
            x->_parent = m_parent;
            m->_right_sib = nullptr;
            m->_parent = nullptr;
            if (x == m_left_child)
                m->_left_child = x_right_sib;
            else
                x_left_sib->_right_sib = x_right_sib;
        }

        // Make m the new rightmost child of t
        m->_parent = t;
        if (!t->_left_child)
            t->_left_child = m;
        else {
            // Find rightmost child of t
            m_left_sib = t->_left_child;
            while (m_left_sib->_right_sib)
                m_left_sib = m_left_sib->_right_sib;

            // Make rightmost child of t the left sib of m
            m_left_sib->_right_sib = m;
        }
    }

    inline void TreeManip::setLeafNames(const vector<string> & leafnames) const {
        for (auto nd : _tree->_preorder) {
            if (nd->_number < static_cast<int>(leafnames.size())) {
                nd->_name = leafnames[nd->_number];
            }
        }
    }

    inline void TreeManip::buildFromNewick(const string & newick, bool rooted, bool allow_polytomies) {
        _tree.reset();
        _tree = std::make_shared<Tree>();
#if defined(ALWAYS_ROOTED)
        _tree->_is_rooted = true;
#else
        _tree->_is_rooted = rooted;
#endif

        string commentless_newick = newick;
        stripOutNexusComments(commentless_newick);

        _tree->_nleaves = countNewickLeaves(commentless_newick);
        if (_tree->_nleaves == 0)
            throw Xop("Expecting newick tree description to have at least 4 leaves");

#if defined(ALWAYS_ROOTED)
        unsigned max_nodes = 2*_tree->_nleaves;
#else
        unsigned max_nodes = 2*_tree->_nleaves - (_tree->_is_rooted ? 0 : 2);
#endif
        _tree->_nodes.resize(max_nodes);

        // Assign all nodes a default node number that is negative to make it easy to tell if we've failed to set it.
        // Leaves will replace this number with the number equivalent of their name. Internal nodes will replace
        // this number with a number higher than any leaf.
        for (unsigned i = 0; i < max_nodes; ++i) {
            _tree->_nodes[i]._number = -1;
        }

        // This will point to the first tip node encountered so that we can reroot at this node before returning

        try {
            unsigned curr_node_index = 0;
            unsigned num_edge_lengths = 0;
            Node * first_tip = nullptr;
            unsigned curr_leaf = 0;
            set<unsigned> used;
            // Root node
            Node * nd = &_tree->_nodes[curr_node_index];
            _tree->_root = nd;

#if defined(ALWAYS_ROOTED)
            nd = &_tree->_nodes[++curr_node_index];
            nd->_parent = &_tree->_nodes[curr_node_index - 1];
            nd->_parent->_left_child = nd;
#else
            if (_tree->_is_rooted) {
                nd = &_tree->_nodes[++curr_node_index];
                nd->_parent = &_tree->_nodes[curr_node_index - 1];
                nd->_parent->_left_child = nd;
            }
#endif

            // Some flags to keep track of what we did last
            enum {
                Prev_Tok_LParen     = 0x01, // the previous token was a left parenthesis ('(')
                Prev_Tok_RParen     = 0x02, // the previous token was a right parenthesis (')')
                Prev_Tok_Colon      = 0x04, // the previous token was a colon (':')
                Prev_Tok_Comma      = 0x08, // the previous token was a comma (',')
                Prev_Tok_Name       = 0x10, // the previous token was a node name (e.g. '2', 'P._articulata')
                Prev_Tok_EdgeLen    = 0x20  // the previous token was an edge length (e.g. '0.1', '1.7e-3')
                };
            unsigned previous = Prev_Tok_LParen;

            // Some useful flag combinations
            unsigned LParen_Valid = (Prev_Tok_LParen | Prev_Tok_Comma);
            unsigned RParen_Valid = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_EdgeLen);
            unsigned Comma_Valid  = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_EdgeLen);
            unsigned Colon_Valid  = (Prev_Tok_RParen | Prev_Tok_Name);
            unsigned Name_Valid   = (Prev_Tok_RParen | Prev_Tok_LParen | Prev_Tok_Comma);

            // Set to true while reading an edge length
            bool inside_edge_length = false;
            string edge_length_str;
            unsigned edge_length_position = 0;

            // Set to true while reading a node name surrounded by (single) quotes
            bool inside_quoted_name = false;

            // Set to true while reading a node name not surrounded by (single) quotes
            bool inside_unquoted_name = false;

            // Set to start of each node name and used in case of error
            unsigned node_name_position = 0;

            // loop through the characters in newick, building up the tree as we go
            unsigned position_in_string = 0;
            for (auto ch : commentless_newick) {
                position_in_string++;

                if (inside_quoted_name) {
                    if (ch == '\'') {
                        inside_quoted_name = false;
                        node_name_position = 0;
                        if (!nd->_left_child) {
                            extractNodeNumberFromName(nd, used, _tree->_nleaves);
                            curr_leaf++;
                            if (!first_tip)
                                first_tip = nd;
                        }
                        previous = Prev_Tok_Name;
                    }
                    else if (iswspace(ch))
                        nd->_name += ' ';
                    else
                        nd->_name += ch;

                    continue;
                }
                else if (inside_unquoted_name) {
                    if (ch == '(')
                        throw Xop(str(format("Unexpected left parenthesis inside node name at position %d in tree description\n%s") % node_name_position % commentless_newick));

                    if (iswspace(ch) || ch == ':' || ch == ',' || ch == ')') {
                        inside_unquoted_name = false;

                        // Expect node name only after a left paren (child's name), a comma (sib's name) or a right paren (parent's name)
                        if (!(previous & Name_Valid))
                            throw Xop(str(format("Unexpected node name (%s) at position %d in tree description\n%s") % nd->_name % node_name_position % commentless_newick));

                        if (!nd->_left_child) {
                            extractNodeNumberFromName(nd, used, _tree->_nleaves);
                            curr_leaf++;
                            if (!first_tip)
                                first_tip = nd;
                        }

                        previous = Prev_Tok_Name;
                    }
                    else {
                        nd->_name += ch;
                        continue;
                    }
                }
                else if (inside_edge_length) {
                    if (ch == ',' || ch == ')' || ch == ';' || iswspace(ch)) {
                        inside_edge_length = false;
                        edge_length_position = 0;
                        extractEdgeLen(nd, edge_length_str);
                        ++num_edge_lengths;
                        previous = Prev_Tok_EdgeLen;
                    }
                    else {
                        bool valid = (ch =='e' || ch == 'E' || ch =='.' || ch == '-' || ch == '+' || isdigit(ch));
                        if (!valid)
                            throw Xop(str(format("Invalid branch length character (%c) at position %d in tree description\n%s") % ch % position_in_string % commentless_newick));
                        edge_length_str += ch;
                        continue;
                    }
                }

                if (iswspace(ch))
                    continue;

                switch(ch) {
                    case ';':
                        break;

                    case ')':
                        // If nd is the bottommost node, expecting left paren or semicolon, but not right paren
                        if (!nd->_parent)
                            throw Xop(str(format("Too many right parentheses at position %d in tree description\n%s") % position_in_string % commentless_newick));

                        // Expect right paren only after an edge length, a node name, or another right paren
                        if (!(previous & RParen_Valid))
                            throw Xop(str(format("Unexpected right parenthesisat position %d in tree description\n%s") % position_in_string % commentless_newick));

                        // Go down a level
                        nd = nd->_parent;
                        if (!nd->_left_child->_right_sib)
                            throw Xop(str(format("Internal node has only one child at position %d in tree description\n%s") % position_in_string % commentless_newick));
                        previous = Prev_Tok_RParen;
                        break;

                    case ':':
                        // Expect colon only after a node name or another right paren
                        if (!(previous & Colon_Valid))
                            throw Xop(str(format("Unexpected colon at position %d in tree description\n%s") % position_in_string % commentless_newick));
                        previous = Prev_Tok_Colon;
                        break;

                    case ',':
                        // Expect comma only after an edge length, a node name, or a right paren
                        if (!nd->_parent || !(previous & Comma_Valid))
                            throw Xop(str(format("Unexpected comma at position %d in tree description\n%s") % position_in_string % commentless_newick));

#if defined(ALWAYS_ROOTED)
                        // Check for polytomies
                        if (!canHaveSibling(nd, true, allow_polytomies))
                            throw Xop(str(format("Polytomy found in the following tree description but polytomies prohibited:\n%s") % commentless_newick));
#else
                        // Check for polytomies
                        if (!canHaveSibling(nd, rooted, allow_polytomies))
                            throw Xop(str(format("Polytomy found in the following tree description but polytomies prohibited:\n%s") % commentless_newick));
#endif

                        // Create the sibling
                        curr_node_index++;
                        if (curr_node_index == _tree->_nodes.size())
                            throw Xop(str(format("Too many nodes (%d nodes allocated for %d leaves) at position %d in tree description\n%s") % _tree->_nodes.size() % _tree->_nleaves % position_in_string % commentless_newick));
                        nd->_right_sib = &_tree->_nodes[curr_node_index];
                        nd->_right_sib->_parent = nd->_parent;
                        nd = nd->_right_sib;
                        previous = Prev_Tok_Comma;
                        break;

                    case '(':
                        // Expect left paren only after a comma or another left paren
                        if (!(previous & LParen_Valid))
                            throw Xop(str(format("Not expecting left parenthesis at position %d in tree description\n%s") % position_in_string % commentless_newick));

                        // Create a new node above and to the left of the current node
                        assert(!nd->_left_child);
                        curr_node_index++;
                        if (curr_node_index == _tree->_nodes.size())
                            throw Xop(str(format("malformed tree description (more than %d nodes specified)\n%s") % _tree->_nodes.size() % commentless_newick));
                        nd->_left_child = &_tree->_nodes[curr_node_index];
                        nd->_left_child->_parent = nd;
                        nd = nd->_left_child;
                        previous = Prev_Tok_LParen;
                        break;

                    case '\'':
                        // Encountered an apostrophe, which always indicates the start of a
                        // node name (but note that node names do not have to be quoted)

                        // Expect node name only after a left paren (child's name), a comma (sib's name)
                        // or a right paren (parent's name)
                        if (!(previous & Name_Valid))
                            throw Xop(str(format("Not expecting node name at position %d in tree description\n%s") % position_in_string % commentless_newick));

                        // Get the rest of the name
                        nd->_name.clear();

                        inside_quoted_name = true;
                        node_name_position = position_in_string;

                        break;

                    default:
                        // Get here if ch is not "(", ")", ";", ":", ",", or "'"

                        // Expecting either an edge length or an unquoted node name
                        if (previous == Prev_Tok_Colon) {
                            // Edge length expected (e.g. "235", "0.12345", "1.7e-3")
                            inside_edge_length = true;
                            edge_length_position = position_in_string;
                            edge_length_str = ch;
                        }
                        else {
                            // Get the node name
                            nd->_name = ch;

                            inside_unquoted_name = true;
                            node_name_position = position_in_string;
                        }

                }   // end of the switch statement
            }   // loop over characters in newick string

            if (inside_unquoted_name)
                throw Xop(str(format("Tree description ended before end of node name starting at position %d was found\n%s") % node_name_position % commentless_newick));
            if (inside_edge_length)
                throw Xop(str(format("Tree description ended before end of edge length starting at position %d was found\n%s") % edge_length_position % commentless_newick));
            if (inside_quoted_name)
                throw Xop(str(format("Expecting single quote to mark the end of node name at position %d in tree description\n%s") % node_name_position % commentless_newick));

#if defined(ALWAYS_ROOTED)
#else
            if (!_tree->_is_rooted) {
                // Root at leaf whose _number = 0
                rerootAt(0);
            }
#endif

            refreshPreorder();
            refreshLevelorder();
        }
        catch(Xop &) {
            _tree.reset();
            throw;
        }
    }

    inline void TreeManip::setEdgeLength(const Split & s, double edge_length) const {
        // Assumes nodes in the tree have splits set correctly
        // Find split s in the tree and change that edge length
        for (auto nd : _tree->_preorder) {
            if (s == nd->_split) {
                nd->_split.setEdgeLen(edge_length);
                nd->setEdgeLength(edge_length);

                return;
            }
        }
    }

    inline void TreeManip::storeSplits(set<Split> & internal_splits, set<Split> & leaf_splits) const {
        // Start by clearing and resizing all splits
        for (auto & nd : _tree->_nodes) {
            nd._split.resize(_tree->_nleaves);
        }

        // Now do a postorder traversal and add the bit corresponding
        // to the current node in its parent node's split
        for (auto nd : adaptors::reverse(_tree->_preorder)) {
            // Set split's edge length
            nd->_split.setEdgeLen(nd->_edge_length);

            if (nd->_left_child) {
                if (nd != _tree->_root->_left_child) {
                    internal_splits.insert(nd->_split);
                }
                // add this internal node's split to splitset
            }
            else {
                // set the bit corresponding to this leaf node's number
                nd->_split.setBitAt(nd->_number);
                leaf_splits.insert(nd->_split);
            }

            if (nd->_parent) {
                // parent's bits are the union of the bits set in all its children
                nd->_parent->_split.addSplit(nd->_split);
            }
        }
    }

    inline void TreeManip::dropSplit(const Split & s) const {
        // Locate node having split s
        auto it = find_if(_tree->_preorder.begin(), _tree->_preorder.end(),[s](const Node * nd){return nd->_split == s;});
        assert(it != _tree->_preorder.end());
        Node * nd = *it;

        // Make each of nd's children a child of nd->_parent
        Node * par = nd->_parent;
        assert(par);
        stack<Node *> child_pile;
        for (Node * child = nd->_left_child; child; child = child->_right_sib) {
            child_pile.push(child);
        }
        nd->_left_child = nullptr;
        while (!child_pile.empty()) {
            // Remove child from the stack
            Node * child = child_pile.top();
            child_pile.pop();

            // Attach child to par
            child->_right_sib = par->_left_child;
            par->_left_child = child;
            child->_parent = par;
        }

        // Now we can drop nd
        Node * par_child = par->_left_child;
        assert(par_child);
        if (par_child == nd) {
            // nd is its parent's left child
            par->_left_child = nd->_right_sib;
        }
        else {
            // nd is not its parent's left child
            while (par_child->_right_sib != nd) {
                par_child = par_child->_right_sib;
            }
            assert(par_child);
            par_child->_right_sib = nd->_right_sib;
        }

        // Clear everything except node number and name
        nd->_left_child = nullptr;
        nd->_right_sib = nullptr;
        nd->_parent = nullptr;
        nd->_edge_length = Node::_smallest_edge_length;
        nd->_split.clear();

        refreshPreorder();
        refreshLevelorder();
        //debugCheckSplits();
    }

    inline void TreeManip::addSplit(const Split & s) const {
        // Find the first node not currently in use
        Node * newnd = nullptr;
        for (auto & node : _tree->_nodes) {
            Node * nd = &node;
            if (nd->_left_child == nullptr && nd->_right_sib == nullptr && nd->_parent == nullptr) {
                newnd = nd;
                break;
            }
        }

        // If no unused node was found, we will not be able to add any splits
        assert (newnd != nullptr);

        // Find the first node in the postorder sequence that contains s
        for (auto nd : adaptors::reverse(_tree->_preorder)) {
            Split & ndsplit = nd->_split;
            if (s.subsumedIn(ndsplit)) {
                // Add all nodes subsumed in s to a new node and detach them from nd
                stack<Node *> child_pile;
                for (Node * child = nd->_left_child; child; child = child->_right_sib) {
                    child_pile.push(child);
                }
                while (!child_pile.empty()) {
                    // Remove child from the stack
                    Node * child = child_pile.top();
                    child_pile.pop();
                    Split & childsplit = child->_split;
                    if (childsplit.subsumedIn(s)) {
                        // Move child to newnd
                        Node * left_sib = (child_pile.empty() ? nullptr : child_pile.top());
                        if (left_sib) {
                            left_sib->_right_sib = child->_right_sib;
                        }
                        else {
                            nd->_left_child = child->_right_sib;
                        }
                        Node * oldlchild = newnd->_left_child;
                        newnd->_left_child = child;
                        child->_parent = newnd;
                        child->_right_sib = oldlchild;
                    }
                }

                // Make newnd a child of nd
                Node * oldndlchild = nd->_left_child;
                nd->_left_child = newnd;
                newnd->_right_sib = oldndlchild;
                newnd->_parent = nd;
                newnd->_split = s;
                break;
            }
        }

        refreshPreorder();
        refreshLevelorder();
        //debugCheckSplits();
    }

    string TreeManip::nexusTranslateCommand() {
        string translate_command = "  translate\n";
        vector<string> entries;
        for (unsigned i = 0; i < _taxon_names.size(); i++) {
            string with_underscores = std::regex_replace(_taxon_names[i], std::regex(" "), "_");
            entries.push_back(str(format("    %d %s") % (i+1) % with_underscores));
        }
        translate_command += join(entries, ",\n");
        translate_command += "\n  ;\n";
        return translate_command;
    }

    inline void TreeManip::stripComments(string & newick) {
        newick = regex_replace(newick, regex("\\[.+?\\]"), "");
    }

    inline string TreeManip::renumberNewick(const string & newick, map<unsigned, unsigned> taxon_index_map) {
        // Returns a newick tree description in which the taxa (assumed to be numbers) have been
        // translated by taxon_index_map. For example, assume taxon number in newick equals 3, which
        // means the index of this taxon is 3-1 = 2. Suppose taxon_index_map[2] = 0, then renumbered_newick
        // will have 0 + 1 = 1 where the original newick had 3.
        string renumbered_newick;

        // The number of replacements should equal the size of taxon_index_map
        unsigned num_replacements = 0;

        // Construct a regular expression pattern that identifies taxon numbers in newick descriptions
        regex taxon_regex(R"([(,]\s*(\d+)\s*(?=[,):]))");
        //                      [(,] <-- leaves always follow a left parenthesis or a comma
        //                          \s* <-- possibly some whitespace before taxon name
        //                             (\d+) <-- a number
        //                                  \s* <-- possibly some whitespace after name
        //                                     (?=[,):]) <-- lookahead to comma, semicolon, or right paren

        // Create a match object
        std::smatch taxon_match;

        // Create a variable to store the starting iterator because we will need
        // to change this after each match is found
        auto search_start = newick.cbegin();

        // These keep track of which parts of the newick string have already been dealt with
        // (i.e. before prev_pos) and the position of the current match (i.e. curr_pos)
        string::size_type prev_pos = 0;
        string::size_type curr_pos = 0;

        // Search through the part of newick from search_start onward, stopping if there is a match
        while (std::regex_search(search_start, newick.cend(), taxon_match, taxon_regex)) {
            // A match was found, so extract the taxon number
            unsigned taxon_number = 0;
            try {
                taxon_number = stod(taxon_match[1]);
            } catch (const std::invalid_argument &) {
                throw Xop(format("In TreeManip::renumberNewick, could not match taxon number at position %d in newick tree description \"%s\"") % taxon_match.position(1) % newick);
            }

            // Taxon number should follow nexus 1-offset standard
            assert(taxon_number > 0);

            // Get the index of this taxon by subtracting 1
            unsigned taxon_index = taxon_number - 1;

            // Look up the new index in taxon_index_map
            unsigned new_index = taxon_index_map.at(taxon_index);

            // The new taxon number will be the taxon index plus 1
            unsigned new_number = new_index + 1;

            // Get the length of the match (i.e. number of characters in the taxon number)
            string::size_type match_len = taxon_match.length(1);

            // Update curr_pos to be the position, relative to search_atart, of the start of the match
            curr_pos = taxon_match.position(1);

            // Get that part of newick starting at prev_pos and ending just before the current match
            // Note that curr_pos here is the length of the substring
            string prev_text = newick.substr(prev_pos, curr_pos);

            // Append prev_text and the new taxon number to renumbered_newick
            renumbered_newick += str(format("%s%d") % prev_text % new_number);

            // Update prev_oos by adding the number of characters since search_start (curr_oos)
            // and the length of the match (match_len)
            prev_pos += curr_pos + match_len;

            // Update search start so that the next search will start where this one left off
            // Why first? taxon_match.suffix() returns a pair of iterators pointing to the
            // beginning and end of the suffix. We want to begin again at the start of the
            // suffix, hence we let search_start equal suffix().first.
            search_start = taxon_match.suffix().first;

            num_replacements++;
        }

        // Add the part of newick that is beyond the last match
        // Note: not specifying a length means that the substring extends
        // from prev_pos to the end of the string
        string final_text = newick.substr(prev_pos);
        renumbered_newick += final_text;

        if (num_replacements != taxon_index_map.size()) {
            throw Xop(format("In TreeManip::renumberNewick, number of replacements (%d) does not equal size of taxon_index_map (%d)") % num_replacements % taxon_index_map.size());
        }

        return renumbered_newick;
    }

    inline void TreeManip::nniNodeSwap(Node * a, Node * b) {
        // Doesn't assume binary tree
        Node *x = a->_parent;
        assert(x);

        Node * y = b->_parent;
        assert(y);
        assert(x->_parent == y);

        bool a_is_left_child = (a == x->_left_child);
        bool x_is_left_child = (x == y->_left_child);

        if (a_is_left_child && x_is_left_child) {
            // Case 1:
            //     a                        b
            //      \   /                    \   /
            //       \ /                      \ /
            //        x     b            a     x
            //         \   /              \   /
            //          \ /      ==>       \ /
            //           y                  y
            //           |                  |
            //           |                  |
            //           |                  |

            // Detach a from tree
            x->_left_child = a->_right_sib;
            a->_parent = nullptr;
            a->_right_sib = nullptr;

            // Detach b from tree
            Node * ychild = y->_left_child;
            while (ychild && ychild->_right_sib != b) {
                ychild = ychild->_right_sib;
            }
            assert(ychild->_right_sib == b);
            ychild->_right_sib = b->_right_sib;
            b->_parent = nullptr;
        }
        else if (a_is_left_child && !x_is_left_child) {
            // Case 2:
            //         a                  b
            //          \   /              \   /
            //           \ /                \ /
            //      b     x            a     x
            //       \   /              \   /
            //        \ /      ==>       \ /
            //         y                  y
            //         |                  |
            //         |                  |
            //         |                  |

            // Detach a from tree
            x->_left_child = a->_right_sib;
            a->_parent = nullptr;
            a->_right_sib = nullptr;

            // Detach b from tree
            if (b == y->_left_child) {
                y->_left_child = b->_right_sib;
                b->_parent = nullptr;
            }
            else {
                Node * ychild = y->_left_child;
                while (ychild && ychild->_right_sib != b) {
                    ychild = ychild->_right_sib;
                }
                assert(ychild->_right_sib == b);
                ychild->_right_sib = b->_right_sib;
                b->_parent = nullptr;
            }
        }
        else if (!a_is_left_child && x_is_left_child) {
            // Case 3
            //           a                  b
            //      \   /                    \   /
            //       \ /                      \ /
            //        x     b            a     x
            //         \   /              \   /
            //          \ /      ==>       \ /
            //           y                  y
            //           |                  |
            //           |                  |
            //           |                  |

            // Detach a from tree
            Node * xchild = x->_left_child;
            while (xchild && xchild->_right_sib != a) {
                xchild = xchild->_right_sib;
            }
            assert(xchild->_right_sib == a);
            xchild->_right_sib = a->_right_sib;
            a->_parent = nullptr;

            // Detach b from tree
            Node * ychild = y->_left_child;
            while (ychild && ychild->_right_sib != b) {
                ychild = ychild->_right_sib;
            }
            assert(ychild->_right_sib == b);
            ychild->_right_sib = b->_right_sib;
            b->_parent = nullptr;

        }
        else if (!a_is_left_child && !x_is_left_child) {
            // Case 4:
            //               a            b
            //          \   /              \   /
            //           \ /                \ /
            //      b     x            a     x
            //       \   /              \   /
            //        \ /      ==>       \ /
            //         y                  y
            //         |                  |
            //         |                  |
            //         |                  |

            Node * xchild = x->_left_child;
            while (xchild && xchild->_right_sib != a) {
                xchild = xchild->_right_sib;
            }
            assert(xchild->_right_sib == a);
            xchild->_right_sib = a->_right_sib;
            a->_parent = nullptr;

            // Detach b from tree
            if (b == y->_left_child) {
                y->_left_child = b->_right_sib;
                b->_parent = nullptr;
            }
            else {
                Node * ychild = y->_left_child;
                while (ychild && ychild->_right_sib != b) {
                    ychild = ychild->_right_sib;
                }
                assert(ychild->_right_sib == b);
                ychild->_right_sib = b->_right_sib;
                b->_parent = nullptr;
            }
        }
        else {
            throw Xop(format("Unexpected case in TreeManip::nniNodeSwap()"));
        }

        // Reattach a to y
        Node * ylchild = y->_left_child;
        a->_right_sib = ylchild;
        a->_parent = y;
        y->_left_child = a;

        // Reattach b to x
        Node * xlchild = x->_left_child;
        b->_right_sib = xlchild;
        b->_parent = x;
        x->_left_child = b;

        refreshPreorder();
    }
}
