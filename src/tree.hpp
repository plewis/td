#pragma once

#include <memory>
#include <iostream>
#include "node.hpp"

using namespace std;
using namespace boost;

namespace op
    {

    class TreeManip;
    class Likelihood;
    //class Updater;

    class Tree
        {

        friend class TreeManip;
        friend class Likelihood;
        //friend class Updater;

        public:

                                        Tree();
                                        ~Tree();

            bool                        isRooted() const;
            unsigned                    numLeaves() const;
            unsigned                    numNodes() const;

        private:

            void                        clear();

            bool                        _is_rooted;
            Node *                      _root;
            unsigned                    _nleaves;
            Node::PtrVector             _preorder;
            Node::PtrVector             _levelorder;
            Node::Vector                _nodes;

        public:

            typedef std::shared_ptr< Tree > SharedPtr;
        };

    inline Tree::Tree()
        {
        //cout << "Constructing a Tree" << endl;
        clear();
        }

    inline Tree::~Tree()
        {
        //cout << "Destroying a Tree" << endl;
        }

    inline void Tree::clear()
        {
        _is_rooted = false;
        _root = 0;
        _nodes.clear();
        _preorder.clear();
        _levelorder.clear();
        }

    inline bool Tree::isRooted() const
        {
        return _is_rooted;
        }

    inline unsigned Tree::numLeaves() const
        {
        return _nleaves;
        }

    inline unsigned Tree::numNodes() const
        {
        return (unsigned)_nodes.size();
        }

    }
