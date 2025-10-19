#pragma once

using namespace boost;

namespace op {

    class TreeManip;
    class Likelihood;
    //class Updater;

    class Tree {

        friend class TreeManip;
        friend class Likelihood;

        public:

                                        Tree();
                                        ~Tree() = default;

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

    inline Tree::Tree() {
        _is_rooted = false;
        _root = nullptr;
        _nleaves = 0;
    }

    inline void Tree::clear() {
        _is_rooted = false;
        _root = nullptr;
        _nleaves = 0;
        _nodes.clear();
        _preorder.clear();
        _levelorder.clear();
    }

    inline bool Tree::isRooted() const {
        return _is_rooted;
    }

    inline unsigned Tree::numLeaves() const {
        return _nleaves;
    }

    inline unsigned Tree::numNodes() const {
        return static_cast<unsigned>(_nodes.size());
    }

}
