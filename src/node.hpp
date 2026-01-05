#pragma once

namespace op
{

    class Tree;
    class TreeManip;
    class Likelihood;
    //class Updater;

    class Node {
            friend class Tree;
            friend class TreeManip;
            friend class Likelihood;
            //friend class Updater;

        public:
                                        Node();
                                        ~Node() = default;

                    Node *              getParent() const    {return _parent;}
                    Node *              getLeftChild() const {return _left_child;}
                    Node *              getRightSib() const  {return _right_sib;}
                    int                 getNumber() const    {return _number;}
                    string              getName()            {return _name;}
                    Split               getSplit()           {return _split;}

                    double              getEdgeLength() const {return _edge_length;}
                    void                setEdgeLength(double v);

                    void                describeNode(const string & title) const;

            static const double _smallest_edge_length;

            typedef vector<Node>    Vector;
            typedef vector<Node *>  PtrVector;

        private:

            void                clear();

            Node *              _left_child;
            Node *              _right_sib;
            Node *              _parent;
            int                 _number;
            string              _name;
            double              _edge_length;
            Split               _split;
    };

    inline Node::Node() {
        _left_child = nullptr;
        _right_sib = nullptr;
        _parent = nullptr;
        _number = 0;
        _name = "";
        _edge_length = _smallest_edge_length;
        _split.clear();
    }

    inline void Node::clear() {
        _left_child = nullptr;
        _right_sib = nullptr;
        _parent = nullptr;
        _number = 0;
        _name = "";
        _edge_length = _smallest_edge_length;
        _split.clear();
    }

    inline void Node::setEdgeLength(double v) {
        _edge_length = (v < _smallest_edge_length ? _smallest_edge_length : v);
    }

    inline void Node::describeNode(const string & title) const {
        cerr << title << ": number = " << _number << ", name = \"" << _name << "\"" << endl;
        if (_parent)
            cerr << "  _parent: number = " << _parent->_number << ", name = \"" << _parent->_name << "\"" << endl;
        else
            cerr << "  _parent: NULL" << endl;
        if (_left_child)
            cerr << "  _left_child: number = " << _left_child->_number << ", name = \"" << _left_child->_name << "\"" << endl;
        else
            cerr << "  _left_child: NULL" << endl;
        if (_right_sib)
            cerr << "  _right_sib: number = " << _right_sib->_number << ", name = \"" << _right_sib->_name << "\"" << endl;
        else
            cerr << "  _right_sib: NULL" << endl;
    }

}

