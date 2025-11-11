#pragma once

#include <vector>
#include <memory>
#include <set>
#include <map>
#include <climits>
#include <cassert>

using namespace std;
using namespace boost;

namespace op
    {

    class Split
        {
        public:
                                                                Split();
                                                                ~Split() = default;

            bool                                                operator==(const Split & other) const;
            bool                                                operator!=(const Split & other) const;
            bool                                                operator<(const Split & other) const;

            void                                                clear();
            void                                                resize(unsigned nleaves);

            void                                                setBitsPerUnit(unsigned nbits_per_unit);
            unsigned                                            getBitsPerUnit() const;

            typedef unsigned long                               split_unit_t;
            typedef vector<split_unit_t>                        split_t;
            typedef set<Split>                                  treeid_t;
            typedef pair<treeid_t, treeid_t>                    treeid_pair_t;
            typedef pair<Split, Split>                          split_pair_t;
            typedef map< treeid_t, vector<unsigned> >           treemap_t;

            void                                                setEdgeLen(double edgelen);
            double                                              getEdgeLen() const;

            const split_t &                                     getBits() const;

            unsigned                                            getSize() const;
            unsigned                                            getNumBitsSet() const;

            double                                              entropy() const;
            double                                              mutualClusteringInfo(const Split & other) const;

            void                                                bitwiseAnd(const Split & other);
            unsigned                                            findFirstSetBit() const;
            void                                                invertBits();
            bool                                                subsumedIn(const Split & other) const;
            bool                                                compatibleWith(const Split & other) const;

            void                                                setBitAt(unsigned leaf_index);
            void                                                setBitsAt(vector<unsigned> leaves);
            bool                                                isBitSetAt(unsigned leaf_index) const;
            void                                                addSplit(const Split & other);

            string                                              createPatternRepresentation(bool show_edge_length = false) const;

        private:

            double                                              _edgelen;
            split_t                                             _bits;
            unsigned                                            _bits_per_unit;
            unsigned                                            _nleaves;

        public:

            typedef std::shared_ptr< Split >                    SharedPtr;
    };

inline Split::Split() {
    _edgelen = 0.0;
    _nleaves = 0;
    _bits_per_unit = (CHAR_BIT)*sizeof(Split::split_unit_t);
    clear();
    //cout << "Constructing a Split" << endl;
}

inline void Split::clear() {
    _edgelen = 0.0;
    for (auto & u : _bits) {
        u = 0L;
    }
}

inline bool Split::operator==(const Split & other) const {
    return (_bits == other._bits);
}

inline bool Split::operator!=(const Split & other) const {
    return (_bits != other._bits);
}

inline bool Split::operator<(const Split & other) const {
    assert(_bits.size() == other._bits.size());
    return (_bits < other._bits);
}

inline const Split::split_t & Split::getBits() const {
    return _bits;
}

inline unsigned Split::getSize() const {
    return _nleaves;
}

inline double Split::entropy() const {
    assert(_nleaves > 0);
    auto b = static_cast<double>(getNumBitsSet());
    auto n = static_cast<double>(_nleaves);
    double p = b/n;
    double q = 1.0 - p;
    double h = -p*log2(p) - q*log2(q);
    return h;
}

inline double Split::mutualClusteringInfo(const Split & other) const {
    // Example: n = 7
    //   A1     --**---  pA1   = 2/7
    //   B1     **--***  pB1   = 5/7
    //   A2     --****-  pA2   = 4/7
    //   B2     **----*  pB2   = 3/7
    //   A1,A2  --**---  pA1A2 = 2/7  iA1A2 = (2/7)*(log2(2/7) - log2(2/7) - log2(4/7)) = 0.23067283
    //   A1,B2  -------  pA1B2 = 0/7  iA1B2 = (0/7)*(log2(0/7) - log2(2/7) - log2(3/7)) = 0.0
    //   B1,A2  ----**-  pB1A2 = 2/7  iB1A2 = (2/7)*(log2(2/7) - log2(5/7) - log2(4/7)) = -0.14702091
    //   B1,B2  **----*  pB1B2 = 3/7  iB1B2 = (3/7)*(log2(3/7) - log2(5/7) - log2(3/7)) = 0.20804007
    //   Icl = 0.23067283 + 0.0 + (-0.14702091) + 0.20804007 = 0.29169199
    //
    // Confirm using TreeDist R package:
    //   library(TreeDist)
    //   tree1 <- ape::read.tree(text = '(A,B,(C,D),E,F,G);')
    //   tree2 <- ape::read.tree(text = '(A,B,(C,D,E,F),G);')
    //   ClusteringInfoDistance(tree1, tree2, reportMatching = TRUE)
    //   [1] 1.264965
    //   attr(,"matching")
    //   [1] 1
    //   attr(,"matchedSplits")
    //   [1] "C D | A B E F G => C D E F | A B G"
    //   attr(,"matchedScores")
    //   [1] 0.291692
    //   attr(,"pairScores")
    //            [,1]
    //   [1,] 0.291692

    assert(_nleaves > 0);
    assert(_nleaves == other.getSize());

    // Define each of the four possible taxon clusters
    const Split & A1 = *this;
    Split B1 = *this;
    B1.invertBits();

    const Split & A2 = other;
    Split B2 = other;
    B2.invertBits();

    // cerr << "mutualClusteringInfo:" << endl;
    // cerr << "  this:  " << this->createPatternRepresentation() << endl;
    // cerr << "  other: " << other.createPatternRepresentation() << endl;
    // cerr << "  A1: " << A1.createPatternRepresentation() << endl;
    // cerr << "  B1: " << B1.createPatternRepresentation() << endl;
    // cerr << "  A2: " << A2.createPatternRepresentation() << endl;
    // cerr << "  B2: " << B2.createPatternRepresentation() << endl;

    // Create intersections
    Split A1A2 = A1;
    Split A1B2 = A1;
    Split B1A2 = B1;
    Split B1B2 = B1;
    A1A2.bitwiseAnd(A2);
    A1B2.bitwiseAnd(B2);
    B1A2.bitwiseAnd(A2);
    B1B2.bitwiseAnd(B2);

    // cerr << "  A1A2: " << A1A2.createPatternRepresentation() << endl;
    // cerr << "  A1B2: " << A1B2.createPatternRepresentation() << endl;
    // cerr << "  B1A2: " << B1A2.createPatternRepresentation() << endl;
    // cerr << "  B1B2: " << B1B2.createPatternRepresentation() << endl;

    auto n = static_cast<double>(_nleaves);
    double pA1 = static_cast<double>(A1.getNumBitsSet())/n;
    double pB1 = static_cast<double>(B1.getNumBitsSet())/n;
    double pA2 = static_cast<double>(A2.getNumBitsSet())/n;
    double pB2 = static_cast<double>(B2.getNumBitsSet())/n;
    double pA1A2 = static_cast<double>(A1A2.getNumBitsSet())/n;
    double pA1B2 = static_cast<double>(A1B2.getNumBitsSet())/n;
    double pB1A2 = static_cast<double>(B1A2.getNumBitsSet())/n;
    double pB1B2 = static_cast<double>(B1B2.getNumBitsSet())/n;

    // Compute mutual clustering information (MCI)
    double Icl = 0.0;
    if (pA1A2 > 0.0)
        Icl += pA1A2*(log2(pA1A2) - log2(pA1) - log2(pA2));
    if (pA1B2 > 0.0)
        Icl += pA1B2*(log2(pA1B2) - log2(pA1) - log2(pB2));
    if (pB1A2 > 0.0)
        Icl += pB1A2*(log2(pB1A2) - log2(pB1) - log2(pA2));
    if (pB1B2 > 0.0)
        Icl += pB1B2*(log2(pB1B2) - log2(pB1) - log2(pB2));
    return Icl;
}

inline void Split::setBitsPerUnit(unsigned nbits_per_unit) {
    _bits_per_unit = nbits_per_unit;
    _nleaves = 0;
    unsigned nunits = 0;
    _bits.resize(nunits);
    clear();
}

inline unsigned Split::getBitsPerUnit() const {
    return static_cast<unsigned>(_bits_per_unit);
}

inline void Split::resize(unsigned nleaves) {
    _nleaves = nleaves;
    unsigned nunits = 1 + ((nleaves - 1)/_bits_per_unit);
    _bits.resize(nunits);
    clear();
}

inline void Split::setEdgeLen(double edgelen) {
    _edgelen = edgelen;
}

inline double Split::getEdgeLen() const {
    return _edgelen;
}

inline void Split::setBitAt(unsigned leaf_index) {
    assert(leaf_index < _nleaves);
    unsigned unit_index = leaf_index/_bits_per_unit;
    unsigned bit_index = leaf_index - unit_index*_bits_per_unit;
    auto bit_to_set = static_cast<split_unit_t>(1) << bit_index;
    _bits[unit_index] |= bit_to_set;
}

inline void Split::setBitsAt(vector<unsigned> leaves) {
    for (unsigned leaf_index : leaves) {
        setBitAt(leaf_index);
    }
}


inline bool Split::isBitSetAt(unsigned leaf_index) const {
    assert(leaf_index < _nleaves);
    unsigned unit_index = leaf_index/_bits_per_unit;
    unsigned bit_index = leaf_index - unit_index*_bits_per_unit;
    auto bit_to_set = static_cast<split_unit_t>(1) << bit_index;
    return static_cast<bool>(_bits[unit_index] & bit_to_set);
}

inline void Split::addSplit(const Split & other) {
    auto nunits = static_cast<unsigned>(_bits.size());
    assert(nunits == other._bits.size());
    for (unsigned i = 0; i < nunits; ++i) {
        _bits[i] |= other._bits[i];
    }
}

inline unsigned Split::getNumBitsSet() const {
    unsigned n = 0;
    unsigned ntax_added = 0;
    for (unsigned long b : _bits) {
        for (unsigned j = 0; j < _bits_per_unit; ++j) {
            auto bitmask = (static_cast<split_unit_t>(1) << j);
            bool bit_is_set = ((b & bitmask) > static_cast<split_unit_t>(0));
            if (bit_is_set)
                n++;
            if (++ntax_added == _nleaves)
                break;
        }
    }
    return n;
}

inline void Split::bitwiseAnd(const Split & other) {
    const split_t & other_bits = other.getBits();
    assert(_bits.size() == other_bits.size());
    for (unsigned i = 0; i < _bits.size(); ++i) {
        _bits[i] &= other_bits[i];
    }
}

inline unsigned Split::findFirstSetBit() const {
    unsigned bit_index = 0;
    for (auto & u : _bits) {
        if (u > 0) {
            while ((u & (static_cast<split_unit_t>(1) << bit_index)) == 0)
                bit_index++;
            break;
        }
    }
    return bit_index;
}

inline void Split::invertBits() {
    for (auto & u : _bits) {
        u = ~u;
    }
}

#if 1
inline bool Split::compatibleWith(const Split & other) const {
    // Suppose units each contain 4 bits.
    // Example of incompatible splits:
    //       a = |0000|1111|
    //       b = |1111|0110|
    //      ab = |0000|0110| > 0 and equals neither a nor b --> false
    // Example of incompatible splits:
    //       a = |0001|1000|
    //       b = |0001|0100|
    //      ab = |0001|0000| > 0 and equals neither a nor b --> false
    // Example of incompatible splits:
    //       a = |0001|1100|
    //       b = |0011|0000|
    //      ab = |0001|0000| > 0 and equals neither a nor b --> false
    // Example of compatible splits:
    //       a = |0001|1000|
    //       b = |0010|0100|
    //      ab = |0000|0000| equals 0
    // Example of compatible splits:
    //       a = |0001|1100|
    //       b = |0001|0100|
    //      ab = |0001|0100| > 0 and equals b --> true
    // Example of compatible splits:
    //       a = |0001|0000|
    //       b = |1011|0000|
    //      ab = |0001|0000| > 0 and equals a --> true
    const split_t & other_bits = other._bits;
    assert(_bits.size() == other_bits.size());
    bool ab_equals_a = true;
    bool ab_equals_b = true;
    bool a_and_b_zero = true;
    for (unsigned i = 0; i < _bits.size(); ++i) {
        split_unit_t a  = _bits[i];
        split_unit_t b  = other_bits[i];
        split_unit_t a_and_b = (a & b);
        split_unit_t a_or_b = (a | b);
        if (a_or_b) {
            if (a_and_b)
                a_and_b_zero = false;
            if (a_and_b != a)
                ab_equals_a = false;
            if (a_and_b != b)
                ab_equals_b = false;
        }
    }
    bool is_compatible = a_and_b_zero || (ab_equals_a || ab_equals_b);
    return is_compatible;
}
#else
    // This version represents a bug only discovered 2025-10-28
    // If set bits straddle two units, then the separate units can be compatible
    // separately but not jointly!
    // Suppose units each contain 4 bits
    //     a = |0001|1000|
    //     b = |0001|0100|
    // unit 0:
    //          a = |1000|
    //          b = |0100|
    //      a & b = |0000| = false
    //   equals_a = false
    //   equals_b = false
    //   equals_a || equals_b = false
    //   !(equals_a || equals_b) = true
    //   is_compatible = true
    // unit 1:
    //          a = |0001|
    //          b = |0001|
    //      a & b = |0001| = true
    //   equals_a = true
    //   equals_b = true
    //   equals_a || equals_b = true
    //   !(equals_a || equals_b) = false
    //   is_compatible = true
    // So, the splits are in reality not compatible, but this version
    // says they are because each unit is comptible!
    //
    inline bool Split::compatibleWith(const Split & other) const {
    const split_t & other_bits = other._bits;
    assert(_bits.size() == other_bits.size());
    bool is_compatible = true;
    for (unsigned i = 0; i < _bits.size(); ++i) {
        split_unit_t a       = _bits[i];
        split_unit_t b       = other_bits[i];
        split_unit_t a_and_b = (a & b);
        bool equals_a   = (a_and_b == a);
        bool equals_b   = (a_and_b == b);
        if (a_and_b && !(equals_a || equals_b)) {
            // A failure of any unit to be compatible makes the entire split incompatible
            is_compatible = false;
            break;
        }
    }
    return is_compatible;
}
#endif

inline bool Split::subsumedIn(const Split & other) const {
    // Example: this is NOT subsumed in other
    // 8421
    // **-* 13 this <----------+
    // **-- 12 other           | not equal
    // **-- 12 this & other <--+
    //
    // Example: this IS subsumed in other
    // 8421
    // **-- 12 this <----------+
    // **-* 13 other           | equal
    // **-- 12 this & other <--+
    //
    // Example: this is compatible with, but not subsumed in, other
    // 8421
    // **-- 12 this <----------+
    // --**  3 other           | not equal
    // ----  0 this & other <--+

    const split_t & other_bits = other.getBits();
    //temporary!
    if (_bits.size() != other_bits.size()) {
        cerr << "oops" << endl;
    }
    assert(_bits.size() == other_bits.size());
    for (unsigned i = 0; i < _bits.size(); ++i) {
        if ((_bits[i] & other_bits[i]) != _bits[i])
            return false;
    }
    return true;
}

inline string Split::createPatternRepresentation(bool show_edge_length) const {
    string s;
    if (show_edge_length)
        s += str(format("%.9f: ") % _edgelen);
    unsigned ntax_added = 0;
    for (unsigned long b : _bits) {
        for (unsigned j = 0; j < _bits_per_unit; ++j) {
            split_unit_t bitmask = (static_cast<split_unit_t>(1) << j);
            bool bit_is_set = ((b & bitmask) > static_cast<split_unit_t>(0));
            if (bit_is_set)
                s += '*';
            else
                s += '-';
            if (++ntax_added == _nleaves)
                break;
            }
        }
    return s;
    }

}
