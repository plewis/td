#pragma once    

#include <set>
#include <map>
#include <vector>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <boost/format.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include "split.hpp"
#include "tree_manip.hpp"
#include "xtreedist.hpp"

#include "ncl/nxsmultiformat.h"

using namespace std;

extern treedist::OutputManager om;

namespace treedist {

	class TreeSummary {
	public:
		TreeSummary();
		~TreeSummary();

        void                        readTreefile(const string filename, unsigned skip, unsigned ref_tree, bool store_all, bool debug);
        void                        calcDistances(vector<double> & kfvect, vector<unsigned> & rfvect, bool quiet, bool debug);
		typename Tree::SharedPtr    getTree(unsigned index);
		string                      getNewick(unsigned index);
		string                      getTreeName(unsigned index);
		bool                        isRooted(unsigned index);
		void                        clear();

	private:

		vector<string>          _newicks;
		vector<string>          _tree_names;
        vector<bool>            _is_rooted;
        vector<Split::treeid_t> _splitset_vect;
        
        string                  _ref_newick;
        string                  _ref_treename;
        bool                    _ref_is_rooted;
        Split::treeid_t         _ref_splitset;
        
        vector<string>          _taxon_names;
        map<string, unsigned>   _master_taxon_map;
        map<unsigned, unsigned> _taxon_map;
        
        // If _splitset_vect.size() < _nprogress, no progress shown
        // If _splitset_vect.size() >= _nprogress, _nprogress lines will
        //                             be displayed showing progress
        const unsigned          _nprogress;

	public:

		typedef std::shared_ptr< TreeSummary > SharedPtr;
	};

	inline TreeSummary::TreeSummary() : _nprogress(100) {
		//cout << "Constructing a TreeSummary" << endl;
	}

	inline TreeSummary::~TreeSummary() {
		//cout << "Destroying a TreeSummary" << endl;
	}

	inline Tree::SharedPtr TreeSummary::getTree(unsigned index) {
		if (index >= _newicks.size())
			throw XTreeDist("getTree called with index >= number of stored trees");
        assert(_newicks.size() == _is_rooted.size());

		TreeManip tm;

		// build the tree
		tm.buildFromNewick(_newicks[index], _is_rooted[index], false);

		return tm.getTree();
	}

	inline string TreeSummary::getNewick(unsigned index) {
		if (index >= _newicks.size())
			throw XTreeDist("getNewick called with index >= number of stored trees");

		return _newicks[index];
	}

	inline string TreeSummary::getTreeName(unsigned index) {
		if (index >= _tree_names.size())
			throw XTreeDist("getTreeName called with index >= number of stored trees");

		return _tree_names[index];
	}

	inline bool TreeSummary::isRooted(unsigned index) {
		if (index >= _is_rooted.size())
			throw XTreeDist("isRooted called with index >= number of stored trees");

		return _is_rooted[index];
	}

	inline void TreeSummary::clear() {
		_newicks.clear();
		_tree_names.clear();
		_is_rooted.clear();
	}

	inline void TreeSummary::readTreefile(const string filename, unsigned skip, unsigned ref_tree, bool store_all, bool debug) {
        // If store_all is true, then
        //   - split sets from all trees (after skip) will be appended to _splitset_vect
        //   - if ref_tree is > 0, splits from tree with index ref_tree - 1 (of the trees
        //     not skipped) will be saved to _ref_splitset and not to _splitset_vect
        // If store_all is false,
        //   - only the splits set from the reference tree will be saved in _ref_splitset
        //   - _splitsset_vect will not be touched

        // If not storing all trees, then ref_tree must be a valid
        // index of the tree to store
        assert (ref_tree > 0 || store_all);
        
		TreeManip tm;
		Split::treeid_t splitset;

		// See http://phylo.bio.ku.edu/ncldocs/v2.1/funcdocs/index.html for NCL documentation

		//MultiFormatReader nexusReader(-1, NxsReader::WARNINGS_TO_STDERR);
		MultiFormatReader nexusReader(-1, NxsReader::IGNORE_WARNINGS);
  
        // Both of these needed to suppress "storing read block" messages
        // see NxsReader::statusMessage in nxsreader.cpp
        nexusReader.SetAlwaysReportStatusMessages(false);
        nexusReader.SetWarningOutputLevel(NxsReader::SUPPRESS_WARNINGS_LEVEL);
  
		try {
			nexusReader.ReadFilepath(filename.c_str(), MultiFormatReader::NEXUS_FORMAT);
		}
		catch (...) {
			nexusReader.DeleteBlocksFromFactories();
			throw;
		}

		int num_taxa_blocks = nexusReader.GetNumTaxaBlocks();
        if (num_taxa_blocks > 1) {
            throw XTreeDist("This program allows only one taxa block per tree file");
        }
        
        NxsTaxaBlock* taxaBlock = nexusReader.GetTaxaBlock(0);
        string taxaBlockTitle = taxaBlock->GetTitle();
        
        if (_taxon_names.empty()) {
            // This is the first taxa block encountered:
            // - store taxon labels in order in _taxon_names
            // - build _master_taxon_map
            _taxon_map.clear();
            if (debug)
                om.outputConsole("\nProcessing master taxa block:\n");
            for (unsigned t = 0; t < taxaBlock->GetNTax(); t++) {
                string nm = taxaBlock->GetTaxonLabel(t);
                _taxon_names.push_back(nm);
                if (debug)
                    om.outputConsole(format("  taxon %d is \"%s\"\n") % t % nm);
                _master_taxon_map[nm] = t;
                _taxon_map[t] = t;
            }
        }
        else {
            // This is NOT the first taxa block encountered:
            // - build _taxon_map to translate node numbers in this taxa block
            //   to node numbers from the master taxa block
            _taxon_map.clear();
            if (debug)
                om.outputConsole("\nProcessing subsequent taxa block:\n");
            for (unsigned t = 0; t < taxaBlock->GetNTax(); t++) {
                string nm = taxaBlock->GetTaxonLabel(t);
                if (debug)
                    om.outputConsole(format("  taxon %d is \"%s\" --> %d\n") % t % nm % _master_taxon_map.at(nm));
                _taxon_map[t] = _master_taxon_map.at(nm);
            }
        }

        const unsigned num_trees_blocks = nexusReader.GetNumTreesBlocks(taxaBlock);
        if (num_trees_blocks > 1) {
            throw XTreeDist("This program allows only one trees block per tree file");
        }

        const NxsTreesBlock* treesBlock = nexusReader.GetTreesBlock(taxaBlock, 0);
        unsigned ntrees = treesBlock->GetNumTrees();
        if (skip >= ntrees) {
            throw XTreeDist(format("Skipping more trees than are in the file: skip = %d, ntrees = %d") % skip % ntrees);
        }
        
        if (store_all) {
            for (unsigned t = skip; t < ntrees; ++t) {
                unsigned tindex = t - skip;
                const NxsFullTreeDescription& d = treesBlock->GetFullTreeDescription(t);
                
                // store the tree name
                string tree_name = d.GetName();
                
                // Assuming that the newick string has been processed to convert
                // taxon names to 1-offset numbers in the order specified in the
                // corresponding taxa block
                assert(d.IsProcessed());

                // store the newick tree description
                bool is_rooted = d.IsRooted();
                string newick = d.GetNewick();

                // build the tree
                tm.buildFromNewick(newick, is_rooted, false);
                
                if (ref_tree > 0 && tindex == ref_tree - 1) {
                    _ref_newick = newick;
                    _ref_treename = tree_name;
                    _ref_is_rooted = is_rooted;
                    _ref_splitset.clear();
                    tm.storeSplits(_ref_splitset, _taxon_map);
                }
                else {
                    _newicks.push_back(newick);
                    _tree_names.push_back(tree_name);
                    _is_rooted.push_back(is_rooted);

                    // store set of splits
                    splitset.clear();
                    tm.storeSplits(splitset, _taxon_map);
                    
                    _splitset_vect.push_back(splitset);
                }
            } // trees loop
        }
        else {
            // Only store reference tree
            if (ref_tree > ntrees - skip) {
                throw XTreeDist(format("Reference tree is beyond last non-skipped tree in file: skip = %d, ntrees = %d, reftree = %d") % skip % ntrees % ref_tree);
            }
            
            const NxsFullTreeDescription& d = treesBlock->GetFullTreeDescription(skip + ref_tree - 1);
            
            // store the newick tree description
            string tree_name = d.GetName();
            bool is_rooted = d.IsRooted();
            string newick = d.GetNewick();
            _ref_newick = newick;
            _ref_treename = tree_name;
            _ref_is_rooted = is_rooted;

            // build the tree
            tm.buildFromNewick(newick, is_rooted, false);

            // store set of splits
            _ref_splitset.clear();
            tm.storeSplits(_ref_splitset, _taxon_map);
        }

		// No longer any need to store raw data from nexus file
		nexusReader.DeleteBlocksFromFactories();
	}
    
    inline void TreeSummary::calcDistances(vector<double> & kfvect, vector<unsigned> & rfvect, bool quiet, bool debug) {
        kfvect.clear();
        rfvect.clear();

        // Map the splits in _ref_splitset
        map<Split, double> ref_split_map;
        for (auto & s : _ref_splitset) {
            ref_split_map[s] = s.getWeight();
        }
    
        map<Split, double> split_map;
        vector<Split> inboth(_ref_splitset.size());
        vector<Split> inone(_ref_splitset.size());
        
        unsigned n = (unsigned)_splitset_vect.size();

        // Decide whether to show progress
        // n = 11                         n = 10
        // w = 11/10 = 1.1 -> 1           w = 10/5 = 2
        //  i  i*w (i+1) % w == 0         i  i+1 (i+1) % w == 0
        //  0  0.0       *                0   1
        //  1  1.1       *                1   2     *
        //  2  2.2       *                2   3
        //  3  3.3       *                3   4     *
        //  4  4.4       *                4   5
        //  5  5.5       *                5   6     *
        //  6  6.6       *                6   7
        //  7  7.7       *                7   8     *
        //  8  8.8       *                8   9
        //  9  9.9       *                9  10     *
        // 10 11.0
        unsigned w = (n > _nprogress) ? (unsigned)(n/_nprogress) : 1;

        unsigned i = 0;
        for (auto & ss : _splitset_vect) {
            if (!quiet && (i+1) % w == 0) {
                om.outputConsole(format("Tree %d of %d (\"%s\"): ") % (i+1) % n % _tree_names[i]);
            }
            else if (debug) {
                om.outputConsole(format("\nTree %d of %d (\"%s\")\n") % (i+1) % n % _tree_names[i]);
            }
            
            // Map the splits in ss
            split_map.clear();
            for (auto & s : ss) {
                split_map[s] = s.getWeight();
            }
            
            // Get the intersection between _ref_splitset and ss
            auto intersection_end = set_intersection(_ref_splitset.begin(), _ref_splitset.end(), ss.begin(), ss.end(), inboth.begin());
                        
            // Get set symmetric difference
            auto symmdiff_end = set_symmetric_difference(_ref_splitset.begin(), _ref_splitset.end(), ss.begin(), ss.end(), inone.begin());

            double KFdist = 0.0;

            if (debug) {
                om.outputConsole("  Intersection of split sets:\n");
                om.outputConsole(format("  %20s %20s %20s %20s %s\n") % "ref" % "other" % "diff^2" % "KF" % "split");
            }

            for (auto it = inboth.begin(); it != intersection_end; ++it) {
                Split & s = *it;
                KFdist += pow(ref_split_map[s] - split_map[s], 2.0);
                
                if (debug) {
                    double ref = ref_split_map[s];
                    double other = split_map[s];
                    double diffsq = pow(ref - other, 2);
                    om.outputConsole(format("  %20.9f %20.9f %20.9f %20.9f %s\n") % ref % other % diffsq % KFdist % s.createPatternRepresentation());
                }
            }

            if (debug) {
                om.outputConsole("\n  Symmetric difference of split sets:\n");
                om.outputConsole(format("  %20s %20s %20s %20s %s\n") % "ref" % "other" % "diff^2" % "KF" % "split");
            }

            unsigned RFdist = (unsigned)std::distance(inone.begin(), symmdiff_end);
            for (auto it = inone.begin(); it != symmdiff_end; ++it) {
                Split & s = *it;
                KFdist += pow(s.getWeight(), 2.0);

                if (debug) {
                    bool s_in_ref = (ref_split_map.count(s) == 1);
                    bool s_in_other = (split_map.count(s) == 1);
                    assert(s_in_ref || s_in_other);
                    if (s_in_ref) {
                        double ref = ref_split_map[s];
                        double diffsq = pow(ref, 2);
                        om.outputConsole(format("  %20.9f %20s %20.9f %20.9f %s\n") % ref % "0" % diffsq % KFdist % s.createPatternRepresentation());
                    }
                    else {
                        double other = split_map[s];
                        double diffsq = pow(other, 2);
                        om.outputConsole(format("  %20s %20.9f %20.9f %20.9f %s\n") % "0" % other % diffsq % KFdist % s.createPatternRepresentation());
                    }
                }
            }

            if (!quiet && (i+1) % w == 0) {
                om.outputConsole(format("KF = %.9f, RF = %d\n") % KFdist % RFdist);
            }
            
            kfvect.push_back(KFdist);
            rfvect.push_back(RFdist);
            i++;
        }
    }

}
