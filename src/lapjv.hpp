#pragma once

using namespace std;
using namespace boost;

namespace op {
    // Downloaded lap.h and lap.cpp on 2025-09-28 from https://github.com/yongyanghz/LAPJV-algorithm-c
    // Combined by Paul O. Lewis into a single file lapjv.hpp
    /************************************************************************
    *
    *  lap.cpp
       version 1.0 - 4 September 1996
       author: Roy Jonker @ MagicLogic Optimization Inc.
       e-mail: roy_jonker@magiclogic.com

       Code for Linear Assignment Problem, according to

       "A Shortest Augmenting Path Algorithm for Dense and Sparse Linear
        Assignment Problems," Computing 38, 325-340, 1987

       by

       R. Jonker and A. Volgenant, University of Amsterdam.

    *
       CHANGED 2016-05-13 by Yong Yang(yongyanglink@gmail.com) in column reduction part according to
       matlab version of LAPJV algorithm(Copyright (c) 2010, Yi Cao All rights reserved)--
       https://www.mathworks.com/matlabcentral/fileexchange/26836-lapjv-jonker-volgenant-algorithm-for-linear-assignment-problem-v3-0:
    *
    *************************************************************************/

    class LAPJV {
    public:
        LAPJV(unsigned dim);

        void   assignCost(unsigned i, unsigned j, double cost);
        double lap();
        void   debugLAP(string title) const;
        void   getOptimalPairings(vector<unsigned> & optpairings) const;

        unsigned                 _dim;
        vector<int>              _rowsol;
        vector<int>              _colsol;
        vector<double>           _u;
        vector<double>           _v;
        vector< vector<double> > _assigncost;

    private:
        vector<int>              _pred;            // row-predecessor of column in augmenting/alternating path
        vector<int>              _freeunassigned;  // list of unassigned rows.
        vector<int>              _collist;         // list of columns to be scanned in various ways
        vector<int>              _matches;         // counts how many times a row could be assigned
        vector<double>           _d;               // 'cost-distance' in augmenting path calculation

    };

    inline LAPJV::LAPJV(unsigned dim) : _dim(dim) {
        _assigncost.resize(_dim);
        for (unsigned i = 0; i < _dim; i++) {
            _assigncost[i].resize(_dim);
        }
        _rowsol.resize(_dim);
        _colsol.resize(_dim);
        _u.resize(_dim);
        _v.resize(_dim);
        _pred.resize(_dim);
        _freeunassigned.resize(_dim);
        _collist.resize(_dim);
        _matches.resize(_dim);
        _d.resize(_dim);
    }

    inline void LAPJV::assignCost(unsigned i, unsigned j, double cost) {
        assert(i < _dim);
        assert(j < _dim);
        _assigncost[i][j] = cost;
    }

    inline void LAPJV::debugLAP(string title) const {
#if defined(DEBUGGING_LAPJV)
        cerr << "\n############## " << title << " ##############\n" << endl;

        // headers
        cerr << str(format("%12s  ") % " ");
        for (int j = 0; j < _dim; j++) {
            cerr << str(format(" %12d") % j);
        }
        cerr << str(format("   %12s %12s %15s\n") % "u" % "rowsol" % "freeunassigned");

        // dashes
        cerr << str(format("%12s +") %  " ");
        for (int j = 0; j < _dim; j++) {
            cerr << str(format(" %12s") % "------------");
        }
        cerr << str(format(" + %12s %12s %15s\n") % " " % " " % " ");

        // matrix
        for (int i = 0; i < _dim; i++) {
            cerr << str(format("%12d |") % i);
            for (int j = 0; j < _dim; j++) {
                cerr << str(format(" %12.0f") % _assigncost[i][j]);
            }
            cerr << str(format(" | %12.0f %12d %15d\n") % _u[i] % _rowsol[i] % _freeunassigned[i]);
        }

        // dashes
        cerr << str(format("%12s +") % " ");
        for (int j = 0; j < _dim; j++) {
            cerr << str(format(" %12s") % "------------");
        }
        cerr << str(format(" + %12s %12s %15s\n") % " " % " " % " ");

        // v row
        cerr << str(format("%12s  ") % "v");
        for (int j = 0; j < _dim; j++) {
            cerr << str(format(" %12.0f") % _v[j]);
        }
        cerr << endl;

        // colsol row
        cerr << str(format("%12s  ") % "colsol");
        for (int j = 0; j < _dim; j++) {
            cerr << str(format(" %12d") % _colsol[j]);
        }
        cerr << endl;

        // matches row
        cerr << str(format("%12s  ") % "matches");
        for (int j = 0; j < _dim; j++) {
            cerr << str(format(" %12d") % _matches[j]);
        }
        cerr << endl;
#endif
    }

    /*This function is the jv shortest augmenting path algorithm to solve the assignment problem*/
    inline double LAPJV::lap() {

        // input:
        // dim        - problem size
        // assigncost - cost matrix

        // output:
        // rowsol     - column assigned to row in solution
        // colsol     - row assigned to column in solution
        // u          - dual variables, row reduction numbers
        // v          - dual variables, column reduction numbers
        bool unassignedfound;
        int i;
        int imin;
        int numfree = 0;
        int prvnumfree;
        int f;
        int i0;
        int k;
        int freerow;
        int j;
        int j1;
        int j2;
        int endofpath;
        int last;
        int low;
        int up;
        double min;
        double h;
        double umin;
        double usubmin;
        double v2;

        // init how many times a row will be assigned in the column reduction.
        _matches.assign(_dim, 0);

        // COLUMN REDUCTION (reverse order gives better results)
        for (j = _dim; j--;) {
            // find minimum cost over rows.
            min = _assigncost[0][j];
            imin = 0;
            for (i = 1; i < _dim; i++)
                if (_assigncost[i][j] < min) {
                    min = _assigncost[i][j];
                    imin = i;
                }

            // v holds minimum cost in each column
            _v[j] = min;

            if (++_matches[imin] == 1) {
                // init assignment if minimum row assigned for first time.
                _rowsol[imin] = j;
                _colsol[j] = imin;
            } else if (_v[j] < _v[_rowsol[imin]]) {
                int j1 = _rowsol[imin];
                _rowsol[imin] = j;
                _colsol[j] = imin;
                _colsol[j1] = -1;
            } else
                _colsol[j] = -1;  // row already assigned, column not assigned.
        }

        debugLAP("After column reduction");

        // REDUCTION TRANSFER
        for (i = 0; i < _dim; i++) {
            if (_matches[i] == 0) {
                // fill list of unassigned 'free' rows.
                _freeunassigned[numfree++] = i;
            }
            else if (_matches[i] == 1) { // transfer reduction from rows that are assigned once.
                j1 = _rowsol[i];
                min = numeric_limits<double>::max();
                for (j = 0; j < _dim; j++)
                    if (j != j1)
                        if (_assigncost[i][j] - _v[j] < min) min = _assigncost[i][j] - _v[j];
                _v[j1] = _v[j1] - min;
            }
        }

        debugLAP("After reduction transfer");

        //   AUGMENTING ROW REDUCTION
        int loopcnt = 0;  // do-loop to be done twice.
        do {
            loopcnt++;

            //     scan all free rows.
            //     in some cases, a free row may be replaced with another one to be scanned next.
            k = 0;
            prvnumfree = numfree;
            numfree = 0;  // start list of rows still free after augmenting row reduction.
            while (k < prvnumfree) {
                i = _freeunassigned[k];
                k++;

                //       find minimum and second minimum reduced cost over columns.
                umin = _assigncost[i][0] - _v[0];
                j1 = 0;
                usubmin = numeric_limits<double>::max();
                for (j = 1; j < _dim; j++) {
                    h = _assigncost[i][j] - _v[j];
                    if (h < usubmin) {  // POL added {
                        if (h >= umin) {
                            usubmin = h;
                            j2 = j;
                        } else {
                            usubmin = umin;
                            umin = h;
                            j2 = j1;
                            j1 = j;
                        }
                    }   // POL added }
                }

                i0 = _colsol[j1];
                if (umin < usubmin)
                    //         change the reduction of the minimum column to increase the minimum
                        //         reduced cost in the row to the subminimum.
                            _v[j1] = _v[j1] - (usubmin - umin);
                else              // minimum and subminimum equal.
                    if (i0 > -1)  // minimum column j1 is assigned.
                    {
                        //           swap columns j1 and j2, as j2 may be unassigned.
                        j1 = j2;
                        i0 = _colsol[j2];
                    }

                //       (re-)assign i to j1, possibly de-assigning an i0.
                _rowsol[i] = j1;
                _colsol[j1] = i;

                if (i0 > -1) {
                    // POL added {
                    // minimum column j1 assigned earlier.
                    if (umin < usubmin)
                        //           put in current k, and go back to that k.
                            //           continue augmenting path i - j1 with i0.
                                _freeunassigned[--k] = i0;
                    else
                        //           no further augmenting reduction possible.
                            //           store i0 in list of free rows for next phase.
                                _freeunassigned[numfree++] = i0;
                }   // POL added }
            }

        } while (loopcnt < 2);  // repeat once.

        debugLAP("After augmenting row reduction");

        // AUGMENT SOLUTION for each free row.
        for (f = 0; f < numfree; f++) {
            freerow = _freeunassigned[f];  // start row of augmenting path.

            // Dijkstra shortest path algorithm.
            // runs until unassigned column added to shortest path tree.
            for (j = _dim; j--;) {
                _d[j] = _assigncost[freerow][j] - _v[j];
                _pred[j] = freerow;
                _collist[j] = j;  // init column list.
            }

            low = 0;  // columns in 0..low-1 are ready, now none.
            up = 0;   // columns in low..up-1 are to be scanned for current minimum, now none.
            // columns in up..dim-1 are to be considered later to find new minimum,
            // at this stage the list simply contains all columns
            unassignedfound = false;
            do {
                if (up == low)  // no more columns to be scanned for current minimum.
                {
                    last = low - 1;

                    // scan columns for up..dim-1 to find all indices for which new minimum occurs.
                    // store these indices between low..up-1 (increasing up).
                    min = _d[_collist[up++]];
                    for (k = up; k < _dim; k++) {
                        j = _collist[k];
                        h = _d[j];
                        if (h <= min) {
                            if (h < min)  // new minimum.
                            {
                                up = low;  // restart list at index low.
                                min = h;
                            }
                            // new index with same minimum, put on undex up, and extend list.
                            _collist[k] = _collist[up];
                            _collist[up++] = j;
                        }
                    }
                    // check if any of the minimum columns happens to be unassigned.
                    // if so, we have an augmenting path right away.
                    for (k = low; k < up; k++)
                        if (_colsol[_collist[k]] < 0) {
                            endofpath = _collist[k];
                            unassignedfound = true;
                            break;
                        }
                }

                if (!unassignedfound) {
                    // update 'distances' between freerow and all unscanned columns, via next scanned
                    // column.
                    j1 = _collist[low];
                    low++;
                    i = _colsol[j1];
                    h = _assigncost[i][j1] - _v[j1] - min;

                    for (k = up; k < _dim; k++) {
                        j = _collist[k];
                        v2 = _assigncost[i][j] - _v[j] - h;
                        if (v2 < _d[j]) {
                            _pred[j] = i;
                            if (v2 == min) {
                                //POL added { // new column found at same minimum value
                                if (_colsol[j] < 0) {
                                    // if unassigned, shortest augmenting path is complete.
                                    endofpath = j;
                                    unassignedfound = true;
                                    break;
                                }
                                // else add to list to be scanned right away.
                                else {
                                    _collist[k] = _collist[up];
                                    _collist[up++] = j;
                                }
                            }   // POL added }
                            _d[j] = v2;
                        }
                    }
                }
            } while (!unassignedfound);

            // update column prices.
            for (k = last + 1; k--;) {
                j1 = _collist[k];
                _v[j1] = _v[j1] + _d[j1] - min;
            }

            // reset row and column assignments along the alternating path.
            do {
                i = _pred[endofpath];
                _colsol[endofpath] = i;
                j1 = endofpath;
                endofpath = _rowsol[i];
                _rowsol[i] = j1;
            } while (i != freerow);
        }

        debugLAP("After augmenting solution for each free row");

        // calculate optimal cost.
        double lapcost = 0;
        //  for (i = 0; i < dim; i++)
        for (i = _dim; i--;) {
            j = _rowsol[i];
            _u[i] = _assigncost[i][j] - _v[j];
            lapcost += _assigncost[i][j];
        }

        debugLAP("After calculating total cost");

        return lapcost;
    }

    inline void LAPJV::getOptimalPairings(vector<unsigned> & optpairings) const {
        for (unsigned i = 0; i < _dim; ++i) {
            optpairings.push_back(_rowsol[i]);
        }
    }


}