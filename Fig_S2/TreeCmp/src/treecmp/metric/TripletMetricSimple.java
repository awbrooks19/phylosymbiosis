/** This file is part of TreeCmp, a tool for comparing phylogenetic trees
    using the Matching Split distance and other metrics.
    Copyright (C) 2011,  Damian Bogdanowicz

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. */

package treecmp.metric;

import pal.misc.IdGroup;
import pal.tree.Tree;
import pal.tree.TreeUtils;
import treecmp.common.TreeCmpUtils;

public class TripletMetricSimple extends BaseMetric implements Metric {
  

    public double getDistance(Tree t1, Tree t2) {
        return getDistForArbitrary(t1, t2);
    }

   
    /**
     * Calculate triplet distance using simple O(n^3) algorithm
     * enumeration all possible triplets
     * @param t1
     * @param t2
     * @return
     */
    public double getDistForArbitrary(Tree t1, Tree t2) {

        IdGroup id1 = TreeUtils.getLeafIdGroup(t1);
        int[][] nsMatrix1 = TreeCmpUtils.calcNodalSplittedMatrix(t1, null);
        int[][] nsMatrix2 = TreeCmpUtils.calcNodalSplittedMatrix(t2, id1);
        long unResolved_T1 = 0;
        long unResolved_T2 = 0;
        long unResolved_Common = 0;
        long resolved_Common = 0;
        long sum = 0;
        int type1, type2;

        int leafNum = t1.getExternalNodeCount();
        for (int i = 0; i < leafNum; i++) {
            for (int j = i + 1; j < leafNum; j++) {
                for (int k = j + 1; k < leafNum; k++) {
                    type1 = getTripletType(i, j, k, nsMatrix1);
                    type2 = getTripletType(i, j, k, nsMatrix2);
                    if (type1 == -1) {
                        unResolved_T1++;
                    }
                    if (type1 == -1) {
                        unResolved_T2++;
                    }
                    if (type1 == type2) {
                        if (type1 == -1) {
                            unResolved_Common++;
                        } else {
                            resolved_Common++;
                        }
                    }
                    sum++;
                }
            }
        }

        long dist = sum - unResolved_Common - resolved_Common;
        return (double) dist;
    }

    /* Retruns the type of a triplet, i.e. the index for which LCA with the others is
     * the closest to the root. Returns -1 if the triplet is unresolved. For example:
     * - type for (i,(j,k)) = i,
     * - type for (j,(i,k)) = j,
     * - type for (k,(i,j)) = k,
     * - type for (i,j,k) = -1.
     */
    private int getTripletType(int i, int j, int k, int nsMatrix[][]) {

        if (nsMatrix[j][i] > nsMatrix[j][k]) {
            return i;
        } else if (nsMatrix[j][i] < nsMatrix[j][k]) {
            return k;
        } else {
            //two situation are possible: type = j or (-1) - triplet is unresoved
            if (nsMatrix[i][j] > nsMatrix[i][k]) {
                return j;
            }
            //nsMatrix[i][j] < nsMatrix[i][k] could not happen
            //only nsMatrix[i][j] = nsMatrix[i][k] is possible here
            //hence, triplet is unresolved
        }
        return -1;
    }
}
