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
import pal.tree.TreeDistanceMatrix;
import pal.tree.TreeUtils;
import treecmp.common.TreeCmpUtils;

public class NodalL2Metric extends BaseMetric implements Metric {

    //This seems to be faster than old implementation
    public double getDistance(Tree t1, Tree t2) {
        double dist, diff;

        IdGroup id1 = TreeUtils.getLeafIdGroup(t1);
        int[][] nsMatrix1 = TreeCmpUtils.calcNodalSplittedMatrix(t1, null);
        int[][] nsMatrix2 = TreeCmpUtils.calcNodalSplittedMatrix(t2, id1);

        dist = 0.0;
        for (int i = 0; i < id1.getIdCount(); i++) {
            for (int j = i + 1; j < id1.getIdCount(); j++) {
                diff = nsMatrix1[i][j] + nsMatrix1[j][i] - nsMatrix2[i][j] - nsMatrix2[j][i];
                dist += diff * diff;
            }
        }
        return Math.sqrt(dist);
    }

    public double getDistanceOld(Tree t1, Tree t2) {
        double dist, diff;


        IdGroup id1 = TreeUtils.getLeafIdGroup(t1);
        TreeDistanceMatrix tr1 = new TreeDistanceMatrix(t1, true, 0);
        TreeDistanceMatrix tr2 = new TreeDistanceMatrix(t2, id1, true, 0);

        dist = 0.0;
        for (int i = 0; i < id1.getIdCount(); i++) {
            for (int j = i + 1; j < id1.getIdCount(); j++) {
                diff = tr1.getDistance(i, j) - tr2.getDistance(i, j);
                dist += diff * diff;
            }
        }
        return Math.sqrt(dist);
    }
}
