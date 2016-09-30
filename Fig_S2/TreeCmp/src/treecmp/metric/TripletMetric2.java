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

import pal.tree.Node;
import pal.tree.Tree;
import pal.tree.TreeUtils;
import treecmp.common.ClustIntersectInfoMatrix;
import treecmp.common.TreeCmpUtils;



public class TripletMetric2 extends BaseMetric implements Metric {

    public double getDistance(Tree t1, Tree t2) {

        long n = t1.getExternalNodeCount();
        Node [] nodesT1 = TreeCmpUtils.getNodesInPostOrder(t1);
        Node [] nodesT2 = TreeCmpUtils.getNodesInPostOrder(t2);
        ClustIntersectInfoMatrix cIM =TreeCmpUtils.calcClustIntersectMatrix(t1, t2, TreeUtils.getLeafIdGroup(t1));
       // long Rt1 = TreeCmpUtils.calcResolvedTriplets(t1, nodesT1, cIM.cSize1);
        long Rt2 = TreeCmpUtils.calcResolvedTriplets(t2, nodesT2, cIM.cSize2);
        long St1t2 = TreeCmpUtils.calcResolvedAndEqualTriplets(cIM, nodesT1, nodesT2);
        long R1t1t2 = TreeCmpUtils.calcResolvedOnlyInT1(cIM, nodesT1, nodesT2);
       // long npo3 = (n*(n-1)*(n-2))/6;

        long dist = Rt2 - St1t2 + R1t1t2;
       
        return (double)dist;

    }

    
}
