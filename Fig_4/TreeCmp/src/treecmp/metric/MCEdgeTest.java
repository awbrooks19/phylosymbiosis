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
import pal.tree.NodeUtils;
import pal.tree.Tree;
import treecmp.common.TreeCmpUtils;

public class MCEdgeTest extends BaseMetric implements Metric {
  
    public double getDistance(Tree t1, Tree t2) {

        int i, j;
        int metric = 0;
        
        int t1IntNum = t1.getInternalNodeCount();
        int t2IntNum = t2.getInternalNodeCount();
        int n = t1.getExternalNodeCount();
        int c1,h1,c2,h2;
        int removedNum = 0;

        for(i = 0; i<t1IntNum; i++){
            Node n1 = t1.getInternalNode(i);
            if (n1.isRoot())
                continue;
            c1 = NodeUtils.getLeafCount(n1);
            h1 = TreeCmpUtils.getNodeDepth(n1) -1;
            for(j = 0; j<t2IntNum; j++){
                Node n2 = t2.getInternalNode(j);
                if (n2.isRoot())
                    continue;
                //1. |C1| <= n-1 -h2
                //2. |C2| <= n-1- h1,
                c2 = NodeUtils.getLeafCount(n2);
                h2 = TreeCmpUtils.getNodeDepth(n2) -1;
                if( !(c1 <= n-1- h2 &&
                      c2 <= n-1 -h1 ) ){
                    removedNum++;
                }
            }

        }

        return (t1IntNum-1)*(t2IntNum-1) - removedNum;
    }
}
