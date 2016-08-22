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
import pal.tree.Node;
import pal.tree.Tree;
import pal.tree.TreeUtils;
import treecmp.common.ClustIntersectInfoMatrix;
import treecmp.common.LapSolver;
import treecmp.common.TreeCmpUtils;

public class MatchingSpliMetricFree extends BaseMetric implements Metric {

    protected int[] costId2NumT1;
    protected int[] costId2NumT2;
    protected int[] rowsol;
    protected int[] colsol;
    protected short[][] assigncost;
    protected ClustIntersectInfoMatrix cIntM;

    public double getDistance(Tree t1, Tree t2) {

        int metric;
        double metricDouble;
        IdGroup idGroup1 = TreeUtils.getLeafIdGroup(t1);
        IdGroup idGroup2 = TreeUtils.getLeafIdGroup(t2);
        IdGroup idGroup = TreeCmpUtils.mergeIdGroups(idGroup1, idGroup2);
        cIntM = TreeCmpUtils.calcClustIntersectMatrix(t1, t2, idGroup);
        int totL = idGroup.getIdCount();

        int intSize1 = t1.getInternalNodeCount();
        int intSize2 = t2.getInternalNodeCount();
        int extSize1 = t1.getExternalNodeCount();
        int extSize2 = t2.getExternalNodeCount();

        int totSize1 = intSize1 + extSize1;
        int totSize2 = intSize2 + extSize2;

        //we do not want to use root
        int sizeIt = Math.max(totSize1, totSize2);
        int size = sizeIt - 1;

        assigncost = new short[size][size];
        rowsol = new int[size];
        colsol = new int[size];
        int[] u = new int[size];
        int[] v = new int[size];

        if (size <= 0) {
            return 0;
        }

        Node t1Root = t1.getRoot();
        Node t2Root = t2.getRoot();
        
    //    int t1RootNum = t1Root.getNumber();
    //    int t2RootNum = t2Root.getNumber();

        int rootsInterSize = cIntM.getInterSize(t1Root, t2Root);
        
        
        Node[] nodeT1 = TreeCmpUtils.getAllNodes(t1);
        Node[] nodeT2 = TreeCmpUtils.getAllNodes(t2);
        int n1Num, n2Num, intCsize = 0;
        int acInterSize = 0, bcInterSize = 0, adInterSize = 0, bdInterSize = 0;
        int aCsize = 0, cCsize = 0, bCsize = 0, dCsize = 0;
        int r1cInterSize = 0, ar2InterSize = 0;
        int extSum = extSize1 + extSize2;
        int max = 0;
        
        short cost = 0;
        boolean isLeafN1, isLeafN2;
        Node n1 = null, n2 = null;

        int ii, jj;
        ii = 0;
        for (int i = 0; i < sizeIt; i++) {
            if (i < totSize1) {

                n1 = nodeT1[i];
                if (n1.isRoot()) {
                    continue;
                }

                n1Num = n1.getNumber();
                isLeafN1 = n1.isLeaf();
                if (isLeafN1) {
                    aCsize = 1;
                } else {
                    aCsize = cIntM.cSize1[n1Num];
                }
                bCsize = extSize1 - aCsize;
            }
            jj = 0;
            for (int j = 0; j < sizeIt; j++) {
                if (j < totSize2) {
                    n2 = nodeT2[j];
                    if (n2.isRoot()) {
                        continue;
                    }
                    n2Num = n2.getNumber();
                    isLeafN2 = n2.isLeaf();
                    if (isLeafN2) {
                        cCsize = 1;
                    } else {
                        cCsize = cIntM.cSize2[n2Num];
                    }
                    dCsize = extSize2 - cCsize;
                }

                if (i < totSize1 && j < totSize2) {
                    //norma distance
                    
                    r1cInterSize = cIntM.getInterSize(t1Root, n2);
                    ar2InterSize = cIntM.getInterSize(n1, t2Root);
                    
                    acInterSize = cIntM.getInterSize(n1, n2);
                    bcInterSize = r1cInterSize - acInterSize;
                    adInterSize = ar2InterSize - acInterSize;
                    bdInterSize = rootsInterSize - acInterSize - bcInterSize - adInterSize;
                    max = Math.max(acInterSize + bdInterSize, adInterSize + bcInterSize);
                    cost = (short) (extSum - (max << 1));


                } else if (i < totSize1) {
                    cost = (short) Math.min(extSum + aCsize - bCsize, extSum + bCsize - aCsize);
                } else {
                    cost = (short) Math.min(extSum + cCsize - dCsize, extSum + dCsize - cCsize);
                }

                assigncost[ii][jj] = cost;
                jj++;
            }
            ii++;
        }

        metric = LapSolver.lapShort(size, assigncost, rowsol, colsol, u, v);
        metricDouble = metric;
        metricDouble = metricDouble / 2.0;
        return metricDouble;
    }

  
}
