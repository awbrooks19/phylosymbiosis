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

public class MatchingClusterMetricFree extends BaseMetric implements Metric {

    protected int[] costId2NumT1;
    protected int[] costId2NumT2;
    protected int[] rowsol;
    protected int[] colsol;
    protected short[][] assigncost;
    protected ClustIntersectInfoMatrix cIntM;

    public double getDistance(Tree t1, Tree t2) {

        int metric;

        IdGroup idGroup1 = TreeUtils.getLeafIdGroup(t1);
        IdGroup idGroup2 = TreeUtils.getLeafIdGroup(t2);
        IdGroup idGroup = TreeCmpUtils.mergeIdGroups(idGroup1,idGroup2);
        cIntM = TreeCmpUtils.calcClustIntersectMatrix(t1, t2, idGroup);
        
        int intSize1 = t1.getInternalNodeCount();
        int intSize2 = t2.getInternalNodeCount();
        int extSize1 = t1.getExternalNodeCount();
        int extSize2 = t2.getExternalNodeCount();
        
        int totSize1 = intSize1 + extSize1;
        int totSize2 = intSize2 + extSize2;

        int size = Math.max(totSize1, totSize2);

        assigncost = new short[size][size];
        rowsol = new int[size];
        colsol = new int[size];
        int[] u = new int[size];
        int[] v = new int[size];
        
        if (size <= 0) {
            return 0;
        }
        
       Node[] nodeT1 = TreeCmpUtils.getAllNodes(t1);
       Node[] nodeT2 = TreeCmpUtils.getAllNodes(t2);
       int n1Num,n2Num, n1Csize = 0, n2Csize = 0, intCsize = 0;
       short cost = 0;
       boolean isLeafN1, isLeafN2;
       Node n1 = null, n2 = null;
       for (int i = 0; i < size; i++){
           if (i < totSize1 ){

           n1 = nodeT1[i];
           n1Num = n1.getNumber();
           isLeafN1 = n1.isLeaf();
           if (isLeafN1)
               n1Csize = 1;
           else
               n1Csize = cIntM.cSize1[n1Num];
           }
           for (int j = 0; j < size; j++){
              if (j < totSize2 ){
               n2 = nodeT2[j];
                n2Num = n2.getNumber();
                isLeafN2 = n2.isLeaf();
                if (isLeafN2)
                    n2Csize = 1;
                else
                    n2Csize = cIntM.cSize2[n2Num];
               }

                if (i < totSize1 && j < totSize2 ){
                    //norma distance
                    intCsize = cIntM.getInterSize(n1, n2);
                    cost = (short) (n1Csize + n2Csize - (intCsize << 1));

                } else if (i < totSize1 ){
                    cost = (short) n1Csize;
                }else{
                    cost = (short) n2Csize;
                }

                assigncost[i][j] = cost;
            }
       }
/*

        il = 0;
        for (int i = 0; i < sizeIt; i++) {
            t1NodeNum = -1;
            if (i < size1){
                t1Node = t1.getInternalNode(i);
                if (t1Node.isRoot()) {
                    continue;
                }
                t1NodeNum = t1Node.getNumber();
                //there is indetical cluster in T2 skip it
                if (cIntM.eqClustT1[t1NodeNum]) {
                    continue;
                }
            }
            jl = 0;
            for (int j = 0; j < sizeIt; j++) {
                t2NodeNum = -1;
                if (j < size2){
                    t2Node = t2.getInternalNode(j);
                    if (t2Node.isRoot()) {
                        continue;
                    }
                        t2NodeNum = t2Node.getNumber();
                //there is indetical cluster in T1 skip it
                    if (cIntM.eqClustT2[t2NodeNum]) {
                        continue;
                    }
                }
                if (i < size1 && j < size2) {
                    assigncost[il][jl] = (short) (cIntM.cSize1[t1NodeNum] + cIntM.cSize2[t2NodeNum] - (cIntM.intCladeSize[t1NodeNum][t2NodeNum] << 1));

                } else if (i >= size1 && j < size2) {
                    assigncost[il][jl] = cIntM.cSize2[t2NodeNum];
                } else {
                    assigncost[il][jl] = cIntM.cSize1[t1NodeNum];
                }
                jl++;
            }
            il++;
        }
*/
        metric = LapSolver.lapShort(size, assigncost, rowsol, colsol, u, v);
        return metric;
    }

  
}
