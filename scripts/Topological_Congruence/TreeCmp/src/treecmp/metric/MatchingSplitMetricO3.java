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
import treecmp.common.AlignInfo;
import treecmp.common.ClustIntersectInfoMatrix;
import treecmp.common.ClustIntersectInfoMatrix.ClustPair;
import treecmp.common.IntNodePair;
import treecmp.common.LapSolver;
import treecmp.common.TreeCmpUtils;
import treecmp.config.IOSettings;

public class MatchingSplitMetricO3 extends BaseMetric implements Metric {
    
    protected int[] costId2NumT1;
    protected int[] costId2NumT2;
    protected int[] rowsol;
    protected int[] colsol;
    protected short[][] assigncost;
    protected ClustIntersectInfoMatrix cIntM;

    public double getDistance(Tree t1, Tree t2) {

        int metric, t1NodeNum, t2NodeNum, il, jl, x1, x2;
        Node t1Node, t2Node;
        short n = (short) t1.getExternalNodeCount();
        IdGroup idGroup = TreeUtils.getLeafIdGroup(t1);
        cIntM = TreeCmpUtils.calcClustIntersectMatrix(t1, t2, idGroup);

        int size1 = t1.getInternalNodeCount();
        int size2 = t2.getInternalNodeCount();
        int eqClustSize = cIntM.eqClustList.size();

        int size = Math.max(size1 - eqClustSize, size2 - eqClustSize);
        int sizeIt = Math.max(size1, size2);
     
        assigncost = new short[size][size];
        rowsol = new int[size];
        colsol = new int[size];
        int[] u = new int[size];
        int[] v = new int[size];

        if (size <= 0) {
            return 0;
        }

        // used for alignement generation
        if (IOSettings.getIOSettings().isGenAlignments()) {
            //start of initialization of alignemnt helper tabels
            costId2NumT1 = new int[sizeIt];
            costId2NumT2 = new int[sizeIt];

            //store id map for alignemnt T1
            int ii = 0;
            for (int i = 0; i < size1; i++) {
                t1Node = t1.getInternalNode(i);
                if (t1Node.isRoot()) {
                    continue;
                }
                t1NodeNum = t1Node.getNumber();
                //there is indetical cluster in T2 skip it
                if (cIntM.eqClustT1[t1NodeNum]) {
                    continue;
                }
                costId2NumT1[ii] = t1NodeNum;
                ii++;
            }
            for (int i = ii; i < sizeIt; i++) {
                // -1 means unparied
                costId2NumT1[i] = -1;
            }

            ii = 0;
            //store id map for alignemnt T2
            for (int i = 0; i < size2; i++) {
                t2Node = t2.getInternalNode(i);
                if (t2Node.isRoot()) {
                    continue;
                }
                t2NodeNum = t2Node.getNumber();
                //there is indetical cluster in T2 skip it
                if (cIntM.eqClustT2[t2NodeNum]) {
                    continue;
                }
                costId2NumT2[ii] = t2NodeNum;
                ii++;
            }

            for (int i = ii; i < sizeIt; i++) {
                // -1 means unparied
                costId2NumT2[i] = -1;
            }
            //end of initialization of alignemnt helper tabels
        }

        il = 0;
        for (int i = 0; i < sizeIt; i++) {
            t1NodeNum = -1;
            if (i < size1) {

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
                if (j < size2) {
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
                    x1 = cIntM.cSize1[t1NodeNum] + cIntM.cSize2[t2NodeNum] - (cIntM.intCladeSize[t1NodeNum][t2NodeNum] << 1);
                    x2 = n - x1;
                    assigncost[il][jl] = (short) Math.min(x1, x2);

                } else if (i >= size1 && j < size2) {
                    assigncost[il][jl] = (short) Math.min(n - cIntM.cSize2[t2NodeNum], cIntM.cSize2[t2NodeNum]);
                } else {
                    assigncost[il][jl] = (short) Math.min(n - cIntM.cSize1[t1NodeNum], cIntM.cSize1[t1NodeNum]);
                }
                jl++;
            }
            il++;
        }

        metric = LapSolver.lapShort(size, assigncost, rowsol, colsol, u, v);
        return metric;
    }

 @Override
    public AlignInfo getAlignment() {

        Tree t1 = cIntM.getT1();
        Tree t2 = cIntM.getT2();

        int leafSize = t1.getExternalNodeCount();
        int size1 = t1.getInternalNodeCount();
        int size2 = t2.getInternalNodeCount();
        AlignInfo alignInfo = new AlignInfo();
        int j, cost;
        int size = Math.max(size1, size2);

        int sizeWithoutRoot = size -1;
        IntNodePair[] aln = new IntNodePair[sizeWithoutRoot];
         //store id map for alignemnt T1
        int alnNum = 0;
        for (ClustPair cp: cIntM.eqClustList){
            //skip root clusters
            if (cIntM.cSize1[cp.t1IntId] == leafSize)
                 continue;
             aln[alnNum] = new IntNodePair();
             aln[alnNum].t1_node = cp.t1IntId;
             aln[alnNum].t2_node = cp.t2IntId;
             aln[alnNum].cost = 0;
             alnNum++;
        }

        int totalCost = 0;
        for (int i = 0; i<rowsol.length; i++){
            j = rowsol[i];
            cost =  assigncost[i][j];
            totalCost += cost;
            aln[alnNum] = new IntNodePair();
             // -1 means unparied
             aln[alnNum].t1_node = costId2NumT1[i];
             aln[alnNum].t2_node = costId2NumT2[j];
             aln[alnNum].cost = cost;
             alnNum++;
        }

        alignInfo.setAln(aln);
        alignInfo.setUseClusters(false);
        alignInfo.setT1(t1);
        alignInfo.setT2(t2);
        alignInfo.setTotalCost(totalCost);
        return alignInfo;
    }



}
