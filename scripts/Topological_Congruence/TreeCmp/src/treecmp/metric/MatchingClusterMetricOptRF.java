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

import treecmp.common.ClusterDist;
import java.util.BitSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.ListIterator;
import pal.misc.IdGroup;
import pal.tree.Tree;
import pal.tree.TreeUtils;
import treecmp.common.LapSolver;

public class MatchingClusterMetricOptRF extends BaseMetric implements Metric {
  
    public double getDistance(Tree t1, Tree t2) {

        int i, j;
        int metric;
        BitSet[] c1,c2;
        BitSet bs_temp;

        IdGroup idGroup = TreeUtils.getLeafIdGroup(t1);
        BitSet[] c1_temp = ClusterDist.RootedTree2BitSetArray(t1, idGroup);
        BitSet[] c2_temp = ClusterDist.RootedTree2BitSetArray(t2, idGroup);

        int n = t1.getExternalNodeCount();

        int size1_temp = c1_temp.length;
        int size2_temp = c2_temp.length;

        if(size1_temp <= size2_temp){
            c1 = c1_temp;
            c2 = c2_temp;
        }else{
            c2 = c1_temp;
            c1 = c2_temp;
        }
        //s1 is the smaller one

        LinkedList<BitSet> ll1 = new LinkedList<BitSet>();
        for (i=0; i< c1.length;i++){
            ll1.add(c1[i]);
        }

        int hashSetSize=(4*(c2.length+1))/3;
        LinkedHashSet<BitSet> lhs2=new LinkedHashSet<BitSet>(hashSetSize);
        for (i=0; i<c2.length;i++){
            lhs2.add(c2[i]);
        }
        ListIterator<BitSet> it1 = ll1.listIterator();

        while(it1.hasNext()){
            bs_temp=it1.next();
            if (lhs2.remove(bs_temp)){
                it1.remove();
            }
       }

        int size1 = ll1.size();
        int size2 = lhs2.size();


        int size =Math.max(size1, size2);
        if(size<=0)
            return 0;


        short[][] assigncost = new short[size][size];
        int[] rowsol = new int[size];
        int[] colsol = new int[size];
        int[] u = new int[size];
        int[] v = new int[size];

        Iterator<BitSet> it2;
        BitSet bs1,bs2;
        if (size1 > size2) {
            i=0;
            it1 = ll1.listIterator();
            while(it1.hasNext()){
                bs1=it1.next();
                it2=lhs2.iterator();
                j=0;
                while(it2.hasNext()){
                    bs2=it2.next();
                    assigncost[i][j] = (short)ClusterDist.getDistXorBit(bs1, bs2);
                    j++;
                }
                 for (j = size2; j < size1; j++) {
                    assigncost[i][j] = (short)ClusterDist.getDistToOAsMinBit(bs1);
                }
                i++;
            }

        } else {
            i=0;
            it2 = lhs2.iterator();
            while(it2.hasNext()){
                bs2=it2.next();
                it1 = ll1.listIterator();
                j=0;
                while(it1.hasNext()){
                    bs1=it1.next();
                    assigncost[i][j] = (short)ClusterDist.getDistXorBit(bs2, bs1);
                    j++;
                }
                for (j = size1; j < size2; j++) {
                    assigncost[i][j] = (short)ClusterDist.getDistToOAsMinBit(bs2);
                }
               i++;
            }
        }

        //metric = LapSolver.lap(size, assigncost, rowsol, colsol, u, v);
        metric = LapSolver.lapShort(size, assigncost, rowsol, colsol, u, v);
        return metric;
    }
}
