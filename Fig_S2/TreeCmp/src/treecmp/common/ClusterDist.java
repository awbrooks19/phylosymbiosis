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

package treecmp.common;

import java.util.BitSet;
import pal.misc.IdGroup;
import pal.tree.Node;
import pal.tree.Tree;

public class ClusterDist {

    public ClusterDist() {
    }

    public static int clusterXor(boolean[] clade_t1, boolean[] clade_t2)
    {
        int n=clade_t1.length;
        int neq=0;

        for(int i=0;i<n;i++) {
            if(clade_t1[i]!=clade_t2[i]) neq++;
        }

        return  neq;
    }
    
    public static BitSet[] RootedTree2BitSetArray(Tree t, IdGroup idGroup) {
        int N = t.getInternalNodeCount();
        int n = t.getExternalNodeCount();
        Node node;
        int j = 0;
        BitSet[] bsA = new BitSet[N - 1];

        for (int i = 0; i < N; i++) {
            node = t.getInternalNode(i);
            if (node.isRoot()) {
                continue;
            }
            bsA[j] = new BitSet(n);
            markNode(idGroup, node, bsA[j]);
            j++;
        }

        return bsA;
    }
     
    static void markNode(IdGroup idGroup, Node node, BitSet cluster) {
        if (node.isLeaf()) {
            String name = node.getIdentifier().getName();
            int index = idGroup.whichIdNumber(name);

            if (index < 0) {
                throw new IllegalArgumentException("INCOMPATIBLE IDENTIFIER (" + name + ")");
            }

            cluster.set(index);
        } else {
            for (int i = 0; i < node.getChildCount(); i++) {
                markNode(idGroup, node.getChild(i), cluster);
            }
        }
    }

    public static int getDistXorBit(BitSet cluster1, BitSet cluster2) {

        BitSet temp = (BitSet) cluster1.clone();
        temp.xor(cluster2);
        int d = temp.cardinality();
        return d;
    }
    public static int getAndBit(BitSet cluster1, BitSet cluster2) {

        BitSet temp = (BitSet) cluster1.clone();
        temp.and(cluster2);
        int d = temp.cardinality();
        return d;
    }

    public static int getDistToOAsMinBit(BitSet cluster) {

        int t = cluster.cardinality();
        return t;

    }
}
