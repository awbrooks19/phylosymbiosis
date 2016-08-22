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
import pal.tree.NodeUtils;
import pal.tree.SplitSystem;
import pal.tree.Tree;

public class SplitDist {

    /** Creates a new instance of SplitDist */
    public SplitDist() {
    }

    //dist([A1|A2],[B1|B2])=0.5*min((A1 xor B1)+(A2 xor B2),(A1 xor B2)+(A2 xor B1))
    public static double getDist1(boolean[] split1, boolean[] split2) {
        int n = split1.length;
        int eq = 0, neq = 0, q;
        double metric = 0.0;

        for (int i = 0; i < n; i++) {
            if (split1[i] == split2[i]) {
                eq++;
            } else {
                neq++;
            }
        }

        if (eq > neq) {
            q = eq;
        } else {
            q = neq;
        }

        return (double) (n - q);
    }

    public static int getDist1Int(boolean[] split1, boolean[] split2) {
        int n = split1.length;
        int eq = 0, neq = 0, q;

        for (int i = 0; i < n; i++) {
            if (split1[i] == split2[i]) {
                eq++;
            } else {
                neq++;
            }
        }

        if (eq > neq) {
            q = eq;
        } else {
            q = neq;
        }

        return (n - q);
    }


    public static double getDistToO_1(boolean[] split) {
        int n = split.length;

        return Math.ceil(Math.floor(n / 2.0) / 2.0);

    }

    public static double getDistToO_2(boolean[] split) {
        int n = split.length;
        int s_true = 0;

        for (int i = 0; i < n; i++) {
            if (split[i] == true) {
                s_true++;
            }
        }

        return Math.min(s_true, n - s_true);

    }

    public static int getMinSize(boolean[] split) {

        int n = split.length;
        int s_true = 0;

        for (int i = 0; i < n; i++) {
            if (split[i] == true) {
                s_true++;
            }
        }
        return Math.min(s_true, n - s_true);
    }

    public static int getMaxSize(boolean[] split) {
        int n = split.length;
        int max = n - getMinSize(split);

        return max;
    }

    public static int getDist1Bit(BitSet split1, BitSet split2, int n) {
        BitSet temp = (BitSet) split1.clone();

        temp.xor(split2);

        int neq = temp.cardinality();
        int eq = n - neq;
        int d;

        if (neq < eq) {
            d = neq;
        } else {
            d = eq;
        }
        return d;
    }

    public static int getDistToOAsMinBit(BitSet split, int n) {

        int t = split.cardinality();
        int f = n - t;
        int d;

        if (t < f) {
            d = t;
        } else {
            d = f;
        }

        return d;
    }

    public static int getDistToOAsN4Bit(BitSet split, int n) {
        return (int) Math.ceil(Math.floor(n / 2.0) / 2.0);
    }

    public static BitSet[] SplitSystem2BitSetArray(SplitSystem s) {
        int N = s.getSplitCount();

        BitSet[] bsA = new BitSet[N];
        int n = s.getLabelCount();
        boolean[] split;

        for (int i = 0; i < N; i++) {
            bsA[i] = new BitSet(n);
            split = s.getSplit(i);
            for (int j = 0; j < n; j++) {
                if (split[j] == true) {
                    bsA[i].set(j);
                }
            }
        }

        return bsA;
    }
    public static BitSet[] getSplits(Tree t, IdGroup idGroup) {

        int n = t.getExternalNodeCount();
        int internal = t.getInternalNodeCount();
        int size = internal - 1;
        BitSet[] splits = new BitSet[size];

        Node curNode = t.getExternalNode(0);

        int ind = 0;
        Node child;
        int childCount = 0;
        int childInd = 0;
        int leafId;
        int i;

        while (!curNode.isRoot()) {
            if (!curNode.isLeaf()) {
                BitSet bs = new BitSet(n);
                ind = curNode.getNumber();
                childCount = curNode.getChildCount();
                for (i = 0; i < childCount; i++) {
                    child = curNode.getChild(i);
                    if (child.isLeaf()) {
                        leafId = idGroup.whichIdNumber(child.getIdentifier().getName());
                        bs.set(leafId);
                    } else {
                        childInd = child.getNumber();
                        bs.or(splits[childInd]);
                    }
                }
                splits[ind] = bs;
            }
            curNode = NodeUtils.postorderSuccessor(curNode);
        }

        // standardize split (i.e. first index is alway true)
	for(i=0;i<size;i++){
            if (splits[i].get(0) == false)
                splits[i].flip(0, n);
	}

        return splits;
    }
}
