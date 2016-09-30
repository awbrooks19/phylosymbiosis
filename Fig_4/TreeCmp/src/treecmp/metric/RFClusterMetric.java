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
import java.util.HashSet;
import pal.misc.IdGroup;
import pal.tree.*;

public class RFClusterMetric extends BaseMetric implements Metric{

    public static double getRFClusterMetric(Tree t1, Tree t2) {

        IdGroup idGroup = TreeUtils.getLeafIdGroup(t1);

        BitSet[] bs1 = ClusterDist.RootedTree2BitSetArray(t1, idGroup);
        BitSet[] bs2 = ClusterDist.RootedTree2BitSetArray(t2, idGroup);

        int size1 = bs1.length;
        int size2 = bs2.length;
        int hashSetSize=(4*(size1+1))/3;

        HashSet<BitSet> hs1=new HashSet<BitSet>(hashSetSize);

        for(int i=0;i<size1;i++){
            hs1.add(bs1[i]);
        }

        int common=0;
        for(int i=0;i<size2;i++){
            if (hs1.contains(bs2[i]))
                common++;
        }

        double dist=((double)size1+(double)size2)*0.5-(double)common;
        return dist;

    }

    

    public double getDistance(Tree t1, Tree t2) {

        return RFClusterMetric.getRFClusterMetric(t1, t2);

    }
}
