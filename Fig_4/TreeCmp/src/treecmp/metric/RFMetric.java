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
import treecmp.common.SplitDist;
import java.util.BitSet;
import java.util.HashSet;
import pal.misc.IdGroup;
import pal.tree.*;

/**
 *
 * @author Damian
 */
public class RFMetric extends BaseMetric implements Metric{

    public static double getRFDistance(Tree t1, Tree t2) {

        int n = t1.getExternalNodeCount();
        if (n <= 3)
            return 0;

        IdGroup idGroup = TreeUtils.getLeafIdGroup(t1);
        BitSet[] s_t1=SplitDist.getSplits(t1, idGroup);
        BitSet[] s_t2=SplitDist.getSplits(t2, idGroup);
        int N1=s_t1.length;
        int N2=s_t2.length;
        int hashSetSize=(4*(N1+1))/3;

        HashSet<BitSet> s_t1_hs=new HashSet<BitSet>(hashSetSize);

        int i;
        for(i=0;i<N1;i++){
            s_t1_hs.add(s_t1[i]);
        }

        int common=0;
        for(i=0;i<N2;i++){
             if (s_t1_hs.contains(s_t2[i])){
                common++;
            }
        }

        double dist=((double)N1+(double)N2)*0.5-(double)common;
        return dist;

    }

    

    public double getDistance(Tree t1, Tree t2) {

        return RFMetric.getRFDistance(t1, t2);

    }
}
