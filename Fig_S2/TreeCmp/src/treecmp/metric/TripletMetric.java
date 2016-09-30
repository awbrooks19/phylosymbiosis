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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import pal.misc.IdGroup;
import pal.tree.Tree;
import pal.tree.TreeUtils;
import treecmp.common.TreeCmpUtils;

public class TripletMetric extends BaseMetric implements Metric {
    private TripletMetric2 tt2;
  public TripletMetric(){
      super();
      tt2 = new TripletMetric2();
  }

    public double getDistance(Tree t1, Tree t2) {

        if (TreeCmpUtils.isBinary(t1, true) && TreeCmpUtils.isBinary(t2, true)) {
            return getDistForBinary(t1, t2);
        }
        //run distance for arbitrary tree in O(n^2) time
        return tt2.getDistance(t1, t2);

    }

    public double getDistForBinary(Tree t1, Tree t2) {
        IdGroup id1 = TreeUtils.getLeafIdGroup(t1);
        int[][] lcaMatrix1 = TreeCmpUtils.calcLcaMatrix(t1, null);
        int[][] lcaMatrix2 = TreeCmpUtils.calcLcaMatrix(t2, id1);
        int n = lcaMatrix1.length;
        long n_l = (long) n;
        long val_l;
        long commonT = 0;

        for (int i = 0; i < n; i++) {
            List<Integer> numList = getPatternNum(i, lcaMatrix1, lcaMatrix2);
            for (Integer val : numList) {
                val_l = (long) val;
                commonT += val_l * (val_l - 1) / 2;
            }
        }

        long dist = n_l * (n_l - 1) * (n_l - 2) / 6 - commonT;
        return (double) dist;
    }

    class Pattern {

        public int a;
        public int b;

        public Pattern(int a, int b) {
            this.a = a;
            this.b = b;
        }

        @Override
        public int hashCode() {
            return 31 * a + b;

        }

        @Override
        public boolean equals(Object obj) {
            if (obj == this) {
                return true;
            }
            if (obj instanceof Pattern) {
                Pattern test = (Pattern) obj;
                if (a == test.a && b == test.b) {
                    return true;
                }
            }
            return false;
        }
    }

    public List<Integer> getPatternNum(int x, int[][] a, int[][] b) {

        int n = a.length;
        int mapSize = (4 * n) / 3;

        List<Integer> mList = new ArrayList<Integer>();

        Map<Pattern, Integer> patternMap = new HashMap<Pattern, Integer>(mapSize);
        for (int i = 0; i < n; i++) {
            if (i == x) {
                continue;
            }
            Pattern p = new Pattern(a[x][i], b[x][i]);
            Integer num = patternMap.get(p);
            if (num == null) {
                patternMap.put(p, new Integer(1));
            } else {
                patternMap.put(p, num.intValue() + 1);
            }

        }
        for (Integer val : patternMap.values()) {
            if (val >= 2) {
                mList.add(val);
            }
        }
        return mList;
    }

    
}
