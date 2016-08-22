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


import pal.tree.Tree;
import treecmp.common.AlignInfo;
import treecmp.config.IOSettings;

public class MatchingClusterMetric extends BaseMetric implements Metric {

    private MatchingClusterMetricO3 mc03;
    private MatchingClusterMetricOptRF mcRF;
    private MatchingClusterMetricFree mcFree;

    public MatchingClusterMetric() {
        super();
        mc03 = new MatchingClusterMetricO3();
        mcRF = new MatchingClusterMetricOptRF();
        mcFree = new MatchingClusterMetricFree();
    }

    public double getDistance(Tree t1, Tree t2) {

        if (IOSettings.getIOSettings().isOptMsMcByRf()) {
            return mcRF.getDistance(t1, t2);
        }
        if (IOSettings.getIOSettings().isUseMsMcFreeLeafSet()) {
            return mcFree.getDistance(t1, t2);
        }

        return mc03.getDistance(t1, t2);
    }

    @Override
    public AlignInfo getAlignment() {
        if (IOSettings.getIOSettings().isOptMsMcByRf()) {
            return null;
        }
        if (IOSettings.getIOSettings().isUseMsMcFreeLeafSet()) {
            return null;
        }

        return mc03.getAlignment();
    }
}
