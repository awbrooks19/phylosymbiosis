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
import pal.tree.*;
import pal.io.*;
import treecmp.qt.*;

/**
 *
 * @author Damian
 */
public class QuartetMetricDouble extends BaseMetric implements Metric {

    public QuartetMetricDouble() {
    }

    public static double getQuartetDistance(pal.tree.Tree tree1, pal.tree.Tree tree2) {


        OutputTarget tree1OT = OutputTarget.openString();
        OutputTarget tree2OT = OutputTarget.openString();

        TreeUtils.printNH(tree1, tree1OT, false, false);
        TreeUtils.printNH(tree2, tree2OT, false, false);

        String tree1Newick = tree1OT.getString();
        String tree2Newick = tree2OT.getString();

        tree1OT.close();
        tree2OT.close();

        Distance d = new GeneralN2DQDistDoubleShort();
        double dist = -1.0;
        try {
            treecmp.qt.Tree tree_tt1 = new treecmp.qt.Tree(tree1Newick);
            treecmp.qt.Tree tree_tt2 = new treecmp.qt.Tree(tree2Newick);

            DistResult dr = d.getMeasures(tree_tt1, tree_tt2);
            dist = (double) (dr.qdist() + dr.q1() + dr.q2());
        } catch (ParseException ex) {
            ex.printStackTrace();
        }
        return dist;

    }

    public double getDistance(pal.tree.Tree t1, pal.tree.Tree t2) {

        return QuartetMetricDouble.getQuartetDistance(t1, t2);
    }
}
