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

import pal.tree.Tree;

public class AlignInfo {

    private IntNodePair[] aln;
    private Tree t1;
    private Tree t2;
    private int totalCost;
    private boolean useClusters;

    public Tree getT1() {
        return t1;
    }

    public void setT1(Tree t1) {
        this.t1 = t1;
    }

    public Tree getT2() {
        return t2;
    }

    public void setT2(Tree t2) {
        this.t2 = t2;
    }

    public IntNodePair[] getAln() {
        return aln;
    }

    public void setAln(IntNodePair[] aln) {
        this.aln = aln;
    }

    public boolean isUseClusters() {
        return useClusters;
    }

    public void setUseClusters(boolean useClusters) {
        this.useClusters = useClusters;
    }

    public double getTotalCost() {
        return totalCost;
    }

    public void setTotalCost(int totalCost) {
        this.totalCost = totalCost;
    }
}
