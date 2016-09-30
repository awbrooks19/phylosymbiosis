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

import java.text.MessageFormat;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Set;
import pal.tree.Node;
import pal.tree.NodeUtils;
import pal.tree.Tree;
import treecmp.config.IOSettings;
import treecmp.io.ResultWriter;

public class AlignWriter {

    public final static String TREE = "Tree: ";
    public final static String TREE1 = "Tree1: ";
    public final static String TREE2 = "Tree2: ";
    public final static String REF_TREE = "Referenced tree";
    public final static String S_SPACE = " ";
    public final static String NEW_LINE = System.getProperty("line.separator");
    public final static String ALIGN = "Alignemnt: ";
    public final static String ALIGN_COST = "Alignemnt cost (distance): {0}";
    public final static String ALIGN_SEP = "-----------------";
    public final static String EQUAL_CLUSTERS = "Equal clusters in both trees:";
    public final static String ALIGN_CLUSTERS = "Aligned clusters:";
    public final static String EQUAL_SPLITS = "Equal splits in both trees:";
    public final static String ALIGN_SPLITS = "Aligned splits:";
    public final static String PAIR_COST = "Cost: {0}";
    public final static String MATCH_CLUST = "{0}. {1} <-> {2}";
    public final static String MATCH_EMPTY = "EMPTY";
 
    private String outFileName;
    private boolean isGenAlignments;
    private ResultWriter[] rw;

    public AlignWriter() {
        this.isGenAlignments = IOSettings.getIOSettings().isGenAlignments();
        this.outFileName = IOSettings.getIOSettings().getOutputFile();
    }

    public void initFiles(StatCalculator[] stats) {
        if (!isGenAlignments) {
            return;
        }

        rw = new ResultWriter[stats.length];

        for (int i = 0; i < stats.length; i++) {
            if (stats[i].getAlnFileSuffix() != null) {
                String name = outFileName + stats[i].getAlnFileSuffix();
                rw[i] = new ResultWriter();
                rw[i].setFileName(name);
                rw[i].isWriteToFile(true);
                rw[i].init();
            }
        }
    }

    public void closeFiles(StatCalculator[] stats) {
        if (!isGenAlignments) {
            return;
        }
        for (int i = 0; i < rw.length; i++) {
            if (rw[i] != null) {
                rw[i].close();
            }
        }
    }

    public void writeAlignments(int rowNum, int t1, int t2, StatCalculator[] stats) {
        if (!isGenAlignments) {
            return;
        }

        for (int i = 0; i < stats.length; i++) {
            AlignInfo alignInfo = stats[i].getAlignment();
            if (alignInfo != null) {
                writeAlignContent(rowNum, t1, t2, alignInfo, rw[i]);
            }
        }
    }

    private void writeAlignContent(int rowNum, int t1, int t2, AlignInfo alignInfo, ResultWriter rw) {

        if (alignInfo == null)
            return;

        Tree tree1 = alignInfo.getT1();
        Tree tree2 = alignInfo.getT2();


        StringBuilder sb = new StringBuilder();
        sb.append(ALIGN);
        sb.append(rowNum);
        sb.append(NEW_LINE);

        sb.append(TREE1);
        sb.append(t1);
        sb.append(NEW_LINE);

        sb.append(TREE2);
        if (t2 != -1) {
            sb.append(t2);
        } else {
            sb.append(REF_TREE);
        }
        sb.append(NEW_LINE);

        String aln_cost = MessageFormat.format(ALIGN_COST, alignInfo.getTotalCost());
        sb.append(aln_cost);
        sb.append(NEW_LINE);

        if (alignInfo.isUseClusters()) {
            sb.append(ALIGN_CLUSTERS);
        } else {
            sb.append(ALIGN_SPLITS);
        }
        sb.append(NEW_LINE);

        rw.setText(sb.toString());
        rw.write_pure();

        IntNodePair[] aln = alignInfo.getAln();

        int num;
        String t1_part, t2_part, aln_row, cost_row;
        for (int i = 0; i < aln.length; i++) {
            sb = new StringBuilder();
            num = i + 1;
            IntNodePair np = aln[i];

            if (np.t1_node != -1) {
                t1_part = getClusterString(np.t1_node, tree1, alignInfo.isUseClusters());
            } else {
                t1_part = MATCH_EMPTY;
            }

            if (np.t2_node != -1) {
                t2_part = getClusterString(np.t2_node, tree2, alignInfo.isUseClusters());
            } else {
                t2_part = MATCH_EMPTY;
            }

            aln_row = MessageFormat.format(MATCH_CLUST, num, t1_part, t2_part);
            sb.append(aln_row);
            sb.append(NEW_LINE);

            cost_row = MessageFormat.format(PAIR_COST, np.cost);
            sb.append(cost_row);
            sb.append(NEW_LINE);

            rw.setText(sb.toString());
            rw.write_pure();
        }
        rw.setText(ALIGN_SEP);
        rw.write();

    }

    private String getClusterString(int num, Tree t, boolean printAsCluster) {

        Node n = t.getInternalNode(num);
        int leafNum = t.getExternalNodeCount();
        Set<String> nameSet = new HashSet<String>((leafNum*4)/3);
        for (int i =0; i<leafNum;i++){
            String name = t.getExternalNode(i).getIdentifier().getName();
                nameSet.add(name);
        }
        
        Node[] leaves = NodeUtils.getExternalNodes(n);
        Arrays.sort(leaves, new LeafNodeComparator());
        String id;
        StringBuilder sb = new StringBuilder();
        sb.append('(');
        id = leaves[0].getIdentifier().getName();
        sb.append(id);
        nameSet.remove(id);
        for (int i = 1; i < leaves.length; i++) {
            id = leaves[i].getIdentifier().getName();
            sb.append(", ");
            sb.append(id);
            nameSet.remove(id);
        }
        if (printAsCluster) {
            sb.append(')');
        } else {
           // sb.append("|...)");
            sb.append("| ");

            String nameTab[] = nameSet.toArray(new String[0]);
            Arrays.sort(nameTab);
            
            id = nameTab[0];
            sb.append(id);
            for (int i=1; i<nameTab.length;i++) {
                id = nameTab[i];
                sb.append(", ");
                sb.append(id);
            }
             sb.append(')');
        }

        return sb.toString();
    }

    private String getHeader(StatCalculator stat) {

        return "";

    }
}
class LeafNodeComparator implements Comparator<Node>{

    public int compare(Node n1, Node n2){
        String id1 = n1.getIdentifier().getName();
        String id2 = n2.getIdentifier().getName();
        return id1.compareTo(id2);
    }

}
