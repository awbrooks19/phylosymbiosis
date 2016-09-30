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

package treecmp.command;

import java.util.ArrayList;
import treecmp.common.ProgressIndicator;
import treecmp.common.StatCalculator;
import treecmp.common.SummaryStatCalculator;
import pal.tree.Tree;
import treecmp.common.AlignWriter;
import treecmp.common.ReportUtils;
import treecmp.common.TreeCmpException;
import treecmp.io.ResultWriter;
import treecmp.io.TreeReader;
import treecmp.config.ActiveMetricsSet;
import treecmp.config.IOSettings;
import treecmp.metric.Metric;

public class RunMCommand extends Command {

    public RunMCommand(int paramNumber, String name) {
        super(paramNumber, name);
    }

    @Override
    public void run() throws TreeCmpException{
        super.run();       
        out.init();
        reader.open();
        
        matrixCompareExecute(reader, out);
        
        reader.close();
        out.close();
    }

    private void matrixCompareEx(TreeReader reader, ResultWriter out, StatCalculator[] metrics) throws TreeCmpException {

        pal.tree.Tree tree1,tree2 ;
        ArrayList<Tree> tree_vec = new ArrayList<Tree>();
        int k,counter,maxIt;
        double val;
        String row;

        int mSize = metrics.length;

        //initialize summary stat calculators
        SummaryStatCalculator[] sStatCalc=new SummaryStatCalculator[mSize];
        for(int i=0;i<mSize;i++){
            sStatCalc[i]=new SummaryStatCalculator(metrics[i]);
        }

        String head = ReportUtils.getHeaderRow(metrics);
        out.setText(head);
        out.write();

        AlignWriter aw = new AlignWriter();
        aw.initFiles(metrics);

        while ((tree1 = reader.readNextTree()) != null) {
            tree_vec.add(tree1);
        }

        int N = tree_vec.size();
        counter = 1;
        maxIt = N*(N-1)/2;
        ProgressIndicator progress = new ProgressIndicator();
        progress.setMaxVal(maxIt);
        progress.setPrintInterval(600);
        progress.setPrintPercentInterval(5.0);
        progress.init();

        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                tree1 = tree_vec.get(i);
                tree2 = tree_vec.get(j);                
                for(k=0; k<metrics.length; k++){
                    val = metrics[k].getDistance(tree1, tree2);
                    sStatCalc[k].insertValue(val);
                }   
                row = ReportUtils.getResultRow(counter, i + 1, j + 1, metrics);
                out.setText(row);
                out.write();

                aw.writeAlignments(counter, i + 1, j + 1, metrics);

                progress.displayProgress(counter);
                counter++;
            }
        }
        aw.closeFiles(metrics);
        
        SummaryStatCalculator.printSummary(out, sStatCalc);
    }  
    
    public void matrixCompareExecute(TreeReader reader, ResultWriter out ) throws TreeCmpException {

        Metric[] metrics = ActiveMetricsSet.getActiveMetricsSet().getActiveMetricsTable();
        StatCalculator[] statsMetrics=new StatCalculator[metrics.length];

        for(int i=0; i<metrics.length; i++){
            statsMetrics[i] = new StatCalculator(metrics[i]);

            if(IOSettings.getIOSettings().isCalcCorrelation())//temprary set statcalc to hold valuse
                statsMetrics[i].setRecordValues(true);
        }
        matrixCompareEx(reader, out, statsMetrics);
    }
}
