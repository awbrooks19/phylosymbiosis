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

import java.util.Locale;
import treecmp.config.IOSettings;


public class ReportUtils {
    private final static int ROW_PRECISION = 4;
    public final static String ROW_DATA_FORMAT="%1$."+ROW_PRECISION+"f";
    private final static IOSettings ioSet = IOSettings.getIOSettings();
    private final static String sep = ioSet.getSSep();
    private final static boolean pruneTrees = ioSet.isPruneTrees();
    private final static boolean randomComparison = ioSet.isRandomComparison();

    public final static String NUM_COLUMN = "No";
    public final static String T1_COLUMN = "Tree1";
    public final static String T2_COLUMN = "Tree2";
    public final static String T_COLUMN = "Tree";
    public final static String T1_TAXA= "Tree1_taxa";
    public final static String T2_TAXA= "Tree2_taxa";
    public final static String T_TAXA= "Tree_taxa";
    public final static String RT_TAXA= "RefTree_taxa";
    public final static String COMMON_TAXA= "Common_taxa";
    public final static String YULE_FRAC= "_toYuleAvg";
    public final static String UNIF_FRAC= "_toUnifAvg";
    public final static String NA_FRAC= "N/A";
    //it t2 ==-1 do not print t2
    public static String getResultRow(int rowNum, int t1, int t2, StatCalculator[] stats){

        StringBuilder sb = new StringBuilder();
        double dist, distToYuleAvg, distToUnifAvg;
        String distStr;

        //we do not need num in ref mode - tree number has the same valuse
        if (t2 != -1){
            sb.append(rowNum);
            sb.append(sep);
        }
        sb.append(t1);
        
        //used in reference tree mode
        if (t2 != -1){
            sb.append(sep);
            sb.append(t2);
        }
        //to do if prune enabled
        if (pruneTrees && stats.length > 0){
            sb.append(sep);
            sb.append(stats[0].getT1TaxaNum());
            sb.append(sep);
            sb.append(stats[0].getT2TaxaNum());
            sb.append(sep);
            sb.append(stats[0].getCommonTaxaNum());
        }

        for (int i=0; i< stats.length; i++){
             sb.append(sep);
             
             dist = stats[i].getLastDist();
             distStr = String.format(Locale.US,ROW_DATA_FORMAT, dist);
             sb.append(distStr);
             if (randomComparison){
                 distToYuleAvg = stats[i].getLastDistToYuleAvg();
                 if (distToYuleAvg != Double.NEGATIVE_INFINITY)
                    distStr = String.format(Locale.US,ROW_DATA_FORMAT, distToYuleAvg);
                 else
                     distStr = NA_FRAC;
                 sb.append(sep);
                 sb.append(distStr);

                 distToUnifAvg = stats[i].getLastDistToUnifAvg();
                 if (distToUnifAvg != Double.NEGATIVE_INFINITY)
                    distStr = String.format(Locale.US,ROW_DATA_FORMAT, distToUnifAvg);
                 else
                     distStr = NA_FRAC;
                 sb.append(sep);
                 sb.append(distStr);           
             }
        }

        //sb.append("\n");
        return sb.toString();
    }

     public static String getHeaderRow(StatCalculator[] stats){
         
         return getHeaderRow(stats, false);
     }
     public static String getHeaderRow(StatCalculator[] stats, boolean ifReefTreeMode){
        StringBuilder sb = new StringBuilder();

        String metricName;

         if (!ifReefTreeMode){
            sb.append(NUM_COLUMN);
            sb.append(sep);
            sb.append(T1_COLUMN);
        }else{
            //in ref tree mode Tree header is enough
            sb.append(T_COLUMN);
        }
       
        //used in reference tree mode
        if (!ifReefTreeMode){
            sb.append(sep);
            sb.append(T2_COLUMN);
        }
        //to do if prune enabled
        if (pruneTrees && stats.length > 0) {
            sb.append(sep);

            if (!ifReefTreeMode) {
                sb.append(T1_TAXA);
            } else {
                sb.append(T_TAXA);
            }

            sb.append(sep);

            if (!ifReefTreeMode) {
                sb.append(T2_TAXA);
            } else {
                sb.append(RT_TAXA);
            }

            sb.append(sep);
            sb.append(COMMON_TAXA);
        }

        for (int i = 0; i < stats.length; i++) {
            sb.append(sep);
            metricName = stats[i].getName();
            sb.append(metricName);
            if (randomComparison){
                sb.append(sep);             
                sb.append(metricName);
                sb.append(YULE_FRAC);
                sb.append(sep);
                sb.append(metricName);
                sb.append(UNIF_FRAC);
            }
        }
        //sb.append("\n");
        return sb.toString();
    }


}
