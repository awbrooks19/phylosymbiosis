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

public class ProgressIndicator {

    private final static int ROW_PRECISION=2;
    private long startTime;
    private int maxVal;

    private long lastPrintTime;

    /**
     * max number of seconds without reporting of the progress
     */
    private int printInterval;
    private long printIntervalMilis;
    private int lastPrintVal;
    private double printPercentInterval;

    private String dataFormat;

    public ProgressIndicator() {
        maxVal=0;
        lastPrintTime=0;
        lastPrintVal=0;
        printInterval=0;
        printPercentInterval=0;
        startTime=0;


    }

    public double getPrintPercentInterval() {
        return printPercentInterval;
    }

    public void setPrintPercentInterval(double printPercentInterval) {
        this.printPercentInterval = printPercentInterval;
    }




    public int getMaxVal() {
        return maxVal;
    }

    public void setMaxVal(int maxVal) {
        this.maxVal = maxVal;
    }

    public int getPrintInterval() {
        return printInterval;
    }

    public void setPrintInterval(int printInterval) {
        this.printInterval = printInterval;
        this.printIntervalMilis=printInterval*1000;
    }


   public void init()
    {
       long now = System.currentTimeMillis();
       this.dataFormat="%1$."+ROW_PRECISION+"f";
       startTime=now;

        if(lastPrintTime==0)
        {
            lastPrintTime=now;
            printStatus("Start of calculation...please wait...");

            String msg=String.format(Locale.US, this.dataFormat, 0.0);
            printStatus(msg+"% completed...");
         }

    }

    public void displayProgress(int currentVal)
    {


       long now = System.currentTimeMillis();

        String msg="";

        double prog=0;

        long timeDiff=now-lastPrintTime;
        double ratio=((currentVal-lastPrintVal)/(double)maxVal)*100.0;

        if(currentVal>=maxVal)
        {

            msg=String.format(Locale.US, this.dataFormat, 100.0);
            printStatus(msg+"% completed.");
            printStatus("End of calculation.");
            printStatus("Total calculation time: "+(now-startTime)+" ms.");

        }else if(timeDiff>=this.printIntervalMilis ||ratio>=printPercentInterval)
        {
            this.lastPrintVal=currentVal;
            this.lastPrintTime=now;
            prog=(currentVal/(double)maxVal)*100.0;
            msg=String.format(Locale.US, this.dataFormat, prog);
            printStatus(msg+"% completed...");

        }

    }
    private void printStatus(String status)
    {

        System.out.println(TimeDate.now()+": "+status);
    }


}


