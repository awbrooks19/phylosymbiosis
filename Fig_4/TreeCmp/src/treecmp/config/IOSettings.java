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

package treecmp.config;


public class IOSettings {

    private static IOSettings IOConf;
    private String inputFile;
    private String outputFile;
    private String sSep;
    private boolean pruneTrees;
    private boolean randomComparison;
    private boolean optMsMcByRf;
    private boolean genAlignments;
    private boolean useMsMcFreeLeafSet;;
    private int iStep;
    //defalut false
    private boolean calcCorrelation;
    private boolean genSummary;

    public boolean isUseMsMcFreeLeafSet() {
        return useMsMcFreeLeafSet;
    }

    public void setUseMsMcFreeLeafSet(boolean useMsMcFreeLeafSet) {
        this.useMsMcFreeLeafSet = useMsMcFreeLeafSet;
    }

    public boolean isGenSummary() {
        return genSummary;
    }

    public void setGenSummary(boolean genSummary) {
        this.genSummary = genSummary;
    }

    public boolean isRandomComparison() {
        return randomComparison;
    }

    public void setRandomComparison(boolean randomComparison) {
        this.randomComparison = randomComparison;
    }

    public boolean isCalcCorrelation() {
        return calcCorrelation;
    }

    public void setCalcCorrelation(boolean calcCorrelation) {
        this.calcCorrelation = calcCorrelation;
    }

    public String getSSep() {
        return sSep;
    }

    public void setSSep(String sSep) {
        this.sSep = sSep;
    }

    public int getIStep() {
        return iStep;
    }

    public void setIStep(int iStep) {
        this.iStep = iStep;
    }

    public String getInputFile() {
        return inputFile;
    }

    public void setInputFile(String inputFile) {
        this.inputFile = inputFile;
    }

    public String getOutputFile() {
        return outputFile;
    }

    public void setOutputFile(String outputFile) {
        this.outputFile = outputFile;
    }
    

     protected IOSettings()
     {
         inputFile = null;
         outputFile = null;
         iStep = 1;
         calcCorrelation = false;
         pruneTrees = false;
         randomComparison = false;
         optMsMcByRf = false;
         genAlignments = false;
         genSummary = false;


     }
     public static IOSettings getIOSettings()
    {
        if(IOConf==null)
        {
            IOConf=new IOSettings();
        }
        return IOConf;
    }

    public boolean isPruneTrees() {
        return pruneTrees;
    }

    public void setPruneTrees(boolean pruneTrees) {
        this.pruneTrees = pruneTrees;
    }

    public boolean isGenAlignments() {
        return genAlignments;
    }

    public void setGenAlignments(boolean genAlignments) {
        this.genAlignments = genAlignments;
    }

    public boolean isOptMsMcByRf() {
        return optMsMcByRf;
    }

    public void setOptMsMcByRf(boolean optMsMcByRf) {
        this.optMsMcByRf = optMsMcByRf;
    }
}


