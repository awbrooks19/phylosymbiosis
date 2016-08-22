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

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import pal.tree.Tree;
import treecmp.common.AlignInfo;
import treecmp.config.ConfigSettings;
import treecmp.statdata.IMetircDistrbHolder;
import treecmp.statdata.MetircDistrbHolder;
import treecmp.statdata.MetricDistribution;

public abstract class BaseMetric implements Metric{

    protected String name;
    protected String commandLineName;
    protected String description;
    protected String unifomFileName;
    protected String yuleFileName;
    protected String alnFileSuffix;

    protected IMetircDistrbHolder unifomRandData;
    protected IMetircDistrbHolder yuleRandData;

    protected boolean rooted;
    protected boolean diffLeafSets;

    public boolean isDiffLeafSets() {
        return diffLeafSets;
    }

    public void setDiffLeafSets(boolean diffLeafSets) {
        this.diffLeafSets = diffLeafSets;
    }

    public boolean isRooted() {
        return rooted;
    }

    public void setRooted(boolean rooted) {
        this.rooted = rooted;
    }

    public String getAlnFileSuffix() {
        return alnFileSuffix;
    }

    public void setAlnFileSuffix(String alnFileSuffix) {
        this.alnFileSuffix = alnFileSuffix;
    }
    
    public String getUnifomFileName() {
        return unifomFileName;
    }

    public void setUnifomFileName(String unifomFileName) {
        this.unifomFileName = unifomFileName;
    }

    public String getYuleFileName() {
        return yuleFileName;
    }

    public void setYuleFileName(String yuleFileName) {
        this.yuleFileName = yuleFileName;
    }


    public IMetircDistrbHolder getUnifomRandData() {
        return unifomRandData;
    }


    public IMetircDistrbHolder getYuleRandData() {
        return yuleRandData;
    }


    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }
    
    public BaseMetric() {
        this.rooted = false;
        this.diffLeafSets = false;
    }

    public String getCommandLineName() {
        return commandLineName;
    }

    public void setCommandLineName(String commandLineName) {
        this.commandLineName = commandLineName;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public abstract double getDistance(Tree t1, Tree t2) ;
    public  AlignInfo getAlignment(){
        return null;
    }
    
    private IMetircDistrbHolder parseData(String dataDir, String dataFileName){
        
        MetircDistrbHolder mdh = new MetircDistrbHolder();

        String fullPath = dataDir + "/" + dataFileName;

        try {
        FileReader fr = new FileReader(fullPath);
        BufferedReader br = new BufferedReader(fr);
        int cnt = 0;
        String line;
        while ((line = br.readLine()) != null){
          cnt++;
          if(cnt > 1){
              MetricDistribution md = new MetricDistribution();
              md.readData(line);
              mdh.insertDistribution(md);
          }
        }
        br.close();
      }
      catch (IOException e) {
          Logger.getLogger(ConfigSettings.class.getName()).log(Level.SEVERE, "Error while reading data file:" + fullPath, e);
          mdh = null;
      }
        return mdh;
    }

    public void initData(){
        ConfigSettings config = ConfigSettings.getConfig();
        String dataDir = config.getDataDir();
        if (unifomRandData == null){
            if (unifomFileName != null){
                unifomRandData = parseData(dataDir, unifomFileName);
            }         
        }
        
        if (yuleRandData == null){
            if (yuleFileName != null){
                yuleRandData = parseData(dataDir, yuleFileName);
            }
        }
    }
      
}
