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

package treecmp;

import java.io.UnsupportedEncodingException;
import java.net.URLDecoder;
import java.util.logging.Level;
import java.util.logging.Logger;
import treecmp.common.TreeCmpException;
import treecmp.io.ResultWriter;
import treecmp.io.TreeReader;
import treecmp.command.Command;
import treecmp.commandline.CommandLineParser;
import treecmp.common.TimeDate;
import treecmp.config.ActiveMetricsSet;
import treecmp.config.ConfigSettings;
import treecmp.config.IOSettings;
import treecmp.config.PersistentInfo;
import treecmp.metric.Metric;

public class Main {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
      
       String runtimePathTemp = Main.class.getProtectionDomain().getCodeSource().getLocation().getPath();
       if (runtimePathTemp.indexOf('+') != -1){
           System.out.println("Unsupported character in TreeCmp path!\n"
                   + "Please move TreeCmp application to a location that does not contain plus (+) character.");
           return;
       }

       String runtimePath = runtimePathTemp;
        try {
            runtimePath = URLDecoder.decode(runtimePathTemp, "UTF-8");
        } catch (UnsupportedEncodingException ex) {
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
        }

 
       String conf="";
       String dataDir="";

       String version = Main.class.getPackage().getImplementationVersion();

       if(version==null){
           conf = runtimePath + "../" + PersistentInfo.configFile;
           dataDir = runtimePath + "../" + PersistentInfo.dataPath;

       }else{
           String tempPath = runtimePath.substring(0,runtimePath.lastIndexOf("/")+1);

           conf = tempPath + PersistentInfo.configFile;
           dataDir = tempPath + PersistentInfo.dataPath;
       }

        ConfigSettings.initConfig(conf, dataDir);
        Command cmd=CommandLineParser.run(args);

        if(cmd!=null)
        {
            IOSettings settings =IOSettings.getIOSettings();
            //init data if needed
            if (settings.isRandomComparison()){
                System.out.println(TimeDate.now()+": Start of reading random statistics.");
                for (Metric m: ActiveMetricsSet.getActiveMetricsSet().getActiveMetrics()){
                    m.initData();
                }
                System.out.println(TimeDate.now()+": End of reading random statistics.");
            }

            TreeReader reader = new TreeReader(settings.getInputFile());
            //scanning all content of the input file
            
            if(reader.open()==-1){
                //an error occured during reading the input file
                return;
            }
                
            System.out.println(TimeDate.now()+": Start of scanning input file: " + settings.getInputFile());
            int numberOfTrees = reader.scan();
            reader.close();
            System.out.println(TimeDate.now()+": End of scanning input file: " + settings.getInputFile());
            System.out.println(TimeDate.now()+": "+numberOfTrees+" valid trees found in file: "+settings.getInputFile());

            reader.setStep(settings.getIStep());
            cmd.setReader(reader);

            ResultWriter out = new ResultWriter();
            out.isWriteToFile(true);
            out.setFileName(settings.getOutputFile());
            cmd.setOut(out);
            try {
                cmd.run();
            } catch (TreeCmpException ex) {
                System.out.println(ex.getError());
            }

        }
    }

}
