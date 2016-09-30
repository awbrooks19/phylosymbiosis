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

import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import org.w3c.dom.*;
import org.xml.sax.SAXException;
import treecmp.metric.BaseMetric;


public class ConfigSettings {

    private static ConfigSettings config;
    private String configFile;
    private String dataDir;

    protected ConfigSettings(String configFile, String dataDir) {
        this.configFile = configFile;
        this.dataDir = dataDir;
    }

    public static ConfigSettings getConfig() {        
        return config;
    }

    public String getConfigFile() {
        return configFile;
    }

    public String getDataDir() {
        return dataDir;
    }
    
     public static void initConfig(String configFile, String dataDir) {
        if (config == null) {
            config = new ConfigSettings(configFile, dataDir);
            config.readConfigFromFile();
        }
    }


    private void readConfigFromFile() {

        try {

            File xmlFile=new File(configFile);
            DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
            // Use the factory to create a builder
            DocumentBuilder builder = factory.newDocumentBuilder();
            Document doc = builder.parse(xmlFile);
            String className = "";
            String metricName = "";
            String commandLineName = "";
            String metricDesc="";

            String uniformFileName = "";
            String yuleFileName = "";
            String alnFileSuffix = "";
            String rooted = "";
            String diff_leaves = "";
            /**
             * Update defined metric set
             * 
             */
            DefinedMetricsSet DMset = DefinedMetricsSet.getDefinedMetricsSet();

            NodeList list = doc.getElementsByTagName("metric");
            for (int i = 0; i < list.getLength(); i++) {
                // Get element
                Element element = (Element) list.item(i);
                //System.out.println(getTextValue(element, "class"));
                className = getTextValue(element, "class");
                metricName = getTextValue(element, "name");
                commandLineName = getTextValue(element, "command_name");
                metricDesc=getTextValue(element, "description");
                uniformFileName = getTextValue(element, "unif_data");
                yuleFileName = getTextValue(element, "yule_data");
                alnFileSuffix = getTextValue(element, "aln_file_suffix");
                rooted = getTextValue(element, "rooted");
                diff_leaves = getTextValue(element, "diff_leaves");

                if(className!=null) {
                    Class cl = Class.forName(className);
                    //Metric m=(Metric) cl.newInstance();
                    BaseMetric m=(BaseMetric) cl.newInstance();

                    m.setName(metricName);
                    m.setCommandLineName(commandLineName);
                    m.setDescription(metricDesc);
                    m.setUnifomFileName(uniformFileName);
                    m.setYuleFileName(yuleFileName);
                    m.setAlnFileSuffix(alnFileSuffix);
                    if (rooted !=null )
                        if (rooted.equals("true"))
                            m.setRooted(true);
                    if (diff_leaves !=null )
                        if (diff_leaves.equals("true"))
                            m.setDiffLeafSets(true);

                    DMset.addMetric(m);
                }
            }

            //parse statistic section

            list = doc.getElementsByTagName("reporting");
            Element element = (Element) list.item(0);
            String sSep=getTextValue(element, "filed_separator");
            IOSettings IOs=IOSettings.getIOSettings();

            if(sSep.compareTo("tab")==0) {
                IOs.setSSep("\t");
            } else {
                IOs.setSSep(sSep);
            }      

        } catch (SAXException ex) {
            Logger.getLogger(ConfigSettings.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(ConfigSettings.class.getName()).log(Level.SEVERE, null, ex);
        } catch (Exception ex) {
            Logger.getLogger(ConfigSettings.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    /**
     * I take a xml element and the tag name, look for the tag and get
     * the text content
     * i.e for <employee><name>John</name></employee> xml snippet if
     * the Element points to employee node and tagName is 'name' I will return John
     */
    private String getTextValue(Element ele, String tagName) {
        String textVal = null;
        NodeList nl = ele.getElementsByTagName(tagName);
        if (nl != null && nl.getLength() > 0) {
            Element el = (Element) nl.item(0);
            textVal = el.getFirstChild().getNodeValue();
            textVal=textVal.trim();
        }
        return textVal;
    }

}
