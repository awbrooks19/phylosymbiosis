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

package treecmp.commandline;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.ListIterator;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import treecmp.config.*;
import treecmp.metric.Metric;
import treecmp.command.*;

public class CommandLineParser {

    private final static String S_DESC = "- Overlapping pair comparison mode. Every two neighboring trees are compared";
    private final static String W_DESC = "- Window comparison mode. Every two trees within a window are compared.";
    private final static String W_ARG = "size";
    private final static String M_DESC = "- Matrix comparison mode. Every two trees in the input file are compared.";
    private final static String R_DESC = "- Single tree to all tree mode. Each tree in the input file is compared to the single referenced tree.";
    private final static String R_ARG = "refTreeFile";
    private final static String I_DESC = "- Input file.";
    private final static String I_ARG = "inputFile";
    private final static String O_DESC = "- Output file.";
    private final static String O_ARG = "outputFile";
    private final static String P_DESC = "- Prune compared trees if needed (trees can have different leaf sets).";
    private final static String SS_DESC = "- Report normalized distances.";
    private final static String II_DESC = "- Include summary section in the output file.";
    private final static String A_DESC = "- Generate alignment files (only for MS and MC metrics). Cannot be used with -O option.";
    private final static String OO_DESC = "- Use MS/MC metrics optimized for similar trees. Cannot be used with -A option.";
   // private final static String F_DESC = "- Use MS/MC metrics for trees with free leaf set. ";
    private final static String CMD_ERROR = "Error. There is a problem with parsing the command line. See the usage below.\n";


    private final static String D_DESC = "- Allow to specify distances (from 1 up to 8):\n"+
                                        "Metrics for unrooted trees:\n" +
                                        " ms - the Matching Split metric,\n"+
                                        " rf - the Robinson-Foulds metric,\n"+
                                        " pd - the Path Difference metric,\n"+
                                        " qt - the Quartet metric,\n"+
                                        "Metrics for rooted trees:\n" +
                                        " mc - the Matching Cluster metric,\n"+
                                        " rc - the Robinson-Foulds metric based on clusters,\n" +
                                        " ns - the Nodal Splitted metric with L2 norm,\n"+
                                        " tt - the Triples metric.\n"+
                                        "Example: -d ms rf\n";

    private final static String D_ARG = "metrics";
    private final static String OPTS_HEADER = "Active options:\n";
    private final static String OPTS_TYPE = "Type of the analysis: ";
    private final static String OPTS_METRICS = "Metrics:\n";
    private final static String OPTS_INPUT = "Input file: ";
    private final static String OPTS_OUTPUT = "Output file: ";
    private final static String OPTS_CUSTOM = "Additional options:\n";


    //private final static String HEADER = "";
    //private final static String FOOTER = "ssd";
    private final static String CMD_LINE_SYNTAX = "java -jar TreeCmp.jar -s|-w <size>|-m|-r <refTreeFile>"
            +" -d <metrics> -i <inputFile> -o <outputFile> [-N] [-P] [-I] [-A|-O]\n"
            + "Options order is important.";

    public static Command run(String args[]) {
        Command cmd = null;
        DefinedMetricsSet DMSet = DefinedMetricsSet.getDefinedMetricsSet();
        List<Metric> DMetrics = DMSet.getDefinedMetrics();

        Option oS = new Option("s", S_DESC);
        Option oW = new Option("w", W_DESC);
        oW.setArgName(W_ARG);
        oW.setArgs(1);
        Option oM = new Option("m", M_DESC);
        Option oR = new Option("r", R_DESC);
        oR.setArgName(R_ARG);
        oR.setArgs(1);

        OptionGroup cmdOpts = new OptionGroup();
        cmdOpts.addOption(oS);
        cmdOpts.addOption(oW);
        cmdOpts.addOption(oM);
        cmdOpts.addOption(oR);
        
        cmdOpts.setRequired(true);
        //set metric option
        Option oD = new Option("d", D_DESC);
        oD.setArgName(D_ARG);
        oD.setValueSeparator(' ');
        oD.setArgs(DMetrics.size());
        oD.setRequired(true);
        
        Option oI = new Option("i", I_DESC);
        oI.setArgName(I_ARG);
        oI.setArgs(1);
        oI.setRequired(true);

        Option oO = new Option("o", O_DESC);
        oO.setArgs(1);
        oO.setArgName(O_ARG);
        oO.setRequired(true);

        Option oP = new Option("P", P_DESC);
        Option oSS = new Option("N", SS_DESC);
        Option oII = new Option("I", II_DESC);

        Option oOO = new Option("O", OO_DESC);
      //  Option oF = new Option("F", F_DESC);
        Option oA = new Option("A", A_DESC);
        OptionGroup customMOpts = new OptionGroup();
        customMOpts.addOption(oOO);
        customMOpts.addOption(oA);
      //  customMOpts.addOption(oF);
        
        Options opts = new Options();

        opts.addOptionGroup(cmdOpts);
        opts.addOption(oD);
        opts.addOption(oI);
        opts.addOption(oO);
        opts.addOption(oP);
        opts.addOption(oSS);
        opts.addOption(oII);
        opts.addOptionGroup(customMOpts);


        //getting version from manifest file
        String version = CommandLineParser.class.getPackage().getImplementationVersion();
        if (version == null) {
            version = "";
        }
        //String HEADER="TreeCmp version "+version;
        //String FOOTER = "TreeCmp version "+version;
        String FOOTER = "";
        String HEADER=" ";
        String APP_NAME="TreeCmp version "+version+"\n";
        GnuParser parser = new GnuParser();
        HelpFormatter formatter = new HelpFormatter();
        formatter.setOptionComparator(new OptOrder());

        System.out.println(APP_NAME);
        if(args.length==0){
            formatter.printHelp(CMD_LINE_SYNTAX, HEADER,opts,FOOTER, false);
            return null;                    
        }



        try {
            //parser.checkRequiredOptions();

            CommandLine commandLine = parser.parse(opts, args);
            if (commandLine != null) {
                // process these values
                //set IO settings
                String inputFileName = (String) commandLine.getOptionValue(oI.getOpt());
                String outputFileName = (String) commandLine.getOptionValue(oO.getOpt());
                String [] metrics= commandLine.getOptionValues(oD.getOpt());


                if(inputFileName == null){
                    System.out.println("Error: input file not specified!");
                    formatter.printHelp(CMD_LINE_SYNTAX, HEADER,opts,FOOTER, false);

                    return null;
                }
                if(outputFileName == null){
                    System.out.println("Error: output file not specified!");
                    formatter.printHelp(CMD_LINE_SYNTAX, HEADER,opts,FOOTER, false);
                    return null;
                }

                //commandLine.
                IOSettings IOset = IOSettings.getIOSettings();
                IOset.setInputFile(inputFileName);
                IOset.setOutputFile(outputFileName);

                //custom additinal options
                ArrayList<Option> custOpts = new ArrayList<Option>();


                if (commandLine.hasOption(oP.getOpt())) {
                    IOset.setPruneTrees(true);
                    custOpts.add(oP);
                }
                if (commandLine.hasOption(oSS.getOpt())) {
                    IOset.setRandomComparison(true);
                    custOpts.add(oSS);
                }
                if (commandLine.hasOption(oA.getOpt())) {
                    IOset.setGenAlignments(true);
                    custOpts.add(oA);
                }
                if (commandLine.hasOption(oOO.getOpt())) {
                    IOset.setOptMsMcByRf(true);
                    custOpts.add(oOO);
                }
                if (commandLine.hasOption(oII.getOpt())) {
                    IOset.setGenSummary(true);
                    custOpts.add(oII);
                }
                /*
                if (commandLine.hasOption(oF.getOpt())) {
                    IOset.setUseMsMcFreeLeafSet(true);
                    //additioanly set prune trees
                    IOset.setPruneTrees(true);
                    custOpts.add(oF);
                }
                */
                Collections.sort(custOpts, new OptOrder());
                /*
                if(commandLine.hasOption(oStep))
                {
                String sStep=(String)commandLine.getValue(oStep);
                int iStep=Integer.parseInt(sStep);
                IOset.setIStep(iStep);
                }
                 */
                //set active metrics
                ActiveMetricsSet AMSet = ActiveMetricsSet.getActiveMetricsSet();
                DMSet = DefinedMetricsSet.getDefinedMetricsSet();
                DMetrics = DMSet.getDefinedMetrics();

                for(int i=0;i<metrics.length;i++){

                    ListIterator<Metric> itDM = DMetrics.listIterator();
                    Metric found=null;
                    while (itDM.hasNext()) {
                        Metric m = itDM.next();
                        if(m.getCommandLineName().equals(metrics[i])){
                            found = m;
                        }
                    }
                    if (found != null){
                        AMSet.addMetric(found);
                    }else{
                        System.out.print("Error: ");
                        System.out.println("Metric: "+metrics[i]+" is unknown\n.");
                        formatter.printHelp(CMD_LINE_SYNTAX, HEADER,opts,FOOTER, false);
                        return null;
                    }

                }

                //set active command
                String analysisType="";

                if (commandLine.hasOption(oW.getOpt())) {
                    String sWindowSize = (String) commandLine.getOptionValue(oW.getOpt());
                    //String sWindowSize2 = (String) commandLine.getOptionValue();
                    int iWindowSize = Integer.parseInt(sWindowSize);
                    cmd = new RunWCommand(1, "-w", iWindowSize);
                    analysisType="window comparison mode (-w) with window size: "+iWindowSize;
                }else if (commandLine.hasOption(oM.getOpt())) {
                    cmd = new RunMCommand(0, "-m");
                    analysisType="matrix comparison mode (-m)";
                }else if (commandLine.hasOption(oS.getOpt())) {
                    cmd = new RunSCommand(0, "-s");
                    analysisType="overlapping pair comparison mode (-s)";
                }else if (commandLine.hasOption(oR.getOpt())) {
                    String sRefTreeFile = (String) commandLine.getOptionValue(oR.getOpt());
                    cmd = new RunRCommand(0, "-r",sRefTreeFile);
                    analysisType=" one-to-all comparison mode (-r)";
                }else{
                    System.out.println("Error: type of the analysis not specified correctly!");
                    formatter.printHelp(CMD_LINE_SYNTAX, HEADER,opts,FOOTER, false);
                    return null;
                }

                printOptionsInEffect(analysisType,AMSet,inputFileName,outputFileName, custOpts);

                return cmd;
            } else {
                //Error during parsing command line
                return null;
            }
        } catch (ParseException ex) {
            System.out.println(CMD_ERROR);
            //System.out.println(ex.getMessage());
           
            formatter.printHelp(CMD_LINE_SYNTAX, HEADER,opts,FOOTER, false);

        } catch (NumberFormatException ex){
              System.out.print("Error: ");
              System.out.println("window size should be an integer.\n");
            formatter.printHelp(CMD_LINE_SYNTAX, HEADER,opts,FOOTER, false);


        }
        return cmd;
    }
    private static void printOptionsInEffect(String analysisType,ActiveMetricsSet AMSet,String inputFileName,String outputFileName, List<Option> custOpts){
        System.out.print(OPTS_HEADER);
        System.out.print(OPTS_TYPE+analysisType+"\n");
        System.out.print(OPTS_METRICS);
        Metric[] metrics=AMSet.getActiveMetricsTable();
        int nr;
        Metric m;
        for(int i=0;i<metrics.length;i++){
            m=metrics[i];
            nr=i+1;
            System.out.print("  "+nr+". "+m.getName()+" ("+m.getCommandLineName()+")\n");
        }
   
        System.out.print(OPTS_INPUT+inputFileName+"\n");
        System.out.print(OPTS_OUTPUT+outputFileName+"\n");
        if (!custOpts.isEmpty()){
            System.out.print(OPTS_CUSTOM);
            for (Option opt: custOpts){
                String optMsg = opt.getOpt() +" " + opt.getDescription()+"\n";
                 System.out.print(optMsg);
            }
        }
        System.out.print("-----\n");


    }
}

class OptOrder implements Comparator {

    private LinkedHashMap<String, Integer> order = new LinkedHashMap<String, Integer>();

    public OptOrder() {
        order.put("s", new Integer(1));
        order.put("w", new Integer(2));
        order.put("m", new Integer(3));
        order.put("r", new Integer(4));
        order.put("d", new Integer(5));
        order.put("i", new Integer(6));
        order.put("o", new Integer(7));
        order.put("N", new Integer(8));
        order.put("P", new Integer(9));
        order.put("I", new Integer(10));
        order.put("A", new Integer(11));
        order.put("O", new Integer(12));
       // order.put("F", new Integer(13));
    }

    public int compare(Object o1, Object o2) {
        Option opt1 = (Option) o1;
        Option opt2 = (Option) o2;
        Integer n1 = (Integer) order.get(opt1.getOpt());
        Integer n2 = (Integer) order.get(opt2.getOpt());
        if (n1 != null || n2 != null) {
            return n1 - n2;
        } else {
            return 0;
        }
    }
}
