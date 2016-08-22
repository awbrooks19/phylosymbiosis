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

import treecmp.common.TreeCmpException;
import treecmp.io.ResultWriter;
import treecmp.io.TreeReader;



public class Command {

    int paramNumber;
    TreeReader reader;
    ResultWriter out;
    String name;
    int param;
    String errorMsg;
    String args[];


    public String[] getArgs() {
        return args;
    }
    
    public int getParamNumber() {
        return paramNumber;
    }

    public String getName() {
        return name;
    }

    public Command(int paramNumber, String name) {
        this.paramNumber = paramNumber;
        this.name = name;
    }

    public void run() throws TreeCmpException {
        //this method need to be overidden

    }
    

    public boolean init(String[] args) {
        
        this.args=args;
        this.errorMsg = "";

        if (paramNumber > 0 && args.length <= paramNumber) {
            errorMsg = "Required parametr not specified!";
            return false;

        } else if (args.length <= paramNumber + 1) {
            errorMsg = "Infut file not specified!";
            return false;

        } else if (args.length <= paramNumber + 2) {

            if (paramNumber > 0) {
                try {

                    int temp = Integer.parseInt(args[paramNumber]);
                    param = temp;
                } catch (NumberFormatException ex) {

                    errorMsg = "Incorrect window size! Window size set to 2.";
                    return false;
                }
            }

            reader = new TreeReader(args[paramNumber + 1]);

            //no output file specified
            out = new ResultWriter();
            out.isWriteToFile(false);
        } else if (args.length <= paramNumber + 3) {
            if (paramNumber > 0) {
                try {

                    int temp = Integer.parseInt(args[paramNumber]);
                    param = temp;
                } catch (NumberFormatException ex) {

                    errorMsg = "Incorrect window size! Window size set to 2.";
                    return false;
                }
            }

            reader = new TreeReader(args[paramNumber + 1]);

            //output file specified
            out = new ResultWriter();
            out.isWriteToFile(true);
            out.setFileName(args[paramNumber + 2]);

        }


        return true;

    }


    public void setOut(ResultWriter out) {
        this.out = out;
    }

    public void setParam(int param) {
        this.param = param;
    }

    public void setReader(TreeReader reader) {
        this.reader = reader;
    }


    


    public ResultWriter getOut() {
        return out;
    }

    public int getParam() {
        return param;
    }

    public TreeReader getReader() {
        return reader;
    }

    public String getErrorMsg() {
        return errorMsg;
    }

   


}
