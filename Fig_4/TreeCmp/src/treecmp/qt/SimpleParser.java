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

/** This code comes from the QuartetDist application, see
    Chris Christiansen, Thomas Mailund, Christian NS Pedersen,
    Martin Randers and Martin Stig Stissing, "Fast calculation of the quartet distance between trees of arbitrary
    degrees", Algorithms for Molecular Biology, 1:16, 2006.
*/


package treecmp.qt;
/**A simple newick parser.*/
public class SimpleParser {
  private char[] newick;
  private int position, startpos;

  /**Parses the given string and creates a node representing the tree
   *@param newickstr the string to parse
   *@return the node representing the tree
   *@exception ParseException if the string is not in newick format
   */
  public Node parse(String newickstr) throws ParseException {
    newickstr = newickstr.trim();
    if (!newickstr.endsWith(";"))
      newickstr += ";"; //parsing is easier if the string always end with a ';'
    newick = newickstr.toCharArray();
    position = 0;
    Node tree = parseNode();
    checkEnd();
    return tree;
  }

  /**Parses a node*/
  private Node parseNode() throws ParseException {
    if (newick[position] != '(')
      throw new ParseException("Input not in newick format. Missing '(' at character "+position+
			       "in the input string, counting from the first '('.");
    Node current = new InnerNode();
    position++;
    
    while (position < newick.length) {
      if (newick[position] == ')') {
	return current;
      }
      else if (newick[position] == '(') {
	Node tmp = parseNode();
	tmp.addNeighbour(current);
	current.addNeighbour(tmp);
      }
      else if (newick[position] == '\'') {
	Node tmp = parseQuotedLabel();
	tmp.addNeighbour(current);
	current.addNeighbour(tmp);
      }
      else if (newick[position] == '[')
	parseComment();
      else if (newick[position] == ':') 
      	parseBranchLength();
      else if (Character.isWhitespace(newick[position]) || newick[position] == ',');
      else if (newick[position] == ';') {
	if (position != newick.length - 1)
	  throw new ParseException("Input not in newick format. Data appears after ';'.");
      }
      else {
	Node tmp = parseUnQuotedLabel();
	tmp.addNeighbour(current);
	current.addNeighbour(tmp);
      }
      position++;
    }
    return current;
  }

  /**Parses a leaf with a quoted label*/
  private Leaf parseQuotedLabel() throws ParseException { 
    StringBuffer leafname = new StringBuffer();
    while (newick[++position] != '\'') {
      leafname.append(newick[position]);
      if (newick[position+1] == '\'' && newick[position+2] == '\'') {
	leafname.append("''");
	position+=2;
      }
    }
    return new Leaf(leafname.toString());
  }

  /**Parses a leaf with an unquoted label*/
  private Leaf parseUnQuotedLabel() throws ParseException {
    StringBuffer leafname = new StringBuffer();
    while (",):[".indexOf(newick[position]) < 0) {
      if ("(]';".indexOf(newick[position]) >= 0) {
	System.out.println(("(]';".indexOf(newick[position])));
	throw new ParseException("Input not in newick format. Unquoted labels may not "+
				 "contain parentheses, square brackets, "+
				 "single_quotes or semicolons");
      }
      leafname.append(newick[position++]);
    }
    String tmp = leafname.toString().trim();
    if (tmp.indexOf(' ') >= 0 || tmp.indexOf('\t') >= 0 || tmp.indexOf('\n') >= 0)
      throw new ParseException("Input not in newick format. Unquoted labels may not "+
			       "contain blanks.");
    position--; //last read char was not part of the unquoted label
    return new Leaf(tmp);
  }

  /**Parses a comment*/
  private void parseComment() throws ParseException {
    startpos = position;
    while (newick[++position] != ']') {
      if (position == newick.length-1)
	throw new ParseException("Unclosed comment started at "+startpos+".");
    }
    position++;
  }

  /**Parses a branch length*/
  private void parseBranchLength() throws ParseException {
    StringBuffer length = new StringBuffer();
    while (",)[;".indexOf(newick[++position]) < 0)
      length.append(newick[position]);
    try {
      Double.parseDouble(length.toString());
      position--;
    }
    catch (NumberFormatException nfe) {
      throw new ParseException("Branch length '"+length.toString()+"' is not a number.");
    }
  }

  /**Checks the end of the string*/
  private void checkEnd() throws ParseException {
    while (++position < newick.length) {
      if (!Character.isWhitespace(newick[position]))
	if (newick[position] == ':')
	  parseBranchLength();
	else if (newick[position] == ';' && position != newick.length - 1)
	  throw new ParseException("Input not in newick format. Data appears after ';'.");
    }
  }
}
