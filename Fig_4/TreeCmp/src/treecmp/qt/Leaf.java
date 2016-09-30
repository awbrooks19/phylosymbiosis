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
import java.util.*;
/**Class representing a leaf in a tree*/
public class Leaf extends Node implements Comparable {
  private String name; //the unique name of this leaf
  private Edge e; //the edge leading from this leaf to the tree
  private int id; //the id of this leaf
  
  /**Constructs a new leaf with the given name
     @param name the name
  */
  public Leaf(String name) {
    this.name = name;
  }
  
  /**Specified in Node*/
  protected Edge addNeighbour(Node neighbour) {
    e = new Edge(this, neighbour);
    return e;
  }  

  /**Specified in Node*/
  public void accumulate(Edge caller, LinkedList edgeacc, LinkedList leafacc, LinkedList innernodeacc, BoolHolder binary) {
    if (caller == null)
      e.to.accumulate(null, edgeacc, leafacc, innernodeacc, binary);
    else {
      leafacc.add(this);
      edgeacc.add(e);
      caller.backedge = e;
      e.backedge = caller;
    }
          //OLD
      // leafacc.add(this);
//     if (e == null) //if there is no edge leading to/from this leaf, it is the tree consisting
//       return;      //of one leaf and nothing else.

//     edgeacc.add(e);
//     if (caller == null) //recurse if initial call to this method
//       e.to.accumulate(e, edgeacc, leafacc, innernodeacc, binary);
//     else {
//       caller.backedge = e;
//       e.backedge = caller;
//     }
  }

  /**Specified in Node*/
  public void findLeaves(Edge caller, LinkedList leafacc) {
    leafacc.add(this);
    if (e == null) //if there is no edge leading to/from this leaf, it is the tree consisting
      return;      //of one leaf and nothing else.

    if (caller == null) //recurse if initial call to this method
      e.to.findLeaves(e, leafacc);
  }

  
  /**Specified in Node*/
  public LinkedList getPathTo(Leaf end, Edge from) {
    LinkedList path = new LinkedList();
    if (end.equals(this) && from != null) { //If the path ends here (and does not start here), just add the last edge
      path.add(from); //add the edge pointing to this node
    }
    if (from == null) { //If the path starts here, search
      path = e.to.getPathTo(end, e);
    }
    return path;
  }
  
  /**Specified in Node*/
  public ListIterator getEdges() {
    LinkedList res = new LinkedList();
    res.add(e);
    return res.listIterator();
  }
  
  /**Specified in Node*/
  protected Node copy(Node caller) {
    if (caller == null && e != null) {
      if (e.to instanceof Leaf) {
	Leaf l1 = new Leaf(name);
	Leaf l2 = new Leaf(((Leaf)e.to).name);
	Edge e1, e2;
	e1 = l1.addNeighbour(l2);
	e2 = l2.addNeighbour(l1);
	e1.backedge = e2;
	e2.backedge = e1;
	return l1;
      }
      return e.to.copy(null);
    }
    return new Leaf(name);
  }

  /**Removes this leaf from the tree. If the neighbour of this leaf
     had degree 3 before the removal, it is collapsed. Returns a node
     that is in the new tree.
  */
  protected Node delete() {
    if (e == null)
      return null;
    return e.to.removeNeighbour(this);
  }
  
  /**Specified in Node*/
  protected Node removeNeighbour(Node neighbour) {
    return this;
  }

  /**Specified in Node*/
  protected String getString(Node caller) {
    String tmpname1 = name;
    //The name might need to be quoted
    for (int i=0; i<tmpname1.length(); i++) {
      if (Character.isWhitespace(tmpname1.charAt(i)) ||
	  i+1 < tmpname1.length() && tmpname1.charAt(i) == '\'' && tmpname1.charAt(i+1) == '\'') {
	tmpname1 = "'"+name+"'";
	break;
      }
    }

    if (e != null && caller == null) {
      if (e.to instanceof InnerNode)
	return e.to.getString(null);
      else {
      	String tmpname2 = ((Leaf)e.to).name;
	for (int i=0; i<tmpname2.length(); i++) {
	  if (Character.isWhitespace(tmpname2.charAt(i)) ||
	      i+1 < tmpname2.length() && tmpname2.charAt(i) == '\'' && tmpname2.charAt(i+1) == '\'') {
	    tmpname2 = "'"+tmpname2+"'";
	    break;
	  }
	}
	return "(L"+tmpname1+",L"+tmpname2+")";
      }
    }
    return "L"+tmpname1;
  }

  /**Specified in Node*/
  protected void getSingleDirEdges(Node caller, Collection c) {
    if (e != null && caller == null) {
      if (e.to instanceof InnerNode)
	e.to.getSingleDirEdges(null,c);
      else
	c.add(e);
    }
  }

  /**Specified in Node*/
  public int getNumEdges() {
    if (e != null)
      return 1;
    else
      return 0;
  }
  
  /**Compares this leaf to the given object, which is assumed to be a
     leaf. A ClassCastException will be thrown if this is not the
     case. Comparison is based on the names of the leaves, and just uses
     String's compareTo.
     @param o the leaf to compare this leaf to
  */
  public int compareTo(Object o) {
    return name.compareTo(((Leaf)o).name);
  }

  /**Returns if this leaf is equal to the given object. Returns true
     only of the given object is a leaf, and have the same name as this
     leaf.
     @param o the object to check for equality
  */
  public boolean equals(Object o) {
    return o instanceof Leaf && name.equals(((Leaf)o).name);
  }
 
  /**Returns the name of this leaf
     @return the name
  */
  public String toString() {
    return name;
  }

  /**Returns the id of this leaf
     @return the leaf
  */
  public int getId() {
    return id;
  }

  /**Sets the id of this leaf
     @param id the id
  */
  protected void setId(int id) {
    this.id = id;
  }

  protected Node makeBinary(Collection newedges, Edge[] edgemapping, Node caller) {
    if (caller == null && e != null) {
      if (e.to instanceof Leaf) {
	Leaf l1 = new Leaf(name);
	Leaf l2 = new Leaf(((Leaf)e.to).name);
	Edge e1, e2;
	e1 = l1.addNeighbour(l2);
	e2 = l2.addNeighbour(l1);
	e1.backedge = e2;
	e2.backedge = e1;
	edgemapping[e.getId()] = e1;
	edgemapping[e.getBackEdge().getId()] = e2;
	return l1;
      }
      return e.to.makeBinary(newedges, edgemapping, null);
    }
    return new Leaf(name);
  }

  
}
