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

/** This is a modified and optimized version of the original code*/


package treecmp.qt;
import java.util.*;
import java.io.*;
import java.io.StringReader;
/**Class representing a general (not necessarily binary) unrooted tree*/
public class Tree {
  private Node tree; //handle to the tree
  private Leaf[] leaves; //the leaves in this tree, sorted by their names
  private InnerNode[] inodes; //the inner nodes in this tree
  private Edge[] edges; //the directed edges in this tree
  private Edge[] onediredges; //half of the directed edges in this tree, no backedges
  private boolean isbinary;
  private boolean isinitialized = false;
  private int mindeg = Integer.MAX_VALUE, maxdeg = -1;
  private double avgdeg = -1;
  
  /**Constructs a new tree represented by the given nodse (including it's
     neighbors, if any)
     @param tree the node representing the tree
  */
  public Tree(Node tree) {
    this.tree = tree;
    init(tree);
  }

  public Tree(String newickstring) throws ParseException {
    this.tree = new SimpleParser().parse(newickstring);
    init(tree);
  }

  /**Calculates the intersection sizes of all induced rooted subtrees
     of this tree and the given tree. Returns them in an array of arrays
     @param other the tree to calculate intersection sizes against
     @return the intersection sizes
  */
 //  public int[][] calcInterSizes(Tree other) {
//     //initialize table for data
//     int[][] intersections = new int[edges.length][other.edges.length];
//     for (int i = 0; i < edges.length; i++)
//       Arrays.fill(intersections[i], Integer.MIN_VALUE);

//     //fill the table
//     for (int i = 0; i < edges.length; i++)
//       for (int j = 0; j < other.edges.length; j++)
// 	if (intersections[i][j] == Integer.MIN_VALUE)
// 	  edges[i].intSize(other.edges[j], intersections, true);

//     return intersections;
//   }

   /**Calculates the intersection sizes of all induced rooted subtrees
     of this tree and the given tree. Returns them in an array of arrays
     @param other the tree to calculate intersection sizes against
     @return the intersection sizes
  */
  public int[][] calcInterSizes(Tree other) {
    if (tree instanceof Leaf || other.tree instanceof Leaf)
      throw new RuntimeException("Currently only works if the node-handle in each Tree is an inner node");
    
    //initialize table for data
    int[][] intersections = new int[edges.length][other.edges.length];
    for (int i = 0; i < edges.length; i++)
      Arrays.fill(intersections[i], Integer.MIN_VALUE);

    //Fill for all 'down-down' edges
    for (int i = 0; i < onediredges.length; i++)
      for (int j = 0; j < other.onediredges.length; j++)
	if (intersections[onediredges[i].getId()][other.onediredges[j].getId()] == Integer.MIN_VALUE)
	  onediredges[i].intSize(other.onediredges[j], intersections, true);
    
    for (int i = 0; i < onediredges.length; i++)
      for (int j = 0; j < other.onediredges.length; j++) {
	//Fill for all 'up-up' edges
	intersections[onediredges[i].getBackEdgeId()][other.onediredges[j].getBackEdgeId()] =
	  leaves.length -
	  onediredges[i].getSubtreeSize() -
	  other.onediredges[j].getSubtreeSize() +
	  intersections[onediredges[i].getId()][other.onediredges[j].getId()]; 
	//Fill for all 'down-up' edges
      	intersections[onediredges[i].getId()][other.onediredges[j].getBackEdgeId()] =
	  onediredges[i].getSubtreeSize() -
	  intersections[onediredges[i].getId()][other.onediredges[j].getId()]; 
	//Fill for all 'up-down' edges
	intersections[onediredges[i].getBackEdgeId()][other.onediredges[j].getId()] =
	  other.onediredges[j].getSubtreeSize() -
	  intersections[onediredges[i].getId()][other.onediredges[j].getId()]; 
      }

    return intersections;    
  }
  //modified by Damian Bogdanowicz
  public short[][] calcInterSizesShort(Tree other) {
    if (tree instanceof Leaf || other.tree instanceof Leaf)
      throw new RuntimeException("Currently only works if the node-handle in each Tree is an inner node");
    
    //initialize table for data
    short[][] intersections = new short[edges.length][other.edges.length];
    for (int i = 0; i < edges.length; i++)
      Arrays.fill(intersections[i], Short.MIN_VALUE);

    //Fill for all 'down-down' edges
    for (int i = 0; i < onediredges.length; i++)
      for (int j = 0; j < other.onediredges.length; j++)
	if (intersections[onediredges[i].getId()][other.onediredges[j].getId()] == Short.MIN_VALUE)
	  onediredges[i].intSize(other.onediredges[j], intersections, true);
    
    for (int i = 0; i < onediredges.length; i++)
      for (int j = 0; j < other.onediredges.length; j++) {
	//Fill for all 'up-up' edges
	intersections[onediredges[i].getBackEdgeId()][other.onediredges[j].getBackEdgeId()] =
	    (short) ( leaves.length - onediredges[i].getSubtreeSize() - other.onediredges[j].getSubtreeSize() + intersections[onediredges[i].getId()][other.onediredges[j].getId()]);
	//Fill for all 'down-up' edges
      	intersections[onediredges[i].getId()][other.onediredges[j].getBackEdgeId()] =
	  (short)(onediredges[i].getSubtreeSize() -
	  intersections[onediredges[i].getId()][other.onediredges[j].getId()]);
	//Fill for all 'up-down' edges
	intersections[onediredges[i].getBackEdgeId()][other.onediredges[j].getId()] =
	  (short)(other.onediredges[j].getSubtreeSize() -
	  intersections[onediredges[i].getId()][other.onediredges[j].getId()]);
      }

    return intersections;    
  }


  /**Returns whether this tree is a binary tree, i.e. whether all
     internal nodes have degree three
     @return whether this tree is binary
  */
  public boolean isBinary() {
    return isbinary;
  }

  /**Makes a copy of this tree*/
  public Tree copy() {
    return new Tree(tree.copy(null));
  }
  
  /**Returns a new tree that is the same as this one, except that the
     leaves at the given positions have been removed. The tree is
     collapsed, so that it does not contain any internal nodes of
     degree 2.
     @param pos the positions of the leaves to delete, given as a list of integers
  */
  public Tree deleteLeaves(LinkedList pos) {
    //make new copy
    Tree newtree = new Tree(tree.copy(null));
      
    //delete the leaves
    Iterator it = pos.iterator();
    while (it.hasNext())
      newtree.tree = newtree.leaves[((Integer)it.next()).intValue()].delete();
    
    //recompute leaf and edge arrays in new tree - return it
    newtree.init(newtree.tree);
    return newtree.copy(); //There is an error somewhere that requires us to make a copy
  }

  /**Returns a string representation of this tree. This representation
     is in standard Newick format
     @return the string representation in Newick format
  */
  public String toString() {
    return tree.getString(null)+";";
  }

  /**Return the number of leaves in this tree
    @return the number of leaves in this tree*/
  public int numLeaves() {
    return leaves.length;
  }

  /**Return the leaves in this tree
    @return the leaves in this tree*/
  public Leaf[] getLeaves() {
    return leaves;
  }

  /**Return the edges in this tree
    @return the edges in this tree*/
  public Edge[] getEdges() {
    return edges;
  }

  public int numEdges() {
    return edges.length;
  }

  /**Return the inner nodes in this tree
   @return the inner nodes in this tree*/
  public InnerNode[] getInnerNodes() {
    return inodes;
  }

  /**Return the number of inner nodes in this tree
   @return the number of inner nodes in this tree*/
  public int numInnerNodes() {
    return inodes.length;
  }

  /**Return half of the edges in this tree such that none of the returned edges are backedges of another
    @return the edges*/
  public Edge[] getOneDirEdges() {
    return onediredges;
  }

  public boolean isInitialized() {
    return isinitialized;
  }
  
  
  /**Does initialization stuff for this tree
     @param tree the tree
  */
  private void init(Node tree) {
    this.tree = tree;
    
    //accumulate leaves and edges, sort leaves
    LinkedList edgeacc = new LinkedList();
    LinkedList leafacc = new LinkedList();
    LinkedList inodeacc = new LinkedList();
    BoolHolder binary = new BoolHolder();
    tree.accumulate(null, edgeacc, leafacc, inodeacc, binary);
    
    leaves = (Leaf[])leafacc.toArray(new Leaf[0]);
    inodes = (InnerNode[])inodeacc.toArray(new InnerNode[0]);
    edges = (Edge[])edgeacc.toArray(new Edge[0]);
    isbinary = binary.bool;

    for (int i = 0; i < edges.length; i++) {
      edges[i].id = i;
      edges[i].setSize();
    }

    Arrays.sort(leaves);
    for (int i = 0; i < leaves.length; i++)
      leaves[i].setId(i);

    for (int i = 0; i < inodes.length; i++) {
      inodes[i].placement = i;
    }
    
    LinkedList singlediredges = new LinkedList();
    tree.getSingleDirEdges(null, singlediredges);
    onediredges = (Edge[])singlediredges.toArray(new Edge[0]);
    isinitialized = true;
  }

  /**Returns the minimal degree of all the nodes in the tree
   *@return the minimum degree
   * */
  public int getMinDegree() {
    if (mindeg == Integer.MAX_VALUE)
      for (int i=0; i<inodes.length; i++)
	mindeg = Math.min(mindeg, inodes[i].getNumEdges());
    return mindeg;
  }

  /**Returns the maximal degree of all the nodes in the tree
   *@return the maximal degree
   * */
  public int getMaxDegree() {
    if (maxdeg < 0)
      for (int i=0; i<inodes.length; i++)
	maxdeg = Math.max(maxdeg, inodes[i].getNumEdges());
    return maxdeg;
  }

  /**Returns the average degree of all the nodes in the tree
   *@return the average degree
   * */
  public double getAvgDegree() {
    if (avgdeg < 0) {
      int sum = 0;
      for (int i=0; i<inodes.length; i++)
	sum += inodes[i].getNumEdges();
      avgdeg = (double) sum / inodes.length;
    }
    return avgdeg;
  }
  
  /**Make a binary version of this tree. While doing this collect all new edges and create
      a mapping from the id of the edges in the original tree to the corresponding new edge.
      This will make it possible to identify every edge in the binary tree that corresponds to
      an edge in the original tree.
      @param newedges the collection in which to collect the added edges
      @param edgemapping the mapping
      @return the new binary tree*/
  public Tree makeBinary(Collection newedges, Edge[] edgemapping) {
    return new Tree(tree.makeBinary(newedges, edgemapping, null));
  }
  public Tree makeBinary(Collection newedges) {
    return new Tree(tree.makeBinary(newedges, new Edge[edges.length], null));
  }
  public Tree makeBinary() {
    return new Tree(tree.makeBinary(new LinkedList(), new Edge[edges.length], null));
  }
  public Tree makeBinary(Edge[] edgemapping) {
    return new Tree(tree.makeBinary(new LinkedList(), edgemapping, null));
  }
   
  
}
