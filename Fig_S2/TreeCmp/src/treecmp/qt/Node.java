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
/**Class representing a general node*/
public abstract class Node {
  
  /**Accumulates all the directed edges and leaves in the tree this
     node belongs to. The edges, leaves and inner nodes are
     accumulated in three given lists. The edge from the calling node is
     also given, and should be null in the initial call. A Boolean is
     threaded to carry information whether this tree is binary
     (i.e. all internal nodes have degree three)
     @param caller the calling node, null if it is the initial call
     @param edgeacc the linked list to accumulate the edges in
     @param leafacc the linked list to accumulate the leaves in
     @param innernodeacc the linked list to accumulate the inner nodes in
     @param binary holds whether the tree is binary
  */
  public abstract void accumulate(Edge caller, LinkedList edgeacc, LinkedList leafacc,
				  LinkedList innernodeacc, BoolHolder binary);

  /**Finds all leaves reachable from this node, except the ones that can
    be reached by going opposite the caller edge
    @param caller the caller edge
    @param leafacc the linked list to save the leaves in.*/
  public abstract void findLeaves(Edge caller, LinkedList leafacc);

  
  /**Finds a path from this node to the given leaf by searching in all directions except
     the direction in which the given node lies
     @param end the leaf we want to find a path to
     @param from the edge that lies in the directions in which we don't want to search.
     <code>from</code> would typically be a node on the path so far.
     @return a linked list of the edges on the path. If the list is empty, the leaf
     (<code>end</code>) cannot be found.*/
  public abstract LinkedList getPathTo(Leaf end, Edge from);

  
  /**Adds a new neighbour to this node
     @param neighbour the neighbour to add
     @return the edge resulting from the addition.
  */
  protected abstract Edge addNeighbour(Node neighbour);
  
  /**Removes a neighbour from this node. Returns a node that is in the
     new tree (in case current node ius collapsed)
     @param neighbour the neighbour to remove
     @return a node in the tree
  */
  protected abstract Node removeNeighbour(Node neighbour);

  /**Returns an iterator of the edges leading from this node
     @return the edge-iterator
  */
  public abstract ListIterator getEdges();
  
  /**Returns a copy of the tree this node belongs to. The calling node
     is given, to build the new tree. caller should be null in the initial call.
     @param caller the calling node, null if it is the initial call
  */
  protected abstract Node copy(Node caller);
  
  /**Returns the number of edges leading from this node
     @return the number of edges
  */
  public abstract int getNumEdges();

  /**Returns a string representation of the tree this node
   represents. The calling node is given, and should be null in the
   initial call.
   @param caller the calling node, null if it is the initial call
   @return the string representation
  */
  protected abstract String getString(Node caller);

  /** Accumulate all edges in the tree represented by this node in the given collection.
      None of the edges will be backedges of any of the other edges in the collection.
      @param caller the calling node, null if it is the initial call
      @param c the collection to accumulate the edges in.
  */
  protected abstract void getSingleDirEdges(Node caller, Collection c);  

  /** Make a binary version of this tree, collecting all new edges created in the process and creating
      a mapping from the placement of the edges in the original tree to the corresponding new edge.
      @param newedges a collection for collecting the new edges
      @param edgemapping the mapping
      @param caller the calling node (nescessary for not making infinite recursion)
      @return the binary version of this tree
  */
  protected abstract Node makeBinary(Collection newedges, Edge[] edgemapping, Node caller);

}
