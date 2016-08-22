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
/**Class representing an inner node in a tree. The node can have an
   arbitrary number of neighbours*/
public class InnerNode extends Node {
  protected int placement;
  private LinkedList edges; //the edges leading to the neighbours
  
  /**Constructs a new inner node without any neighbours to other nodes*/
  public InnerNode() {
    edges = new LinkedList();
  }
  
  /**Adds a new neighbour to this node
     @param neighbour the neighbour to add
     @return the edge resulting from the addition
  */
  protected Edge addNeighbour(Node neighbour) {
    Edge ret = new Edge(this, neighbour);
    edges.add(ret);
    return ret;
  }
  
  /**Specified in Node*/
  protected Node removeNeighbour(Node neighbour) {
    Iterator it = edges.iterator();
    Edge next;
    while(it.hasNext()) {
      next = (Edge)it.next();
      if (next.to == neighbour) {
	it.remove();
	break;
      }
    }
    
    //collapse this node if it has degree 2 after removal of neighbour
    if (edges.size() == 2) {
      Edge e1 = (Edge)edges.getFirst();
      Edge e2 = (Edge)edges.getLast();
      Edge ne1, ne2;
            
      ne1 = e1.to.addNeighbour(e2.to);
      ne2 = e2.to.addNeighbour(e1.to);
      ne1.backedge = ne2;
      ne2.backedge = ne1;
      e1.to.removeNeighbour(this);
      return e2.to.removeNeighbour(this);
    }
    return this;
  }

  /**Specified in Node*/
  public LinkedList getPathTo(Leaf end, Edge from) {
    Iterator it = getEdges();
    Edge next;
    LinkedList path = new LinkedList();
    while (it.hasNext()) { //look through all outgoing edges
      next = (Edge)it.next();
      if (next.backedge != from) { //ignore the outgoing edge pointed back to where we came from
	path = next.to.getPathTo(end, next);
	if (path.size() > 0) { //if the size of the path in one directoin is positive, this is 
	  path.addFirst(from); //the way to go.
	  return path;
	}
      }
    }
    return path;
  }

  
  /**Specified in Node*/
  public void accumulate(Edge caller, LinkedList edgeacc, LinkedList leafacc,
			 LinkedList innernodeacc, BoolHolder binary) {
    Iterator it = edges.iterator();
    Edge next;
    binary.bool = binary.bool && edges.size() == 3;
    innernodeacc.add(this);
    while (it.hasNext()) {
      next = (Edge)it.next();
      edgeacc.add(next);
      if (caller == null || next.to != caller.from) { //recurse on all neighbours that are not caller
	next.to.accumulate(next, edgeacc, leafacc, innernodeacc, binary);
      }
      else { //this edge is the only edge that points back to where we came from
	next.backedge = caller;
	caller.backedge = next;
      }
    }
  }

  /**Specified in Node*/
  public void findLeaves(Edge caller, LinkedList leafacc) {
    Iterator it = edges.iterator();
    Edge next;
    while (it.hasNext()) {
      next = (Edge)it.next();
      if (caller == null || next.to != caller.from) { //recurse on all neighbours that are not caller
	next.to.findLeaves(next, leafacc);
      }
    }
  }

  
  /**Specified in Node*/
  protected Node copy(Node caller) {
    Iterator it = edges.iterator();
    Edge next, e1, e2;
    InnerNode thiscopy = new InnerNode();
    Node tmp;
    while (it.hasNext()) { //Invoke recursively
      next = (Edge)it.next();
      if (next.to != caller) {
	tmp = next.to.copy(this);
	e1 = thiscopy.addNeighbour(tmp);
	e2 = tmp.addNeighbour(thiscopy);
	e1.backedge = e2;
	e2.backedge = e1;
      }
    }
    return thiscopy;
  }

  
  /**Specified in Node*/
  protected String getString(Node caller) {
    String result = "(";
    Iterator it = edges.iterator();
    Edge next;
    while(it.hasNext()) {
      next = (Edge)it.next();
      if (next.to != caller) 
	result += next.to.getString(this) + ",";
    }
    return result.substring(0, result.length() - 1) + ")";
  }

  /**Specified in Node*/
  protected void getSingleDirEdges(Node caller, Collection c) {
    Iterator it = edges.iterator();
    Edge next;
    while(it.hasNext()) {
      next = (Edge)it.next();
      if (next.to != caller) {
	c.add(next);
	next.to.getSingleDirEdges(this,c);
      }      
    }
  }

  /**Specified in Node*/
  public ListIterator getEdges() {
    return edges.listIterator();
  }

  public Edge[] getEdgesAsArray() {
    return (Edge[])edges.toArray(new Edge[0]);
  }

  /**Specified in Node*/
  public int getNumEdges() {
    return edges.size();
  }

  public String toString() {
    return super.toString().substring(super.toString().lastIndexOf('@') + 1);
  }

  public int getPlacement() {
    return placement;
  }
  


  /**Specified in Node*/
  protected Node makeBinary(Collection newedges, Edge[] edgemapping, Node caller) {
    //O(n) each edge is visited 1 time in the first part
    //In the second part an upper bound of the total number
    //of runs through the while-loop is O(n) (number of edges).
    Iterator it = edges.iterator();
    Edge e1, e2, newedge1, newedge2, next;
    InnerNode thiscopy = new InnerNode();
    Node tmp;

    //first part (recursion)                                    
    while (it.hasNext()) { //Invoke recursively    //Like the copy-method
      next = (Edge)it.next();
      if (next.to != caller) {
	tmp = next.to.makeBinary(newedges, edgemapping, this);
	e1 = thiscopy.addNeighbour(tmp);
	e2 = tmp.addNeighbour(thiscopy);
	e1.backedge = e2;
	e2.backedge = e1;
	edgemapping[next.getId()] = e1;
	edgemapping[next.getBackEdge().getId()] = e2;
      }
    }

    //second part (fixing degree with a while-loop)
    while ((thiscopy.edges.size() > 2 && caller != null) || thiscopy.edges.size() > 3) { //not binary
      //System.out.println(thiscopy+" has "+thiscopy.edges.size()+" edges");
      it = thiscopy.edges.iterator();
      e1 = (Edge)it.next();  //1               //Save two outgoing edges (e1 and e2) for binarification
      if (e1.to == caller)
	e1 = (Edge)it.next(); //1 or 2
      it.remove(); //Remove this edge from the edge list, see why below **
      e2 = (Edge)it.next(); //2 or 3
      if (e2.to == caller)
	e2 = (Edge)it.next(); //3  //since edges.size() >=4 I can get three elements out without problems
      it.remove(); //Remove this edge from the edge list, see why below **
      

      //e1 and e2 are two outgoing edges that does not lead back to the caller

      InnerNode newnode = new InnerNode(); //Make a new node to attach e1 and e2 to

      newedge1 = thiscopy.addNeighbour(newnode); //attach this node to the new node
      newedge2 = newnode.addNeighbour(thiscopy); //and vice versa
      newedge1.backedge = newedge2;
      newedge2.backedge = newedge1;
      newedges.add(newedge1); //Since these edges are new edges, add them to the collection
      newedges.add(newedge2); 

      
      //Now instead of using removeNeighbour and that stuff, we use a 'hack'. We change the endpoint of the
      //edges between 'e1' and 'thiscopy' to 'e1' and 'newnode', and 'e2' and 'thiscopy' to 'e2' and 'newnode' respectively.
      //This way any new edge already added to the newedges-collection is still valid.

      e1.from = newnode; //change the endpoints
      e2.from = newnode;
      e1.getBackEdge().to = newnode; //
      e2.getBackEdge().to = newnode;
      
      //now we need to remove e1 and e2 from this nodes edgelist and add them to newnode's edge list, but we
      //already did half of that above **
      newnode.edges.add(e1);
      newnode.edges.add(e2); //the other half well done.
      
      //For each run of this while-loop one neighbour is added and two are removed, hence it must terminate.
    }
  
    return thiscopy;
  }
}
