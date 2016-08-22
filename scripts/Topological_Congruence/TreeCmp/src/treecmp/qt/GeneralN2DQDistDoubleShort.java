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

public class GeneralN2DQDistDoubleShort extends Distance{
/**Calculates the quartet distance between two general unrooted
     trees. The computation takes O(n^2d) time for trees that are
     well resolved but up to quartic time depending on how unresulved the tree is.
     @param t1 the first tree
     @param t2 the second tree
     @return the quartet distance
  */
  protected long calcDistance(Tree gt1, Tree gt2) {
    if (!gt1.isInitialized() || !gt2.isInitialized())
      throw new RuntimeException("Both trees must be initialized to use this algorithm");

    short[][] intersect = gt1.calcInterSizesShort(gt2);

    double shared = 0, diff = 0;
    InnerNode[] inodes1 = gt1.getInnerNodes();
    InnerNode[] inodes2 = gt2.getInnerNodes();

    for (int i = 0; i < inodes1.length; i++) { //O(#inner nodes in t1)
      Edge[] outedges1 = inodes1[i].getEdgesAsArray();

      for (int j = 0; j < inodes2.length; j++) { //O(#inner nodes in t2)
	Edge[] outedges2 = inodes2[j].getEdgesAsArray();

	double[] c2_in_t1_allout_t2 = new double[inodes1[i].getNumEdges()];
	double[] c2_in_t2_allout_t1 = new double[inodes2[j].getNumEdges()];

	double[] c2_out_t1_allout_t2 = new double[inodes1[i].getNumEdges()];
	double[] c2_out_t2_allout_t1 = new double[inodes2[j].getNumEdges()];

	double c2_allout_t1_allout_t2 = 0;


	double[] fixed_t1_allout_t2 = new double[inodes1[i].getNumEdges()];
	double[] fixed_t2_allout_t1 = new double[inodes2[j].getNumEdges()];

	//Preprocess edges leading from the pair of inner nodes
	for (int pos1 = 0; pos1 < outedges1.length; pos1++) { //O(degree of current node in t1)
	  for (int pos2 = 0; pos2 < outedges2.length; pos2++) { //O(degree of current node in t2)
	    double out_i_out_j = choose2(intersect[outedges1[pos1].getId()][outedges2[pos2].getId()]);
	    double in_i_out_j = choose2(intersect[outedges1[pos1].getBackEdgeId()][outedges2[pos2].getId()]);
	    double out_i_in_j = choose2(intersect[outedges1[pos1].getId()][outedges2[pos2].getBackEdgeId()]);

	    c2_in_t1_allout_t2[pos1] += in_i_out_j;
	    c2_in_t2_allout_t1[pos2] += out_i_in_j;

	    c2_out_t1_allout_t2[pos1] += out_i_out_j;
	    c2_out_t2_allout_t1[pos2] += out_i_out_j;

	    c2_allout_t1_allout_t2 += out_i_out_j;

	    fixed_t1_allout_t2[pos1] +=
	      (double)intersect[outedges1[pos1].getBackEdgeId()][outedges2[pos2].getId()] *
	      (double)intersect[outedges1[pos1].getId()][outedges2[pos2].getId()];

	    fixed_t2_allout_t1[pos2] +=
	      (double)intersect[outedges1[pos1].getId()][outedges2[pos2].getBackEdgeId()] *
	      (double)intersect[outedges1[pos1].getId()][outedges2[pos2].getId()];
	  }
	}

	//More preprocessing, for diff butterflys
	double[][] diff_innersums = new double[inodes1[i].getNumEdges()][inodes1[i].getNumEdges()];

	for (int pos1 = 0; pos1 < outedges1.length; pos1++) { //O(degree of current node in t1)
	  for (int row = 0; row < outedges1.length; row++) { //O(degree of current node in t1)
	    for (int col = 0; col < outedges2.length; col++) { //O(degree of current node in t2)
	      if (pos1 != row) {
		diff_innersums[pos1][row] +=
		  (double)intersect[outedges1[row].getId()][outedges2[col].getId()] *
		  (double)intersect[outedges1[pos1].getId()][outedges2[col].getId()];
	      }
	    }
	  }
	}

	//Now do the actual processing of edges leading to the nodes
	for (int pos1 = 0; pos1 < outedges1.length; pos1++) { //O(degree of current node in t1)
	  for (int pos2 = 0; pos2 < outedges2.length; pos2++) { //O(degree of current node in t2)
	    shared +=
	      choose2(intersect[outedges1[pos1].getId()][outedges2[pos2].getId()]) *
	      (choose2(intersect[outedges1[pos1].getBackEdgeId()][outedges2[pos2].getBackEdgeId()]) -
	       (c2_in_t2_allout_t1[pos2] - choose2(intersect[outedges1[pos1].getId()][outedges2[pos2].getBackEdgeId()])) -
	       (c2_in_t1_allout_t2[pos1] - choose2(intersect[outedges1[pos1].getBackEdgeId()][outedges2[pos2].getId()])) +
	       c2_allout_t1_allout_t2 - c2_out_t1_allout_t2[pos1] - c2_out_t2_allout_t1[pos2] +
	       choose2(intersect[outedges1[pos1].getId()][outedges2[pos2].getId()]));

	    diff +=
	      (double)intersect[outedges1[pos1].getId()][outedges2[pos2].getId()] *
	      ((double)intersect[outedges1[pos1].getBackEdgeId()][outedges2[pos2].getBackEdgeId()] *
	      (double)intersect[outedges1[pos1].getBackEdgeId()][outedges2[pos2].getId()] *
	      (double)intersect[outedges1[pos1].getId()][outedges2[pos2].getBackEdgeId()] -
	      (double)intersect[outedges1[pos1].getBackEdgeId()][outedges2[pos2].getId()] *
	      (fixed_t1_allout_t2[pos1] -
	       (double)intersect[outedges1[pos1].getId()][outedges2[pos2].getId()] *
	       (double)intersect[outedges1[pos1].getBackEdgeId()][outedges2[pos2].getId()]) -
	      (double)intersect[outedges1[pos1].getId()][outedges2[pos2].getBackEdgeId()] *
	      (fixed_t2_allout_t1[pos2] -
	       (double)intersect[outedges1[pos1].getId()][outedges2[pos2].getBackEdgeId()] *
	       (double)intersect[outedges1[pos1].getId()][outedges2[pos2].getId()]));

	    //Add diff butterflys that was deducted twice
	    for (int pos3 = 0; pos3 < outedges1.length; pos3++) { //O(degree of current node in t1) {
	      if (pos1 != pos3) {
		diff +=
		  (double)intersect[outedges1[pos1].getId()][outedges2[pos2].getId()] *
		  (double)intersect[outedges1[pos3].getId()][outedges2[pos2].getId()] *
		  (diff_innersums[pos1][pos3] -
		   (double)intersect[outedges1[pos3].getId()][outedges2[pos2].getId()] *
		   (double)intersect[outedges1[pos1].getId()][outedges2[pos2].getId()]);
	      }
	    }
	  }
	}
      }
    }

    //All shared quartets are counted twice - divide with 2 to get actual number
    shared /= 2;

    //All different quartets are counted four times - divide by 4 to get actual number
    diff /= 4;

    //Calculate BQs O(nd)
    double bq1 = calcBQ(gt1);
    double bq2 = calcBQ(gt2);

    //Use the found values to calculate and return the qdist
    double dist=bq1 + bq2 - 2 * shared - diff;
    Double distD=new Double(dist);
    return distD.longValue();
  }

  /**Computes the number of ways to select two elements from a set of
     size n (also known as 'n choose 2')
     @param n the size of the set
     @return n choose 2
  */
  private double choose2(double n) {
    return (n * (n - 1)) / 2;
  }

  /**Convenience method, returns the two outgoing from the node the
     edge points to that does not lead to where the edge comes from*/
  private Edge[] otherEdges(Edge e) {
    Iterator it = e.pointsTo().getEdges();
    LinkedList ll = new LinkedList();
    while (it.hasNext()) {
      Edge next = (Edge)it.next();
      if (e.getBackEdge() != next)
	ll.add(next);
    }
    if (ll.size() == 0)
      return null; //we are looking at a leaf
    return (Edge[])ll.toArray(new Edge[0]);
  }

  /**Calculates the number of butterfly quartets in the given tree
     @param t the tree
     @return the number of butterflys in the tree
  */
  private double calcBQ(Tree t) {
    Edge[] edges = t.getEdges();
    double res = 0;
    for (int i = 0; i < edges.length; i++) {
      Edge[] other = otherEdges(edges[i]);
      if (other == null)
	continue;
      double sum = 0;
      for (int j = 0; j < other.length; j++)
	sum += choose2(other[j].getSubtreeSize());
      res +=
	choose2(edges[i].getBackEdge().getSubtreeSize()) *
	(choose2(edges[i].getSubtreeSize()) - sum);
    }
    return res / 2;
  }
}
