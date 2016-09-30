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
/**An abstract class for algorithms calculating distances between
   trees. Distance methods need only implement the calcDistance
   method, and can assume that the trees contain the same leaves when
   this method is called.
*/
public abstract class Distance {
  
  /**Finds the distance between the two given trees. The distance is
     calculated on the subset of the trees that contain only the
     intersection of the leaves. Basically just thins out the trees
     ofg excess leaves, and then call the calcDistance, which must be
     implemented for specific metrics.
     @param t1 the first tree
     @param t2 the second tree
     @return the distance 
  */
  public long getDistance(Tree t1, Tree t2) {
    if (!t1.isInitialized())
      throw new RuntimeException("Tree t1 has not been initialized");
    if (!t2.isInitialized())
      throw new RuntimeException("Tree t2 has not been initialized");
    try {
      //Check if trees have same leaves;
      LinkedList delete1 = new LinkedList();
      LinkedList delete2 = new LinkedList();
      Tree[] pruned = pruneTrees(t1, t2, delete1, delete2);

      pruned[0] = pruned[0];
      pruned[1] = pruned[1];
      
      if (delete1.isEmpty() && delete2.isEmpty())
	//same leaf set, just do an ordinary calculation
	return calcDistance(t1, t2);

      //Different leaf sets.
      long qdisttop = calcDistance(pruned[0], pruned[1]);        //qdist on pruned
      long qsimtop  = choose4(pruned[0].numLeaves()) - qdisttop; //agree on pruned
      //qdisttop is the number of quartets existing in both trees that have
      //different topology
      //qsimtop  is the number of quartets existing in both trees that have same
      //topology
      return ((choose4(t1.numLeaves()) + choose4(t2.numLeaves())) -
	      2 * qsimtop - qdisttop);
    }
    catch (PruneException pe) {
      //PruneException is thrown when one of the pruned versions of
      //the trees will have less than 4 leaves left. Hence there will
      //be no quartets in common between the trees.
      return choose4(t1.numLeaves()) + choose4(t2.numLeaves());
    }
  }

  public DistResult getMeasures(Tree t1, Tree t2) {
    //Fixme count and prune trees
    //for putting indexes of leaves to delete
    long qdist = 0;
    long qsim = 0;
    int n = 0;
    LinkedList delete1 = new LinkedList();
    LinkedList delete2 = new LinkedList();
    try {
      Tree[] pruned = pruneTrees(t1, t2, delete1, delete2);
      n = pruned[0].numLeaves();
      qdist = calcDistance(pruned[0], pruned[1]);
      qsim  = choose4(n) - qdist;
    }
    catch (PruneException pe) { System.out.println("Pruneexception\n"+pe.getMessage());}
    return new DistResult(qdist, qsim,
			  choose4(n+delete1.size())-choose4(n),
			  choose4(n+delete2.size())-choose4(n));
  }

  private long choose4(long n) {
    return (n * (n - 1) * (n - 2) * (n - 3)) / 24;
  }

  public class PruneException extends Exception {
    public PruneException(String msg) {
      super(msg);
    }
  }

  
  /**Prunes the given trees, so that they contain the same leaves, and
     returns the resulting trees. If leaves are deleted, new trees
     will be created to preserve the old ones. The two resulting trees
     are returned in an array of size 2
     @param t1 the first tree
     @param t2 the second tree
     @return the two resulting trees, in an array of size 2
  */
  public Tree[] pruneTrees(Tree t1, Tree t2, LinkedList delete1, LinkedList delete2)
    throws PruneException {
    Leaf[] leaves1= t1.getLeaves();
    Leaf[] leaves2= t2.getLeaves();

    //run through arrays of leaves, and identify surplus leaves
    int i = 0, j = 0;
    int comp;
    while(i < leaves1.length && j < leaves2.length) {
      comp = leaves1[i].compareTo(leaves2[j]);
      if (comp == 0) { //leaves equal, all ok
	i++;
	j++;	
      }
      else if (comp < 0) { //the leaf in the first tree is not in the second
	delete1.add(new Integer(i++));
      }
      else { //comp > 0 - ie. the leaf in the second tree is not in the first
	delete2.add(new Integer(j++));
      }
    }

    
    //If any leaves left in one of the arrays, they must be deleted
    while(i < leaves1.length)
      delete1.add(new Integer(i++));
    while(j < leaves2.length)
      delete2.add(new Integer(j++));

    if (delete1.size() >= t1.numLeaves()-3 || delete2.size() >= t2.numLeaves()-3)
      throw new PruneException("One of the trees will have no quartets left.");
    //delete leaves, if needed
    if (!delete1.isEmpty())
      t1 = t1.deleteLeaves(delete1);

    if (!delete2.isEmpty())
      t2 = t2.deleteLeaves(delete2);

    
    return new Tree[] {t1, t2};
  }
  
			  
  /**Calculates the distance between the two given trees
     @param t1 the first tree
     @param t2 the second tree
     @return the distance     
  */
  protected abstract long calcDistance(Tree t1, Tree t2);

}
