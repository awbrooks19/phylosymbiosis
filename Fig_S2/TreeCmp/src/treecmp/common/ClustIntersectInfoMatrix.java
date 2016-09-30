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

package treecmp.common;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import pal.misc.IdGroup;
import pal.tree.Node;
import pal.tree.Tree;
import pal.tree.TreeUtils;


public class ClustIntersectInfoMatrix {

    protected int intT1Num;
    protected int extT1Num;
    protected int intT2Num;
    protected int extT2Num;

    protected Tree t1;
    protected Tree t2;
    protected IdGroup idGroup;

    protected int[] alias1;
    protected int[] alias2;

    public short [][] intCladeSize;
    public short [] cSize1;
    public short [] cSize2;
    public BitSet [] intT1toLeafT2;
    public BitSet [] intT2toLeafT1;
    public List<ClustPair> eqClustList;
    public boolean [] eqClustT1;
    public boolean [] eqClustT2;

    public ClustIntersectInfoMatrix(Tree t1, Tree t2, IdGroup idGroup) {
           this.t1 = t1;
           this.t2 = t2;
           this.idGroup = idGroup;
    }
    
    public void init(){
        intT1Num = t1.getInternalNodeCount();
        extT1Num = t1.getExternalNodeCount();

        intT2Num = t2.getInternalNodeCount();
        extT2Num = t2.getExternalNodeCount();

        alias1 = TreeUtils.mapExternalIdentifiers(idGroup, t1);
        alias2 = TreeUtils.mapExternalIdentifiers(idGroup, t2);

        intCladeSize = new short[intT1Num][intT2Num];
        intT1toLeafT2 = new BitSet[intT1Num];
        intT2toLeafT1 = new BitSet[intT2Num];

        for(int i = 0; i<intT1Num; i++)
            intT1toLeafT2[i] = new BitSet(extT2Num);
         for(int i = 0; i<intT2Num; i++)
            intT2toLeafT1[i] = new BitSet(extT1Num);

        cSize1 = new short[intT1Num];
        cSize2 = new short[intT2Num];
        
        int minIntNodeNum = Math.min(intT1Num, intT2Num);

        eqClustList = new ArrayList<ClustPair>(minIntNodeNum);
        eqClustT1 = new boolean[intT1Num];
        eqClustT2 = new boolean[intT2Num];

    }

    public int getExtT1Num() {
        return extT1Num;
    }

    public int getExtT2Num() {
        return extT2Num;
    }

    public int getIntT1Num() {
        return intT1Num;
    }

    public int getIntT2Num() {
        return intT2Num;
    }

    public int[] getAlias1() {
        return alias1;
    }

    public int[] getAlias2() {
        return alias2;
    }

    public IdGroup getIdGroup() {
        return idGroup;
    }

    public Tree getT1() {
        return t1;
    }

    public Tree getT2() {
        return t2;
    }
     
    public short getT1Ext_T2Ext(int t1ExtId, int t2ExtId){

        if (alias1[t1ExtId] == alias2[t2ExtId])
            return 1;
        else
            return 0;
    }

    public short getT1Int_T2Ext(int t1IntId, int t2ExtId){

        if (intT1toLeafT2[t1IntId].get(alias2[t2ExtId]))
            return 1;
        else
            return 0;
    }
    
    public short getT1Ext_T2Int(int t1ExtId, int t2IntId){
    
        if (intT2toLeafT1[t2IntId].get(alias1[t1ExtId]))
            return 1;
        else
            return 0;
    }
    
    public short getT1Int_T2Int(int t1IntId, int t2IntId){
    
        return intCladeSize[t1IntId][t2IntId];
    }


    public void setT1Int_T2Ext(int t1IntId, int t2ExtId, short intSize) {

        if (intSize == 1) {
            intT1toLeafT2[t1IntId].set(alias2[t2ExtId]);
        }
    }

    public void setT1Ext_T2Int(int t1ExtId, int t2IntId, short intSize) {

        if (intSize == 1) {
            intT2toLeafT1[t2IntId].set(alias1[t1ExtId]);
        }
    }

    public void setT1Int_T2Int(int t1IntId, int t2IntId, short intSize){

        intCladeSize[t1IntId][t2IntId] = intSize;

        //check if these are the same clusters
        //assume cSize1 and cSize2 have been computed
        if (intSize == cSize1[t1IntId] && intSize == cSize2[t2IntId] ){
            ClustPair cp = new ClustPair();
            cp.t1IntId = t1IntId;
            cp.t2IntId = t2IntId;
            eqClustList.add(cp);
            //mark that the clusters have identical equivalent clusters in the other tree
            eqClustT1[t1IntId] = true;
            eqClustT2[t2IntId] = true;
        }
    }


    public short getInterSize(Node n1, Node n2){

        int n1Num = n1.getNumber();
        int n2Num = n2.getNumber();
        boolean n1Leaf = n1.isLeaf();
        boolean n2Leaf = n2.isLeaf();

        if ((!n1Leaf) && (!n2Leaf))
            return getT1Int_T2Int(n1Num,n2Num);
        
        if (n1Leaf && n2Leaf)
            return getT1Ext_T2Ext(n1Num,n2Num);

         if (n1Leaf && (!n2Leaf))
            return getT1Ext_T2Int(n1Num,n2Num);

         return getT1Int_T2Ext(n1Num,n2Num);
    }

     public short getSizeT1(Node n){
         if (n.isLeaf())
             return 1;
         else
             return cSize1[n.getNumber()];

     }

     public short getSizeT2(Node n){
            if (n.isLeaf())
             return 1;
         else
             return cSize2[n.getNumber()];
     }


public class ClustPair{
    public int t1IntId;
    public int t2IntId;
}

}
