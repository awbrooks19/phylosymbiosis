/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package treecmp.statdata;

import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author Damian
 */
public class MetircDistrbHolder implements IMetircDistrbHolder{


    //protected IMetricDistribution distributions[];
    private Map<Integer,IMetricDistribution> distMap;
    private int minLeafNum;
    private int maxLeafNum;

    public MetircDistrbHolder(){
        distMap  = new HashMap<Integer,IMetricDistribution>(1000);
        minLeafNum = Integer.MAX_VALUE;
        maxLeafNum = -1;
    }

    public IMetricDistribution getDistribution(int n) {
        return distMap.get(n);
    }

    public int getMaxLeafNum(int n) {
        return maxLeafNum;
    }

    public int getMinLeafNum(int n) {
        return minLeafNum;
    }

    public void insertDistribution(IMetricDistribution distrb){
        int leafNum = distrb.getLeafNum();
        distMap.put(leafNum, distrb);
        if(leafNum < minLeafNum)
            minLeafNum = leafNum;
        if(leafNum > maxLeafNum)
            maxLeafNum = leafNum;
    }
}
