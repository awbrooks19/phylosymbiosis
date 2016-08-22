/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package treecmp.statdata;

/**
 *
 * @author Damian
 */
public interface IMetricDistribution {

     public int getLeafNum();

     public double getAvg();
     public double getStd();

     public double getMin();
     public double getMax();

     public double get02Quantile();
     public double get05Quantile();
     public double get10Quantile();
     public double get20Quantile();
     public double get30Quantile();
     public double get40Quantile();
     public double get50Quantile();
     public double get60Quantile();
     public double get70Quantile();
     public double get80Quantile();
     public double get90Quantile();
     public double get95Quantile();
     public double get97Quantile();

}
