/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package treecmp.statdata;

/**
 *
 * @author Damian
 */
public interface IMetircDistrbHolder {
           IMetricDistribution getDistribution(int n);
           int getMaxLeafNum(int n);
           int getMinLeafNum(int n);
}
