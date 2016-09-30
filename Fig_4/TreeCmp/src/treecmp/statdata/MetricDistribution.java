/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package treecmp.statdata;

import java.util.regex.Pattern;

/**
 *
 * @author Damian
 */
public class MetricDistribution implements IMetricDistribution {

    private final static String SPLIT_REGEX = "\t";
    private final static Pattern p = Pattern.compile(SPLIT_REGEX);
    private final static int QUANTILE_NUM = 13;

    private int n;
    private double avg;
    private double std;
    private double min;
    private double max;

    private double quantile[];

    public MetricDistribution(){
        this.n = -1;
        this.avg = -1;
        this.std = -1;
        this.min = -1;
        this.max = -1;
        this.quantile = null;
        quantile = new double[QUANTILE_NUM];
    }

    public void readData(String dataRow){

       String[] fileds = p.split(dataRow);
       n = Integer.parseInt(fileds[0]);
       avg = Double.parseDouble(fileds[1]);
       std = Double.parseDouble(fileds[2]);
       min = Double.parseDouble(fileds[3]);
       max = Double.parseDouble(fileds[4]);

       int index = 5;
       for(int i = 0; i< QUANTILE_NUM; i++){
           quantile[i] = Double.parseDouble(fileds[index]);
           index++;
       }
    }

    public int getLeafNum() {
        return n;
    }

    public double getAvg() {
         return avg;
    }

    public double getStd() {
        return std;
    }

    public double getMin() {
        return min;
    }

    public double getMax() {
        return max;
    }

    public double get02Quantile() {
        return quantile[0];
    }

    public double get05Quantile() {
        return quantile[1];
    }

    public double get10Quantile() {
        return quantile[2];
    }

    public double get20Quantile() {
        return quantile[3];
    }

    public double get30Quantile() {
        return quantile[4];
    }

    public double get40Quantile() {
        return quantile[5];
    }

    public double get50Quantile() {
        return quantile[6];
    }

    public double get60Quantile() {
        return quantile[7];
    }

    public double get70Quantile() {
        return quantile[8];
    }

    public double get80Quantile() {
        return quantile[9];
    }

    public double get90Quantile() {
        return quantile[10];
    }

    public double get95Quantile() {
        return quantile[11];
    }

    public double get97Quantile() {
        return quantile[12];
    }

}
