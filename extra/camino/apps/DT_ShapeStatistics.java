package apps;
import data.*;
import misc.*;
import tools.*;

import java.io.*;
import java.util.logging.*;

/**
 * Computes indices of tensor shape. These include standard indices: l1, l2, l3, tr, md, rd, fa, ra, 2dfa.
 * Also linearity, planarity, and isotropy as defined by Westin et al Proc MICCAI 441-452 (1999).
 *
 * @author Philip Cook
 * @version $Id$
 */
public class DT_ShapeStatistics extends Executable{
    
    public DT_ShapeStatistics(String[] args){
        super(args);
    }
    
    /**
     * logging object
     */
    private static Logger logger= Logger.getLogger("camino.apps.DT_ShapeStatistics");
    
   
    
    int maxComponents;
    DataSource input;
    String statistic; // was final
    
    //String stat = "";
    String stat;
    
    public void initDefaultVals(){
        stat = " ";
    }
    
    
    public void initOptions(String[] args){
        CL_Initializer.maxTensorComponents = 1;
        CL_Initializer.inputDataType = "double";
        
        CL_Initializer.CL_init(args);
        
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-stat")) {
                stat = args[i+1];
                CL_Initializer.markAsParsed(i,2);
            }
        }
        
        CL_Initializer.checkParsing(args);
    }
    
    
    public void initVariables(){
        maxComponents = CL_Initializer.maxTensorComponents;
        input = ExternalDataSource.getDataSource(CL_Initializer.inputFile, 12 * maxComponents, CL_Initializer.inputDataType);
        statistic = stat;
    }
    
    
    public void execute(OutputManager om){
        while (input.more()) {
            
            double[] eig = input.nextVoxel();

            double[] output = new double[maxComponents];

            for (int i = 0; i < maxComponents; i++) {
        
                output[i] = computeStat(eig[12*i], eig[12*i+4], eig[12*i+8], statistic);
            
            }
            
            om.output(output);
        }
        
        
        om.close();
    }
    
    
    public static final double computeStat(double l1, double l2, double l3, String statistic) {
        
        double answer = 0.0;
        
        if (statistic.equals("cl")) {
            answer = (l1 - l2) / l1;
        }
        else if (statistic.equals("cp")) {
            answer = (l2 - l3) / l1;
        }
        else if (statistic.equals("cs")) {
            answer = l3 / l1;
        }
        else if (statistic.equals("l1")) {
            answer = l1;
        }
        else if (statistic.equals("l2")) {
            answer = l2;
        }
        else if (statistic.equals("l3")) {
            answer = l3;
        }
        else if (statistic.equals("tr")) {
            answer = l1 + l2 + l3;
        }
        else if (statistic.equals("md")) {
            answer = (l1 + l2 + l3) / 3.0;
        }
        else if (statistic.equals("rd")) {
            answer = (l2 + l3) / 2.0;
        }
        else if (statistic.equals("fa")) {
            
            double lmean = (l1 + l2 + l3) / 3.0;
            
            if (lmean > 0.0) {
                
                answer =
                Math.sqrt(1.5 * ((l1 - lmean) * (l1 - lmean) +
                (l2 - lmean) * (l2 - lmean) +
                (l3 - lmean) * (l3 - lmean)) /
                (l1 * l1 + l2 * l2 + l3 * l3));
            }
            else {
                answer = 0.0;
            }
            
        }
        else if (statistic.equals("ra")) {
            
            double lmean = (l1 + l2 + l3) / 3.0;
            
            if (lmean > 0.0) {
                
                answer =
                Math.sqrt( 1.0 / 3.0 * ((l1 - lmean) * (l1 - lmean) +
                (l2 - lmean) * (l2 - lmean) +
                (l3 - lmean) * (l3 - lmean)) ) / lmean;
                
            }
            else {
                answer = 0.0;
            }
            
        }
        else if (statistic.equals("2dfa")) {
            double lmean = (l2 + l3) / 2.0;
            
            if (lmean > 0.0) {
                
                answer =
                Math.sqrt(2.0 * ((l2 - lmean) * (l2 - lmean) +
                (l3 - lmean) * (l3 - lmean)) /
                (l2 * l2 + l3 * l3));
            }
            else {
                answer = 0.0;
            }
        }
        else {
            throw new LoggedException("unrecognized statistic " + statistic);
        }
        
        return answer;
        
    }
    
    
    
}
