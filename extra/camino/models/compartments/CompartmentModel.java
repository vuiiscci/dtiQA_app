package models.compartments;

import java.util.logging.Logger;

import tools.CL_Initializer;

import data.DataSource;
import data.DataSourceException;
import data.DataSynthesizer;
import numerics.MTRandom;
import numerics.RealMatrix;
import imaging.DW_Scheme;

import misc.LoggedException;
import models.ParametricModel;


/**
 * generalised compartment model for generating data from a weighted
 * sum of two or more compartments.
 * 
 * @author matt
 *
 */
public class CompartmentModel extends ParametricModel implements DataSource{

	/** logging object */
	private final Logger logger= Logger.getLogger(this.getClass().getName());
	
	/** scheme for data synthesis */
	private DW_Scheme scheme= null;
	
    /** array of individual compartments */
    private final ParametricModel[] compartment;
    
    /** array of individual compartment names*/
    private final String[] compartmentnames;
        
    /** arrays for passing model parameters to individual compartments */
    private final double[][] compParams;
    
    /** linearised array for data synthesis (initial params for fitting) */
    private final double[] linParams;
    
    /** table of indices within parameter array corresponding to each compartment */
    private final int[][] index;
    
    /** number of voxels synthesised so far (ignored for fitting) */
    private int voxel=0;
    
    /** constant defining column of start indices in index table */
    private final int START=0;
    
    /** constant defining column of end indices in index table */
    private final int END=1;
    
    /** RNG for adding noise to measurements */
	private MTRandom twister = new MTRandom(CL_Initializer.seed);

    
    /**
     * constructor. needs names of compartments, (initial) model parameters 
     * for each compartment and a set of volume fractions. The volume fractions
     * array should be one element shorter than the array of compartment names
     * because the final compartment size will be assumed to be the difference 
     * between the sum of the sizes and 1. 
     * 
     * If the given vol fracs sum to 1 or greater an error will be generated.
     * 
     * @param comps names of compartments
     * @param params array of vol fracs plus model params for each compartment
     */
    public CompartmentModel(String[] comps, double[] params){
        
        super(CompartmentModel.getNumParamsFor(comps));
        
        // instantiate array space for compartments
        compartment= new ParametricModel[comps.length];
        compParams= new double[comps.length][];
        
        // instantiate compartmentnames
        compartmentnames = comps;
        
        // remember linearised params array
        this.linParams= params;
        
        // instantiate individual compartments
        for(int i=0; i<comps.length; i++){
            // make params array
            compartment[i]= CompartmentFactory.getCompartment(comps[i]);
        }
        
        // assemble indices of entries in the global params array where
        // params to individual compartments are stored
        int numComps= compartment.length;
        
        index= new int[numComps][2];
        
        int count= numComps+1; // first param is 1+numcomps
        for(int i=0; i<numComps; i++){
            
            CompartmentType type= CompartmentType.getCompartmentType(comps[i]);
            
            index[i][START]=count;
            count+=type.numParams;
            index[i][END]= count;
            
            // initialise arrays for parameter passing
            compParams[i]= new double[type.numParams];
        }
    }
    
    /**
     * constructor with scheme for use in data synthesis.
     * 
     * @param comps string specifying compartments
     * @param params linearised model params array
     * @param scheme scheme file to use
     */
    public CompartmentModel(String[] comps, double[] params, DW_Scheme scheme){
    	this(comps, params);
    	
    	this.scheme= scheme;
    }    
    
    private static int getNumParamsFor(String[] comps){
        
        int total = comps.length+1;
        
        for(int i=0; i<comps.length; i++){
            total+= CompartmentType.getCompartmentType(comps[i]).numParams;
        }
        
        return total;
    }
    
    
    /**
     * sum up the signals from each compartment
     */
    public RealMatrix getSignals(double[] modParams, DW_Scheme scheme){
        
        RealMatrix retVals= new RealMatrix(scheme.numMeasurements(), 1);

	double[] normmixpars = new double[compartment.length];
	double nmpsum = 0.0;
        for(int i=0; i<compartment.length-1; i++){
            // volume fraction is the (i+1)th entry.
            nmpsum += modParams[i+1];
            normmixpars[i] = modParams[i+1];
	}
        // Set the last one to ensure they sum to one and check.
        normmixpars[compartment.length-1] = 1-nmpsum;
        if(normmixpars[compartment.length-1]<0 || normmixpars[compartment.length-1]>1) {
            // Actually this happens alot just to compute derivatives
            // numerically to shut off warning.
            //logger.warning("Compartment mixing fractions outside normal bounds: " + normmixpars[compartment.length-1]);
        }
            

        for(int i=0; i<compartment.length; i++){
            
            double[] compParams= new double[compartment[i].numParams()];
            
            for(int j=index[i][START]; j< index[i][END]; j++){
                compParams[j-index[i][START]]= modParams[j];
            }
            
            RealMatrix compSignals= compartment[i].getSignals(compParams, scheme);
            compSignals= compSignals.scalarMult(normmixpars[i]);
            
            retVals= retVals.add(compSignals);
            
        }
        
        // multiply by unweighted measurement
        retVals = retVals.scalarMult(modParams[0]);  
        
        //To print resulting signals on the screen for debugging!
        //for(int i=0; i<retVals.rows(); i++){
        	//System.err.println(retVals.entry(i, 0));
        //}
        
        
        return retVals;
        
    }
    
    /**
     * sum up the signals from each compartment
     */
    public double getSignal(double[] modParams, DW_Scheme scheme, int line){
        
        double retVal= 0.0;

	double[] normmixpars = new double[compartment.length];
	double nmpsum = 0.0;
        for(int i=0; i<compartment.length-1; i++){
            // volume fraction is the (i+1)th entry.
            nmpsum += modParams[i+1];
            normmixpars[i] = modParams[i+1];
	}
        // Set the last one to ensure they sum to one and check.
        normmixpars[compartment.length-1] = 1-nmpsum;
        if(normmixpars[compartment.length-1]<0 || normmixpars[compartment.length-1]>1) {
            // Actually this happens alot just to compute derivatives
            // numerically to shut off warning.
            //logger.warning("Compartment mixing fractions outside normal bounds: " + normmixpars[compartment.length-1]);
        }
            
        for(int i=0; i<compartment.length; i++){
            
            for(int j=index[i][START]; j< index[i][END]; j++){
                compParams[i][j-index[i][START]]= modParams[j];
            }
            
            double compSignal= compartment[i].getSignal(compParams[i], scheme, line);
            
            compSignal*=normmixpars[i];
            
            retVal+=compSignal;
            
        }
        
        return modParams[0]*retVal;
        
    }

    
    /**
     * Get Jacobian for each compartment separately.
     *
     * This could be more efficient often if combined with
     * getSignals, as it has to call getSignals and often
     * computations are repeated.  Future work...
     */
    public RealMatrix getJacobian(final double[] modParams, final DW_Scheme scheme) {
        RealMatrix J = new RealMatrix(scheme.numMeasurements(),modParams.length);
        //double eps= 1e-6;
        //s0
        double s0 = modParams[0];
        RealMatrix s0J = this.getSignals(modParams, scheme);
        for(int r=0; r < scheme.numMeasurements();r++){
        	J.setEntry(r, 0, s0J.entry(r, 0)/s0);
        }
        
        // Determine the mixing parameters for each compartment.
	double[] normmixpars = new double[compartment.length];
	double nmpsum = 0.0;
        for(int i=0; i<compartment.length; i++){
            // volume fraction is the (i+1)th entry.
            nmpsum += modParams[i+1];
	}
        for(int i=0; i<compartment.length; i++){
	    normmixpars[i] = modParams[i+1]/nmpsum;
	}
            
        // Extract the parameters for each compartment
        for(int i=0; i<compartment.length; i++){
            for(int j=index[i][START]; j< index[i][END]; j++){
                compParams[i][j-index[i][START]]= modParams[j];
            }
        }
            

        // Need signals from last compartment, as it contributes
        // to the derivative wrt each earlier mixing parameter, so
        // do it first.
        RealMatrix lastCompJ = compartment[compartment.length-1].getJacobian(compParams[compartment.length-1], scheme);
            
        for(int j=index[compartment.length-1][START]; j< index[compartment.length-1][END]; j++){
            int c_comp = j-index[compartment.length-1][START];
            int c_all = j;
            for(int r=0; r < scheme.numMeasurements();r++){
                //if (f==0){
                //	J.setEntry(r, c_all, eps);
                //}
                //else{
                J.setEntry(r, c_all, (s0*normmixpars[compartment.length-1])*lastCompJ.entry(r, c_comp));
                //}  
            }
        }
        RealMatrix lastCompSignals = compartment[compartment.length-1].getSignals(compParams[compartment.length-1], scheme);


        // Now do the remaining compartments.
        for(int i=0; i<compartment.length-1; i++){
            
            double f = normmixpars[i];
            
            RealMatrix compJ = compartment[i].getJacobian(compParams[i], scheme);
            
            for(int j=index[i][START]; j< index[i][END]; j++){
                int c_comp = j-index[i][START];
                int c_all = j;
                for(int r=0; r < scheme.numMeasurements();r++){
                    //if (f==0){
                    //	J.setEntry(r, c_all, eps);
                    //}
                    //else{
                    J.setEntry(r, c_all, (s0*f)*compJ.entry(r, c_comp));
                    //}  
                 }
                
            }
            
            // volume fraction is the (i+1)th
            RealMatrix fJ = compartment[i].getSignals(compParams[i], scheme);
            
            for(int r=0; r < scheme.numMeasurements();r++){
            	J.setEntry(r, i+1, s0*(fJ.entry(r, 0)-lastCompSignals.entry(r, 0)));
            }  
            
        }
        
        return J;
        
    }


	/** have we done all the voxels we've been asked for? */
	public boolean more() {
		
		return (voxel<CL_Initializer.numVoxels);
	}

	/**
	 * interface for data synthesis. uses the paramters and scheme specified
	 * at Construction to synthesise a voxel of data using the current compartments
	 * 
	 * @return a voxel's worth of data (a "voxeth" in 19th-century parlance)
	 */
	public double[] nextVoxel() throws DataSourceException {
		
		// get data using model fitting interface
		RealMatrix signals= getSignals(this.linParams, scheme);
		
		// convert to rawer format
		double[] signalArray= new double[scheme.numMeasurements()];
		
		for(int i=0; i<signalArray.length; i++){
			signalArray[i]= signals.entry(i, 0);
		}
		
		// increment voxel counter
		voxel++;
		
		// if an SNR is specified, add noise to the outputs
		if(CL_Initializer.SNR>0.0){
			boolean gaussianNoise=false;
			double[] noisyVoxel;
    		
			
			// check noise type
			if(CL_Initializer.noiseType.toLowerCase().equals("gaussian")){
				gaussianNoise=true;
			}
			
			// add noise. noise sigma is 1/SNR (unwieghted signals are 1.0)
    		if(gaussianNoise) {
    			noisyVoxel=DataSynthesizer.addGaussianNoise(signalArray, 1.0/CL_Initializer.SNR, twister);
    		}
    		else {
    			noisyVoxel=DataSynthesizer.addNoise(signalArray, 1.0/CL_Initializer.SNR, twister);
    		}
    		        		
		
			return noisyVoxel;
		}
		
		
		// return the data
		return signalArray;
	}

	/**
	 * @return the compartmentnames
	 */
	public String[] getCompartmentnames() {
		return compartmentnames;
	}
	

}
