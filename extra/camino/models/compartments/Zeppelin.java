package models.compartments;

import numerics.RealMatrix;
import numerics.Vector3D;
import misc.DT;
import misc.LoggedException;
import models.ParametricModel;
import imaging.DW_Scheme;
import imaging.StejskalTannerScheme;


public class Zeppelin extends ParametricModel {

    /** constructor. needs array of params. 
     * in this case, an array of 4 parameters, the diffusivity and the angles theta and phi
     *  and perpendicular diffusivity.
     * @param params params array
     */
    public Zeppelin(){
        super(CompartmentType.ZEPPELIN.numParams);        
    }
    

    /**
     * generates signals from this compartment for each line
     * in the scheme given.
     * 
     * @param scheme scan specifics
     */
    public RealMatrix getSignals(double[] params, DW_Scheme scheme) {
        
        RealMatrix signals= new RealMatrix(scheme.numMeasurements(),1);
    
        for(int i=0; i<scheme.numMeasurements(); i++){
            signals.setEntry(i, 0, getSignal(params, scheme, i));
        }
        
        return signals;
    }

    public double getSignal(double[] params, DW_Scheme  scheme, int i){
        
	double theta = params[1];
        double phi = params[2];
        
        double sinT = Math.sin(theta);
        double cosT = Math.cos(theta);
        double sinP = Math.sin(phi);
        double cosP = Math.cos(phi);
        
        // fibre orientation
        double[] n = new double[3];
        n[0] = cosP * sinT;
        n[1] = sinP * sinT;
        n[2] = cosT;
        
        /*double[] q= scheme.getQ(i);
      //double t= scheme.getDiffusionTime(i);
        double t= scheme.getDELTA(i)-scheme.getDelta(i)/3;
        
        // Diffusion tensor, cylindrically symmetric
        //constructing the compartments of the tensor
        double dcomp1 = (params[0] - params[3]);
        double[] dnnT = new double[6];

        dnnT[0] = dcomp1 * n[0] * n[0];
        dnnT[1] = dcomp1 * n[0] * n[1];
        dnnT[2] = dcomp1 * n[0] * n[2];
        dnnT[3] = dcomp1 * n[1] * n[1];
        dnnT[4] = dcomp1 * n[1] * n[2];
        dnnT[5] = dcomp1 * n[2] * n[2];

        double[] e1 = new double[] { 1.0, 0.0, 0.0 };

        double[] dperpI = new double[6];

        dperpI[0] = params[3] * e1[0];
        dperpI[1] = params[3] * e1[1];
        dperpI[2] = params[3] * e1[1];
        dperpI[3] = params[3] * e1[0];
        dperpI[4] = params[3] * e1[1];
        dperpI[5] = params[3] * e1[0];

        double[] Dcomparts = new double[6];

        Dcomparts[0] = dnnT[0] + dperpI[0];
        Dcomparts[1] = dnnT[1] + dperpI[1];
        Dcomparts[2] = dnnT[2] + dperpI[2];
        Dcomparts[3] = dnnT[3] + dperpI[3];
        Dcomparts[4] = dnnT[4] + dperpI[4];
        Dcomparts[5] = dnnT[5] + dperpI[5];

        //the diffusion tensor
        DT d = new DT(Dcomparts[0], Dcomparts[1], Dcomparts[2], Dcomparts[3],
                Dcomparts[4], Dcomparts[5]);
        
        d.sortedEigenSystem();
        
        
        
        
        double qDq = d.contractBy(q);
         double zeppelin= Math.exp(-t * qDq);
*/

	 double b = scheme.getB_Value(i);
	 double[] gd = scheme.getG_Dir(i);
            
	 //dot product of the fibre orientation n and the wavenumber q
	 double dotqn= n[0]*gd[0] + n[1]*gd[1] + n[2]*gd[2];

	 double diffPar = params[0];
	 double diffPerp = params[3];
	 double zepSignal = Math.exp(-b*((diffPar-diffPerp)*dotqn*dotqn + diffPerp));

        return zepSignal;        
    }
    

    public RealMatrix getJacobian(final double[] modParams, final DW_Scheme scheme) {
        RealMatrix jac = new RealMatrix(scheme.numMeasurements(), modParams.length);
        
        double parDiff = modParams[0];
        double theta = modParams[1];
        double phi = modParams[2];
        double perpDiff = modParams[3];
        
        double sinT = Math.sin(theta);
        double cosT = Math.cos(theta);
        double sinP = Math.sin(phi);
        double cosP = Math.cos(phi);
        
        // fibre orientation
        double[] n = new double[3];
        n[0] = cosP * sinT;
        n[1] = sinP * sinT;
        n[2] = cosT;
        
        for(int i=0; i < scheme.numMeasurements();i++) {

            double b = scheme.getB_Value(i);
            double[] gd = scheme.getG_Dir(i);
            
            //dot product of the fibre orientation n and the wavenumber q
            double dotqn= n[0]*gd[0] + n[1]*gd[1] + n[2]*gd[2];
        
	    double zepSignal = Math.exp(-b*((parDiff-perpDiff)*dotqn*dotqn + perpDiff));

            // deriv wrt parallel diffusivity
            jac.setEntry(i, 0, -b*dotqn*dotqn*zepSignal);

            // deriv wrt theta
            double cfac = -2*b*(parDiff - perpDiff)*dotqn*zepSignal;
            jac.setEntry(i, 1, cfac*(gd[0]*cosT*cosP + gd[1]*cosT*sinP - gd[2]*sinT));

            // deriv wrt phi
            jac.setEntry(i, 2, cfac*(gd[1]*sinT*cosP - gd[0]*sinT*sinP));

	    // deriv wrt perp diffusivity
            jac.setEntry(i, 3, -b*(1-dotqn*dotqn)*zepSignal);
        } 

        return jac;
    }


    /*    public RealMatrix getJacobian(final double[] modParams, final DW_Scheme scheme) {
    	RealMatrix jac = new RealMatrix(scheme.numMeasurements(), modParams.length);
    	
    	//diff
    	RealMatrix diffJ = super.getJacobian(modParams, scheme, 0);
        
        for(int r=0; r < scheme.numMeasurements();r++){
        	jac.setEntry(r, 0, diffJ.entry(r, 0));
        }
        
        
        
        //theta
        double theta = modParams[1];
        
        if (theta < 1e-6)
        {
        	modParams[1] = 1e-6;
        }

        
    	RealMatrix thetaJ = super.getJacobian(modParams, scheme, 1);
        
        for(int r=0; r < scheme.numMeasurements();r++){
        	jac.setEntry(r, 1, thetaJ.entry(r, 1));
        } 
        
        //phi
        RealMatrix phiJ = super.getJacobian(modParams, scheme, 2);
        
        for(int r=0; r < scheme.numMeasurements();r++){
        	jac.setEntry(r, 2, theta > 1e-6 ? phiJ.entry(r, 2) : 1);
        } 
        
    	//diffperp
    	RealMatrix diffprpJ = super.getJacobian(modParams, scheme, 3);
        
        for(int r=0; r < scheme.numMeasurements();r++){
        	jac.setEntry(r, 3, diffprpJ.entry(r, 3));
        } 
    	
    	return jac;
    	
	}*/
    
    
}
