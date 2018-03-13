package models.compartments;

import java.util.Vector;

import numerics.RealMatrix;
import numerics.Vector3D;
import misc.DT;
import misc.LoggedException;
import models.ParametricModel;
import imaging.DW_Scheme;
import imaging.StejskalTannerScheme;


public class Tensor extends ParametricModel {
    
    /** constructor. needs array of params. 
     * in this case, an array of 6 parameters, the diffusivity and the angles theta, phi
     * that determine the fibre orientation, perpendicular diffusivity1, perpendicular diffusivity2 and the angle alpha
     * aplha is the angle by which the perpendicular fibre components are rotated by.
     * @param params params array
     */
    public Tensor(){
        super(CompartmentType.TENSOR.numParams);
    }
    

    /**
     * generates signals from this compartment for each line
     * in the scheme given.
     * 
     * @param scheme scan specifics
     */
    public RealMatrix getSignals(double[] params, DW_Scheme rawScheme) {
        
        StejskalTannerScheme scheme;
        
        try{
            scheme=(StejskalTannerScheme)rawScheme;
        }
        catch(ClassCastException cce){
            throw new LoggedException("Tensor compartment model assumes Stejskal Tanner schemes but was passed something else");
        }
        
        RealMatrix signals= new RealMatrix(scheme.numMeasurements(), 1);
    
        for(int i=0; i<scheme.numMeasurements(); i++){
            signals.setEntry(i, 0, getSignal(params, scheme, i));
        }
        
        return signals;
    }

    
    
    /**
     * generates signals from this compartment for each line
     * in the scheme given.
     * 
     * @param scheme scan specifics
     */
    public double getSignal(double[] params, DW_Scheme rawScheme, int i) {
        
        StejskalTannerScheme scheme;
        
        try{
            scheme=(StejskalTannerScheme)rawScheme;
        }
        catch(ClassCastException cce){
            throw new LoggedException("Tensor compartment model assumes Stejskal Tanner schemes but was passed something else");
        }
        
        return getSignal(params, scheme, i);
    }
    
    
    private double getSignal(double[] params, StejskalTannerScheme scheme, int i){
        
        double theta = params[1];
        double phi = params[2];
        double alpha = params[5];
        
        double sinT = Math.sin(theta);
        double cosT = Math.cos(theta);
        double sinP = Math.sin(phi);
        double cosP = Math.cos(phi);
        
        // fibre orientation
        double[] n = new double[3];
        n[0] = cosP * sinT;
        n[1] = sinP * sinT;
        n[2] = cosT;
        
        Vector3D n1=new Vector3D(n);
                
        //rotate n1 by 90 degrees in xz plane
        Vector3D n2_zrot=Vector3D.vectorFromSPC(1, theta+(Math.PI/2), phi);
        
        //rotate around n1 by alpha
        Vector3D n2 = n2_zrot.rotate(alpha, n1);
        
        Vector3D n3 = n1.cross(n2);
     
        /** Diffusion tensor, nons symmetric */
       
        
//        double[] dnnT1 = new double[6];
//
//        dnnT1[0] = params[0] * n1.x * n1.x;
//        dnnT1[1] = params[0] * n1.x * n1.y;
//        dnnT1[2] = params[0] * n1.x * n1.z;
//        dnnT1[3] = params[0] * n1.y * n1.y;
//        dnnT1[4] = params[0] * n1.y * n1.z;
//        dnnT1[5] = params[0] * n1.z * n1.z;
//
//        double[] dnnT2 = new double[6];
//
//        dnnT2[0] = params[3] * n2.x * n2.x;
//        dnnT2[1] = params[3] * n2.x * n2.y;
//        dnnT2[2] = params[3] * n2.x * n2.z;
//        dnnT2[3] = params[3] * n2.y * n2.y;
//        dnnT2[4] = params[3] * n2.y * n2.z;
//        dnnT2[5] = params[3] * n2.z * n2.z;
//
//        double[] dnnT3 = new double[6];
//
//        dnnT3[0] = params[4] * n3.x * n3.x;
//        dnnT3[1] = params[4] * n3.x * n3.y;
//        dnnT3[2] = params[4] * n3.x * n3.z;
//        dnnT3[3] = params[4] * n3.y * n3.y;
//        dnnT3[4] = params[4] * n3.y * n3.z;
//        dnnT3[5] = params[4] * n3.z * n3.z;
//
//        double[] Dcomparts = new double[6];
//
//        Dcomparts[0] = dnnT1[0] + dnnT2[0] + dnnT3[0];
//        Dcomparts[1] = dnnT1[1] + dnnT2[1] + dnnT3[1];
//        Dcomparts[2] = dnnT1[2] + dnnT2[2] + dnnT3[2];
//        Dcomparts[3] = dnnT1[3] + dnnT2[3] + dnnT3[3];
//        Dcomparts[4] = dnnT1[4] + dnnT2[4] + dnnT3[4];
//        Dcomparts[5] = dnnT1[5] + dnnT2[5] + dnnT3[5];
//
//        // the diffusion tensor
//        DT d = new DT(Dcomparts[0], Dcomparts[1], Dcomparts[2], Dcomparts[3],
//                Dcomparts[4], Dcomparts[5]);
//
//        /** Contract the diffusion tensor by a vector */
//        double qDq = d.contractBy(q);
//
//        return Math.exp(-t * qDq);        
//       
        double b = scheme.getB_Value(i);
   	 double[] gd = scheme.getG_Dir(i);
               
   	 //dot product of the fibre orientations n1,n2,n3 and the wavenumber q
   	 double dotqn1= n1.x*gd[0] + n1.y*gd[1] + n1.z*gd[2];
   	double dotqn2= n2.x*gd[0] + n2.y*gd[1] + n2.z*gd[2];
   	double dotqn3= n3.x*gd[0] + n3.y*gd[1] + n3.z*gd[2];

   	 double diffPar = params[0];
   	 double diffPerp1 = params[3];
   	 double diffPerp2 = params[4];
   	 
	 double tensorSignal = Math.exp(-b*(diffPar*dotqn1*dotqn1+ diffPerp1*dotqn2*dotqn2+ diffPerp2*dotqn3*dotqn3));

     return tensorSignal;   
        
    }
    public RealMatrix getJacobian(final double[] modParams, final DW_Scheme scheme) {
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
        
    	//diffperp1
    	RealMatrix diffprpJ = super.getJacobian(modParams, scheme, 3);
        
        for(int r=0; r < scheme.numMeasurements();r++){
        	jac.setEntry(r, 3, diffprpJ.entry(r, 3));
        } 
        
      //diffperp2
    	RealMatrix diffprp2J = super.getJacobian(modParams, scheme, 4);
        
        for(int r=0; r < scheme.numMeasurements();r++){
        	jac.setEntry(r, 4, diffprp2J.entry(r, 4));
        } 
        
        //alpha
        RealMatrix alphaJ = super.getJacobian(modParams, scheme, 5);
        
        for(int r=0; r < scheme.numMeasurements();r++){
        	jac.setEntry(r, 5, alphaJ.entry(r, 5));
        } 
    	
    	return jac;
    	
    }
}
