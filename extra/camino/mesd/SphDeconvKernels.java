package mesd;

import java.util.logging.Logger;

import tools.CL_Initializer;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Container class for deconvolution kernels.
 * 
 * <dt>Description:
 * 
 * <dd>Implements various deconvolution kernels for diffusion MRI. All kernels
 * are accessed through the single public interface function kernel. The type of
 * kernel and its parameter settings are specified in the argument params. The
 * first element of params specifies the type of kernel:
 * 
 * params[0] = 1.0 specifies the PASMRI deconvolution kernel
 * 
 * Any other setting specifies the SPIKE deconvolution kernel.
 * 
 * </dl>
 * 
 * $Id$
 * 
 * @author Danny Alexander
 *  
 */
public class SphDeconvKernels {

    private static Logger logger = Logger.getLogger("mesd.SphDeconvKernels");
    
    /**
     * Deconvolution kernel for perfectly anisotropic diffusion tensor.
     * 
     * @param q
     *            Normalized wavenumber or gradient direction.
     * 
     * @param x
     *            The normalizes displacement or fiber direction.
     * 
     * @param params
     *            Specification of the kernel.
     * 
     * @return The kernel evaluated at q and x.
     */
    public static double kernel(double[] q, double[] x, double[] params) {

        double kernelValue;
        
        // 1 is the identifier for the standard PAS kernel.
        if (params[0] == 1.0) kernelValue= pasKernel(q, x, params[1]);

        else if (params[0] == 2.0) kernelValue= tensorKernel(q, x, params[1]);
        
        // Default is spike deconvolution kernel.
        else kernelValue= spikeKernel(q, x, params[1]);

        if(Double.isNaN(kernelValue)){
            logger.warning("kernel return value with q=("+q[0]+","+q[1]+","+q[2]+") x="
                +x[0]+","+x[1]+","+x[2]+") is "+kernelValue);
        }
        
        return kernelValue;
    }

    /**
     * The spherical wave kernel used in basic PASMRI.
     * 
     * @param q
     *            Normalized wavenumber or gradient direction.
     * 
     * @param x
     *            The normalizes displacement or fiber direction.
     * 
     * @param r
     *            The regularization parameter r.
     * 
     * @return cos(r q.x).
     */
    protected static double pasKernel(double[] q, double[] x, double r) {

        double qdotx = q[0] * x[0] + q[1] * x[1] + q[2] * x[2];

        return Math.cos(r * qdotx);
    }

    /**
     * Kernel assuming perfectly anisotropic diffusion in fibres.
     * 
     * @param q
     *            Normalized wavenumber or gradient direction.
     * 
     * @param x
     *            The normalized displacement or fiber direction.
     * 
     * @param bd
     *            t |q|^2 d, where t is the diffusion time, |q| is the radial
     *            wavenumber and d is the diffusivity in the fibre direction.
     * 
     * @return exp(bd (q.x)^2).
     */
    protected static double spikeKernel(double[] q, double[] x, double bd) {

        double qdotx = q[0] * x[0] + q[1] * x[1] + q[2] * x[2];

        return Math.exp(-bd * qdotx * qdotx);
    }
    
    
    /** 
     * TODO tensor-like kernel for tensor-distribution deconvolution
     * 
     * @param q
     * 	        Normalised wwavenumber gradient direction
     * 
     * @param x 
     * 	        Normalised displacement vector or fibre direction
     * 
     * @param t
     *          Diffusion time
     * 
     */
    protected static double tensorKernel(double[] q, double[] x, double t){
       
        double[][] Dx=getDmatrix(x);
        
        double[] Dxq= new double[3];
        
        for(int i=0; i<3; i++){
            Dxq[i]=Dx[i][0]*q[0] + Dx[i][1]*q[1] + Dx[i][2]*q[2];
        }
        
        double qTDxq = q[0]*Dxq[0] + q[1]*Dxq[1] + q[2]*Dxq[2];
        
        return Math.exp(-t*qTDxq);
        
    }
    
    private static double[][] getDmatrix(double[] x){
        
        double[][] Rx=getRmatrix(x);
        
        double[][] D= new double[3][3];
        
        //double[] l= CL_Initializer.imPars.getTensorResponseEigenVals();
        double[] l = new double[] {0.9, 0.9, 0.1};
        
        double[][] lR=new double[3][3];
        
        // construct eigenVals x Rx 
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                lR[i][j]=l[i]*Rx[i][j];
            }
        }
        
        // left-multiply lR by transpose of Rx
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                D[i][j] = 0.0;
                for(int k=0; k<3; k++){
                    D[i][j]+=lR[i][k]*Rx[j][k];
                }
            }
        }

        return D;
        
    }
    
    
    /** constructs the matrix that rotates the z-axis onto
     * 	the given direction vector x from euler angle representation 
     *  of rotation
     * 
     * @param x 
     *           direction vector to rotate onto 
     * 
     *  @return 
     * 	         3x3 rotation matrix
     */
    private static double[][] getRmatrix(double[] x){
        
        // catch special xases where x is parallel to z axis
        /*if((x[0]==0.0)&&(x[1]==0.0)){
            double[][] R= new double[3][3];
            
            R[0][0]=1.0;
            R[1][1]=1.0;
            if(Math.abs(x[2]-1.0)<=1E-5){
                R[2][2]=1.0;
            }
            else if(Math.abs(x[2]+1.0)<=1E-5){
                R[2][2]=-1.0;
            }
            else{
                throw new RuntimeException("x passed to tensor kernel vector not normalised: x=("+
                        x[0]+","+x[1]+","+x[2]+")");
            }
            
            return R;
            
        }*/
        
        double psi=Math.atan2(x[1], x[0]);
        double cosPsi= Math.cos(psi);
        double sinPsi= Math.sin(psi);
        
        double cosTheta= x[2];
        double sinTheta= Math.sqrt(1.0-cosTheta*cosTheta);
        
        double[][] R=new double[3][3];
        
        R[0][0] =  cosPsi;             R[0][1] =  cosTheta*sinPsi;             R[0][2] = sinTheta*sinPsi;
        R[1][0] = -sinPsi;             R[1][1] =  cosTheta*cosPsi;             R[1][2] = sinTheta*cosPsi;
        R[2][0] =  0.0;                R[2][1] =  0.0;                         R[2][2] = cosTheta;
        
        
        return R;
        
    }
 
    public static void main(String[] args){
        
        
        double[] x= new double[3];
        double[] z= new double[3];
        
        x[0]=1.0/Math.sqrt(2.0);   z[0]=0.0;
        x[1]=1.0/Math.sqrt(2.0);   z[1]=0.0;
        x[2]=0.0;                  z[2]=1.0;
        
        double[][] R=getRmatrix(x);
        
        double[] y= new double[3];
        
        for(int i=0; i<3; i++){
            y[i]=R[i][0]*z[0] + R[i][1]*z[1] + R[i][2]*z[2];
        }
        
        System.err.println("y=("+y[0]+","+y[1]+","+y[2]+")");
        
        
    }
    
}

	