package imaging;

/**
 * This interface contains all the methods required to apply the spatial gradient 
 * deviations assoiated with the Human Connectome Project data. It is currently implemented only 
 * by B_Vector schemes and therefore the methods will not work if the wrong type of scheme is 
 * provided. 
 *
 * @author Tara Ganepola
 * @version $Id$ 
 *
 */

public interface HCPScheme {


    /**
     * Gets the i-th gradient direction in a DW_Scheme.
     *
     * @return gDir[i] the i-th gradient direction.
     */
    public double[] getG_Dir(int i);

    /**
     * Gets the b-values for a scheme such as B_VectorScheme.
     *
     * @return b-value array.
     */
    public double[] getB_Values();

    /**
     * Get ith b-value for a scheme such as B_VectorScheme.
     *
     * @return ith b-value array.
     */
    public double getB_Value(int i);


    /**
     * Creates a copy of the scheme.
     *
     * @return a new scheme that is identical to this one.
     */
    public HCPScheme copyScheme();

	/**
    * Resets the values of gDirs, bValues to equal those in source.
    * 
    * @param source a master copy of the original scheme which you want to reset to.
    */
    public void resetScheme(HCPScheme source);

    /**
    * Modifies the gradients and bvals of the scheme per voxel according to 
    * the Human Connectome Project gradiant deviation information. 
    *
    * @param gradAdj specifies the gradient adjustment values.
    */
    public void modifyScheme(double[] gradAdj, HCPScheme source);

}



