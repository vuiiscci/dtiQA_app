package inverters;

import imaging.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Fits the diffusion tensor.
 * 
 * <dt>Description:
 * 
 * <dd>General class for diffusion inversions that fit the single diffusion
 * tensor.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public abstract class DT_Inversion extends DiffusionInversion {


    /**
     * The number of data items per voxel in the output format
     * from this inversion.
     */
    public static final int ITEMSPERVOX = 8;


    /**
     * Returns an indexed DT_Inversion.
     * 
     * @param index
     *            The index for the required inverter. 1 means linear, 
     *            2 means non-linear positive definite, 4 means
     *            non-linear unconstrained, 7 means weighted linear least squares.
     * 
     * @return The inverter.
     */
    public static DT_Inversion getIndexedDT_Inversion(ModelIndex index, DW_Scheme imPars) {
       
	if (index == ModelIndex.LDT || index == ModelIndex.LDT_ALIAS) {
            return new LinearDT_Inversion(imPars);
        }
	else if (index == ModelIndex.LDT_WTD) {
            return new WeightedLinearDT_Inversion(imPars);
        }
	else { // some form of nonlinear inversion
            return new NonLinearDT_Inversion(imPars, index);
        }
    }


    public int itemsPerVoxel() {
        return ITEMSPERVOX;
    }

}
