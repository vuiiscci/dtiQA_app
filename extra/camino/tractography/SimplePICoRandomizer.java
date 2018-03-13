package tractography;

import numerics.*;

/**
 *
 * Provides samples from the fibre orientation PDF. There may be multiple fibres
 * in each voxel, and hence multiple PDFs. 
 *
 * @version $Id $
 * @author  Philip Cook
 * 
 * 
 */
public abstract class SimplePICoRandomizer implements PICoRandomizer {

    private final AxialDistribution[][][][] pdfs;


    /**
     * Construct by telling the data dimensions.
     *
     */
    protected SimplePICoRandomizer(int xDataDim, int yDataDim, int zDataDim) {
	pdfs = new AxialDistribution[xDataDim][yDataDim][zDataDim][];
    }
    
     
    public final Vector3D[] getRandomizedPDs(int i, int j, int k) {
	
	if (pdfs[i][j][k] == null) {
	    pdfs[i][j][k] = getPDFs(i,j,k);
	}
	
	Vector3D[] pds = new Vector3D[pdfs[i][j][k].length];
	
	for (int p = 0; p < pdfs[i][j][k].length; p++) {
	    pds[p] = pdfs[i][j][k][p].nextVector();
	}
	
	return pds;
	
    }


   
    public final double pdf(int i, int j, int k, int pdIndex, Vector3D v) {
	if (pdfs[i][j][k] == null) {
	    pdfs[i][j][k] = getPDFs(i,j,k);
	}
	return pdfs[i][j][k][pdIndex].pdf(v);
    }


    /**
     * Get the PDF(s) for this voxel. No caching is necessary, since this will be called 
     * once by <code>SimplePICoRandomizer</code>.
     */
    protected abstract AxialDistribution[] getPDFs(int i, int j, int k);



}
