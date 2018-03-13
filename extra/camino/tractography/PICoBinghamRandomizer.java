package tractography;

import data.*;
import misc.DT;
import numerics.*;
import java.util.Random;

/**
 *
 * Provides samples from the Bingham fibre PDF. 
 * There may be multiple fibres in each voxel, and hence multiple PDFs. This class 
 * will currently only work with data containing one or two tensors in each voxel.
 *
 * @version $Id $
 * @author  Philip Cook
 * 
 * 
 */
public class PICoBinghamRandomizer extends SimplePICoRandomizer {

    private final PICoTractographyImage image;
    private final Random ran;
   

    /**
     * @param image the image to randomize. 
     *
     */
    protected PICoBinghamRandomizer(PICoTractographyImage im, Random r) {
	super(im.xDataDim(), im.yDataDim(), im.zDataDim());
	image = im;
	ran = r;
    }


    
    protected AxialDistribution[] getPDFs(int i, int j, int k) {
	double[] params = image.getPICoPDFParams(i,j,k);

	Vector3D[] vecs = image.getEigenvectors(i,j,k);

	AxialDistribution[] pdfs = new AxialDistribution[vecs.length / 3];

        for (int p = 0; p < pdfs.length; p++) {
	    try {

                Vector3D[] pdfVecs = new Vector3D[3];

                pdfVecs[0] = vecs[3*p];
                pdfVecs[1] = vecs[3*p+1];
                pdfVecs[2] = vecs[3*p+2];

		pdfs[p] = BinghamDistribution.getBinghamDistribution(pdfVecs, new double[] 
		    {params[2*p], params[2*p + 1]}, ran);
	    }
	    catch (ConvergenceException e) {

		double k1 = params[2*p];
		double k2 = params[2*p+1];

		// hack alert - catch case where k1 < k2 = 0 throws ConvergenceException
		if (k1 < k2 && k2 == 0.0) {
		    pdfs[p] = new WatsonDistribution(vecs[3*p+2], k1, ran);
		}
		else {

		    // nothing to do here but quit
		    
		    throw new RuntimeException("Can't get distribution for pdf " + i + ", " + j + ", " + k + " with params " + params[2*p] + " " + params[2*p + 1]);
		}
	    }

	}

	return pdfs;
    }


    

}
