package tractography;

import data.*;
import misc.DT;
import numerics.*;
import java.util.Random;

/**
 *
 * Provides samples from the Watson fibre (PDF). 
 * There may be multiple fibres in each voxel, and hence multiple PDFs. This class 
 * will currently only work with data containing one or two tensors in each voxel.
 *
 * @version $Id$
 * @author  Philip Cook
 * 
 * 
 */
public class PICoWatsonRandomizer extends SimplePICoRandomizer {

    protected final PICoTractographyImage image;
    protected final Random ran;


    /**
     * @param image the image to randomize.
     *
     */
    protected PICoWatsonRandomizer(PICoTractographyImage im, Random r) {
	super(im.xDataDim(), im.yDataDim(), im.zDataDim());
	image = im;
	ran = r;
	
    }


    protected AxialDistribution[] getPDFs(int i, int j, int k) {
	double[] params = image.getPICoPDFParams(i,j,k);

	Vector3D[] vectors = image.getEigenvectors(i,j,k);

	AxialDistribution[] pdfs = new WatsonDistribution[vectors.length / 3];

	for (int n = 0; n < pdfs.length; n++) {
            if (params[n] >= 0.0) {
                pdfs[n] = new WatsonDistribution(vectors[3*n], params[n], ran);
            }
            else {
                pdfs[n] = new WatsonDistribution(vectors[3*n+2], params[n], ran);
            }
	}

	return pdfs;
    }


    


}
