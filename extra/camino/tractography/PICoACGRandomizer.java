package tractography;

import data.*;
import misc.DT;
import numerics.*;
import java.util.Random;

/**
 *
 * Provides samples from the ACG fibre orientation PDF. There may be multiple fibres
 * in each voxel, and hence multiple PDFs. 
 *
 * @version $Id $
 * @author  Philip Cook
 * 
 * 
 */
public class PICoACGRandomizer extends SimplePICoRandomizer {

    private final PICoTractographyImage image;
    private final Random ran;

   
    
    protected PICoACGRandomizer(PICoTractographyImage im, Random r) {
	super(im.xDataDim(), im.yDataDim(), im.zDataDim());
	image = im;
	ran = r;
    }



    protected AxialDistribution[] getPDFs(int i, int j, int k) {

	double[] params = image.getPICoPDFParams(i,j,k);

	Vector3D[] vecs = image.getEigenvectors(i,j,k);

	ACG_Distribution[] pdfs = new ACG_Distribution[vecs.length / 3];

        for (int p = 0; p < pdfs.length; p++) {
            
            Vector3D[] pdfVecs = new Vector3D[3];
            
            pdfVecs[0] = vecs[3*p];
            pdfVecs[1] = vecs[3*p+1];
            pdfVecs[2] = vecs[3*p+2];
            
	    pdfs[p] = new ACG_Distribution(pdfVecs, new double[] 
                {params[p*0], params[p*3 + 1], params[p*3 + 2]}, ran);
	}

	return pdfs;
    }



}
