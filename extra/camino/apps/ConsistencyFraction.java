package apps;

import java.util.*;
import java.util.logging.Logger;

import tools.CL_Initializer;
import numerics.*;
import inverters.ModelIndex;
import data.*;


/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd> Computes the consistency fraction from multiple-fibre reconstruction
 * outputs. The consistency fraction is the fraction of trials in which the 
 * peaks computed by <code>SphFuncPD_Stats</code> are consistent with a known 
 * set of expected peak directions.
 * 
 * <dt>Description:
 * 
 * <dd>Reads in the output of <code>SphFuncPD_Stats</code> and computes the consistency
 * fraction with a specified test function.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 * @see apps.SphFuncPD_Stats
 */
public class ConsistencyFraction extends MultiFibreReconStats {

    public static void main(String[] args) {

        ConsistencyFraction is = new ConsistencyFraction(args);
    }

    /**
     * Runs the main program given the command line arguments.
     * 
     * @param args
     *            The command line arguments.
     */
    protected ConsistencyFraction(String[] args) {
        super(args);

    }

    /**
     * Computes the consistency fraction, which is the only statistic.
     * 
     * @param props
     *            The array of properties.
     * 
     * @param inversion
     *            The inversion index.
     * 
     * @return The list of statistics: {C}.
     */
    protected double[] computeStats(double[][] props, ModelIndex inversion) {

        // Sort the corresponding components.
        int firstDirIndex = SphFuncPD_Stats.GLOBALSTATS;
        int propsPerComponent = SphFuncPD_Stats.STATSPERPD;
        orderComponentsByDirection(props, CL_Initializer.numPDsIO, propsPerComponent,
                firstDirIndex);

        double[] stats = new double[1];

        // Construct the test function.
        if (CL_Initializer.rotationIndex > 0) {
            StandardTestFunctions.setTransformation(Rotations.randomRotMat(new Random(
                    CL_Initializer.rotationIndex)));
        }
        if(CL_Initializer.dt2rotangle != 0.0) {
            StandardTestFunctions.setDT2RotationAngle(CL_Initializer.dt2rotangle);
        }
        ModelPDF p = StandardTestFunctions.getFunction(CL_Initializer.testFunction);
        double[][] pds = p.getPDs();

        // Now loop over all voxels.
        int numOK = 0;
        for (int i = 0; i < props.length; i++) {

            // Innocent until proven guilty.
            boolean ok = true;

            // If the number of directions is incorrect, guilty.
            if ((int) props[i][2] != pds.length) {
                ok = false;
            }
            else {

                for (int k = 0; k < props[i][2]; k++) {

                    // Compute the index of this direction in
                    // the props array.
                    int col = firstDirIndex + k * propsPerComponent;

                    // Find the maximum absolute dot product.
                    double maxDP = -1.0;
                    for (int j = 0; j < pds.length; j++) {
                        double adp = Math.abs(props[i][col] * pds[j][0]
                                + props[i][col + 1] * pds[j][1] + props[i][col + 2]
                                * pds[j][2]);
                        maxDP = (adp > maxDP) ? adp : maxDP;
                    }

                    //System.err.println(maxDP);

                    // If the direction is not within a threshold of a
                    // true direction, the reconstruction is not consistent.
                    if (maxDP < dpThresh) {
                        ok = false;
                    }
               }
            }

            if (ok) {
                numOK++;
            }

        }

        stats[0] = (double) numOK / (double) props.length;

        return stats;
    }

}
