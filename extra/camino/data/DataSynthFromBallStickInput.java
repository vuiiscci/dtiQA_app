package data;

import imaging.*;
import misc.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Class for synthesizing diffusion weighted data from ball and
 * stick model parameters read from an input stream.
 * 
 * <dt>Description:
 * 
 * <dd>This data source provides synthetic data from each consecutive
 * ball and stick model read from an input stream. The data is
 * synthesized by emulating a specified imaging sequence.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @see inverters.DT_Inversion
 * @version $Id: DataSynthFromDT_Input.java,v 1.4 2005/08/18 10:59:37 ucacmgh
 *          Exp $
 *  
 */
public class DataSynthFromBallStickInput extends DataSynthFromInput {

    public DataSynthFromBallStickInput(String filename, String inputDataType,
            DW_Scheme ip, double s) {
        DATAITEMSPERVOXEL = 7;
        init(filename, inputDataType, ip, s, 0);
    }


    public DataSynthFromBallStickInput(String filename, String inputDataType,
            DW_Scheme ip, double s, int seed) {
        DATAITEMSPERVOXEL = 7;
        init(filename, inputDataType, ip, s, seed);
    }


    protected ModelPDF getNextModel(double[] modelData) {

        double diffusivity = modelData[2];
        double vFrac = modelData[3];
        double[] ori = new double[3];
        for(int i=0; i<3; i++) {
            ori[i] = modelData[4+i];
        }

        return new BallStick(diffusivity, vFrac, ori);

    }

}
