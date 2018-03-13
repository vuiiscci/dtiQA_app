package tractography;

import numerics.*;
import java.util.List;
import java.util.ArrayList;
import misc.LoggedException;

public class HistBin {
    
    /**
     * contents
     */
    private double eigenvalue1;
    private double eigenvalue2;
    private List <Vector3D> fibreOrientationsList;

    public HistBin(double e1, double e2) {
	eigenvalue1 = e1;
	eigenvalue2 = e2;
	fibreOrientationsList = new ArrayList<Vector3D>();
    }

    public double [] getEigs() {

	double [] eigs = new double [2];
	eigs[0] = eigenvalue1;
	eigs[1] = eigenvalue2;

	return eigs;
    }

    public double [] getLogEigs() {

	double [] eigs = new double [2];
	eigs[0] = Math.log(eigenvalue1);
	eigs[1] = Math.log(eigenvalue2);

	return eigs;
    }
    
    public double getTrace() {
	
	return (eigenvalue1 + eigenvalue2);
    }

    public double getLogTrace() {
	
	return Math.log(eigenvalue1 + eigenvalue2);
    }

    public Vector3D [] getDirsList() {

	if(fibreOrientationsList.size()>0) {
	    Vector3D [] fibreDirsArray = new Vector3D [fibreOrientationsList.size()];
	    fibreOrientationsList.toArray(fibreDirsArray);
	    return fibreDirsArray;
	}
	else
	    return null;
    }

    public int getNumDirs() {
	return fibreOrientationsList.size();
    }

    public void addDirToList(Vector3D direction) {
	if(!fibreOrientationsList.add(direction)){
	    throw new LoggedException("could not add fibre-orientation estimate to histgram!");
	}
    }
}
