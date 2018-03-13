/* SyntheticScan.java created on 28-Nov-2005
 * (simulation)
 * 
 * author: Matt Hall (m.hall@cs.ucl.ac.uk)
 * 
 */
package simulation.measurement;

import simulation.dynamics.Walker;
import simulation.measurement.ScanFactory.ScanType;


/**
 *  Camino fibre reconstruction and tracking toolkit
 * 
 * SyntheticScan (simulation)
 * 
 * the object that constructs the final readings from the
 * simulation.
 * 
 * 
 *
 * @author Matt Hall (m.hall@cs.ucl.ac.uk)
 *
 */
public interface SyntheticScan {

    
    
    /**
     * gets the signals from the summations of the displacements in
     * the walkers.
     * 
     * @return array of signals, zeroth is unweighted, the rest 
     *         correspond to each gradient direction.
     * 
     */
    public double[] getSignals();        
            
    /**
     * gets the signals from the summations of the displacements in
     * the walkers.
     * 
     * @return array of signals, zeroth is unweighted, the rest 
     *         correspond to each gradient direction.
     * 
     */
    public double[] getCompartmentalSignals(boolean intra);        
    
    /**
     * returns the phase shift at a given location at a given time
     * 
     * @param r position in space (walker location)
     * @param t current time
     * @param gradient direction index
     * @param tLast the last time the walker queried the scan
     * 
     * @return dPhi (phase shift due to gradient)
     */
    public double getPhaseShift(Walker walker, double t, int dir, double tLast);
    

    
    //public double getMagnetisationChange(Walker walker, double t, int dir, double tLast);
    
    /**
     * returns number of gradient directions (including unweighted)
     * 
     * @return N+M
     */
    public int getNumMeasurements();
    
    /**
     * updates the scan variables in each timestep
     * 
     * @param t the time
     */
    public void update(int t);    

    /**
     * returns the scan type
     */
    public ScanType getScanType();
    
}
