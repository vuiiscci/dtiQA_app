/* PercCellularLattice.java created on 08-Dec-2005
 * (simulation)
 * 
 * author: Matt Hall (m.hall@cs.ucl.ac.uk)
 * 
 */
package simulation.geometry.substrates;

import java.util.logging.Logger;

import simulation.SimulationParams;

import numerics.MTRandom;

/**
 *  Camino fibre reconstruction and tracking toolkit
 * 
 * PercCellularLattice (simulation)
 * 
 * Initialises a cellular lattice occupied like a using a
 * percolation-like technique. Each site on the lattice
 * is occupied with a fixed probability p_perc, or unoccupied
 * with prob 1-p_perc.
 * 
 * This is a simple way of including the effects of disorder 
 * in the substrate. SUnstrate geometry varies from free (p_perc=0)
 * to isotropically occupied (p_perc=1.0) with the usual transition
 * in largest cluster size, but this is not so important to us at
 * in this context. 
 *
 * @author Matt Hall (m.hall@cs.ucl.ac.uk)
 *
 */
public class PercCellularLattice extends CellularLattice {

    /** logging object */
     private Logger logger= Logger.getLogger(this.getClass().getName());
     
     /** the percolation probability */
     private final double p_perc;
    
    /**
     * @param l
     * @param L
     * @param p_perc
     */
    public PercCellularLattice(SimulationParams simParams) {
        super(simParams);
        
        this.p_perc=SimulationParams.sim_p_perc;
        
        initLattice();
    }

    /** initialises a lattice of cells with each location
     *  being occupied with propbability p_perc or unoccupied
     *  with probability (1-p_perc)
     * @see simulation.geometry.substrates.CellularLattice#initLattice()
     */
    public void initLattice() {
        
        MTRandom percTwister= new MTRandom((274512654<<32)|1734692047);
        
        for(int i=0; i<occupiedLength; i++){
            double p=percTwister.nextDouble();
            
            if(p<p_perc){
                occupied[i]=true;
            }
            else{
                occupied[i]=false;
            }
        }

    }

    public static void main(String[] args) {
    }
}
