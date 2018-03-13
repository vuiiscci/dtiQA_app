/* StripedCellularLattice.java created on 28-Nov-2005
 * (simulation)
 * 
 * author: Matt Hall (m.hall@cs.ucl.ac.uk)
 * 
 */
package simulation.geometry.substrates;

import simulation.SimulationParams;


/**
 *  Camino fibre reconstruction and tracking toolkit
 * 
 * StripedCellularLattice (simulation)
 * 
 * @author Matt Hall (m.hall@cs.ucl.ac.uk)
 *
 */
public class StripedCellularLattice extends CellularLattice {

    private final int stripeThickness;
    
    
    /**
     * initialises stripes of occupied cells of given thickness (given
     * in number of cells) stripes will be parallel to the second 
     * coordinate direction.
     * 
     * @param l cell size
     * @param L lattice size
     * @param stripethickness thickness of stripes (and gaps)
     */
    public StripedCellularLattice(SimulationParams simParams) {
        super(simParams);
        this.stripeThickness= SimulationParams.sim_stripethickness;
        
        initLattice();
    }

    /** initialises D dimensional lattice in stripes of proscribed thickness
     * 
     * @see simulation.geometry.substrates.CellularLattice#initLattice()
     */
    public void initLattice() {
        int i,j;

        for(i=0; i<occupiedLength; i++){
          int sum=0;
          int Lpower=1;

          for(j=0; j<D; j++){

            int n_j= (i/Lpower)%L;
            
            // ignore the second coord -- this ensures
            // stripes instead of checkers and means
            // that we still get alternation in the
            // first coord so 1D case is still ok.
            if(j!=1){
      	sum+=n_j;
            }

            Lpower*=L;
          }

          if(((sum/stripeThickness)%2)==1){
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
