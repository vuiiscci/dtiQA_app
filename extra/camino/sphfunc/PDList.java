package sphfunc;

import java.util.Vector;

/**
 * <dl>
 * <dt>Purpose: Stores a list of PDs.
 * <dd>
 * 
 * <dt>Description:
 * <dd>This class contains a list of principal directions that can be added to
 * or retrieved from. When a new direction is added a test is made to see if it
 * is close to an existing entry in the list and if so no additional entry is
 * made.
 * 
 * </dl>
 * 
 * @author Danny Alexander $Id: PDList.java,v 1.3 2005/08/18 11:16:06 ucacmgh
 *         Exp $
 *  
 */
public class PDList {

    private Vector<PD> pds;

    private double dpThresh;

    private double TINY = 1.0E-6;

    /**
     * Constructs a PDList object given a threshold for the difference of dot
     * products from 1 below which PDs are considered the same.
     */
    public PDList(double dpt) {
	pds = new Vector<PD>();
        dpThresh = dpt;
    }

    /**
     * Constructs a PDList object with a default threshold of zero, which means
     * that all added PDs remain in the list.
     */
    public PDList() {
        pds = new Vector<PD>();
        dpThresh = TINY;
    }

    /**
     * Adds the new PD to the list provided the list does not already contain an
     * entry close to the new PD.
     */
    public void addPD(PD p) {
        boolean add = true;
        for (int i = 0; i < pds.size(); i++) {
            add &= notCloseTo(p, (PD) pds.elementAt(i));
        }
        if (add) {
            // Add the new one so that they are ordered by their
            // strength.
            int addPos = 0;
            while (addPos < pds.size()
                    && p.getProp() < ((PD) pds.elementAt(addPos)).getProp()) {
                addPos = addPos + 1;
            }
            pds.insertElementAt(p, addPos);
        }
    }

    /**
     * Tests to see if the PD is already in the list to within the specified
     * threshold.
     */
    public boolean isIn(PD p) {
        boolean notIn = true;
        for (int i = 0; i < pds.size(); i++) {
            notIn &= notCloseTo(p, pds.elementAt(i));
        }

        return !notIn;
    }

    /**
     * Returns the i-th PD in the list.
     */
    public PD getPD(int i) {
        return pds.elementAt(i);
    }

    /**
     * Returns the number of principal directions stored.
     */
    public int getNoPDs() {
        return pds.size();
    }

    /**
     * Prunes from the list of PDs all of those whose strength is less than the
     * specified fraction of the sum of the strengths of all the PDs in the
     * list. The process is repeated until it has no effect on the list of PDs.
     */
    public void prune(double thresh) {
        boolean done = false;
        while (!done) {
            done = true;
            double sum = getPropSum();
            for (int i = pds.size() - 1; i >= 0; i--) {
                if (pds.elementAt(i).getProp() < thresh * sum) {
                    pds.removeElementAt(i);
                    done = false;
                }
            }
        }
    }

    /**
     * Prunes from the list of PDs by removing all those with value less than
     * the threshold.
     */
    public void pruneByValue(double thresh) {
        for (int i = pds.size() - 1; i >= 0; i--) {
            if (((PD) pds.elementAt(i)).getProp() < thresh) {
                pds.removeElementAt(i);
            }
        }
    }

    /**
     * Returns the sum of the strengths of all the pds in the list.
     */
    public double getPropSum() {
        double sum = 0.0;
        for (int i = 0; i < pds.size(); i++) {
            sum += pds.elementAt(i).getProp();
        }

        return sum;
    }

    /**
     * Compares this PDList to another to see if they are the same to within the
     * specified tolerance on the dot product between unit vectors.
     * 
     * @param pd
     *            The other PDList to compare.
     * 
     * @param thresh
     *            The threshold on the dot product between two unit vectors
     *            above which they are considered equivalent.
     * 
     * @return The result of the test.
     */
    public boolean equivalent(PDList pd, double thresh) {

        boolean b = (getNoPDs() == pd.getNoPDs());

        int i = 0;
        while (b && i < getNoPDs()) {
            b &= pd.isIn(getPD(i));
            i++;
        }

        return b;
    }

    /**
     * Tests two PDs to see if their dot product is within a threshold of unity.
     */
    private boolean notCloseTo(PD p1, PD p2) {

        //Compute dot product.
        double dp = p1.getPDX() * p2.getPDX() + p1.getPDY() * p2.getPDY() + p1.getPDZ()
                * p2.getPDZ();

        boolean t = (1.0 - Math.abs(dp) > dpThresh);

        // 	if(!t) {
        // 	    System.err.println(p1);
        // 	    System.err.println(p2);
        // 	    System.err.println(dp);
        // 	    System.err.println(dpThresh + "\n");
        // 	}

        //If they are the same, the dot product is one.
        return (1.0 - Math.abs(dp) > dpThresh);
    }

}
