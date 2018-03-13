package apps;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import javax.swing.event.*;

import data.*;
import numerics.*;
import tools.*;

/**
 *
 * Simple tool to visualize principal directions of tensor data. Does nothing but extend 
 * PD_OrientationViewer.
 * 
 * @author Philip Cook
 * @version $Id $
 *  
 */
public class TensorOrientationViewer extends PD_OrientationViewer {

   
    /**
     * Creates and displays the viewer.
     */
    public TensorOrientationViewer(RGB_ScalarImage im, double scalarThresh) {
        
        super(im, scalarThresh);
	
        display();

    }



}

