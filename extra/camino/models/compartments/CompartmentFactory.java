package models.compartments;

import misc.LoggedException;
import models.ParametricModel;

/**
 * holds factory method for generating individual compartments
 * for the compartment model
 * 
 * @author laura
 *
 */
public class CompartmentFactory {

    public static ParametricModel getCompartment(String compartmentName){

        if(compartmentName.equalsIgnoreCase("ball")){

            return new Ball();
        }
        else if(compartmentName.equalsIgnoreCase("stick")){

            return new Stick();
        }
        else if(compartmentName.equalsIgnoreCase("zeppelin")){

            return new Zeppelin();
        }
        else if(compartmentName.equalsIgnoreCase("tensor")){

            return new Tensor();
        }
        else if(compartmentName.equalsIgnoreCase("cylinderGPD")){

            return new CylinderGPD();
        }
        else if(compartmentName.equalsIgnoreCase("sphereGPD")){

            return new SphereGPD();
        }
        else if(compartmentName.equalsIgnoreCase("dot")){

            return new Dot();
        }
        else if(compartmentName.equalsIgnoreCase("astrosticks")){

            return new Astrosticks();
        }
        else if(compartmentName.equalsIgnoreCase("astrocylinders")){

            return new Astrocylinders();
        }
        else if(compartmentName.equalsIgnoreCase("gammadistribradiicylinders")){

            return new GDRCylinders();
        }
        else{
            throw new LoggedException("unrecognised compartment type name '"+compartmentName+"'");
        }   

    }


}
