package fitters;

import misc.LoggedException;
import models.compartments.CompartmentType;

/**
 * list of models for fitting. these are used to specify CompartmentModel types
 * during model fitting.
 * 
 * 
 * @author Matt (m.hall@cs.ucl.ac.uk)
 * 
 */
public enum FitModel {

	BALLSTICK, BALLCYLINDER, ZEPPELINSTICK, ZEPPELINCYLINDER, TENSORSTICK, TENSORCYLINDER, 
	BALLSTICKDOT, BALLSTICKASTROSTICKS, BALLSTICKASTROCYLINDERS, ZEPPELINSTICKDOT, 
	ZEPPELINSTICKASTROSTICKS, ZEPPELINSTICKASTROCYLINDERS, TENSORSTICKDOT, TENSORSTICKASTROSTICKS, 
	TENSORSTICKASTROCYLINDERS, BALLCYLINDERDOT, ZEPPELINCYLINDERDOT, TENSORCYLINDERDOT, 
	ZEPPELINCYLINDERASTROSTICKS, TENSORCYLINDERASTROSTICKS, BALLCYLINDERASTROSTICKS, 
	BALLCYLINDERASTROCYLINDERS, ZEPPELINCYLINDERASTROCYLINDERS, TENSORCYLINDERASTROCYLINDERS,
	BALLGDRCYLINDERS, BIZEPPELIN,ZEPPELINGDRCYLINDERS,TENSORGDRCYLINDERS,BALLGDRCYLINDERSDOT,
	ZEPPELINGDRCYLINDERSDOT,TENSORGDRCYLINDERSDOT, BALLGDRCYLINDERSASTROSTICKS,BALLGDRCYLINDERSASTROCYLINDERS,
	ZEPPELINGDRCYLINDERSASTROSTICKS,ZEPPELINGDRCYLINDERSASTROCYLINDERS,
	TENSORGDRCYLINDERSASTROSTICKS, TENSORGDRCYLINDERSASTROCYLINDERS,  
	  BALLSTICKSPHERE, ZEPPELINSTICKSPHERE,TENSORSTICKSPHERE,
	BALLCYLINDERSPHERE,ZEPPELINCYLINDERSPHERE,  
	TENSORCYLINDERSPHERE, BALLGDRCYLINDERSSPHERE,ZEPPELINGDRCYLINDERSSPHERE, TENSORGDRCYLINDERSSPHERE,
            ZEPPELINSTICKDIRECT, ZEPPELINSTICKTORT, ZEPPELINCYLINDERDIRECT, ZEPPELINCYLINDERTORT, 
	    ZEPPELINCYLINDERDOTDIRECT, ZEPPELINCYLINDERDOTCSFDIRECT, MMWMDBASIC, MMWMDINVIVO, MMWMDFIXED, MMWMDFIXEDNOCSF, MMWMDINVIVONOCSF, VERDICTCOLORECTAL;


	/**
	 * returns the compartments array used to initialise the compartment model
	 * for fitting from commandline string.
	 * 
	 * @param fitModelString
	 *            which model to fit
	 * @return compartment type array for fitting model
	 */
	public static final String[] getCompartmentList(String fitModelString) {

		FitModel fitModel = getFitModel(fitModelString);

		if (fitModel == BALLSTICK) {
			return new String[] { CompartmentType.STICK.toString(),
					CompartmentType.BALL.toString() };
		} else if (fitModel == BALLCYLINDER) {
			return new String[] { CompartmentType.CYLINDERGPD.toString(),
					CompartmentType.BALL.toString() };
		}  else if (fitModel == BIZEPPELIN) {
			return new String[] { CompartmentType.ZEPPELIN.toString(),
					CompartmentType.ZEPPELIN.toString() };
		}else if (fitModel == ZEPPELINSTICK) {
			return new String[] { CompartmentType.STICK.toString(),
					CompartmentType.ZEPPELIN.toString() };
		}else if (fitModel == ZEPPELINSTICKDIRECT) {
			return new String[] { CompartmentType.STICK.toString(),
					CompartmentType.ZEPPELIN.toString() };
		}else if (fitModel == ZEPPELINSTICKTORT) {
			return new String[] { CompartmentType.STICK.toString(),
					CompartmentType.ZEPPELIN.toString() };
		} else if (fitModel == ZEPPELINCYLINDER) {
			return new String[] { CompartmentType.CYLINDERGPD.toString(),
					CompartmentType.ZEPPELIN.toString() };
		} else if (fitModel == ZEPPELINCYLINDERDIRECT) {
			return new String[] { CompartmentType.CYLINDERGPD.toString(),
					CompartmentType.ZEPPELIN.toString() };
		} else if (fitModel == ZEPPELINCYLINDERTORT) {
			return new String[] { CompartmentType.CYLINDERGPD.toString(),
					CompartmentType.ZEPPELIN.toString() };
		} else if (fitModel == TENSORSTICK) {
			return new String[] { CompartmentType.STICK.toString(),
					CompartmentType.TENSOR.toString() };
		} else if (fitModel == TENSORCYLINDER) {
			return new String[] { CompartmentType.CYLINDERGPD.toString(),
					CompartmentType.TENSOR.toString() };
		} else if (fitModel == BALLSTICKDOT) {
			return new String[] { CompartmentType.STICK.toString(),
					CompartmentType.BALL.toString(),
					CompartmentType.DOT.toString() };
		} else if (fitModel == BALLSTICKASTROSTICKS) {
			return new String[] { CompartmentType.STICK.toString(),
					CompartmentType.BALL.toString(),
					CompartmentType.ASTROSTICKS.toString() };
		} else if (fitModel == BALLSTICKASTROCYLINDERS) {
			return new String[] { CompartmentType.STICK.toString(),
					CompartmentType.BALL.toString(),
					CompartmentType.ASTROCYLINDERS.toString() };
		} else if (fitModel == BALLSTICKSPHERE) {
			return new String[] { CompartmentType.STICK.toString(),
					CompartmentType.BALL.toString(),
					CompartmentType.SPHEREGPD.toString() };
		}else if (fitModel == ZEPPELINSTICKDOT) {
			return new String[] { CompartmentType.STICK.toString(),
					CompartmentType.ZEPPELIN.toString(),
					CompartmentType.DOT.toString() };
		} else if (fitModel == ZEPPELINSTICKASTROSTICKS) {
			return new String[] { CompartmentType.STICK.toString(),
					CompartmentType.ZEPPELIN.toString(),
					CompartmentType.ASTROSTICKS.toString() };
		} else if (fitModel == ZEPPELINSTICKASTROCYLINDERS) {
			return new String[] { CompartmentType.STICK.toString(),
					CompartmentType.ZEPPELIN.toString(),
					CompartmentType.ASTROCYLINDERS.toString() };
		} else if (fitModel == TENSORSTICKDOT) {
			return new String[] { CompartmentType.STICK.toString(),
					CompartmentType.TENSOR.toString(),
					CompartmentType.DOT.toString() };
		} else if (fitModel == TENSORSTICKSPHERE) {
			return new String[] { CompartmentType.STICK.toString(),
					CompartmentType.TENSOR.toString(),
					CompartmentType.SPHEREGPD.toString() };
		} else if (fitModel == TENSORSTICKASTROSTICKS) {
			return new String[] { CompartmentType.STICK.toString(),
					CompartmentType.TENSOR.toString(),
					CompartmentType.ASTROSTICKS.toString() };
		} else if (fitModel == TENSORSTICKASTROCYLINDERS) {
			return new String[] { CompartmentType.STICK.toString(),
					CompartmentType.TENSOR.toString(),
					CompartmentType.ASTROCYLINDERS.toString() };
		} else if (fitModel == BALLCYLINDERDOT) {
			return new String[] { CompartmentType.CYLINDERGPD.toString(),
					CompartmentType.BALL.toString(),
					CompartmentType.DOT.toString() };
		} else if (fitModel == BALLCYLINDERSPHERE) {
			return new String[] { CompartmentType.CYLINDERGPD.toString(),
					CompartmentType.BALL.toString(),
					CompartmentType.SPHEREGPD.toString() };
		} else if (fitModel == ZEPPELINCYLINDERDOT) {
			return new String[] { CompartmentType.CYLINDERGPD.toString(),
					CompartmentType.ZEPPELIN.toString(),
					CompartmentType.DOT.toString() };
		} else if (fitModel == ZEPPELINCYLINDERDOTDIRECT) {
			return new String[] { CompartmentType.CYLINDERGPD.toString(),
					CompartmentType.ZEPPELIN.toString(),
					CompartmentType.DOT.toString() };
		} else if (fitModel == ZEPPELINCYLINDERDOTCSFDIRECT) {
			return new String[] { CompartmentType.CYLINDERGPD.toString(),
					CompartmentType.ZEPPELIN.toString(),
					CompartmentType.DOT.toString(),
					CompartmentType.BALL.toString() };
		} else if (fitModel == ZEPPELINCYLINDERSPHERE) {
			return new String[] { CompartmentType.CYLINDERGPD.toString(),
					CompartmentType.ZEPPELIN.toString(),
					CompartmentType.SPHEREGPD.toString() };
		}else if (fitModel == TENSORCYLINDERDOT) {
			return new String[] { CompartmentType.CYLINDERGPD.toString(),
					CompartmentType.TENSOR.toString(),
					CompartmentType.DOT.toString() };
		}else if (fitModel == TENSORCYLINDERSPHERE) {
			return new String[] { CompartmentType.CYLINDERGPD.toString(),
					CompartmentType.TENSOR.toString(),
					CompartmentType.SPHEREGPD.toString() };
		} else if (fitModel == BALLCYLINDERASTROSTICKS) {
			return new String[] { CompartmentType.CYLINDERGPD.toString(),
					CompartmentType.BALL.toString(),
					CompartmentType.ASTROSTICKS.toString() };
		} else if (fitModel == ZEPPELINCYLINDERASTROSTICKS) {
			return new String[] { CompartmentType.CYLINDERGPD.toString(),
					CompartmentType.ZEPPELIN.toString(),
					CompartmentType.ASTROSTICKS.toString() };
		} else if (fitModel == TENSORCYLINDERASTROSTICKS) {
			return new String[] { CompartmentType.CYLINDERGPD.toString(),
					CompartmentType.TENSOR.toString(),
					CompartmentType.ASTROSTICKS.toString() };
		} else if (fitModel == TENSORCYLINDERASTROCYLINDERS) {
			return new String[] { CompartmentType.CYLINDERGPD.toString(),
					CompartmentType.TENSOR.toString(),
					CompartmentType.ASTROCYLINDERS.toString() };
		} else if (fitModel == BALLCYLINDERASTROCYLINDERS) {
			return new String[] { CompartmentType.CYLINDERGPD.toString(),
					CompartmentType.BALL.toString(),
					CompartmentType.ASTROCYLINDERS.toString() };
		} else if (fitModel == ZEPPELINCYLINDERASTROCYLINDERS) {
			return new String[] { CompartmentType.CYLINDERGPD.toString(),
					CompartmentType.ZEPPELIN.toString(),
					CompartmentType.ASTROCYLINDERS.toString() };
		} else if (fitModel == BALLGDRCYLINDERS) {
			return new String[] { CompartmentType.GAMMADISTRIBRADIICYLINDERS.toString(),
					CompartmentType.BALL.toString() };
		} else if (fitModel == BALLGDRCYLINDERSASTROSTICKS) {
			return new String[] { CompartmentType.GAMMADISTRIBRADIICYLINDERS.toString(),
					CompartmentType.BALL.toString(),CompartmentType.ASTROSTICKS.toString() };
		}else if (fitModel == BALLGDRCYLINDERSSPHERE) {
			return new String[] { CompartmentType.GAMMADISTRIBRADIICYLINDERS.toString(),
					CompartmentType.BALL.toString(),CompartmentType.SPHEREGPD.toString() };
		}
		else if (fitModel == BALLGDRCYLINDERSASTROCYLINDERS) {
			return new String[] { CompartmentType.GAMMADISTRIBRADIICYLINDERS.toString(),
					CompartmentType.BALL.toString(),CompartmentType.ASTROCYLINDERS.toString() };
		}else if (fitModel == BALLGDRCYLINDERSDOT) {
			return new String[] { CompartmentType.GAMMADISTRIBRADIICYLINDERS.toString(),
					CompartmentType.BALL.toString(),CompartmentType.DOT.toString() };
		}
		else if (fitModel == ZEPPELINGDRCYLINDERSASTROSTICKS) {
			return new String[] { CompartmentType.GAMMADISTRIBRADIICYLINDERS.toString(),
					CompartmentType.ZEPPELIN.toString(),CompartmentType.ASTROSTICKS.toString() };
		}else if (fitModel == ZEPPELINGDRCYLINDERSSPHERE) {
			return new String[] { CompartmentType.GAMMADISTRIBRADIICYLINDERS.toString(),
					CompartmentType.ZEPPELIN.toString(),CompartmentType.SPHEREGPD.toString() };
		}
		else if (fitModel == ZEPPELINGDRCYLINDERSASTROCYLINDERS) {
			return new String[] { CompartmentType.GAMMADISTRIBRADIICYLINDERS.toString(),
					CompartmentType.ZEPPELIN.toString(),CompartmentType.ASTROCYLINDERS.toString() };
		}
		else if (fitModel == TENSORGDRCYLINDERSASTROSTICKS) {
			return new String[] { CompartmentType.GAMMADISTRIBRADIICYLINDERS.toString(),
					CompartmentType.TENSOR.toString(),CompartmentType.ASTROSTICKS.toString() };
		}
		else if (fitModel == TENSORGDRCYLINDERSSPHERE) {
			return new String[] { CompartmentType.GAMMADISTRIBRADIICYLINDERS.toString(),
					CompartmentType.TENSOR.toString(),CompartmentType.SPHEREGPD.toString() };
		}
		else if (fitModel == TENSORGDRCYLINDERSASTROCYLINDERS) {
			return new String[] { CompartmentType.GAMMADISTRIBRADIICYLINDERS.toString(),
					CompartmentType.TENSOR.toString(),CompartmentType.ASTROCYLINDERS.toString() };
		}		
		else if (fitModel == ZEPPELINGDRCYLINDERS) {
			return new String[] { CompartmentType.GAMMADISTRIBRADIICYLINDERS.toString(),
					CompartmentType.ZEPPELIN.toString() };
		} else if (fitModel == ZEPPELINGDRCYLINDERSDOT) {
			return new String[] { CompartmentType.GAMMADISTRIBRADIICYLINDERS.toString(),
					CompartmentType.ZEPPELIN.toString(),CompartmentType.DOT.toString() };
		} else if (fitModel == TENSORGDRCYLINDERS) {
			return new String[] { CompartmentType.GAMMADISTRIBRADIICYLINDERS.toString(),
					CompartmentType.TENSOR.toString() };
		} 
		 else if (fitModel == TENSORGDRCYLINDERSDOT) {
				return new String[] { CompartmentType.GAMMADISTRIBRADIICYLINDERS.toString(),
						CompartmentType.TENSOR.toString(),CompartmentType.DOT.toString() };
                 }  else if (fitModel == MMWMDBASIC) {
			return new String[] { CompartmentType.CYLINDERGPD.toString(),
					      CompartmentType.ZEPPELIN.toString() };
                 }  else if (fitModel == MMWMDINVIVO) {
			return new String[] { CompartmentType.CYLINDERGPD.toString(),
					CompartmentType.ZEPPELIN.toString(),
					CompartmentType.BALL.toString() };
                 }  else if (fitModel == MMWMDFIXED) {
			return new String[] { CompartmentType.CYLINDERGPD.toString(),
					CompartmentType.ZEPPELIN.toString(),
                                              CompartmentType.DOT.toString(),
					CompartmentType.BALL.toString() };
                 }  else if (fitModel == MMWMDFIXEDNOCSF) {
			return new String[] { CompartmentType.CYLINDERGPD.toString(),
					CompartmentType.ZEPPELIN.toString(),
                                              CompartmentType.DOT.toString() };
                 }
		else if (fitModel == MMWMDINVIVONOCSF)
		{
			return new String[] { CompartmentType.CYLINDERGPD.toString(),
					CompartmentType.ZEPPELIN.toString()};
		}
		else if (fitModel == VERDICTCOLORECTAL)
		{
			return new String[] { CompartmentType.STICK.toString(),
					CompartmentType.BALL.toString(), CompartmentType.SPHEREGPD.toString()};
		}
		else {
			throw new LoggedException("Unknown fitting model " + fitModel);
		}

	}

	/**
	 * returns dummy params array for the compartment model used in fitting.
	 * these params are completely ignored by the fitting procedure but are
	 * required at initialisation.
	 * 
	 * @param fitModelString
	 * @return dummy params array
	 */
	public static final double[] getCompartmentModelParams(String fitModelString) {

		FitModel fitModel = getFitModel(fitModelString);

		if (fitModel == BALLSTICK) {
			return new double[] { 2, 0.5, 0.5, 1.0, 0.0, 1.0, 1.0 };
		} else if (fitModel == BALLCYLINDER) {
			return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0 };
		} else if (fitModel == ZEPPELINSTICK) {
			return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
					0.5 };
		} else if (fitModel == ZEPPELINSTICKDIRECT) {
			return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
					0.5 };
		} else if (fitModel == ZEPPELINSTICKTORT) {
			return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
					0.5 };
		} else if (fitModel == ZEPPELINCYLINDER) {
			return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
					0.5, 1.0 };
		} else if (fitModel == ZEPPELINCYLINDERDIRECT) {
			return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
					0.5, 1.0 };
		} else if (fitModel == ZEPPELINCYLINDERTORT) {
			return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
					0.5, 1.0 };
		} else if (fitModel == TENSORSTICK) {
			return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
					0.5, 0.5, 0.0 };
		} else if (fitModel == TENSORCYLINDER) {
			return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
					0.5, 0.5, 0.5, 0.0 };
		} else if (fitModel == BALLSTICKDOT) {
			return new double[] { 2, 0.5, 0.5, 1.0, 0.0, 1.0, 1.0, 0.5 };
		} else if (fitModel == BALLSTICKASTROSTICKS) {
			return new double[] { 2, 0.5, 0.5, 1.0, 0.0, 1.0, 1.0, 0.5, 0.5 };
		} else if (fitModel == BALLSTICKASTROCYLINDERS) {
			return new double[] { 2, 0.5, 0.5, 1.0, 0.0, 1.0, 1.0, 0.5, 0.5,
					0.5 };
		} 
		else if (fitModel == BALLSTICKSPHERE) {
			return new double[] { 2, 0.5, 0.5, 1.0, 0.0, 1.0, 1.0, 0.5, 0.5,
					0.5 };
		}else if (fitModel == ZEPPELINSTICKDOT) {
			return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
					0.5 };
		} else if (fitModel == ZEPPELINSTICKASTROSTICKS) {
			return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
					0.5, 0.5 };
		} else if (fitModel == ZEPPELINSTICKASTROCYLINDERS) {
			return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
					0.5, 0.5, 0.5 };
		}  else if (fitModel == ZEPPELINSTICKSPHERE) {
			return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
					0.5, 0.5, 0.5 };
		}else if (fitModel == TENSORSTICKDOT) {
			return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
					0.5, 1, 1, 1 };
		} else if (fitModel == BALLCYLINDERDOT) {
			return new double[] { 1, 0.5, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0 };
		} else if (fitModel == ZEPPELINCYLINDERDOT) {
			return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
					0.5, 1.0, 0.5 };
		} else if (fitModel == ZEPPELINCYLINDERDOTDIRECT) {
			return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
					0.5, 1.0, 0.5 };
		} else if (fitModel == ZEPPELINCYLINDERDOTCSFDIRECT) {
                    return new double[] { 1, 0.5, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
					0.5, 1.0, 0.5, 0.5 };
		} else if (fitModel == TENSORCYLINDERDOT) {
			return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
					0.5, 1.0, 0.5, 0.5, 0.5 };
		} else if (fitModel == BALLCYLINDERASTROSTICKS) {
			return new double[] { 1, 0.5, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0,
					1.0 };
		} else if (fitModel == BALLCYLINDERASTROCYLINDERS) {
			return new double[] { 1, 0.5, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0,
					1.0, 0.5 };
		} 
		else if (fitModel == ZEPPELINCYLINDERASTROSTICKS) {
			return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
					0.5, 0.5, 0.5 };
		} 
		else if (fitModel == ZEPPELINCYLINDERSPHERE) {
			return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
					0.5, 0.5, 0.5 ,0.5,0.5};
		}else if (fitModel == TENSORSTICKASTROSTICKS) {
			return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
					0.5, 1.0, 0.5, 0.5 };
		} else if (fitModel == TENSORSTICKASTROCYLINDERS) {
			return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
					0.5, 1.0, 0.5, 0.5, 0.5 };
		} else if (fitModel == TENSORSTICKSPHERE) {
			return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
					0.5, 1.0, 0.5, 0.5, 0.5 };
		} else if (fitModel == TENSORCYLINDERASTROSTICKS) {
			return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
					0.5, 1.0, 0.5, 0.5, 0.5, 0.5 };
		} else if (fitModel == TENSORCYLINDERASTROCYLINDERS) {
			return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
					0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5 };
		} else if (fitModel == TENSORCYLINDERSPHERE) {
			return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
					0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5,0.5 };
		}else if (fitModel == ZEPPELINCYLINDERASTROCYLINDERS) {
			return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
					0.5, 0.5, 0.5, 0.5 };
		} else if (fitModel == MMWMDBASIC) {
                    return new double[] { 1, 0.5, 0.25, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0};
		} else if (fitModel == MMWMDINVIVO) {
                    return new double[] { 1, 0.5, 0.25, 0.25, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
		} else if (fitModel == MMWMDFIXED) {
                    return new double[] { 1, 0.5, 0.25, 0.25, 0.25, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
		}
		else if (fitModel == MMWMDFIXEDNOCSF) {
                    return new double[] { 1, 0.5, 0.25, 0.25, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
		}
		else if (fitModel == MMWMDINVIVONOCSF) {
                    return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.5, 1.0 };
		}
		 else if (fitModel == BALLGDRCYLINDERS) {
				return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0 };
			}
		 else if (fitModel == BALLGDRCYLINDERSASTROSTICKS) {
				return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0 };
			}
		 else if (fitModel == BALLGDRCYLINDERSASTROCYLINDERS) {
				return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0 };
			}
		 else if (fitModel == BALLGDRCYLINDERSSPHERE) {
				return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,0.5 };
			}
		 else if (fitModel == ZEPPELINGDRCYLINDERSASTROSTICKS) {
				return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0 };
			}
		 else if (fitModel == ZEPPELINGDRCYLINDERSASTROCYLINDERS) {
				return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0 };
			}
		 else if (fitModel == ZEPPELINGDRCYLINDERSSPHERE) {
				return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,0.5 };
			}
		 else if (fitModel == TENSORGDRCYLINDERSASTROSTICKS) {
				return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0 };
			}
		 else if (fitModel == TENSORGDRCYLINDERSASTROCYLINDERS) {
				return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0 };
			}
		 else if (fitModel == BALLGDRCYLINDERSDOT) {
				return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,1 };
			}
		 else if (fitModel == BIZEPPELIN) {
				return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,1,1 };
			}
		 else if (fitModel == ZEPPELINGDRCYLINDERS) {
				return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,1,1,1 };
			}
		 else if (fitModel == ZEPPELINGDRCYLINDERSDOT) {
				return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,1,1,1,1 };
			}
		 else if (fitModel == TENSORGDRCYLINDERS) {
				return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,1,1,1,1,1 };
			}
		 else if (fitModel == TENSORGDRCYLINDERSDOT) {
				return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,1,1,1,1,1,1 };
			}
		 else if (fitModel == VERDICTCOLORECTAL) {
				return new double[] { 1, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 1, 1 };
			}

		else {
			throw new LoggedException("Unknown fitting model " + fitModel);
		}

	}

	/**
	 * get the model name from the commandline string
	 * 
	 * @param fitModelString
	 *            which model do we want?
	 * @return the enum object matching the string
	 */
	public static final FitModel getFitModel(String fitModelString) {

		for (FitModel ft : FitModel.values()) {
			if (ft.toString().equalsIgnoreCase(fitModelString)) {
				return ft;
			}
		}

		throw new LoggedException("unrecognised fit model name "
				+ fitModelString);
	}

}
