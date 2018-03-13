package fitters;

import misc.LoggedException;

public enum FitAlgorithm {
	
	LM, MULTIRUNLM,MCMC;
	
	public static final FitAlgorithm getFitAlgorithm(String fitAlgString){
		
		for(FitAlgorithm fa: FitAlgorithm.values()){
			if(fa.toString().equalsIgnoreCase(fitAlgString)){
				return fa;
			}
		}
		
		throw new LoggedException("unknown fitting algorithm "+fitAlgString);
	}

}
