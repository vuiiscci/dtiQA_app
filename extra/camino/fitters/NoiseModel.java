package fitters;

import misc.LoggedException;

public enum NoiseModel {
	
	GAUSSIAN,
        OFFGAUSS,
	RICIAN;
	
	public static final NoiseModel getNoiseModel(String noiseModelString){
		
		for(NoiseModel nm: NoiseModel.values()){
			if(nm.toString().equalsIgnoreCase(noiseModelString)){
				return nm;
			}
		}
		
		throw new LoggedException("unknown noise model "+noiseModelString);
	}

}
