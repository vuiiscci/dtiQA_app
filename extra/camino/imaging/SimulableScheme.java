package imaging;


public interface SimulableScheme {



    /**
     * Gets the average gradient strength in Tesla in effect during specific time from t to 
     * tLast during a measurement.
     * 
     * @param i the index of the measurement.
     * @param t the elapsed time from the start of the particular measurement sequence, in seconds.
     *
     * @return if the gradient is on, G * duration of currrent timestep inside gradient, usually Gdt.
     */
    public double[] getGradImpulse(int i, double t, double tLast);


    /**
     * gets the duration of the scan in seconds. In most cases this will be the 
     * time of the end of the final gradient pulse, but this is different
     * 
     *  @return duration in seconds
     */
    public double getDuration();
    
    
}
