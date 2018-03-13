package simulation.geometry.elements;

import java.util.ArrayList;

import misc.LoggedException;
import simulation.DiffusionSimulation;
import simulation.SimulationParams;
import simulation.dynamics.exceptions.TooDamnCloseException;
import simulation.geometry.elements.CylinderFactory.CylType;


public class NestedCylinder extends BasicCylinder {
	
	/** dimensionality of space */
	private final int D= DiffusionSimulation.D;
	
	/** outer cylinder - contains the nested ones */
	private Cylinder outer;
	
	/** array of 19 inner cylinders */
	private Cylinder[] nested;
	
	/** is this a bottom-level cylinder? */
	private final boolean isBasic;
	
	/** geoetric packing constant -- assumes optimal 19 cylinder packing */
	private final double c= 1.0/(1.0 + Math.sqrt(2.0) + Math.sqrt(6.0));
	
	/** 
	 *  scaling constant used to make cylinders a little smaller than optimal.
	 *  this ensures that there's a little space around them
	 */
	private final double d=0.999;
	
	public NestedCylinder(double[] P, double r, double p, int depth){
		super(P, r, p);				// invoke superconstructor
		
		isBasic= (depth==0);
		
		outer= new BasicCylinder(P, r, p);
		
		/* if depth is zero then there's nothing more to do,
		 * otherwise we need to add the sub-cylinders 
		 */
		if(depth<0){
			throw new LoggedException("cylinder nesting depth given as "+depth+". depth cannot be negative.");
		}
		
		
		if(depth>0){
			nested= new NestedCylinder[19];
			

			double[] Pnew= new double[D];
			
			
			// small cylinder radius: outer radius times packing constant
			double R= r*c;
			
			// zeroth cylinder shares same centre as outer
			nested[0]= new NestedCylinder(P, R*d, p, depth-1);
			
			
			// the next six are at the vertices of a hexagon optimally surrounding the central one
			for(int i=0; i<6; i++){
				
				double theta= i*Math.PI/3;
				
				// new position offset from centre
				Pnew[0]= P[0]+ 2*R*Math.cos(theta);
				Pnew[1]= P[1]+ 2*R*Math.sin(theta);
				Pnew[2]= P[2];
				
				// construct new cylinder
				nested[i+1]= new NestedCylinder(Pnew, R*d, p, depth-1);
			}
			
			// the remaining 12 cylinders are offset from each of the hexagonal 6
			for(int i=0; i<6; i++){						// upper offset
				
				double theta1= i*Math.PI/3;
				double theta2= (2*i+1)*Math.PI/6;
				
				Pnew[0]= P[0]+2*R*(Math.cos(theta1)+Math.cos(theta2));
				Pnew[1]= P[1]+2*R*(Math.sin(theta1)+Math.sin(theta2));
				Pnew[2]= P[2];
				
				nested[7+i]= new NestedCylinder(Pnew, R*d, p, depth-1);
			}

			for(int i=0; i<6; i++){						// lower offset
				
				double theta1= i*Math.PI/3;
				double theta2= (2*i-1)*Math.PI/6;
				
				Pnew[0]= P[0]+2*R*(Math.cos(theta1)+Math.cos(theta2));
				Pnew[1]= P[1]+2*R*(Math.sin(theta1)+Math.sin(theta2));
				Pnew[2]= P[2];
				
				nested[13+i]= new NestedCylinder(Pnew, R*d, p, depth-1);
			}
		}
	}
	
	
	// methods that need to be overridden
	/*public boolean crosses(double[] rawPos, double[] rawStep,
			double[] normal, double[] d, boolean skipCurrent, double origLength, 
			double[] intDist, boolean[] in, double[] p, double walkerRad)
			throws TooDamnCloseException {
		
		final double[][] normals= new double[20][D];
		//final double[][] intPoints= new double[20][D];
		final double[][] ds= new double[20][1];
		final double[][] intDists= new double[20][1];
		final double[][] ps= new double[20][1];
		final boolean[][] ins= new boolean[20][1];
		
		int crossings=0;
		
		for(int i=0; i<20; i++){
			for(int j=0; j<D; j++){
				normals[i][j]= normal[j];
			}
			
			ds[i][0]=d[0];
			intDists[i][0]=intDist[0];
			ps[i][0]=p[0];
			ins[i][0]=in[0];
		}
		
		
		if(outer.crosses(rawPos,  rawStep, normal, d, skipCurrent, origLength, intDist, in, p, walkerRad))
		
	}*/
	
	
	/**
	 * overridden get triangles method which gets triangles from
	 * outer and all nested triangles recursively
	 */
	public ArrayList<Triangle> getTriangles(){
		if(isBasic){
			return super.getTriangles();
		}
		
		ArrayList<Triangle> tris= new ArrayList<Triangle>();
		
		tris.addAll(outer.getTriangles());
		
		for(int i=0; i<nested.length; i++){
			tris.addAll(nested[i].getTriangles());
		}
		
		return tris;
	}
	
	/**
	 * overrdides the intracellular check for the cylinder.
	 * This works as follows: only locations inside the bottom-level
	 * cylinders are considered intracellular. If we're inside the outer
	 * cylinder but outside the inner nested cylinders, this is 
	 * extracellular.
	 * 
	 * (TODO: this may be incompatible with the intracellular vol frac
	 * calculation in the containign substrate)
	 * 
	 * @param R location to check.
	 * @return true if intracellular.
	 */
	public boolean inside(double[] R){
		
		// bottom-level cylinders are covered by the superclass
		if(isBasic){
			return super.inside(R);
		}
		
		// if we're inside the outer, check the inners
		if(outer.inside(R)){
			boolean inside= false;
			for(int i=0; i<nested.length; i++){
				if(nested[i].inside(R)){
					return true;
				}
			}
		}
		
		// if we're not inside any of the inners, return false
		return false;
	}

	
	public ArrayList<Cylinder> allCylinders(){
		
		ArrayList<Cylinder> allCyls= new ArrayList<Cylinder>();
		
		allCyls.add(outer);
		if(!isBasic){			
			for(int i=0; i<nested.length; i++){
				allCyls.addAll(((NestedCylinder)nested[i]).allCylinders());
			}
		}
		
		return allCyls;
	}
	
	
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
