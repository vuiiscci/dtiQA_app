package simulation.geometry.elements;

import simulation.DiffusionSimulation;

public class Line {

	/** dimensionality of space */
	private final int D= DiffusionSimulation.D;
	
	/** x-coord of starting point of line */
	private final double lowerX;
	
	/** y-coord of starting point of line */
	private final double lowerY;
	
	/** x-coord of end point of line */
	private final double upperX;
	
	/** y-coord of starting point of line */
	private final double upperY;
	
	/** length of line */
	private final double len;
	
	/** normal to line */
	private final double[] normal= new double[D];
	
	
	public Line(double lowerX, double lowerY, double upperX, double upperY){
		
		this.lowerX=lowerX;
		this.lowerY=lowerY;
		this.upperX=upperX;
		this.upperY=upperY;
		
		this.len= Math.sqrt((upperX-lowerX)*(upperX-lowerX) + 
				(upperY-lowerY)*(upperY-lowerY));
		
		
		normal[0]=-(upperY-lowerY);
		normal[1]=(upperX-lowerX);

		double normLen=Math.sqrt(normal[0]*normal[0] + normal[1]*normal[1]);
		
		normal[0]/=normLen;
		normal[1]/=normLen;
		
		if(D>2){
			normal[2]=0.0;
		}
		
	}
	
	
	public boolean crossedBy(double[] walkerPos, double[] step, double[] normal, double[] d, boolean skipCurrent, double origLength, double[] intDist){

		double[] t=getIntersectionParams(walkerPos[0], walkerPos[1], walkerPos[0]+step[0], walkerPos[1]+step[1]);
		
		double stepLen=0.0;
		
		for(int i=0; i<2; i++){
			stepLen+=step[i]*step[i];
		}
		stepLen=Math.sqrt(stepLen);
		
		if(stepLen/origLength<=1E-14){
			return false;
		}
		
		// check is intersection on step is in step range
		if(t[0]>stepLen){
			return false;
		}
		if(t[0]<0.0){
			return false;
		}
		
		// calculate interesection distance from origin dot normal 
		double intersectionDotNormal=0.0;
		
		// the absolute distance needs to be normalised
		// by the length of the projection into the plane
		// of the cross section
		double arcLen= t[0]/stepLen;
		//double arcLen2= t[1]/this.len;
		
		double[] intPoint= new double[2];
		
		for(int i=0; i<2; i++){
			intPoint[i]=walkerPos[i]+arcLen*step[i];
			intersectionDotNormal+=(walkerPos[i]+arcLen*step[i])*this.normal[i];
		}
		
		// check if walker is sitting on current barrier
		if(t[0]<1E-14*origLength){
			if(skipCurrent){
				// check for actual equality, as this would be set here 
				// from same calculation as this.
				if(intersectionDotNormal==d[0]){
					return false;
				}
			}
		}
		
		// check if the chord intersection param is in the chord's range
		if(t[1]<0.0){
			return false;
		}
		
		if(t[1]>len){
			return false;
		}
		
		// if we've got past all that then there is an intersection
		for(int i=0; i<2; i++){
			normal[i]=this.normal[i];
		}
		normal[2]=0.0;
		
		d[0]=intersectionDotNormal;
		
		intDist[0]=t[0]/stepLen;
		
		return true;
	}
	
	/**
	 * access to the normal to the line
	 * 
	 * @return new copy of line normal
	 */
	public double[] normal(){
		return new double[]{normal[0], normal[1], normal[2]};
	}
	
	/** 
	 * returns the midpoint of the line in new array
	 * 
	 * @return {(upperX+lowerX)/2, (upperY+lowerY)/2}
	 */
	public double[] midpoint(){
		
		return new double[]{(upperX+lowerX)/2.0, (upperY+lowerY)/2.0};
	}
	
	private double[] getIntersectionParams(double xmin, double ymin, double xmax, double ymax){
		
		// lower end of given line
		double x1=xmin;
		double y1=ymin;
		
		// upper end of given line
		double x2=xmax;
		double y2=ymax;
		
		// lower end of this line
		double x3=lowerX;
		double y3=lowerY;
		
		// upper end of this line
		double x4=upperX;
		double y4=upperY;
		
		// elements of numerator matrix
		double a=det(x1, y1, 
				     x2, y2);
		double b=x1-x2;
		double c=det(x3, y3, 
				     x4, y4);
		double d=x3-x4;
		
		// numerator for x coord 
		double num=det(a, b, c, d);
		
		// demoninator matrix determinant
		double denom=det(x1-x2, y1-y2, 
					     x3-x4, y3-y4);
		
		// quotient x
		double x= num/denom;
		
		// redo numerator matrix elements for y coord
		b=y1-y2;
		d=y3-y4;
		
		// recalc numerator
		num=det(a,b,
				c,d);
		
		// denominator unchanged, get y coord
		double y=num/denom;
		
		
		// assemble vector from lower coord of given line to interesection point
		double[] intVec= new double[]{x-xmin, y-ymin};

		// get the dot product of the step and the displacement vector to 
		// check if the intersection is actually on the lines in question
		// or somewhere in the negative arclength region
		double dp=intVec[0]*(xmax-xmin) + intVec[1]*(ymax-ymin);
		double sgn;
		if(dp>0){
			sgn=+1.0;
		}
		else{
			sgn=-1.0;
		}
		
		// length is arclength along given vector
		double t_that=sgn*Math.sqrt(intVec[0]*intVec[0] + intVec[1]*intVec[1]);
		
		// same again, but vector from lower coord on this curve to interesection
		intVec[0] = x-lowerX;
		intVec[1] = y-lowerY;

		// get the dot product of the step and the displacement vector to 
		// check if the intersection is actually on the lines in question
		// or somewhere in the negative arclength region
		dp=intVec[0]*(upperX-lowerX) + intVec[1]*(upperY-lowerY);

		if(dp>0){
			sgn=+1.0;
		}
		else{
			sgn=-1.0;
		}

		
		// arclength along current curve
		double t_this=sgn*Math.sqrt(intVec[0]*intVec[0] + intVec[1]*intVec[1]);
		
		// return both arclengths
		return new double[]{t_that, t_this};
		
	}

	/**
	 * returns the determinant of a 2x2 matrix given by 4 values
	 * @param a top-right
	 * @param b top-left
	 * @param c bottom-right
	 * @param d bottom-left
	 * 
	 * @return ad-bc
	 */
	private double det(double a, double b, double c, double d){
		
		return a*d-b*c;
	}
	
	
	/**
	 * local entrypoint. tests crossing and non-crossing of line and steps
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		
		// construct line
		double lowerX=1.0;
		double lowerY=0.0;
		
		double upperX=1.0;
		double upperY=1.0;
		
		Line line = new Line(lowerX, lowerY, upperX, upperY);
		
		// construct walker and step
		double[] walkerPos= new double[]{0.0, 0.5, 0.0};
		double[] step= new double[]{2.0, 0.0, 0.0};
		double[] normal = new double[walkerPos.length];
		double[] d= new double[1];
		double[] intDist = new double[3];
		boolean crosses;
		
/*		// test intersection
		System.err.println("intersecting lines");
		crosses= line.crossedBy(walkerPos, step, normal, d, false, 2.0, intDist);
		System.err.println("crosses is "+crosses);
		System.err.println("normal= ("+normal[0]+","+normal[1]+","+normal[2]+")");
		System.err.println("d = "+d[0]);
		System.err.println("intDist= "+intDist[0]);
		
		step[0]=0.4;
		System.err.println("\nstep too short");
		crosses= line.crossedBy(walkerPos, step, normal, d, false, 0.4, intDist);
		System.err.println("crosses is "+crosses);
		System.err.println("normal= ("+normal[0]+","+normal[1]+","+normal[2]+")");
		System.err.println("d = "+d[0]);
		System.err.println("intDist= "+intDist[0]);
		
		walkerPos[0]=1.1;
		System.err.println("\nstep is away from line");
		crosses=line.crossedBy(walkerPos, step, normal, d, false, 0.4, intDist);
		System.err.println("crosses is "+crosses);
		System.err.println("normal= ("+normal[0]+","+normal[1]+","+normal[2]+")");
		System.err.println("d = "+d[0]);
		System.err.println("intDist= "+intDist[0]);
*/		
		
		//problem case from debugging
/*		walkerPos= new double[]{1.534991029659386E-5, 1.4898035067115133E-5, 0.0};
		step= new double[]{1.6227999661478954E-5-1.534991029659386E-5, 1.4412789631539843E-5-1.4898035067115133E-5, 0.0};
		line= new Line(1.6E-5, 1.5E-5, 1.580901699437496E-5, 1.5587785252292472E-5);
*/
		walkerPos= new double[]{1.534991029659386, 1.4898035067115133, 0.0};
		step= new double[]{1.6227999661478954-1.534991029659386, 1.4412789631539843-1.4898035067115133, 0.0};
		line= new Line(1.6, 1.5, 1.580901699437496, 1.5587785252292472);

		
		double origLen=0.0;
		for(int i=0; i<step.length; i++){
			origLen+=step[i]*step[i];
		}
		origLen= Math.sqrt(origLen);
		
		
		System.err.println("\nproblem case\n\n");
		crosses=line.crossedBy(walkerPos, step, normal, d, false, origLen, intDist);
		System.err.println("walkerPos =("+walkerPos[0]+","+walkerPos[1]+")");
		System.err.println("end pos =("+(walkerPos[0]+step[0])+","+(walkerPos[1]+step[1])+")");
		System.err.println("line lower = ("+line.lowerX+","+line.lowerY+")");
		System.err.println("line upper = ("+line.upperX+","+line.upperY+")");
		System.err.println("crosses is "+crosses);
		System.err.println("normal = ("+normal[0]+","+normal[1]+","+normal[2]+")");
		System.err.println("d = "+d[0]);
		System.err.println("intDist= ("+intDist[0]+","+intDist[1]+","+intDist[2]+")");

	
		walkerPos= new double[]{1.53, 1.48, 0.0};
		step= new double[]{0.09, -0.04, 0.0};
		line= new Line(1.6, 1.5, 1.58, 1.56);

		
		origLen=0.0;
		for(int i=0; i<step.length; i++){
			origLen+=step[i]*step[i];
		}
		origLen= Math.sqrt(origLen);
		
		
		System.err.println("\ntest case for awkwardness\n\n");
		crosses=line.crossedBy(walkerPos, step, normal, d, false, origLen, intDist);
		System.err.println("walkerPos =("+walkerPos[0]+","+walkerPos[1]+")");
		System.err.println("end pos =("+(walkerPos[0]+step[0])+","+(walkerPos[1]+step[1])+")");
		System.err.println("line lower = ("+line.lowerX+","+line.lowerY+")");
		System.err.println("line upper = ("+line.upperX+","+line.upperY+")");
		System.err.println("crosses is "+crosses);
		System.err.println("normal = ("+normal[0]+","+normal[1]+","+normal[2]+")");
		System.err.println("d = "+d[0]);
		System.err.println("intDist= ("+intDist[0]+","+intDist[1]+","+intDist[2]+")");

	}

}
