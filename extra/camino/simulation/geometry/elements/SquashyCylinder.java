package simulation.geometry.elements;

import imaging.RectGradSteTanScheme;
import imaging.SimulableScheme;

import java.io.BufferedWriter;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.logging.Logger;

import simulation.DiffusionSimulation;
import simulation.SimulationParams;
import simulation.dynamics.Walker;
import simulation.dynamics.StepGeneratorFactory.StepType;
import simulation.geometry.substrates.Substrate;
import simulation.geometry.substrates.SubstrateFactory;
import tools.CL_Initializer;

import misc.LoggedException;
import numerics.Vector3D;

/**
 * A squashy cylinder is one that deforms when it abutts another
 * cylinder. It is assumed that the other cylinder is squashy as 
 * well so that the two cylinders form a straight line barrier
 * between then along the Chord defined by the points of intersection
 * on the circumference. It is also assumed that the two cylinders 
 * have parallel axes.
 * 
 * 
 * @author matt (m.hall@cs.ucl.ac.uk)
 * 
 */
public class SquashyCylinder extends BasicCylinder {

	/** logging object */
	private Logger logger=Logger.getLogger(this.getClass().getName());
	
	/** dimensionality of space */
	private final int D=DiffusionSimulation.D;
	
	
	/** max angle of intersection before cylinder stops expanding */
	private static final double MaxIntersectionAngle = 0.7*Math.PI;
	
	/** array of ranges of Chords */
	private ArrayList<AnglePair> anglePairs= new ArrayList<AnglePair>();
	
	/** array of Chords and lengths */
	public ArrayList<Chord> chords= new ArrayList<Chord>();
	
	/** array of cylinders we intersect with */
	private ArrayList<SquashyCylinder> bumps= new ArrayList<SquashyCylinder>();
	
	/** hashmap relating anglePairs to chords */
	private HashMap<AnglePair, Chord> chordsForAp= 
							new HashMap<AnglePair, Chord>();
	
	/** hashmap relating chords to anglepairs */
	private HashMap<Chord, AnglePair> apsForChord=
							new HashMap<Chord, AnglePair>();
	
	/** hashmap relating cylinders to anglePairs */
	private HashMap<SquashyCylinder, AnglePair> apsForCylinders= 
							new HashMap<SquashyCylinder, AnglePair>();
	
	/** hashmap relating chords to intersecting cylinders */
	private HashMap<Chord, SquashyCylinder> cylsForChords=new HashMap<Chord, SquashyCylinder>(); 
	
	/** flag that says whether cylinder is still expanding or not */
	private boolean expanding = true;
	
	/** index of initial positon of walker before reflection loop */
	private double[] initialPos= new double[D];
	
	/** records whether the initial walker position is inside or outside the cylinder */
	private boolean initiallyIn;

	/** tag for ignoring chord intersections */
	private int skipChord =-1;
	
	// position and radius of the cylinder
	private final double[] P;
	private final double r;

	/** clones of this cylinder */
	private ArrayList<SquashyCylinder> myClones=null;
	
	/** membrane permeability */
	private final double p;
	
	/**
	 * instantiates a cylinder with no squishes parallel to z-axis
	 * 
	 * @param P position of cylinder
	 * @param r radius of cylinder
	 */
	public SquashyCylinder(double[] P, double r, double p) {
		super(P, r, p);
		
		this.P=super.getPosition();
		this.r=super.getRadius();
		this.p=p;
		
	}

	public SquashyCylinder(double[] V, double[] P, double r, double p){
	    super(V, P, r, p);
	    
	    logger.warning("arbitrarily oriented squashy cylinders are not yet implemented!");
	    
	    this.P=P;
	    this.r=r;
	    this.p=p;
	}	    
	    
	
	
	/**
	 * checks if this cylinder abuts another
	 * 
	 *  @param that the other cylinder
	 *  @return true if abutting, otherwise false
	 */
	public boolean abutts(SquashyCylinder that){
		
		double dist = this.getDistanceFrom(that.getPosition());
		
		if(dist<=(this.getRadius()+that.getRadius())){
			return true;
		}
		else{
			return false;
		}
	}

	/**
	 * checks if cylinder is expanding or not
	 * 
	 * @return true or false
	 */
	public final boolean isExpanding(){
		return expanding;
	}
	
	
	/**
	 * stops this cylinder and all clones expanding due to pressure 
	 * from another one.
	 */
	public final void stopExpanding(){
		expanding = false;
		
		// if there are clones, stop them expanding as well
		if(myClones!=null){
			for(Iterator<SquashyCylinder> cloneIt= myClones.iterator(); cloneIt.hasNext(); ){
				SquashyCylinder clone= cloneIt.next();
			
				clone.cloneHasStopped();
			}
		}
	}
	
	/**
	 * a stop instruction from a clone. This must be a separate method
	 * from the stopExpanding() call in order to prevent an infinite 
	 * loop.
	 */
	public final void cloneHasStopped(){
		expanding=false;
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
	
	
	/** performs circle-line intersection for the 2d problem
	 *  in the plane of the cross-section of the cylinder.
	 *  
	 *  because the arclengths are different in the projection 
	 *  into 
	 *  
	 * @param walkerPos
	 * @param step
	 * 
	 * @return
	 */
	private final double[] get2dIntersectionRoots(double[] walkerPos, double[] step){
		
		//position of cylinder
		double[] P= getPosition();
		
		// position of walker in cylinder coords
		double x1 = walkerPos[0]-P[0];
		double y1 = walkerPos[1]-P[1];
		
		// position of end of step in cylinder coords
		double x2 = walkerPos[0]+step[0]-P[0];
		double y2 = walkerPos[1]+step[1]-P[1];
		
		// differences
		double dx= x2-x1;
		double dy= y2-y1;
		
		// length of difference vector (step length)
		double dr= Math.sqrt(dx*dx+dy*dy);
		
		// determinant of step matrix
		double D= x1*y2-x2*y1;

		// radius of cylinder
		double r = getRadius();
		
		// solve quadratic
		double surd;
		
		double discriminant=r*r*dr*dr-D*D;
		
		if(discriminant>=0.0){
			surd= Math.sqrt(discriminant);
		}
		else{
			return null;
		}
		
		double sgn;
		
		if(dy<=0){
			sgn=-1;
		}
		else{
			sgn=1;
		}
		
		double xSurdTerm = sgn*dx*surd;
		double ySurdTerm = Math.abs(dy)*surd;
		
		// construct solution pair
		double X1 = (D*dy + xSurdTerm)/(dr*dr);
		double Y1 = (-D*dx + ySurdTerm)/(dr*dr);
		
		double X2 = (D*dy - xSurdTerm)/(dr*dr);
		double Y2 = (-D*dx - ySurdTerm)/(dr*dr);
		
		
		// construct pair of arclengths
		double[] roots= new double[2];
		
		roots[0] = Math.sqrt((X1-x1)*(X1-x1) + (Y1-y1)*(Y1-y1));
		roots[1] = Math.sqrt((X2-x2)*(X2-x2) + (Y2-y2)*(Y2-y2));

		// normalise to length of step
		double stepLen=0.0;
		for(int i=0; i<step.length; i++){
			stepLen+=step[i]*step[i];
		}
		stepLen=Math.sqrt(stepLen);
		
		roots[0]/=stepLen;
		roots[1]/=stepLen;
		
		return roots;
	}
	
	
	
	
	/**
	 *  checks if a given step crosses this cylinder, taking account
	 *  of all current chords where the cylinder abutts others
	 *  
	 *  @param walkerPos the position of the walker making the step
	 *  @param step the step vector
	 *  
	 *  @return true if crossing, otherwise false
	 */
	public final boolean crosses(double[] walkerPos, 
			double[] step, double[] normal, double[] d, boolean skipCurrent, 
			double origLength, double[] intDist, boolean[] in, double[] p){
				
		// start and end points of current step projected into plane in axis coords
		double[] w1= new double[]{walkerPos[0]-P[0], walkerPos[1]-P[1], 0.0};
		double[] w2= new double[]{walkerPos[0]-P[0]+step[0], 
								  walkerPos[1]-P[1]+step[1], 
								  0.0};
		
		double[] newPos= new double[D];
		

		
		// distance of start and end points from cylinder axis
		double d1=Math.sqrt(w1[0]*w1[0] + w1[1]*w1[1]);
		double d2=Math.sqrt(w2[0]*w2[0] + w2[1]*w2[1]);

		
		// if skipCurrent is false, this is an inital call
		// before the reflection loop, so we need to update 
		// the initial position of walker and the initially
		// in flag for future skip checks
		if(!skipCurrent){
			initiallyIn=inside(w1);
			in[0]=initiallyIn;
		}

		
		// if start and end points and both inside the cylinder, 
		// we're definitely not crossing the circumference
		if((d1<=r)&&(d2<=r)){
			// we're indside the circle, just check the chords and go

			int chordNum=0;
			for(Iterator<Chord> chordIt=chords.iterator(); chordIt.hasNext(); chordNum++){
				
				Chord chord= (Chord)chordIt.next();

				if(skipCurrent){
					if(chordNum==skipChord){
						continue;
					}
				}
				
				// chord will handle its own geometry
				if(chord.crossedBy(w1, step, normal, d, skipCurrent, origLength, P, intDist)){
					skipChord=chordNum;
					
					SquashyCylinder other=cylsForChords.get(chord);
					other.ignoreIntersectionWith(this);
					p[0]= this.p;
					return true;
				}
				
			}
			
			// if no chord intersections, there's no crossing
			skipChord=-1;
			return false;
		}
		
		
		// otherwise there potentially is an intersection
		// get the intersection parameters
		double[] root=super.getIntersectionRoots(walkerPos, step);
		
		// no real roots = no intersections
		if(root==null){
			return false;
		}
		
		// both roots less than zero = no intersections
		if((root[0]<0.0) && (root[1]<0.0)){
			return false;
		}

		// if first root is out of range but second is in, swap them
		if(!((root[0]>=0.0)&&(root[0]<=1.0))){
			if((root[1]>=0.0)&&(root[1]<=1.0)){
				double temp=root[0];
				
				root[0]=root[1];
				root[1]=temp;				
			}
		}
		
		
		// if both roots in range, make sure lowest one
		// is first so that correct intersection is used
		if((root[0]>=0.0)&&(root[0]<=1.0)){
			if((root[1]>=0.0)&&(root[1]<=1.0)){
				if(root[0]>root[1]){
					double temp=root[0];
					
					root[0]=root[1];
					root[1]=temp;
				}
			}
		}
				
		// calculate step length
		double stepLen=0.0;
		for(int i=0; i<2; i++){
			stepLen+=step[i]*step[i];
			newPos[i]=walkerPos[i]+step[i];
		}
		stepLen=Math.sqrt(stepLen);
				
		// calculate intersection angles
		double[] thetas= new double[2];
		int numThetas=0;
		
		// check angle of interrsection point around circumference
		// of cylinder cross-section circle
		for(int i=0; i<root.length; i++){
			double[] intPoint=new double[]{walkerPos[0]+root[i]*step[0], 
										   walkerPos[1]+root[i]*step[1], 
										   walkerPos[2]+root[i]*step[2]};
			double[] cylPos=getPosition();
			double theta= Math.atan2(intPoint[1]-cylPos[1], intPoint[0]-cylPos[0]);

			thetas[numThetas++] = mapToCircle(theta);	
			//thetas[numThetas++] = Math.atan2(intPoint[1]-cylPos[1], intPoint[0]-cylPos[0]);
		}

		boolean[] inChordRange={false, false};
		
		// check if intersection is in any of the chord intervals
		for(int i=0; i<numThetas; i++){
			double theta = thetas[i];
		
			
			int chordNum=0;
			for(Iterator<AnglePair> apInt= anglePairs.iterator(); apInt.hasNext(); chordNum++){
				AnglePair ap= (AnglePair)apInt.next();
				
				if(ap.contains(theta)){
					inChordRange[i]= true;
					Chord chord=chordsForAp.get(ap);
					
					if(skipCurrent){
						if(chordNum==skipChord){
							continue;
						}
					}

					
					if(chord.crossedBy(w1, step, normal, d, skipCurrent, origLength, P, intDist)){
						skipChord=chordNum;

						SquashyCylinder other=cylsForChords.get(chord);
						other.ignoreIntersectionWith(this);
						return true;
					}	
				}
			}
			// if the first root is in [0,1) and doesn't intersect a chord
			// we should leave the loop and ignor the second root.
			if(i==0){
				if(inChordRange[i]==false){
					break;
				}
			}
		}

		// if we're not in a chord interval, check intersection
		// with circumference.
		skipChord=-1;
		
		for(int i=0; i<root.length; i++){			
			if(!inChordRange[i]){
								
				if((root[i]>0.0)&&(root[i]<=1.0)){
					
					// set normal and distance from origin in normal direction
					double[] intPoint= new double[]{walkerPos[0]+root[i]*step[0]-P[0], 
													walkerPos[1]+root[i]*step[1]-P[1], 
													walkerPos[2]+root[i]*step[2]-P[2]}; 
					
					double theta=Math.atan2(intPoint[1], intPoint[0]);

					double newNormal[] = new double[normal.length];
					
					double newD;
					
					// set normal using radial unit vector
					newNormal[0]=Math.cos(theta);
					newNormal[1]=Math.sin(theta);
					newNormal[2]=0.0;
					
					newD=0.0;
					for(int j=0; j<2; j++){
						newD+=(intPoint[j]+P[j])*newNormal[j];
					}
					
					
					// check if we need to skip this intersection
					if(skipCurrent){
						
						initiallyIn=in[0];
						
						boolean skipIt=checkSkipping(newPos);
						
						// if we're skipping, jump to the next root in the i loop
						if(skipIt){
							if(skipChord!=-1){
								continue;
							}
						}
					}
				
					// if we've got here then we need to return 
					// the distance and normal to the amending
					// routines
					d[0]=newD;
					intDist[0]=0.0;
					for(int j=0; j<normal.length; j++){
						normal[j]=newNormal[j];
					}
					intDist[0]=root[i];
					
					
					
					// if skipCurrent is false, this is an inital call
					// before the reflection loop, so we need to update 
					// the initial position of walker and the initially
					// in flag for future skip checks
					if(!skipCurrent){
						for(int j=0; j<D; j++){
							initialPos[j]=w1[j];
						}
						if(inside(walkerPos)){
							initiallyIn=true;
							in[0]=true;
						}
						else{
							initiallyIn=false;
							in[0]=false;
						}
					}
					p[0]=this.p;
					return true;
				}
			}
		}
		
		
		return false;
	}
	

	/**
	 * instructs the current cylinder to ignore the chord associated with the 
	 * given cylinder. This prevents the "catapult" effect due to chords 
	 * being defined more thasn once.
	 * 
	 * @param cyl
	 */
	protected final void ignoreIntersectionWith(SquashyCylinder cyl){
		AnglePair ap=apsForCylinders.get(cyl);
		Chord chord=chordsForAp.get(ap);
		int chordIndex= chords.indexOf(chord);
		
		skipChord=chordIndex;
	}

	/**
	 * checks if the current intersection should be ignored.
	 * The safest way to do this is to use the fact that the
	 * cylinder is circular and to check the initially inside
	 * flag, which is true if the walker starts inside the 
	 * cylinder, and compare it to whether the end of the step 
	 * is inside or outside.
	 * 
	 * Given that skipping is always true when this method is
	 * called, and the fact that the cylinder is circular it means
	 * that the current intersection should be ignored if the
	 * flags are equal (start and finish inside or start and 
	 * finish outside) because in that case we're sitting on the
	 * circumference when this method is called.
	 * 
	 * simple routine, very complex reasoning!
	 * 
	 *  @param newPos the new position after the step is made
	 */
	private final boolean checkSkipping(double[] newPos){
		
		boolean willBeIn = inside(newPos);
		
		if(initiallyIn==willBeIn){
			return true;
		}
		else{
			return false;
		}
	}
	
	/**
	 *  gets points of intersection on radii of abutting cylinders
	 *  
	 *  @param that the other cylinder
	 *  
	 *  @return new AnglePair object of the intersection
	 */
	private AnglePair getIntersesctionPair(SquashyCylinder that){
		
		// solve the problem along an x-axis connecting the cylinders
		// with the current cylinder at the origin
		double d = this.getDistanceFrom(that.getPosition());
		double R = this.getRadius();
		double r = that.getRadius();
		
		// do the geometry (from mathworld, circle-circle intersection article)
		double x = (d*d - r*r + R*R)/(2.0*d);
		double y = Math.sqrt(R*R-x*x);

		double lower = Math.atan2(y, x);
		double upper = Math.atan2(-y, x);
		
		if(lower>upper){
			double tmp=lower;
			lower=upper;
			upper=tmp;
		}

		// now rotate the angles through an offset angle, which allows
		// for the fact that the axis we've just solved on is in general
		// at an angle to the x-axis that defines the angle system on the
		// circumference of the circle
		double[] thisPos= this.getPosition();
		double[] thatPos= that.getPosition();
		double[] dispVec= new double[]{thatPos[0]-thisPos[0], thatPos[1]-thisPos[1]};
		
		
		//System.err.println("\t\tdispVec= ("+dispVec[0]+", "+dispVec[1]+")");
		// offset is the angle of the displacement vector 
		double offset=Math.atan2(dispVec[1], dispVec[0]);
		
		// add offset to lower and upper
		lower+=offset;
		upper+=offset;

		lower=mapToCircle(lower);
		upper=mapToCircle(upper);
		
		if(lower>upper){
			double tmp=lower;
			lower=upper;
			upper=tmp;
		}

		
		// construct anglePair
		return new AnglePair(lower, upper, R, that);
	}
	
	
        /** 
	 * maps an angle into the range 0 to 2 pi
	 * 
	 * @param theta angle
	 * @return number in [0, 2pi)
	 */
	private final double mapToCircle(double theta){
		
		final double twoPi= 2.0*Math.PI;
		
		double windingNum= Math.floor(theta/twoPi);
		double offset=twoPi*windingNum;
		
		double retval=  theta-offset;
		
		return retval;
	}



	/**
	 * gets a new chord for a given anglePair. includes
	 * checking against other chords already present and
	 * update of the tmin and tmax of each chord (ahem)
	 * accordingly.
	 * 
	 * @param anglePair the pair of angles defining the intersection
	 * 
	 * @return a new Chord with normal, tmin and tmacx
	 */
	private Chord getChord(AnglePair anglePair){
	
		// construct the chord from the endpoints
		double[] lowerCoord= anglePair.getLowerCoords();
		double[] upperCoord= anglePair.getUpperCoords();
		
		double xmin= lowerCoord[0];
		double ymin= lowerCoord[1];
		
		double xmax= upperCoord[0];
		double ymax= upperCoord[1];
		
		double len=Math.sqrt((xmax-xmin)*(xmax-xmin) + (ymax-ymin)*(ymax-ymin));
		
		double tmin=0.0;
		double tmax=len;
	
		Chord newChord= new Chord(anglePair, tmin, tmax);
		
		// now amend end points of chord to take account of
		// overlaps with other chords
		int chordCounter=0;
		for(Iterator<AnglePair> it= anglePairs.iterator(); it.hasNext(); chordCounter++){
			AnglePair ap=(AnglePair)it.next();
			Chord chord= (Chord)chordsForAp.get(ap);
			
			
			if(ap.contains(anglePair.lower)){
				// amend lower end of cord param range
				double[] t=chord.getIntersectionParams(xmin, ymin, xmax, ymax);
								
				newChord.tmin=t[0];
				chord.tmax=t[1];
				
				if(tmin<0.0){
					String errMess= new String("negative intersection parameter between overlapping chords. t= "+tmin);
					logger.severe(errMess);
					throw new RuntimeException(errMess);
				}
			}
			
			if(ap.contains(anglePair.upper)){
				// amend upper end of cord
				double[] t=chord.getIntersectionParams(xmin, ymin, xmax, ymax);				
				newChord.tmax=t[0];
				chord.tmin=t[1];
				
				if(tmax>len){
					String errMess= new String("intersection parameter between overlapping chords greater than cord length. t="+tmax+", len= "+len);
					logger.severe(errMess);
					throw new RuntimeException(errMess);
				}
			}
			
		}
		
		
		return newChord;
	}
	
	private void removeZeroLengthChords(){
		
		Iterator<Chord> chIt= chords.iterator();
		
		ArrayList<Chord> toRemove= new ArrayList<Chord>();
		
		// assemble current list of current 
		// zero- or negative-length chords
		while(chIt.hasNext()){
			Chord chord= (Chord)chIt.next();
			
			if(chord.tmax-chord.tmin<=1E-14){
				toRemove.add(chord);
			}
		}
		
		if(toRemove.size()>0){
			chIt=toRemove.iterator();
			
			while(chIt.hasNext()){
				Chord chord= (Chord)chIt.next();
				AnglePair ap= apsForChord.get(chord);
				
				// remove chord and anglePair
				// (we'll leave the associations)
				chords.remove(chord);
				anglePairs.remove(ap);
			}
		}
	}
	
	
	/** adds a new intersection with the given cylinder 
	 * 
	 * @param cylinder the cylinder to intersect with
	 */
	public void addIntersectionWith(SquashyCylinder cylinder){
		
		AnglePair anglePair = getIntersesctionPair(cylinder);
		Chord chord = getChord(anglePair);
		
		// if the overlap is too great, stop them both exapnding
		if(Math.abs(anglePair.upper-anglePair.lower)>=MaxIntersectionAngle){
			this.stopExpanding();
			cylinder.stopExpanding();
		}
		
		bumps.add(cylinder);
		anglePairs.add(anglePair);
		chords.add(chord);
		chordsForAp.put(anglePair, chord);
		apsForChord.put(chord, anglePair);
		apsForCylinders.put(cylinder, anglePair);
		cylsForChords.put(chord, cylinder);
		
		removeZeroLengthChords();
	}
	
	/**
	 * checks if the current cylinder already has an intersection with 
	 * the given cylinder
	 * 
	 * @param otherCylinder
	 * 
	 * @return true or false
	 * 
	 */
	public final boolean hasIntersectionWith(SquashyCylinder otherCylinder){
		
		return bumps.contains(otherCylinder);
	}
	
	private AnglePair testGetAnglePair(double[] lowerCoord, double[] upperCoord){
		return new AnglePair(lowerCoord, upperCoord);
	}
	
	
	/**
	 * outputs coords around the circumference of the cross-section 
	 * (chords not included). the file should be open before being
	 * passed to this method, and will not be closed afterwards
	 * 
	 * @param out file to write to
	 * @throws IOException if something goes wrong with the file
	 */
	public final void drawCrossSection(BufferedWriter out) throws IOException{
		
		// number of divisions of circumference
		int MAX_PTS=100;
		
		// position of cylinder
		double[] P=getPosition();
		
		// radius of cylinder
		double r= getRadius();
		
		for(int i=0; i<=MAX_PTS; i++){
			double theta= 2.0*Math.PI*(double)i/(double)MAX_PTS;
			
			double x=P[0]+r*Math.cos(theta);
			double y=P[1]+r*Math.sin(theta);
			
			out.write(x+","+y+"\n");
		}
		
		int count=0;
		
		
		out.flush();
		
	}
	
	/**
	 * add a clone to the list.
	 * 
	 * @param clone
	 */
	public final void addToMyClones(SquashyCylinder clone){
		
		if(myClones==null){
			myClones= new ArrayList<SquashyCylinder>();
		}
		
		myClones.add(clone);
	}
	
	/**
	 * before this is called, the orginal cylinder knows about its clones
	 * but the clones don't know about each other, or about the original.
	 * This method tells the clones about the original and about the other 
	 * clones so that the graph of clonal relationships is fully connected.
	 * 
	 */
	public final void propagateCloneRelationships(){
		
		if(myClones!=null){
			for(Iterator<SquashyCylinder> cloneIt= myClones.iterator(); cloneIt.hasNext(); ){
			
				SquashyCylinder clone= cloneIt.next();
			
				// it's a clone of this one...
				clone.addToMyClones(this);
				
				/// ... and all the other clones except itself
				for(Iterator<SquashyCylinder> otherCloneIt= myClones.iterator(); otherCloneIt.hasNext(); ){
					
					SquashyCylinder otherClone= otherCloneIt.next();
					
					if(otherClone!=clone){
						clone.addToMyClones(otherClone);
					}
				}
			}
		}
	}
	
	/**
	 * a little class that is just a pair of angles that define the
	 * range of angles which 
	 * 
	 * @author matt
	 *
	 */
	protected class AnglePair{
		
		private final double lower;
		
		private final double upper;

		private final double[] lowerCoord;
		
		private final double[] upperCoord;

                private final boolean midPointInOtherCyl;
		
		protected AnglePair(double lower, double upper, double r, SquashyCylinder otherCyl){
						
		    /*if(lower>upper){
				double tmp=lower;
				lower=upper;
				upper=tmp;
				}*/
			
			this.upper=upper;
			this.lower=lower;
			
			this.lowerCoord= new double[]{r*Math.cos(lower), r*Math.sin(lower)};
			this.upperCoord= new double[]{r*Math.cos(upper), r*Math.sin(upper)};
			
						// coords to the midpoint on the circumference
			double midAngle= (lower+upper)/2.0;
			double[] midPoint= new double[]{P[0]+r*Math.cos(midAngle), P[1]+r*Math.sin(midAngle), 0.0};
			
			midPointInOtherCyl=otherCyl.inside(midPoint);

		}
		
		
		/**
		 * CONSTRUCTOR FOR TESTING PURPOSES ONLY. ingores angles and allows
		 * endpoint coords to be specified directly. This allows arbitrary cords
		 * to be easily constructed.
		 * 
		 * @param lowerCoord
		 * @param upperCoord
		 */
		private AnglePair(double[] lowerCoord, double[] upperCoord){
			
			this.lower=0.0;
			this.upper=1.0;
			
			this.lowerCoord= new double[2];
			this.upperCoord= new double[2];
			
			for(int i=0; i<2; i++){
				this.lowerCoord[i]= lowerCoord[i];
				this.upperCoord[i]= upperCoord[i];
			}

			midPointInOtherCyl= true;
		}
		
		/**
		 * returns euclidean coords of intersection point associated with 
		 * lower angle of intersection angle pair
		 * 
		 * @return the coords of the the lower interaction pair (same instance)
		 */
		protected double[] getLowerCoords(){
			
			return lowerCoord;
		}
		
		/**
		 * returns euclidean coords of intersection point associated with 
		 * upper angle of intersection angle pair
		 * 
		 * @return the coords of the the upper interaction pair (same instance)
		 */		
		protected double[] getUpperCoords(){
			
			return upperCoord;
		}
		
		/** does this range contain the current angle?
		 * 
		 * @param theta the angle to check
		 * @return true if in range, false otherwise
		 */
		private final boolean contains(double theta){

		    if(midPointInOtherCyl){
			return (theta>=lower)&&(theta<=upper);
		    }
		    else{
			return !((theta>=lower)&&(theta<=upper));
		    }
		}
	}
	
	/**
	 * little class containing the chord
	 * 
	 * @author matt
	 *
	 */
	public class Chord{
		
		private final double[] normal;
		
		public double tmin;
		
		public double tmax;
		
		private final Vector3D P;
		
		private final double[] cordVec;
		
		/**
		 * construct a new Chord between the specified angles
		 * 
		 * TODO: Chord constructor does not construct correct normal of position
		 * 
		 * @param angles the anglePair object for the chord
		 * @param tmin the min length param on the chord
		 * @param tmax the max length param on the chord
		 * @param R radius of parent cylinder
		 */
		protected Chord(AnglePair anglePair, double tmin, double tmax){
			
			double[] lowerCoord= anglePair.getLowerCoords();
			double[] upperCoord= anglePair.getUpperCoords();
			
			double xmin= lowerCoord[0];
			double ymin= lowerCoord[1];
			
			double xmax= upperCoord[0];
			double ymax= upperCoord[1];

			cordVec= new double[]{xmax-xmin, ymax-ymin};
			
			double len=Math.sqrt(cordVec[0]*cordVec[0] + cordVec[1]*cordVec[1]);
			
			cordVec[0]/=len;
			cordVec[1]/=len;
			
			normal=new double[]{cordVec[1], -cordVec[0], 0.0};
			
			P=new Vector3D(xmin, ymin, 0.0);
			
			this.tmin=tmin;
			
			this.tmax=tmax;
			
			if(Double.isNaN(tmin)){
				System.err.println("tmin is NaN in chord "+this);
			}
			if(Double.isNaN(tmax)){
				System.err.println("tmax is NaN in chord "+this);
			}
		}
		
		/**
		 * checks if a step crosses this chord
		 *  
		 * @param walkerPos the position of walker making the step
		 * @param step the step vector
		 * @param normal space top store the normal
		 * @param length of original step
		 * @param cylPos position of parent cylinder
		 * @param array of length one for passing back the fraction of arc to intersection
		 * @return true if crosses, false otherwise
		 */
		private boolean crossedBy(double[] walkerPos, double[] step, double[] normal, double[] d, boolean skipCurrent, double origLength, double[] cylPos, double[] intDist){

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

			double[] intPoint= new double[2];
			
			for(int i=0; i<2; i++){
				intPoint[i]=walkerPos[i]+arcLen*step[i];
				intersectionDotNormal+=(walkerPos[i]+arcLen*step[i]+cylPos[i])*this.normal[i];
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
			if(t[1]<tmin){
				return false;
			}
			
			if(t[1]>tmax){
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
		 * calculates the arclength to the intersection point between two lines
		 * in terms of the arclengths of both lines. This is annoyingly difficult!
		 *  
		 * see mathworld line-line intersection article for explanation of symbols 
		 * and method. @link http://mathoworld.wolfram.com/Line-LineIntersection.html
		 * 
		 * Solution requires evaluating determinanats of 2x2 matrics. Throughout this code
		 * these are denoted as the ordered quadruple {a, b, c, d}, representing the 2x2 
		 * matrix
		 * 
		 * a b
		 * c d
		 * 
		 * with determinant det(a, b, c, d) = ad-bc
		 * 
		 * coords are given by
		 * 
		 * x = det(det(x1, y1, x2, y2), x1-x2, det(x3, y3, x4, y4) x3-x4)/ D
		 * y = det(det(x1, y1, x2, y2), y1-y2, det(x3, y3, x4, y4) y3-y4)/ D
		 * 
		 * where D=det(x1-x2, y1-y2, x3-x4, y3-y4)
		 * 
		 * @return intersection parameters {given line arclenth, this line arc length}
		 */
		protected double[] getIntersectionParams(double xmin, double ymin, double xmax, double ymax){
			
			// lower end of given line
			double x1=xmin;
			double y1=ymin;
			
			// upper end of given line
			double x2=xmax;
			double y2=ymax;
			
			// lower end of this chord
			double x3=P.x+tmin*cordVec[0];
			double y3=P.y+tmin*cordVec[1];
			
			// upper end of this chord
			double x4=P.x+tmax*cordVec[0];
			double y4=P.y+tmax*cordVec[1];
			
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
			double dp=intVec[0]*(x2-x1) + intVec[1]*(y2-y1);
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
			intVec[0] = x-P.x;
			intVec[1] = y-P.y;
			
			// arclength along current curve
			double t_this=Math.sqrt(intVec[0]*intVec[0] + intVec[1]*intVec[1]);
			
			// return both arclengths
			return new double[]{t_that, t_this};
			
		}
		
	}

	

	
	
	/**
	 * tests the geomtery of intersection
	 *
	 */
	private static final void testEqualSizeIntersection(){
		
		double r= 1.0;
		double d=3.0*r/2.0;
		
		double[] P1 = new double[]{0.0, 0.0, 0.0};
		double[] P2 = new double[]{d, 0.0, 0.0};
		
		SquashyCylinder squashy1 = new SquashyCylinder(P1, r, 0.0);
		SquashyCylinder squashy2 = new SquashyCylinder(P2, r, 0.0);

		double dist = squashy1.getDistanceFrom(squashy2.getPosition());

		
		System.err.println("testing overlapping cylinders of equal size:");
		
		System.err.println("\tdist = "+dist+", d="+d);

		AnglePair anglePair=squashy1.getIntersesctionPair(squashy2);

		System.err.println("\tanglePair from squashy1 is ("+anglePair.lower+", "+anglePair.upper+")");
		
		System.err.println("\twhich is ("+r*Math.sin(anglePair.lower)+", "+r*Math.cos(anglePair.lower)
						+") and ("+r*Math.sin(anglePair.upper)+", "+Math.cos(anglePair.upper)+")");
		
		
		System.err.println("\ttesting chord construction");
		Chord chord=squashy1.getChord(anglePair);
	
		double[] chordVec= chord.cordVec;
		double[] chordNormal= chord.normal;
		
		System.err.println("\t\tChord vector is ("+chordVec[0]+", "+chordVec[1]+")");
		System.err.println("\t\tChord normal is ("+chordNormal[0]+", "+chordNormal[1]+")");
	}

	/**
	 * tests the geomtery of intersection
	 *
	 */
	private static final void testUnequalSizeIntersection(){
		
		double r1= 1.0;
		double r2= 0.6;
		double d=3.0/2.0;
		
		double[] P1 = new double[]{0.0, 0.0, 0.0};
		double[] P2 = new double[]{d, 0.0, 0.0};
		
		SquashyCylinder squashy1 = new SquashyCylinder(P1, r1, 0.0);
		SquashyCylinder squashy2 = new SquashyCylinder(P2, r2, 0.0);

		double dist = squashy1.getDistanceFrom(squashy2.getPosition());

		
		System.err.println("testing overlapping cylinders of unequal sizes:");
		
		System.err.println("\tr1 = "+r1+"    r2 = "+r2);
		System.err.println("\tdist = "+dist+", d="+d);

		AnglePair anglePair=squashy1.getIntersesctionPair(squashy2);

		System.err.println("\tanglePair from squashy1 is ("+anglePair.lower+", "+anglePair.upper+")");
		
		double[] lowerCoords=anglePair.getLowerCoords();
		double[] upperCoords=anglePair.getUpperCoords();
		
		System.err.println("\twhich is ("+lowerCoords[0]+", "+lowerCoords[1]
						+") and ("+upperCoords[0]+", "+upperCoords[1]+")");
		
		
		System.err.println("\ttesting chord construction");
		
		Chord chord=squashy1.getChord(anglePair);

		double[] chordVec= chord.cordVec;
		double[] chordNormal= chord.normal;
		
		System.err.println("\t\tChord vector is ("+chordVec[0]+", "+chordVec[1]+")");
		System.err.println("\t\tChord normal is ("+chordNormal[0]+", "+chordNormal[1]+")");
		
	}

	
	private static void testChordIntersection(){
		
		double A=3.0;
		double B=2.0;
		
		double[] lowerCoord1= new double[]{A, B+0.75};
		double[] upperCoord1= new double[]{A, B-0.4};
		
		double[] lowerCoord2= new double[]{A-0.97, B};
		double[] upperCoord2= new double[]{A+0.5, B};
		
		SquashyCylinder shcyl= new SquashyCylinder(new double[]{0.0, 0.0, 0.0}, 1.0, 0.0);
		
		AnglePair anglePair1 = shcyl.testGetAnglePair(lowerCoord1, upperCoord1);
		AnglePair anglePair2 = shcyl.testGetAnglePair(lowerCoord2, upperCoord2);
		
		Chord chord1 = shcyl.getChord(anglePair1);
		Chord chord2 = shcyl.getChord(anglePair2);
		
		double[] t=chord2.getIntersectionParams(lowerCoord1[0], lowerCoord1[1], upperCoord1[0], upperCoord1[1]);
		double[] s=chord1.getIntersectionParams(lowerCoord2[0], lowerCoord2[1], upperCoord2[0], upperCoord2[1]);
		
		System.err.println("t1= "+t[0]+"   t2= "+t[1]);
		System.err.println("s1= "+s[0]+"   s2= "+s[1]);		
		
		
	}
	

	private static void testChordAmendment(){
		
		double r1= 1.0;
		double r2= 0.6;
		double d=3.0/2.0;
		
		double[] P1 = new double[]{0.0, 0.0, 0.0};
		double[] P2 = new double[]{d, 0.0, 0.0};
		double[] P3 = new double[]{d/2.0, r2, 0.0};
		
		SquashyCylinder squashy1 = new SquashyCylinder(P1, r1, 0.0);
		SquashyCylinder squashy2 = new SquashyCylinder(P2, r2, 0.0);
		SquashyCylinder squashy3 = new SquashyCylinder(P3, r2, 0.0);
		
		// add interesection with cylinder2
		squashy1.addIntersectionWith(squashy2);
		
		// compare internal anglePair with direct version
		AnglePair intAnglePair12 = squashy1.apsForCylinders.get(squashy2);
		AnglePair anglePair12=squashy1.getIntersesctionPair(squashy2);

		System.err.println("intAnglepair12= "+intAnglePair12.lower+", "+intAnglePair12.upper);
		System.err.println("anglepair12   = "+anglePair12.lower+", "+anglePair12.upper);
		
		// intersection chord between 1 & 2 without 3
		Chord chord12=squashy1.chordsForAp.get(intAnglePair12);
		
		System.err.println("chord between 1 & 2 without 3 is");
		System.err.println("chord12 vec    = ("+chord12.cordVec[0]+", "+chord12.cordVec[1]+")");
		System.err.println("chord12 normal = ("+chord12.normal[0]+", "+chord12.normal[1]+")");
		System.err.println("chord12 limits: tmin= "+chord12.tmin+" tmax= "+chord12.tmax);
		
		// same for second cylinder
		squashy1.addIntersectionWith(squashy3);
		
		AnglePair intAnglePair13a=squashy1.apsForCylinders.get(squashy3);
		AnglePair anglePair13a=squashy1.getIntersesctionPair(squashy3);
		
		System.err.println("intAnglepair13= "+intAnglePair13a.lower+", "+intAnglePair13a.upper);
		System.err.println("anglepair13   = "+anglePair13a.lower+", "+anglePair13a.upper);

		Chord chord13a= squashy1.chordsForAp.get(intAnglePair13a);

		System.err.println("chords after 3rd cylinder added are:");
		System.err.println("chord13 vec    = ("+chord13a.cordVec[0]+", "+chord13a.cordVec[1]+")");
		System.err.println("chord13 normal = ("+chord13a.normal[0]+", "+chord13a.normal[1]+")");
		System.err.println("chord13 limits: tmin= "+chord13a.tmin+" tmax= "+chord13a.tmax);
	
		System.err.println();
		System.err.println("chord12 vec    = ("+chord12.cordVec[0]+", "+chord12.cordVec[1]+")");
		System.err.println("chord12 normal = ("+chord12.normal[0]+", "+chord12.normal[1]+")");
		System.err.println("chord12 limits: tmin= "+chord12.tmin+" tmax= "+chord12.tmax);

		// now reset the first cylinder and do the same thing in another order
		squashy1= new SquashyCylinder(P1, r1, 0.0);
		
		squashy1.addIntersectionWith(squashy3);

	}

	
	private static void test2DIntersectionRoots(){
	
		double[] P=new double[]{1.0, 1.0, 0.0};
		double r=2.0;
		
		SquashyCylinder cyl= new SquashyCylinder(P, r, 0.0);
		
		double[] pos=new double[]{1.0, 2.9, 0.0};
		double[] step=new double[]{0.0, 0.2, 0.0};
		
		double[] roots=cyl.get2dIntersectionRoots(pos, step);
		
		System.err.println("roots= ("+roots[0]+" "+roots[1]+")");
		
	}
	
	
	private static void testExtraCellularReflection(){
	
		int tmax=10;
		
		double r=1.0;
		double[] P= new double[]{0.0, 0.0, 0.0};
		
		double[] startPos = new double[]{2.0, 0.0, 0.0};
		double[] step = new double[]{-1.5, 0.0, 0.0};
		
		SquashyCylinder cyl = new SquashyCylinder(P, r, 0.0);
		
		double[] root= cyl.getIntersectionRoots(startPos, step);

		System.err.println("root[0]= "+root[0]);
		System.err.println("root[1]= "+root[1]);
		
		
		double[] normal= new double[]{0.0, 0.0, 0.0};
		double[] d= new double[1];
		double[] intDist= new double[1];
		
		boolean[] in= new boolean[]{false}; 
		double[] intP= new double[1];
		
		boolean crosses = cyl.crosses(startPos, step, normal, d, false, 1.5, intDist, in, intP);
		
		System.err.println("crosses = "+crosses);
		System.err.println("normal = ("+normal[0]+", "+normal[1]+", "+normal[2]+")");
		System.err.println("d = "+d[0]);
		System.err.println("intDist = "+intDist[0]);
		
		
		SimulableScheme scheme;
		try{
	           URI uri= DiffusionSimulation.class.getResource("/test/bmx7_ED.scheme1").toURI();
	           
	           String path= uri.getPath();
	           
	           scheme= (SimulableScheme)RectGradSteTanScheme.readScheme(path);
	       }
	       catch(URISyntaxException urise){
	           throw new LoggedException(urise);
	    }
		
		
		CL_Initializer.gamma_k=1.0;
		CL_Initializer.gamma_beta=1.0;
		SimulationParams.sim_cyl_dist_size=1;
		CL_Initializer.numVoxels=1;
		SimulationParams.sim_inflamm_increments=1;
		SimulationParams.sim_L=3.0;
		
		
		
		SimulationParams simParams= new SimulationParams(10, 100, 1.0, 2, 
								SubstrateFactory.SubstrateType.CYL_1_INFLAM, 
								StepType.FIXEDLENGTH, 1.0, scheme);
		
		Substrate substrate = SubstrateFactory.getSubstrate(SubstrateFactory.SubstrateType.CYL_1_INFLAM, simParams);
		
		Walker walker = new Walker(startPos);
		
		double[] toBarrier= new double[3];
		double[] amended= new double[3];
		double[] unamended= new double[3];
		
		substrate.testAmendment(walker, step, normal, d, 1.5, toBarrier, amended, unamended);
		
		System.err.println("startpos= ("+startPos[0]+", "+startPos[1]+", "+startPos[2]+")");
		
		System.err.println("after amendment:");
		System.err.println("normal = ("+normal[0]+", "+normal[1]+", "+normal[2]+")");
		System.err.println("toBarrier =("+toBarrier[0]+", "+toBarrier[1]+", "+toBarrier[2]+")");
		System.err.println("unamended = ("+unamended[0]+", "+unamended[1]+", "+unamended[2]+")");		
		System.err.println("amended = ("+amended[0]+", "+amended[1]+", "+amended[2]+")");
	}

	
	/**
	 * tests the geomtery of intersection
	 *
	 */
	private static final void testNonHorizUnequalIntersection(){
		
		double r1= 1.0;
		double r2= 0.1;
		double d= 1.07;
		
		double[] P1 = new double[]{0.0, 0.0, 0.0};
		double[] P2 = new double[]{0.0, -d, 0.0};
		
		SquashyCylinder squashy1 = new SquashyCylinder(P1, r1, 0.0);
		SquashyCylinder squashy2 = new SquashyCylinder(P2, r2, 0.0);

		double dist = squashy1.getDistanceFrom(squashy2.getPosition());

		System.err.println("testing overlapping cylinders of unequal sizes separated on non-horizontal axis:");
		
		System.err.println("\tr1 = "+r1+"    r2 = "+r2);
		System.err.println("\tdist = "+dist+", d="+d);

		AnglePair anglePair=squashy1.getIntersesctionPair(squashy2);

		System.err.println("\tanglePair from squashy1 is ("+anglePair.lower+", "+anglePair.upper+")");
		
		double[] lowerCoords=anglePair.getLowerCoords();
		double[] upperCoords=anglePair.getUpperCoords();
		
		System.err.println("\twhich is ("+lowerCoords[0]+", "+lowerCoords[1]
						+") and ("+upperCoords[0]+", "+upperCoords[1]+")");
		
		
		System.err.println("\ttesting chord construction");
		
		Chord chord=squashy1.getChord(anglePair);

		double[] chordVec= chord.cordVec;
		double[] chordNormal= chord.normal;
		
		System.err.println("\t\tChord vector is ("+chordVec[0]+", "+chordVec[1]+")");
		System.err.println("\t\tChord normal is ("+chordNormal[0]+", "+chordNormal[1]+")");
		
	}

	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// tests various bits of geometry
		//SquashyCylinder.testEqualSizeIntersection();

		//SquashyCylinder.testUnequalSizeIntersection();
		
		//SquashyCylinder.testChordIntersection();

		//SquashyCylinder.testChordAmendment();
	
		//SquashyCylinder.test2DIntersectionRoots();
	
		//SquashyCylinder.testExtraCellularReflection();
		
		SquashyCylinder.testNonHorizUnequalIntersection();
	}

	

}
