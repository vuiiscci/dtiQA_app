package tractography;


public class Voxel implements Comparable {


    public final int x;
    public final int y;
    public final int z;

    
    public Voxel(int i, int j, int k) {
        x = i;
        y = j;
        z = k;
    }


    public int compareTo(Object o) {
        Voxel v = (Voxel)o;
        
        if (x != v.x) {
            return x - v.x;
        } 
        if (y != v.y) { 
            return y - v.y;
        } 
        
        return z - v.z;
    }


    public boolean equals(Object o) {
        
        if (!(o instanceof Voxel)) {
            return false;
        }
        else {
            Voxel v = (Voxel)o;

            return x == v.x && y == v.y && z == v.z;
        }
        
    }
    

    public String toString() {
	return "Voxel " + x + " " + y + " " + z;
    }

    public int hashCode() {
	return 13 * x + 29 * y + 37 * z;
    }

}
