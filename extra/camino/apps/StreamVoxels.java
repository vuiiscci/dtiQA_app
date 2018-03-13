package apps;

import data.*;
import imaging.*;
import misc.*;
import tools.*;

import java.io.*;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.logging.Logger;


/**
 * Streams data from voxel-order data files. 
 *
 * @author Philip Cook
 * @version $Id$
 *
 */
public class StreamVoxels extends Executable {
    
    public StreamVoxels(String[] args){
        super(args);
    }
    
    private String imageListFile;
    private String imagePrefix;
    private boolean imagePrefixSet;
    private int components;

    /**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.apps.StreamVoxels");
    
    
    public void initDefaultVals(){
        imageListFile = null;
        imagePrefix = "";
        imagePrefixSet = false;
        components = 3;
    }
    
    public void initOptions(String[] args){
        CL_Initializer.inputDataType = "double";
        OutputManager.outputDataType = "double";

        CL_Initializer.CL_init(args);

        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-imagelist")) {
                imageListFile = args[i+1];
                CL_Initializer.markAsParsed(i,2);
            }
            if (args[i].equals("-imageprefix")) {
                imagePrefix = args[i+1];
                CL_Initializer.markAsParsed(i,2);
                imagePrefixSet = true;
            }
            if (args[i].equals("-components")) {
                components = Integer.parseInt(args[i + 1]);
                CL_Initializer.markAsParsed(i);
                CL_Initializer.markAsParsed(i + 1);
            }
        }
        
        CL_Initializer.checkParsing(args);
    }
    
    
    public void initVariables() {
    }
    
    public void execute(OutputManager om){
        
        if (!imagePrefixSet) {
            // default prefix is the path to the imageListFile
            // Cygwin uses / but Windows Java will report \ as separator. Thus we always use /
            String slashie = "/";
            int index = imageListFile.lastIndexOf(slashie);
            if (index > -1) {
                imagePrefix = imageListFile.substring(0, index+1);
            }
        }
            
        // set up 1 meg buffers
        int bufferSize = 1024 * 1024 * 1;
            
        ExternalDataSource.FILEBUFFERSIZE = bufferSize;
        ArrayList<DataSource> streamList = new ArrayList<DataSource>();
        
        DataSource[] streams = null;
        
        int numStreams = 0;
        
        try {
            Scanner imageList = new Scanner(new File(imageListFile));
            imageList.useDelimiter("\r\n|\n");
            while (imageList.hasNext()) {
                String imageFile = imagePrefix.concat(imageList.next().trim());
                streamList.add(new VoxelOrderDataSource(imageFile, components, CL_Initializer.inputDataType));
                numStreams++;
            }
            
            imageList.close();
            streams = new DataSource[numStreams];
            streamList.toArray(streams);
            
            // now output voxels
            double[] nextVoxel = new double[numStreams * components];

            while (streams[0].more()) {
              
                for (int i = 0; i < numStreams; i++) {
                    double[] streamVoxel = streams[i].nextVoxel();
                    System.arraycopy(streamVoxel, 0, nextVoxel, i * components, components);
                }
                om.output(nextVoxel);
            }
        }
        catch (IOException e) {
            throw new LoggedException(e);
        }
        
        
        
        om.close();
        
    
    }
    
    
}
