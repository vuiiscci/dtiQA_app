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
 * Converts Analyze / NIFTI / MHA files to voxel order.
 * Either takes a 4D file (all measurements in single image)
 * or a list of 3D images.
 *
 * @author Philip Cook
 * @version $Id$
 *
 */
public class ImageToVoxel extends Executable{
    
    public ImageToVoxel(String[] args){
        super(args);
    }
    
    private String singleVolumeFile;
    private String imageListFile;
    private String imagePrefix;
    private	boolean imagePrefixSet;
    
    /**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.apps.ImageToVoxel");
    
    /**
     * Output manager
     */
//private static OutputManager om;
    
    public void initDefaultVals(){
        singleVolumeFile = null;
        imageListFile = null;
        imagePrefix = "";
        imagePrefixSet = false;
    }
    
    public void initOptions(String[] args){
        CL_Initializer.inputDataType = "float";
        OutputManager.outputDataType = "float";
        CL_Initializer.CL_init(args);

        if (CL_Initializer.inputFile != null) {
            if (ImageHeader.imageExists(CL_Initializer.inputFile)) {
                singleVolumeFile = CL_Initializer.inputFile;
            }
            else {
                // no volume matches inputFile, assume it's an image list
                imageListFile = CL_Initializer.inputFile;
            }
        }
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-4dimage") || args[i].equals("-4dimageroot")) {
                singleVolumeFile = args[i+1];
                CL_Initializer.markAsParsed(i,2);
            }
            if (args[i].equals("-imagelist")) {
                imageListFile = args[i+1];
                CL_Initializer.markAsParsed(i,2);
            }
            if (args[i].equals("-imageprefix")) {
                imagePrefix = args[i+1];
                CL_Initializer.markAsParsed(i,2);
                imagePrefixSet = true;
            }
            if (args[i].equals("-ignorespmscale")) {
                logger.warning("Option -ignorespmscale is deprecated. Scaling is applied. Modify the " +
                "header to alter scaling");
                CL_Initializer.markAsParsed(i);
            }
        }
        
        CL_Initializer.checkParsing(args);
    }
    
    
    public void initVariables() {
    }
    
    public void execute(OutputManager om){
        if (singleVolumeFile != null) {
            // Construct the data source.
            ImageHeader ih = null;
            
            try {
                ih = ImageHeader.readHeader(singleVolumeFile);
            }
            catch (IOException e) {
                throw new LoggedException("Cannot read image matching " + singleVolumeFile);
            }
            DataSource data = ih.getImageDataSource();
            // Loop over the data
            while (data.more()) {
                om.output(data.nextVoxel());
            }
            // Tidy up.
            om.close();
        }
        else {
            if (!imagePrefixSet) {
                // default prefix is the path to the imageListFile
                // Cygwin uses / but Windows Java will report \ as separator. Thus we always use /
                String slashie = "/";
                int index = imageListFile.lastIndexOf(slashie);
                if (index > -1) {
                    imagePrefix = imageListFile.substring(0, index+1);
                }
            }
            
            // set up 4 meg buffers
            int bufferSize = 1024 * 1024 * 4;
            
            ExternalDataSource.FILEBUFFERSIZE = bufferSize;
            ArrayList<DataSource> streamList = new ArrayList<DataSource>();
            
            DataSource[] streams = null;
            
            int components = 0;
            
            try {
                Scanner imageList = new Scanner(new File(imageListFile));
                imageList.useDelimiter("\r\n|\n");
                while (imageList.hasNext()) {
                    String imageFile = imagePrefix.concat(imageList.next().trim());
                    streamList.add(ImageHeader.readHeader(imageFile).getImageDataSource());
                    components++;
                }
                
                imageList.close();
                streams = new DataSource[components];
                streamList.toArray(streams);
                
                // now output voxels
                while (streams[0].more()) {
                    double[] nextVoxel = new double[components];
                    for (int i = 0; i < components; i++) {
                        nextVoxel[i] = streams[i].nextVoxel()[0];
                    }
                    om.output(nextVoxel);
                }
            }
            catch (IOException e) {
                throw new LoggedException(e);
            }
            
        }
        
        om.close();
        
        
    }
    
    
}
