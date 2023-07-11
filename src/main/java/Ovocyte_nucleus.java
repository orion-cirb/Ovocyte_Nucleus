

import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.Duplicator;
import ij.plugin.PlugIn;
import ij.plugin.filter.Analyzer;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.ImageIcon;
import org.apache.commons.io.FilenameUtils;


/**
 *
 * @author phm
 */
public class Ovocyte_nucleus implements PlugIn {
    File outDir, inDir;
    String imageDir = "";
    String OutDir = "";
    float nucXcent, nucYcent;
    int theta = 10;
    int enlarge = 5;
    private Calibration cal = new Calibration();
    double pixelWidth = 0.1135;
    double time = 500;
    public final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
    
    /**
     * Create out and temp folders
     * @param inDir 
     */
    private void createOuPutDir() {
        OutDir = imageDir + File.separator + "Out" + File.separator;
        outDir = new File(OutDir);
        if (!outDir.isDirectory())
            outDir.mkdirs();
    }
    
    
    private float calculateDistance(float Xc, float Yc, float Xp, float Yp) {
        double dist =  Math.sqrt( Math.pow( (Xc - Xp), 2) + Math.pow( (Yc - Yp), 2));
        return((float)dist);
    }

    /**
     * Compute perimeter
     * return coordinates 
     * @param img
     * @return coord 
     */
    
    /** find centroid of the Nucleus  
    /*
    */
    public void getCentroid(ImagePlus imgBin, int t) {
        ResultsTable table = new ResultsTable();
        Analyzer measure = new Analyzer(imgBin,Measurements.CENTROID,table);
        imgBin.setT(t);
        // measure parameters
        measure.measure();
	nucXcent = (float)table.getValue("X",table.size()-1);
	nucYcent = (float)table.getValue("Y",table.size()-1);
        table.reset();

    }
    
    
    /**
     * Compute centroid distance to border for teta angle
     * @param perimeter, Xc, Yc
     * @return double []
     */
     
    private ArrayList<Float> centerToBorder(ImagePlus img, double Xc, double Yc) {
        ArrayList distance = new ArrayList();
        for (int i = 0; i < 360; i+=theta) {
            double rad = i*(Math.PI/180);            
            for (int p = 0; p < img.getWidth(); p++) {
                int ptX = (int)(cal.getRawX(Xc) + Math.cos(rad)*p);
                int ptY = (int)(cal.getRawY(Yc) - Math.sin(rad)*p);
                if (img.getProcessor().getPixelValue(ptX, ptY) == 0) {
                    distance.add(calculateDistance((float)Xc, (float)Yc, (float)cal.getX(ptX), (float)cal.getY(ptY)));
                    p =  img.getWidth(); 
                }
            }
        }
        return(distance);
    }
    
    private void crop_nucleus(ImagePlus img) {
        ImagePlus imp = new Duplicator().run(img);
        // move to middle stack
        imp.setT(imp.getNFrames()/2);
        // gaussian filter
        IJ.run(imp, "Gaussian Blur...", "sigma=10 slice");
        // Threshold nucleus
        IJ.setAutoThreshold(imp, "Huang dark");
        IJ.run(imp, "Create Selection", "");
        IJ.run(imp, "To Bounding Box", "");
        IJ.run(imp, "Enlarge...", "enlarge="+enlarge);
        Roi rect = imp.getRoi();
        img.setRoi(rect);
        IJ.run(img, "Crop", "");
        imp.flush();
        imp.close();
    }
    
    private void computeDistance (ImagePlus img, String imageName) throws IOException {
        ResultsTable distResults = new ResultsTable();
        distResults.setPrecision(3);
        for (int t = 1; t <= img.getNFrames(); t++) {
           nucXcent = 0;
           nucYcent = 0;
           int angle = 0;
           getCentroid(img, t);
           ArrayList<Float> centerBorderDist = centerToBorder(img, nucXcent, nucYcent);
           String header = "t"+t;
           for (int i = 0; i < centerBorderDist.size(); i++) { 
               distResults.setValue("tetha", i, angle);
               angle += theta;
               distResults.setValue(header, i,centerBorderDist.get(i));
            }     
        }
        // write table
        distResults.saveAs(OutDir+imageName+"_distances.xls");
        distResults.reset();
    }
    
     /**
     * Generate dialog box
     */
    public boolean dialog() {
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 80, 0);
        gd.addImage(icon);
        gd.addDirectoryField("Image folder : ", "");
        gd.addMessage("Parameters", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Theta : ", theta);
        gd.addNumericField("Enlarge before crop (µm)", enlarge);
        gd.addNumericField("XY calibration (µm)  :", pixelWidth);
        gd.addNumericField("Time interval (msec) : ", time);
        gd.showDialog();
        imageDir = gd.getNextString();
        theta = (int)gd.getNextNumber();
        enlarge = (int)gd.getNextNumber();
        cal.pixelWidth = cal.pixelHeight = gd.getNextNumber();
        cal.frameInterval = gd.getNextNumber();
        cal.setUnit("micron");
        if(gd.wasCanceled())
            return(false);
        return(true);
    }
    
   public void run(String arg) {
            if (dialog() != true) {
                 IJ.showStatus("Plugin canceled");
                return;
            }
            inDir = new File(imageDir);
            String [] imageFile = inDir.list();
            if (imageFile == null) {
                IJ.showStatus("No image found ....");
                return;
            }
            createOuPutDir();
            // Open tif files
            for (String f : imageFile) {
                if (f.endsWith(".tif")) {
                    try {
                        String imageName = FilenameUtils.getBaseName(f);
                        String tifFile = inDir + File.separator + f;
                        // write headers
                        ImagePlus img = IJ.openImage(tifFile);
                        img.setCalibration(cal);
                        img.setDimensions(1, 1, img.getNSlices());
                        // crop nucleus
                        crop_nucleus(img);
                        String imgTitle = img.getTitle();
                        // run pure denoise filter
                        
                        IJ.run(img, "PureDenoise ...", "parameters='1 1' estimation='Auto Global' ");
                        while (WindowManager.getWindow("Denoised-"+imgTitle) == null)
                            IJ.wait(10);
                        WindowManager.getWindow("Log").setVisible(false);
                        ImagePlus imgFilter = WindowManager.getCurrentImage();
                        imgFilter.setDimensions(1, 1, img.getNFrames());
                        imgFilter.setCalibration(cal);
                        img.close();
                        img.flush();
                        imgFilter.hide();
                        IJ.run(imgFilter,"8-bit","");
                        IJ.run(imgFilter,"Gaussian Blur...", "sigma=2 stack");
                        // Threshold nucleus
                        Prefs.blackBackground = false;
                        IJ.run(imgFilter, "Convert to Mask", "method=Default background=Dark calculate");
                        //Close arround centroid
                        IJ.run(imgFilter, "Options...", "interations=2 count=1 do=Close stack");
                        // fill hole
                        IJ.run(imgFilter, "Fill Holes", "stack");
                        // stack registration
                        IJ.run(imgFilter, "StackReg", "transformation=[Rigid Body]");
                        FileSaver imgSave = new FileSaver(imgFilter);
                        imgSave.saveAsTiff(OutDir+imageName+"_Mask.tif");
                        // for each time compute distance center to border print all values
                        computeDistance(imgFilter, imageName);
                        imgFilter.close();
                        imgFilter.flush();
                    } catch (IOException ex) {
                        Logger.getLogger(Ovocyte_nucleus.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
            }
            IJ.showStatus("Process done");
       } 
    }

