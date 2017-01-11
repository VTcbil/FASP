
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.awt.Color;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.NewImage;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.Wand;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import ij.process.*;

/**
 * Functional AStrocyte Phenotyping (FASP)
 * 
 * @author Yinxue Wang, Guilai Shi, David J. Miller, Yizhi Wang, Congchao Wang, Gerard Broussard, Yue Wang, Lin Tian, and Guoqiang Yu
 */


public class FASP_ implements PlugIn{
	protected ImagePlus img;
	protected ImageStack stack;
	protected short[][][] imStackArray;

	// Parameters:
	protected double threshActRegZ;
	protected double threshFitRegZ;
	protected int minSize;

	// Basic features
	protected int width;
	protected int height;
	protected int nSlices;
	protected int[] imgDims;
	// For roi analysis
	protected RoiManager roiManager;
	protected RoiManager manager;

	public void run(String arg) {
		try {
			if (showDialog())
				processFASP();
		} catch (IndexOutOfBoundsException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	public boolean showDialog() {
		GaussianLUT calGaussObj = new GaussianLUT(); 
		// Get input parameter
		GenericDialog gd = new GenericDialog("FASP - Parameter Setting");
		gd.addNumericField("Significance level", 0.0100, 4);
		gd.addNumericField("Minimum_allowed_size_of_an_FIU", 30, 0);
		gd.showDialog();
		if (gd.wasCanceled()){
			return false;
		}
		// Threshold for z-score of activity level
		double sigLevel = gd.getNextNumber();
		if(Double.isNaN(sigLevel) || (sigLevel<=0) || (sigLevel>=1)){
			IJ.showMessage("Invalid parameter(s).\n" + "Significance level must be between 0 and 1.");
			return false;
		}
		threshActRegZ = calGaussObj.myQnorm(1.0-sigLevel/10);
		// Threshold for z-score of functional uniqueness
		threshFitRegZ = 0;
		minSize = (int)Math.round(gd.getNextNumber());
		if(Double.isNaN(threshActRegZ) || Double.isNaN(threshFitRegZ) || (minSize==0)){
			IJ.showMessage("Invalid parameter(s).\n" + "Numerical input needed.");
			return false;
		}
		if(minSize<10) minSize = 10;
		return true;
	}


	public void showAbout() {
		IJ.showMessage("Functional AStrocyte Phenotyping (FASP)",
				"Detect and characterize astrocytic functionally independent units (FIUs).\n"
						+ "@ Author: Yinxue Wang (yxwang90@vt.edu)"
				);
	}



	/**
	 * Process an image.
	 * @throws IOException 
	 * @throws IndexOutOfBoundsException 
	 */
	public void processFASP() throws IndexOutOfBoundsException, IOException{
		img = WindowManager.getCurrentImage();
		imgDims = img.getDimensions();
		
		if (img.getBitDepth() != 8 && img.getBitDepth() != 16 || imgDims[2]>1){
			IJ.showMessage("Only 8-bit/16-bit grayscale, single-channel image stack supported.\n");
			return;
		}
		
		nSlices=img.getStackSize();
		if (nSlices == 1) {
			IJ.showMessage("Time-lapse_data_needed.");
			return;
		}

		if (img.getType() != ImagePlus.GRAY8 & img.getType() != ImagePlus.GRAY16){
			IJ.showMessage("Only 8-bit/16-bit grayscale image stack supported.");
			return;
		}
		if (img.getType() == ImagePlus.GRAY16){
			ImageConverter ic = new ImageConverter(img);
			ic.convertToGray8();
			img.updateAndDraw();
		}
		stack=img.getStack();
		height=img.getWidth();
		width=img.getHeight();
		IJ.showStatus(5+"% Finished!") ; 
		IJ.showProgress(5, 100);
		// Transform the image stack object into an 3D array
		imStackArray = new short[height][width][nSlices];
		int mask=0xff;
		int nPixels=width*height;
		int tempx, tempy;
		byte[] pixels;
		short p;

		for(int iSlice=1; iSlice<=nSlices; iSlice++){
			pixels = (byte[])stack.getPixels(iSlice);  // The first frame (slice);
			for (int iPix=0; iPix<nPixels; iPix++){
				tempy = (int)Math.floor((double)iPix/(double)height);
				tempx = iPix - tempy*height;
				p = (short)(mask&pixels[iPix]);
				imStackArray[tempx][tempy][iSlice-1] = p;//height,width,slice
			}
		}


				long startTime=System.nanoTime();     // Testing running time: start

		AstrocyteCaMovie caMovie1 = new AstrocyteCaMovie(imStackArray,threshActRegZ,threshFitRegZ,minSize);

				long endTime=System.nanoTime(); 	  // Testing running time: end
				System.out.println("Overall running time: "+(endTime-startTime)/1e9 +" sec"); 

		if(caMovie1.getNumFIUs()==0){
			IJ.showMessage("No FIU found. Looser thresholds might help to find potential FIUs.");
		}else{
			ROI_Show(caMovie1.getFIUidMap(),caMovie1.getCharXmat(), caMovie1.getNumFIUs());
		}

		IJ.showStatus("Complete!") ; 
		IJ.showProgress(100, 100);

	}


	public void ROI_Show(int[][] fiuIDmap_hw/*fiuIDmap's transform version*/,double[][] charXmat, int max_fiu){
		//Image close
		int[][] fiuIDmapclose = fiuIDmapclose(fiuIDmap_hw,max_fiu);
		//FIU relabel
		int[][] fiuIDmap = ConnectedComponents(fiuIDmapclose);
		//update charXmat
		max_fiu = 0;
		for(int i = 0; i<fiuIDmap.length;i++){
			for(int j = 0; j<fiuIDmap[0].length;j++){
				if(max_fiu<fiuIDmap[i][j])
					max_fiu=fiuIDmap[i][j];
			}
		}
		charXmat = charXupdate(fiuIDmap,fiuIDmapclose,charXmat,max_fiu);
		//properties of FIUs
		long [] roi_size = new long[max_fiu];
		int fiu_width=fiuIDmap[0].length;
		int fiu_height=fiuIDmap.length;

		int chax_width=charXmat[0].length;
		int chax_height=charXmat.length; //height equal to #FIUs

		int[][] roi_pt1 = new int[max_fiu][2]; //larger than enough

		for(int i = 0; i<fiu_height;i++){
			for(int j = 0; j<fiu_width;j++){
				int tmp = fiuIDmap[i][j];
				if(tmp!=0){
					roi_pt1[tmp-1][0] = i;
					roi_pt1[tmp-1][1] = j;
					roi_size[tmp-1] = roi_size[tmp-1]+1;
				}
			}
		}
		//Rebuild charXmat
		charXmat = charXRebuid(fiuIDmap, charXmat, roi_size);
		////FIU show start
		// Generate a new image 
		String imtitle = "FIU MAP";
		ImagePlus imp = NewImage.createByteImage (imtitle, fiu_width, fiu_height, 1,/*stack=1*/ 
				NewImage.FILL_WHITE);
		ImageProcessor impNP = imp.getProcessor(); 
		for(int i = 0; i<fiu_height;i++){
			for(int j = 0; j<fiu_width;j++){
				impNP.putPixel(j,i,fiuIDmap[i][j]);
			}
		}
		ImageStack impstack = imp.getStack();
		ByteProcessor ip = (ByteProcessor)impstack.getProcessor(1).convertToByte(true);
		// Generate roimanager
		RoiManager manager = RoiManager.getInstance();
		if (manager == null)
			manager = new RoiManager();
		int[] roi2fiu = new int[max_fiu];
		int roi_cnt = 0;
		for(int i=0;i<max_fiu;i++){
			if(roi_size[i]<minSize)
				continue;
			roi2fiu[roi_cnt] = i;
			roi_cnt++;
			Wand w = new Wand(ip);
			w.autoOutline(roi_pt1[i][1],roi_pt1[i][0],0.01,Wand.EIGHT_CONNECTED); 
			if (w.npoints>0) { // we have an roi from the wand... 
				Roi roi = new PolygonRoi(w.xpoints, w.ypoints, w.npoints, Roi.TRACED_ROI); 
				imp.setRoi(roi);
				manager.addRoi(roi);
			}
		}

		int numFIU = manager.getCount();
		//generate tables


		//// generate feature table		
		long total_area = 0;//for summary
		ResultsTable Ft_table = new ResultsTable();
		// Generate table of charXmat
		ResultsTable charXmat_table = new ResultsTable();
		double[][] F_charXmat = new double[chax_height][chax_width];//numFIU = chax_height
		double max_Int = 0;
		double min_Int = 100;
		double[] Max_F = new double[chax_height];
		double total_MaxF = 0;
		int FIU_cnt = 0;
		for (int i=0;i<chax_height;i++) {
			if(roi_size[i]<minSize)
				continue;
			double[] meanF = getMeanF(charXmat[i]); //if use charXAvg, do not use this use the next line
			for(int j = 0;j<chax_width;j++){
				if(max_Int<meanF[j])
					max_Int=meanF[j];
				if(min_Int>meanF[j])
					min_Int = meanF[j];
				F_charXmat[i][j] = meanF[j];
			}
			Arrays.sort(meanF);
			Max_F[i] = meanF[meanF.length-1];
			Ft_table.incrementCounter();
			FIU_cnt = FIU_cnt+1;
			Ft_table.addLabel("FIU"+FIU_cnt);
			Ft_table.addValue("area(pixel)",roi_size[i]);
			Ft_table.addValue("max DeltaF/F0",Max_F[i]);

			charXmat_table.incrementCounter();
			charXmat_table.addLabel("FIU"+FIU_cnt);
			for(int k=0; k<chax_width; k++){
				charXmat_table.addValue( "t"+(k+1) ,F_charXmat[i][k]);
			}

			total_MaxF = total_MaxF+Max_F[i];
			total_area = total_area+roi_size[i];
		}
		Ft_table.showRowNumbers(false);
		Ft_table.show("Feature Table");

		charXmat_table.showRowNumbers(false);
		charXmat_table.show("Characteristic Curve Table");

		////generate summary table
		ResultsTable sm_table = new ResultsTable();

		sm_table.incrementCounter();
		sm_table.addValue("Number of FIUs",numFIU);
		sm_table.addValue("Total area(pixel) of FIUs",total_area);
		sm_table.addValue("Average area of FIUs",(float)total_area/numFIU);

		double std_MaxF = 0;
		double std_area = 0;
		for(int i=0;i<numFIU;i++){
			if(roi_size[i]<minSize)
				continue;
			std_MaxF = std_MaxF+(Max_F[i]-total_MaxF/numFIU)*(Max_F[i]-total_MaxF/numFIU);
			std_area = std_area+(roi_size[i]-total_area/numFIU)*(roi_size[i]-total_area/numFIU);

		}
		std_MaxF = Math.sqrt(std_MaxF);
		std_area = Math.sqrt(std_area);
		sm_table.addValue("Std of areas",(float)std_area);
		sm_table.addValue("Average of max DeltaF/F0",(float)total_MaxF/numFIU);
		sm_table.addValue("Std of max DeltaF/F0",(float)std_MaxF);
		sm_table.showRowNumbers(false);
		sm_table.show("Summary Table");

		manager.runCommand("show all with labels");
		manager.setSize(300, 400);
		//listen to roi manager to show the item selected  
		double[] x_data = new double[chax_width];// X-axis
		for(int i=0;i<chax_width;i++)
			x_data[i] = i+1;
		// Select the ROIs of interest
		int former_sl=chax_width+1, now_sl=chax_width+2; //to tell if need to show again
		Plot plot;
		boolean first_plot = true;
		PlotWindow plotWin = null;
		Color c = Color.RED;
		while(true){
			IJ.wait(100);
			int[] sl_idx = manager.getSelectedIndexes();
			if(sl_idx.length>0)
				now_sl = sl_idx[0];
			if(former_sl==now_sl)
				continue;
			for(int i=0;i<sl_idx.length;i++){
				// output FIU curve of FIU : saveRoi[sl_idx[i]]

				int idx_fiu = roi2fiu[sl_idx[i]];
				double[] fiu_int=F_charXmat[idx_fiu];

				PlotWindow.noGridLines = false; // draw grid lines
				plot = new Plot("FIU curve","t","deltaF / F0",x_data,fiu_int);
				plot.setLimits(0, chax_width, min_Int, max_Int);
				plot.setLineWidth(2);
				if(first_plot){
					plot.setColor(c);
					plotWin = plot.show();
					first_plot = false;
				}
				else{
					if(plotWin.isClosed()){
						plot.setColor(c);
						plotWin = plot.show();
					}else{
						plot.setColor(c);
						plotWin.drawPlot(plot);
					}
				}
			}
			former_sl = now_sl;
		}
	}


	public double[][] charXRebuid(int[][] fiuIDmap, double[][] charXmat, long[] roi_size) {
		myBasicMath basicMathObj = new myBasicMath();
		int numFIU = roi_size.length;
		double[][] FIUcurves = new double[numFIU][nSlices];
		for(int i=0;i<nSlices;i++){
			for(int y=0;y<fiuIDmap.length;y++){
				for(int x=0;x<fiuIDmap[0].length;x++){
					if(fiuIDmap[y][x]!=0){
						FIUcurves[fiuIDmap[y][x]-1][i]+=imStackArray[x][y][i];
					}
				}
			}
		}
		for(int i=0;i<numFIU;i++){
			for(int j=0;j<nSlices;j++){
				FIUcurves[i][j] = FIUcurves[i][j]/roi_size[i];
			}
		}
		double[] fiumean = new double[numFIU];
		double[] fiustd = new double[numFIU];
		for(int i=0;i<numFIU;i++){
			fiumean[i] = basicMathObj.sampleMean(FIUcurves[i]);
			fiustd[i] = Math.sqrt(basicMathObj.sampleVar(FIUcurves[i]));
			for(int j=0;j<charXmat[0].length;j++){
				charXmat[i][j] = charXmat[i][j]*fiustd[i]+fiumean[i];
			}
		}
		return charXmat;
	}


	public double[][] charXAvg(int[][] fiuIDmap, double[][] charXmat, long[] roi_size) {
		int numFIU = roi_size.length;
		double[][] FIUcurves = new double[numFIU][nSlices];
		for(int i=0;i<nSlices;i++){
			for(int y=0;y<fiuIDmap.length;y++){
				for(int x=0;x<fiuIDmap[0].length;x++){
					if(fiuIDmap[y][x]!=0){
						FIUcurves[fiuIDmap[y][x]-1][i]+=imStackArray[y][x][i];
					}
				}
			}
		}
		for(int i=0;i<numFIU;i++)
			for(int j=0;j<nSlices;j++)
				FIUcurves[i][j] = FIUcurves[i][j]/roi_size[i];

		return charXmat;
	}


	private double[][] charXupdate(int[][] fiuIDmap, int[][] fiuIDmapclose, double[][] charXmat, int max_fiu) {
		double[][] newcharXmat = new double [max_fiu][charXmat[0].length];
		for(int i=0;i<fiuIDmap.length;i++){
			for(int j=0;j<fiuIDmap[0].length;j++){
				if(fiuIDmap[i][j]!=0){ 
					if(newcharXmat[fiuIDmap[i][j]-1][0]==0.0){
						for(int m = 0;m<newcharXmat[0].length;m++)
							newcharXmat[fiuIDmap[i][j]-1][m] = charXmat[fiuIDmapclose[i][j]-1][m];
					}
				}
			}
		}
		return newcharXmat;
	}


	/**open operation to the ID map with 4-connected
	 * only small region can open to large region**/
	private int[][] fiuIDmapclose(int[][] fiuIDmap_hw/*transform of true ID map*/, int max_fiu) {
		long [] roi_size = new long[max_fiu];

		int[][] fiuIDmap = new int[fiuIDmap_hw[0].length][fiuIDmap_hw.length]; 
		int[][] fiuIDmap_unclosed = new int[fiuIDmap_hw[0].length][fiuIDmap_hw.length]; 
		for(int i = 0; i<fiuIDmap_hw[0].length;i++){
			for(int j = 0; j<fiuIDmap_hw.length;j++){
				fiuIDmap[i][j] = fiuIDmap_hw[j][i];
				fiuIDmap_unclosed[i][j] = fiuIDmap_hw[j][i];
				if(fiuIDmap[i][j]>0)
					roi_size[fiuIDmap[i][j]-1] = roi_size[fiuIDmap[i][j]-1]+1;
			}
		}
		// Close for the output ID image
		for(int i = 0; i<fiuIDmap.length;i++){
			for(int j = 0; j<fiuIDmap[0].length;j++){
				if(fiuIDmap_unclosed[i][j]!=0){
					if(i>0){
						if(fiuIDmap[i-1][j]==0){
							fiuIDmap[i-1][j] = fiuIDmap_unclosed[i][j];
							roi_size[fiuIDmap_unclosed[i][j]-1] = roi_size[fiuIDmap_unclosed[i][j]-1]+1;
						}
						else if(fiuIDmap_unclosed[i-1][j]-1>=0  & roi_size[fiuIDmap[i-1][j]-1]<roi_size[fiuIDmap_unclosed[i][j]-1]){
							fiuIDmap[i-1][j] = fiuIDmap_unclosed[i][j];
							roi_size[fiuIDmap_unclosed[i-1][j]-1] = roi_size[fiuIDmap_unclosed[i-1][j]-1]-1;
						}
					}
					if(j>0){
						if(fiuIDmap[i][j-1]==0){
							fiuIDmap[i][j-1] = fiuIDmap_unclosed[i][j];
							roi_size[fiuIDmap_unclosed[i][j]-1] = roi_size[fiuIDmap_unclosed[i][j]-1]+1;
						}
						else if(fiuIDmap_unclosed[i][j-1]-1>=0  & roi_size[fiuIDmap[i][j-1]-1]<roi_size[fiuIDmap_unclosed[i][j]-1]){
							fiuIDmap[i][j-1] = fiuIDmap_unclosed[i][j];
							roi_size[fiuIDmap_unclosed[i][j-1]-1] = roi_size[fiuIDmap_unclosed[i][j-1]-1]-1;
						}
					}
					if(i<fiuIDmap.length-1){
						if(fiuIDmap[i+1][j]==0){
							fiuIDmap[i+1][j] = fiuIDmap_unclosed[i][j];
							roi_size[fiuIDmap_unclosed[i][j]-1] = roi_size[fiuIDmap_unclosed[i][j]-1]+1;
						}
						else if(fiuIDmap_unclosed[i+1][j]-1>=0 & roi_size[fiuIDmap[i+1][j]-1]<roi_size[fiuIDmap_unclosed[i][j]-1]){
							fiuIDmap[i+1][j] = fiuIDmap_unclosed[i][j];
							roi_size[fiuIDmap_unclosed[i+1][j]-1] = roi_size[fiuIDmap_unclosed[i+1][j]-1]-1;
						}
					}
					if(j<fiuIDmap[0].length-1){
						if(fiuIDmap[i][j+1]==0){
							fiuIDmap[i][j+1] = fiuIDmap_unclosed[i][j];
							roi_size[fiuIDmap_unclosed[i][j]-1] = roi_size[fiuIDmap_unclosed[i][j]-1]+1;
						}
						else if(fiuIDmap_unclosed[i][j+1]-1>=0 & roi_size[fiuIDmap[i][j+1]-1]<roi_size[fiuIDmap_unclosed[i][j]-1]){
							fiuIDmap[i][j+1] = fiuIDmap_unclosed[i][j];
							roi_size[fiuIDmap_unclosed[i][j+1]-1] = roi_size[fiuIDmap_unclosed[i][j+1]-1]-1;
						}
					}
				}

			}
		}
		return fiuIDmap;
	}


	public int [][] ConnectedComponents(int[][] matrix) {
		int[][] labels = new int[matrix.length][matrix[0].length];
		ArrayList<Integer[]> newEdgePixList = new ArrayList<Integer[]>();
		Integer[] curPix = new Integer[2];
		int deltx, delty;
		int neighbX, neighbY;
		int curOldID;
		int curNewID = 0;

		for(int i=0; i<matrix.length; i++){
			for(int j=0; j<matrix[0].length; j++){
				if(labels[i][j]!=0 || matrix[i][j]==0) continue;

				// Once we find an unconsidered pixel
				newEdgePixList.clear();
				curNewID ++;
				curOldID = matrix[i][j];
				newEdgePixList.add(new Integer[]{i,j});

				// Search neighbors layer by layer, until reach the boundary
				while(newEdgePixList.size()>0){
					curPix[0] = newEdgePixList.get(0)[0];
					curPix[1] = newEdgePixList.get(0)[1];
					for(deltx=-1; deltx<=1; deltx++){
						for(delty=-1; delty<=1; delty++){
							if(deltx==0 && delty==0) continue;
							neighbX = curPix[0]+deltx;
							neighbY = curPix[1]+delty;
							if(neighbX<0 || neighbX>=matrix.length || neighbY<0 || neighbY>=matrix[0].length) continue;
							if(matrix[neighbX][neighbY]==curOldID && labels[neighbX][neighbY]==0){
								labels[neighbX][neighbY] = curNewID;
								newEdgePixList.add(new Integer[]{neighbX,neighbY});
							}
						}
					}
					newEdgePixList.remove(0);
				}
			}
		}
		return labels;
	}


	public <T> ArrayList<T> union(ArrayList<T> list1, ArrayList<T> list2) {
		Set<T> set = new HashSet<T>();

		set.addAll(list1);
		set.addAll(list2);

		return new ArrayList<T>(set);
	}


	/**calculate the deltaF/F values**/
	private double[] getMeanF(double[] charXmat_V){
		double[] inFIU = new double [charXmat_V.length];
		for(int i=0;i<inFIU.length;i++)
			inFIU[i] = charXmat_V[i];
		Arrays.sort(inFIU);
		int F0_len = (int) Math.round(0.1*inFIU.length);
		double Fmin = inFIU[0];
		double F0_val = inFIU[F0_len];

		double[] out_values = new double[inFIU.length];
		for(int i=0;i<inFIU.length;i++)
			out_values[i] = (charXmat_V[i]-Fmin)/F0_val;
		return out_values;
	}



	public static void main(String[] args) {
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = FASP_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);

		//		// start ImageJ
		//		new ImageJ();
	}

}
