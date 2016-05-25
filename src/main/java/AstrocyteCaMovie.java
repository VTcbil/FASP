
import java.io.*;
import java.util.ArrayList;
import ij.IJ;


public class AstrocyteCaMovie {
	// Original imaging data
	protected short[][][] imageStack;

	// Features learned by FASP
	protected double[][] neighbCorrMap;
	protected double[][] zNeighbCorrMap;
	protected int[][] fiuIDmap;
	protected double[][] charXmat;
	protected int nFIU;
	protected int nActReg;

	// Parameters:
	protected double threshActPixZ;
	protected double threshActRegZ;
	protected double threshFitRegZ;
	protected int minSize;

	AstrocyteCaMovie(){}

	AstrocyteCaMovie(short[][][] imageStack1, double threshActRegZ1, double threshFitRegZ1, int minSize1) throws IndexOutOfBoundsException, IOException{
		imageStack = imageStack1;
		threshActRegZ = threshActRegZ1;
		threshActPixZ = (double)Math.round(threshActRegZ * 1e3 * (1.0/2.0)) / 1e3;
		threshFitRegZ = threshFitRegZ1;
		minSize = minSize1;

		this.fasp();
	}


	void fasp() throws IndexOutOfBoundsException, IOException{
		// Get the active regions
		ActiveMap actMap = new ActiveMap(imageStack,threshActPixZ,threshActRegZ,minSize,1);
		neighbCorrMap = actMap.corrMap;
		zNeighbCorrMap = actMap.actZmap;
		nActReg = actMap.numConnDomain;


		// Get FIUs from each of the active regions
		boolean[][] actRegPosMap;
		ActiveRegion curActRegion = new ActiveRegion(imageStack);
		int[][] subFIUidMap;
		ArrayList<Double[]> charXlist = new ArrayList<Double[]>();


		fiuIDmap = new int[imageStack.length][imageStack[0].length];

		for(int k=1; k<=nActReg; k++){
			actRegPosMap = actMap.getSingleConnDmMap(k);

			curActRegion.updateActReg(imageStack, zNeighbCorrMap, actRegPosMap, threshFitRegZ, minSize);
			// Record the newly found FIU
			if(curActRegion.nRegFIU > 0){
				subFIUidMap = curActRegion.regFIUidMap;
				for(int m=1; m<=curActRegion.nRegFIU; m++){
					nFIU++;
					for(int x=0; x<fiuIDmap.length; x++){
						for(int y=0; y<fiuIDmap[0].length; y++){
							if(subFIUidMap[x][y]==m){
								fiuIDmap[x][y] = nFIU;
							}
						}
					}
				}
				for(int m=0; m<curActRegion.regCharXlist.size(); m++){
					charXlist.add(curActRegion.regCharXlist.get(m));
				}
			}

			IJ.showStatus(20+(int)(k*78/nActReg)+"% Finished!") ; 
			IJ.showProgress(20+(int)(k*78/nActReg), 100);
		}
		
		if(charXlist.size()>0){
			charXmat = new double[charXlist.size()][charXlist.get(0).length];
			for(int i=0; i<charXlist.size(); i++){
				for(int m=0; m<charXmat[0].length; m++){
					charXmat[i][m] = charXlist.get(i)[m];
				}
			}
		}else{
			charXmat = null;
		}
		
	}

	double[][] getNeighbCorrMap(){
		return neighbCorrMap;
	}

	double[][] getZneighbCorrMap(){
		return zNeighbCorrMap;
	}

	int[][] getFIUidMap(){
		return fiuIDmap;
	}

	double[][] getCharXmat(){
		return charXmat;
	}

	int getNumFIUs(){
		return nFIU;
	}

	int getNumActRegions(){
		return nActReg;
	}

	void setOriMovie(short[][][] imageStack1){
		imageStack = imageStack1;
	}

	void setPara(){
	}

	void showPara(){
		// print parameters on the screen
	}
	
}
