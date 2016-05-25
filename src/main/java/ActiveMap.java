import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

class ActiveMap {
	short[][][] imageStack;
	double[][] corrMap;
	double[][] actZmap;
	double threshActPixZ;
	double threshActRegZ;
	int minSize;
	int corrMethod;  // method of calculating correlation map
	// 1: average of 8 neighbors;
	// 2: maximum of 4 directions
	int numConnDomain = 0;
	int[][] connDmIDmap;



	ActiveMap(short[][][] imageStack1,double threshPixZ1,double threshRegZ1,int minSize1,int method) throws IndexOutOfBoundsException, IOException{
		imageStack = imageStack1;
		threshActPixZ = threshPixZ1;
		threshActRegZ = threshRegZ1;
		minSize = minSize1;

		try{
			corrMap = new double[imageStack.length][imageStack[0].length];
			actZmap = new double[imageStack.length][imageStack[0].length];
			connDmIDmap = new int[imageStack.length][imageStack[0].length];
		} catch(IndexOutOfBoundsException e){
			//			System.out.println("ActiveMap step: imageStack is empty");
			throw e;
		}

		// Calculate neighborhood correlation map and the corresponding z-score map
		corrMethod = method;
		switch(corrMethod){
		case 1:
			setCorrMethod1(); // set corrMap and actZMap
			break;

		case 2:
			setCorrMethod2(); // set corrMap and actZMap
			break;

		default:
			//			System.out.println("Incorrect method ID for calculating correlation map");	
		}

		boolean[][] availableMap = new boolean[actZmap.length][actZmap[0].length];
		for(int i=0; i<actZmap.length; i++){
			for(int j=0; j<actZmap[0].length; j++){
				availableMap[i][j] = true;
			}
		}

		// Binarize the z-score map of neighborhood correlation, find the active regions
		HTregionGrowing neighbZseg = new HTregionGrowing(actZmap,availableMap,threshActPixZ,threshActRegZ,minSize,false);
		connDmIDmap = neighbZseg.connDmIDmap;
		numConnDomain = neighbZseg.nConnDomain;

	}



	private void setCorrMethod1() throws IOException{   // method 1: correlation with the average of 8 neighbors;

		/* *****************     Set corrMap     **************** */
		int NtimePt = imageStack[0][0].length;
		double[] V0 = new double[NtimePt];   // The curve of current pixel (x,y)
		double[] Vnbavg = new double[NtimePt];  // The average curve of 8 neighbors
		PearsonCorr calCorObj = new PearsonCorr();
		Random rNorm = new Random();

		for(int x=1; x<(corrMap.length-1); x++){
			for(int y=1; y<(corrMap[0].length-1); y++){
				for(int k=0; k<NtimePt; k++){
					// The curve of current pixel (x,y)
					V0[k] = (double)imageStack[x][y][k] + rNorm.nextGaussian()*1e-6;
				}
				for(int deltx=-1; deltx<=1; deltx++){
					for(int delty=-1; delty<=1; delty++){
						if(deltx==0 && delty==0) continue;
						for(int k=0; k<NtimePt; k++){
							Vnbavg[k] += ((double)imageStack[x+deltx][y+delty][k] + rNorm.nextGaussian()*1e-6);
						}
					}
				}
				for(int k=0; k<NtimePt; k++){
					// The curve of current pixel (x,y)
					Vnbavg[k] /= 8;
				}

				// Calculate correlation between V0 and its neighbors' average curve
				corrMap[x][y] = calCorObj.calculateCorr(V0, Vnbavg);
			}
		}


		/* ************   Set actZmap (z-score map of neighborhood correlation)   *********** */
		for(int x=1; x<(actZmap.length-1); x++){
			for(int y=1; y<(actZmap[0].length-1); y++){
				// Normalized Fisher transformation
				actZmap[x][y] = (Math.sqrt((double)NtimePt-3.0)/2.0) * Math.log((1.0+corrMap[x][y])/(1.0-corrMap[x][y]));
			}
		}


		// Correction
		double perceZ;
		double perc = 0.1;
		GaussianLUT normLUTObj = new GaussianLUT();
		double targetPerceZ = normLUTObj.myQnorm(perc);
		ArrayList<Double> actZlist = new ArrayList<Double>();
		for(int x=1; x<(actZmap.length-1); x++){
			for(int y=1; y<(actZmap[0].length-1); y++){
				actZlist.add(actZmap[x][y]);
			}
		}
		Collections.sort(actZlist);

		perceZ = actZlist.get((int)Math.round(actZlist.size()*perc));

		for(int x=1; x<(actZmap.length-1); x++){
			for(int y=1; y<(actZmap[0].length-1); y++){
				actZmap[x][y] -= (perceZ - targetPerceZ);
			}
		}

	}




	private void setCorrMethod2(){  // method 2: correlation with the maximum of 4 directions

		/* *****************     Set corrMap     **************** */
		int NtimePt = imageStack[0][0].length;
		double[] V0 = new double[NtimePt];
		double[] V1 = new double[NtimePt];
		double[] V2 = new double[NtimePt];
		double[] V3 = new double[NtimePt];
		double[] V4 = new double[NtimePt];
		double[] curr_Dirs = new double[4];
		PearsonCorr calCorrObj = new PearsonCorr();
		Random rNorm = new Random();

		for(int x=1; x<(corrMap.length-1); x++){
			for(int y=1; y<(corrMap[0].length-1); y++){

				for(int k=0; k<NtimePt; k++){
					// The curve of current pixel (x,y)
					V0[k] = (double)imageStack[x][y][k] + rNorm.nextGaussian()*1e-6;

					// The average curve of direction 1
					V1[k] = (double)(imageStack[x-1][y][k] + imageStack[x+1][y][k]) / 2  + rNorm.nextGaussian()*1e-6;

					// The average curve of direction 2
					V2[k] = (double)(imageStack[x][y-1][k] + imageStack[x][y+1][k]) / 2  + rNorm.nextGaussian()*1e-6;

					// The average curve of direction 3
					V3[k] = (double)(imageStack[x-1][y-1][k] + imageStack[x+1][y+1][k]) / 2 + rNorm.nextGaussian()*1e-6;

					// The average curve of direction 4
					V4[k] = (double)(imageStack[x-1][y+1][k] + imageStack[x+1][y-1][k]) / 2 + rNorm.nextGaussian()*1e-6;
				}

				// Calculate correlation between V0 and the average curve in each direction
				curr_Dirs[0] = calCorrObj.calculateCorr(V0, V1);
				curr_Dirs[1] = calCorrObj.calculateCorr(V0, V2);
				curr_Dirs[2] = calCorrObj.calculateCorr(V0, V3);
				curr_Dirs[3] = calCorrObj.calculateCorr(V0, V4);

				// Find the largest correlation
				corrMap[x][y] = -2;
				for(int di=0; di<4; di++){
					if(curr_Dirs[di] > corrMap[x][y]){
						corrMap[x][y] = curr_Dirs[di];
					}
				}
			}
		}


		/* ************   Set actZmap (z-score map of neighborhood correlation)   *********** */
		for(int x=1; x<(actZmap.length-1); x++){
			for(int y=1; y<(actZmap[0].length-1); y++){
				// Normalized Fisher transformation
				actZmap[x][y] = (Math.sqrt((double)NtimePt-3.0)/2.0) * Math.log((1.0+corrMap[x][y])/(1.0-corrMap[x][y]));

				// Correction for order statistics (Approximation)
				actZmap[x][y] /= Math.sqrt(0.516);
			}
		}


		// Correction
		double medZ;
		double lower,upper;
		ArrayList<Double> actZlist = new ArrayList<Double>();
		for(int x=1; x<(actZmap.length-1); x++){
			for(int y=1; y<(actZmap[0].length-1); y++){
				actZlist.add(actZmap[x][y]);
			}
		}
		Collections.sort(actZlist);
		if (actZlist.size() % 2 == 1){
			medZ = actZlist.get((actZlist.size()+1)/2-1);
		}else{
			lower = actZlist.get(actZlist.size()/2-1);
			upper = actZlist.get(actZlist.size()/2);
			medZ = (lower + upper) / 2.0;
		}
		for(int x=1; x<(actZmap.length-1); x++){
			for(int y=1; y<(actZmap[0].length-1); y++){
				actZmap[x][y] -= medZ;
			}
		}

	}



	boolean[][] getSingleConnDmMap(int domainID){   // Get the position map of one active region 

		boolean[][] singleConnDmMap = new boolean[connDmIDmap.length][connDmIDmap[0].length];
		for(int i=0; i<connDmIDmap.length; i++){
			for(int j=0; j<connDmIDmap[0].length; j++){
				if(connDmIDmap[i][j]==domainID){
					singleConnDmMap[i][j] = true;	
				}else{
					singleConnDmMap[i][j] = false;
				}
			}
		}
		return singleConnDmMap;
	}

}
