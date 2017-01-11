
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

class ActiveRegion {   // A SINGLE active region
	// Original data and static features of an object
	short[][][] imageStack;
	boolean[][] actRegPosMap;
	double[][] corrZmapAll;
	double threshFitPixZ;
	double threshFitRegZ;
	int minSize;

	// Features (to be) learned by FASP
	int[][] regFIUidMap;      // FIUs found in the currently given active region
	ArrayList<Double[]> regCharXlist;
	int nRegFIU;

	private class ModelXm{
		boolean flagValidFIU;
		ArrayList<Integer[]> regPixList;
		Double[] charX;
		double[] betaBank;
		double[] shiftBank;
		double[] sigma2Bank;
		double[][] fitZmap;
		int[] seed;
	}
	private class Update_Info{
		int[] updt_order;
		int[][] valid_pred;
		int[][] map_tbl;
	}
	private class RegPixLstModel{
		double[] charX;
		int[] seed;
		double[] betaBank;
		double[] shiftBank;
		double[] sigma2Bank;
		double[] r;
		ArrayList<Integer[]> rg;
		boolean isConverged;
	}
	private class PixModel{
		int shift;
		double coeff;
		double sigma2;
		double r;
	}

	
	ActiveRegion(short[][][] imageStack1){
		imageStack = imageStack1;
		regFIUidMap = new int[imageStack.length][imageStack[0].length];
		regCharXlist = new ArrayList<Double[]>();
	};

	
	void updateActReg(short[][][] imageStack1, double[][] corrZmapAll1, boolean[][] actRegPosMap1,
			double threshFitRegZ1, int minSize1){
		imageStack = imageStack1;
		corrZmapAll = corrZmapAll1;
		actRegPosMap = actRegPosMap1;     // It's the position map of ONLY ONE active region
		threshFitRegZ = threshFitRegZ1;
		minSize = minSize1;
		regCharXlist.clear();
		nRegFIU = 0;

		for(int i=0; i<regFIUidMap.length; i++){
			for(int j=0; j<regFIUidMap[0].length; j++){
				regFIUidMap[i][j] = -100;
			}
		}

		// Private variables in this function
		ModelXm modelOfXm = new ModelXm();
		boolean flagEmptyAvail = false;
		boolean flagOneFIUfound = false;
		int availCount;
		int seedAvailCount;
		boolean[][] pixAvailMap = new boolean[actRegPosMap.length][actRegPosMap[0].length];  // The pixels not labeled as any FIU's
		boolean[][] seedAvailMap = new boolean[actRegPosMap.length][actRegPosMap[0].length];
		//		The pixels that aren't considered as seed/bad seed yet in finding a certain FIU
		int[] badSeed;
		ArrayList<Integer[]> badSeedList = new ArrayList<Integer[]>();
		ArrayList<Integer[]> badSeedHoleList = new ArrayList<Integer[]>();
		boolean[][] tinyRegCheckedMap = new boolean[actRegPosMap.length][actRegPosMap[0].length];  // The pixels already checked for tiny region 

		int timesDigBadSeedHole = 0;
		int x, y;
		int tempPixCount;


		availCount = 0;
		for(x=0; x<pixAvailMap.length; x++){
			for(y=0; y<pixAvailMap[0].length; y++){
				pixAvailMap[x][y] = actRegPosMap[x][y];
				if(pixAvailMap[x][y]) availCount++;
			}
		}
		if(availCount < minSize){
			flagEmptyAvail = true;
		}

		nRegFIU = 0;

		// Iteration layer2: As long as more than minSize pixels in this active region haven't been labeled as any FIU's
		while(!flagEmptyAvail){
			// Or, as long as some FIU are still not found, we continue to work on the remaining region to find another FIU

			// (Different from R code, here we examine both pixaVailMap and seedAvailMap outside findOneFIU(), namely, findUniqUnit())

			// Reset seedAvailMap as all pixels in the remaining part of this active region
			seedAvailCount = 0;
			for(x=0; x<seedAvailMap.length; x++){
				for(y=0; y<seedAvailMap[0].length; y++){
					seedAvailMap[x][y] = pixAvailMap[x][y];
					if(seedAvailMap[x][y]) seedAvailCount ++;
				}
			}
			flagOneFIUfound = false;

			// Iteration layter3: In finding one FIU, iterate across possible seed pixels, until one FIU is found
			while((!flagOneFIUfound) && (seedAvailCount>5)){
				// TRY ONCE to find the next FIU, given current seedAvailMap
				// Possible to return empty result (because of bad seed)
				boolean zeroornot = false;
				for(int i=0; i<seedAvailMap.length; i++){
					if(zeroornot)
						break;
					for(int j=0; j<seedAvailMap[0].length; j++){
						if(zeroornot)
							break;
						if(seedAvailMap[i][j]){
							zeroornot = true;
						}
					}
				}
				if(!zeroornot)
					zeroornot = false;
				modelOfXm = findOneFIU(pixAvailMap, availCount, seedAvailMap);

				// If modelOfXm is empty, renew the seedAvailMap and TRY AGAIN IN THE NEXT ITERATION
				// Renew seedAvailMap: remove pixels around the bad seed from the seed available map
				// Renew both pixAvailMap and seedAvailMap: remove isolated domains of area less than minSize
				if(!modelOfXm.flagValidFIU){

					// System.out.println("Digging seed available map...\n");
					badSeed = modelOfXm.seed;
					seedAvailMap[badSeed[0]][badSeed[1]] = false;
					seedAvailCount --;
					badSeedList.clear();
					badSeedList.add(new Integer[]{badSeed[0],badSeed[1]});
					timesDigBadSeedHole = (int)Math.ceil(Math.max(Math.sqrt((double)availCount)/6,Math.sqrt((double)minSize)/2));
					if(timesDigBadSeedHole<1) timesDigBadSeedHole = 1;

					for(int i=0; i<timesDigBadSeedHole; i++){
						badSeedHoleList.clear();
						findNeighbors(badSeedHoleList, badSeedList, seedAvailMap);
						if(badSeedHoleList.size()==0) break;
						for(int j=0; j<badSeedHoleList.size(); j++){
							badSeedList.add(badSeedHoleList.get(j));
							seedAvailMap[badSeedHoleList.get(j)[0]][badSeedHoleList.get(j)[1]] = false;
							seedAvailCount --;
						}
					}

				}else{
					flagOneFIUfound = true;     // As long as ANY ONE FIU has been found, break this while loop
				}
			}


			if(flagOneFIUfound){
				//// Record the currently-got FIU
				nRegFIU ++;
				for(int ii=0; ii<modelOfXm.regPixList.size(); ii++){
					x = modelOfXm.regPixList.get(ii)[0];
					y = modelOfXm.regPixList.get(ii)[1];
					regFIUidMap[x][y] = nRegFIU;
				}
				regCharXlist.add(modelOfXm.charX);

				//// Prepare for searching for the next FIU
				// Eliminate the current FIU from available map
				for(int ii=0; ii<modelOfXm.regPixList.size(); ii++){
					x = modelOfXm.regPixList.get(ii)[0];
					y = modelOfXm.regPixList.get(ii)[1];
					pixAvailMap[x][y] = false;
					availCount --;
					seedAvailMap[x][y] = false;
					seedAvailCount --;
				}

				// Remove tiny domains in the remaining available maps
				// System.out.println("Removing tiny domains in available maps...\n");
				for(x=0; x<pixAvailMap.length; x++){
					for(y=0; y<pixAvailMap[0].length; y++){
						tinyRegCheckedMap[x][y] = !seedAvailMap[x][y];
					}
				}
				for(x=0; x<pixAvailMap.length; x++){
					for(y=0; y<pixAvailMap[0].length; y++){
						if(!tinyRegCheckedMap[x][y]){	  // once a unconsidered available pixel is found, search from it
							badSeedList.clear();
							badSeedHoleList.clear();
							badSeedList.add(new Integer[]{x,y});
							tinyRegCheckedMap[x][y] = true;
							tempPixCount = 1;
							findNeighbors(badSeedHoleList, badSeedList, pixAvailMap);
							while(!badSeedHoleList.isEmpty()){   // search until it cannot be enlarged any more
								for(int j=0; j<badSeedHoleList.size(); j++){
									badSeedList.add(badSeedHoleList.get(j));
									tinyRegCheckedMap[badSeedHoleList.get(j)[0]][badSeedHoleList.get(j)[1]] = true;
									tempPixCount ++;
								}
								findNeighbors(badSeedHoleList, badSeedList, pixAvailMap);
							}
							if(tempPixCount < minSize){
								for(int j=0; j<badSeedList.size(); j++){
									pixAvailMap[badSeedList.get(j)[0]][badSeedList.get(j)[1]] = false;
									availCount --;
									seedAvailMap[badSeedList.get(j)[0]][badSeedList.get(j)[1]] = false;
									seedAvailCount --;
								}
							}
						}
					}
				}

			}else if(seedAvailCount<=5){
				flagEmptyAvail = true;
			}

			if(availCount < minSize){
				flagEmptyAvail = true;
			}

		}
	}



	// TRY ONCE to find the next FIU, given current seedAvailMap
	// Possible to return empty result (because of bad seed)
	// (R code findUniqUnit() without update seedAvailMap inside)
	private ModelXm findOneFIU(boolean[][] pixAvailMap, int pixAvailCount, boolean[][] seedAvailMap){  
		// Iteratively learn an Xm and the corresponding parameters,
		// using remaining part of the domain
		ModelXm curModel = new ModelXm();
		int numTimePoints = imageStack[0][0].length;

		double[] charX;
		double[] charX_old = new double[numTimePoints];
		int[] seed_em1 = new int[2];
		int[] tempSeed = new int[2];
		boolean[][] newFIUpos = null; 
		boolean[][] posDiffZpixPos = null;

		//	Generate dataMatCore (curve matrix of the pixels in the given region)
		//  Adding in normal noises to avoid error in calculating correlation
		double[][] dataMatCore = new double[pixAvailCount][numTimePoints];
		ArrayList<Integer[]> curRegPixList = new ArrayList<Integer[]>();
		double[] tempCurve = new double[numTimePoints];
		int pixID = -1;
		double meanOfCurve;
		double maxAbsOfCurve;
		Random rNorm = new Random();
		for(int i=0; i<pixAvailMap.length; i++){
			for(int j=0; j<pixAvailMap[0].length; j++){
				if(pixAvailMap[i][j]){
					curRegPixList.add(new Integer[]{i, j});
					pixID ++;
					meanOfCurve = 0;
					for(int k=0; k<numTimePoints; k++){
						meanOfCurve += imageStack[i][j][k];
					}
					meanOfCurve /= numTimePoints;

					maxAbsOfCurve = 0;
					for(int k=0; k<numTimePoints; k++){
						tempCurve[k] = (imageStack[i][j][k] - meanOfCurve) / 255;
						if(Math.abs(tempCurve[k])>maxAbsOfCurve) maxAbsOfCurve = Math.abs(tempCurve[k]);
					}

					for(int k=0; k<numTimePoints; k++){
						dataMatCore[pixID][k] = tempCurve[k] + maxAbsOfCurve*rNorm.nextGaussian()*0.001;
					}
				}
			}
		}


		//	Generate dataMatNeib to compute neighbor based correlation map
		//  Find one layer of neighbors of current ROI 
		ArrayList<Integer[]> curNeighborList = new ArrayList<Integer[]>();
		boolean[][] neighbAvailMap = new boolean[pixAvailMap.length][pixAvailMap[0].length];
		for(int i=0; i<pixAvailMap.length; i++){
			for(int j=0; j<pixAvailMap[0].length; j++){
				if(pixAvailMap[i][j]){
					neighbAvailMap[i][j] = false;
				}else{
					neighbAvailMap[i][j] = true;
				}
			}
		}
		findNeighbors(curNeighborList, curRegPixList, neighbAvailMap);
		int numNeighb = curNeighborList.size();
		double[][] dataMatNeib = new double[numNeighb][numTimePoints];
		int neighbX, neighbY;
		for(pixID=0; pixID<numNeighb; pixID++){
			neighbX = curNeighborList.get(pixID)[0];
			neighbY = curNeighborList.get(pixID)[1];

			meanOfCurve = 0;
			for(int k=0; k<numTimePoints; k++){
				meanOfCurve += imageStack[neighbX][neighbY][k];
			}
			meanOfCurve /= numTimePoints;

			maxAbsOfCurve = 0;
			for(int k=0; k<numTimePoints; k++){
				tempCurve[k] = (imageStack[neighbX][neighbY][k] - meanOfCurve) / 255;
				if(Math.abs(tempCurve[k])>maxAbsOfCurve) maxAbsOfCurve = Math.abs(tempCurve[k]);
			}

			for(int k=0; k<numTimePoints; k++){
				dataMatNeib[pixID][k] = tempCurve[k] + maxAbsOfCurve*rNorm.nextGaussian()*0.001;
			}
		}

		//  locate the IDs of pixels in dataMatCore in a spatial map (in order to construct dataMatFIU from dataMatCore)
		int[][] idPixInCoreMatMap = new int[imageStack.length][imageStack[0].length];
		for(int i=0; i<idPixInCoreMatMap.length; i++){
			for(int j=0; j<idPixInCoreMatMap[0].length; j++){
				idPixInCoreMatMap[i][j] = -100;
			}
		}
		for(int k=0; k<curRegPixList.size(); k++){
			idPixInCoreMatMap[curRegPixList.get(k)[0]][curRegPixList.get(k)[1]] = k;
		}


		//	EM1 for core pixels
		RegPixLstModel tempModel = em1(dataMatCore, curRegPixList, seedAvailMap);
		charX = tempModel.charX;
		seed_em1[0] = tempModel.seed[0];
		seed_em1[1] = tempModel.seed[1];


		//	EM for neighbor
		tempSeed[0] = curNeighborList.get(0)[0];
		tempSeed[1] = curNeighborList.get(0)[1];
		Update_Info updt_info1 = get_updt_info_simple(curNeighborList,tempSeed);
		RegPixLstModel neibAlpha = calAlpha_local_simple(dataMatNeib, charX, updt_info1);

		//  Correlation of residual map
		double[][] resiMatCore = calResidu(dataMatCore, charX, tempModel.betaBank, tempModel.shiftBank, tempModel.sigma2Bank);
		double[][] resiMatNeib = calResidu(dataMatNeib, charX, neibAlpha.betaBank, neibAlpha.shiftBank, neibAlpha.sigma2Bank);

		//  Calculate z-score map of fitness
		//  methods: charX_residual
		double[][] diffZmap = get_z_map_charX_residual(charX, curRegPixList, curNeighborList,resiMatCore,
				resiMatNeib, tempModel.shiftBank);



		/////////		  Re-fit the model in the non-zero-diffZ regions		/////////
		//		Find regions of positive z-score
		ArrayList<Integer[]> fiuPixList = new ArrayList<Integer[]>();
		double[][] dataMatFIU;
		posDiffZpixPos = new boolean[imageStack.length][imageStack[0].length];
		for(int i=0; i<posDiffZpixPos.length; i++){
			for(int j=0; j<posDiffZpixPos[0].length; j++){
				if(diffZmap[i][j]>=0 & idPixInCoreMatMap[i][j]>=0){
					posDiffZpixPos[i][j] = true;
					fiuPixList.add(new Integer[]{i,j});
				}
			}
		}

		double deltaCharX = 0;
		for(int k=0; k<numTimePoints; k++)
			deltaCharX += (charX_old[k]-charX[k])*(charX_old[k]-charX[k]);

		int nIterRefit = 0;

		while(fiuPixList.size()>10 & deltaCharX>1e-1 & nIterRefit<10){

			nIterRefit++;
			
			for(int k=0; k<numTimePoints; k++)  charX_old[k] = charX[k];

			// EM2 for positive-z area
			dataMatFIU = new double[fiuPixList.size()][];
			for(int jj=0; jj<fiuPixList.size(); jj++){
				neighbX = fiuPixList.get(jj)[0];
				neighbY = fiuPixList.get(jj)[1];
				dataMatFIU[jj] = dataMatCore[idPixInCoreMatMap[neighbX][neighbY]];
			}
			tempModel = em2(charX, dataMatFIU, fiuPixList, posDiffZpixPos);


			// Re-calculate diffZmap based on the newly-fitted curve
			charX = tempModel.charX;

			//		EM for neighbor
			neibAlpha = calAlpha_local_simple(dataMatNeib, charX, updt_info1);

			//  Correlation of residual map
			Update_Info updt_info2 = get_updt_info_simple(curRegPixList,seed_em1);
			RegPixLstModel coreAlpha = calAlpha_local_simple(dataMatCore, charX, updt_info2);

			resiMatCore = calResidu(dataMatCore, charX, coreAlpha.betaBank, coreAlpha.shiftBank, coreAlpha.sigma2Bank);
			resiMatNeib = calResidu(dataMatNeib, charX, neibAlpha.betaBank, neibAlpha.shiftBank, neibAlpha.sigma2Bank);

			//  Calculate z-score map of fitness
			//  methods: charX_residual
			diffZmap = get_z_map_charX_residual(charX, curRegPixList, curNeighborList,resiMatCore,
					resiMatNeib, coreAlpha.shiftBank);


			//	Find regions of positive z-score
			fiuPixList.clear();
			posDiffZpixPos = new boolean[imageStack.length][imageStack[0].length];
			for(int i=0; i<posDiffZpixPos.length; i++){
				for(int j=0; j<posDiffZpixPos[0].length; j++){
					if(diffZmap[i][j]>=0 & idPixInCoreMatMap[i][j]>=0){
						posDiffZpixPos[i][j] = true;
						fiuPixList.add(new Integer[]{i,j});
					}else{
						posDiffZpixPos[i][j] = false;
					}
				}
			}

			deltaCharX = 0;
			for(int k=0; k<numTimePoints; k++)
				deltaCharX += (charX_old[k]-charX[k])*(charX_old[k]-charX[k]);

		}

		
		//  Binarize the z-score map of fitness, find one FIU
		threshFitPixZ = 0; 
		HTregionGrowing regFitZseg = new HTregionGrowing(diffZmap,pixAvailMap,threshFitPixZ,threshFitRegZ,minSize,true);
		// pixAvailMap will not be changed in this function



		/////////		  Re-fit the model in the found FIU		/////////
		fiuPixList.clear();

		//		EM2 for FIU pixels (Only pick out the largest meaningful region)
		if(regFitZseg.nConnDomain!=0){
			newFIUpos = regFitZseg.getMaxSizedConnDmMap();
			for(int i=0; i<newFIUpos.length; i++){
				for(int j=0; j<newFIUpos[0].length; j++){
					if(newFIUpos[i][j])  fiuPixList.add(new Integer[]{i,j});
				}
			}
			dataMatFIU = new double[fiuPixList.size()][];
			for(int jj=0; jj<fiuPixList.size(); jj++){
				neighbX = fiuPixList.get(jj)[0];
				neighbY = fiuPixList.get(jj)[1];
				dataMatFIU[jj] = dataMatCore[idPixInCoreMatMap[neighbX][neighbY]];
			}
			tempModel = em2(charX, dataMatFIU, fiuPixList, newFIUpos);

			curModel.flagValidFIU = true;
			curModel.charX = new Double[numTimePoints];
			for(int k=0; k<numTimePoints; k++) curModel.charX[k] = tempModel.charX[k];
			curModel.regPixList = fiuPixList;
			curModel.betaBank = tempModel.betaBank;
			curModel.shiftBank = tempModel.shiftBank;
			curModel.sigma2Bank = tempModel.sigma2Bank;
			curModel.fitZmap = diffZmap;
			curModel.seed = seed_em1;


		}else{  // not successfully found any FIU
			curModel.flagValidFIU = false;
			curModel.seed = seed_em1;
		}
		return curModel;
	}



	private void findNeighbors(ArrayList<Integer[]> curNeighborList, ArrayList<Integer[]> curRegPixList, boolean[][] neighbAvailMap){
		// Update curNeighborList and neighbAvailMap (by using the pointers)
		Integer[] curPix;
		int neighbX, neighbY;
		boolean[][] consideredNeighbMap = new boolean[neighbAvailMap.length][neighbAvailMap[0].length];
		for(int i=0; i<curRegPixList.size(); i++){
			consideredNeighbMap[curRegPixList.get(i)[0]][curRegPixList.get(i)[1]] = true;
		}

		curNeighborList.clear();

		for(int i=0; i<curRegPixList.size(); i++){
			curPix = curRegPixList.get(i);

			for(int deltx=-1; deltx<=1; deltx++){
				for(int delty=-1; delty<=1; delty++){
					if(deltx==0 && delty==0) continue;
					neighbX = curPix[0]+deltx;
					neighbY = curPix[1]+delty;
					if(neighbX<0 || neighbX>=neighbAvailMap.length || neighbY<0 || neighbY>=neighbAvailMap[0].length) continue;
					if(neighbAvailMap[neighbX][neighbY] && !consideredNeighbMap[neighbX][neighbY]){
						curNeighborList.add(new Integer[]{neighbX, neighbY});
						consideredNeighbMap[neighbX][neighbY] = true;
					}
				}
			}
		}
	}



	private RegPixLstModel em1(double[][] dataMatCore, ArrayList<Integer[]> curRegPixList, boolean[][] seedAvailMap){
		int max_iter = 50;
		RegPixLstModel coreModel = new RegPixLstModel();
		double[][] coreCorrMap = new double[seedAvailMap.length][seedAvailMap[0].length];
		for(int i=0;i<coreCorrMap.length;i++)
			for(int j=0;j<coreCorrMap[0].length;j++)
				coreCorrMap[i][j] = -100;
		double maxCorr = 0;
		int[] seedPixel = new int[2];
		int numTimePoints = imageStack[0][0].length;
		double[] charX;
		double[] oldCharX = new double[numTimePoints];
		RegPixLstModel coreAlpha;
		myBasicMath basicMathObj = new myBasicMath();


		//// Initiation
		for(int i=0; i<seedAvailMap.length; i++){
			for(int j=0; j<seedAvailMap[0].length; j++){
				if(seedAvailMap[i][j]){
					coreCorrMap[i][j] = corrZmapAll[i][j];
					if(coreCorrMap[i][j] > maxCorr){
						maxCorr = coreCorrMap[i][j];
						seedPixel[0] = i;
						seedPixel[1] = j;	
					}
				}else{
					coreCorrMap[i][j] = 0;
				}
			}
		}
		if(seedPixel[0]==0 & seedPixel[1]==0)
			seedPixel[1] = 0;
		charX = getInitCharX(seedPixel);
		for(int k=0; k<numTimePoints; k++) oldCharX[k] = charX[k];

		// Stop criteria
		double delta = 1;
		int loop = 0;
		coreModel.isConverged = true;

		// Initiate updating order
		Update_Info updt_info = get_updt_info_simple(curRegPixList, seedPixel);		

		// Iterative process
		double delta_thr = 1e-3;
		loop = 0;
		while(!Double.isNaN(delta) && delta>delta_thr && loop<max_iter){
			loop ++;

			// 'Adaptive' threshold
			if(loop>20){
				delta_thr = 1e-2;
			}else if(loop>40){
				delta_thr = 1e-1;
			}else{
				delta_thr = 1e-3;
			}

			coreAlpha = calAlpha_local_simple(dataMatCore, charX, updt_info);

			charX = calCharX(dataMatCore,coreAlpha.betaBank,coreAlpha.shiftBank,coreAlpha.sigma2Bank);
			for(int k=0; k<numTimePoints; k++) oldCharX[k] = charX[k] - oldCharX[k];
			delta = Math.sqrt(basicMathObj.sampleVar(oldCharX)) / Math.sqrt(basicMathObj.sampleVar(charX));
			for(int k=0; k<numTimePoints; k++) oldCharX[k] = charX[k];
		}
		if(loop==max_iter) {
			//			System.out.println("charX not converged in EM1\n");
			coreModel.isConverged = false;
		}

		coreAlpha = calAlpha_local_simple(dataMatCore, charX, updt_info);

		coreModel.charX = charX;
		coreModel.seed = seedPixel;
		coreModel.betaBank = coreAlpha.betaBank;
		coreModel.shiftBank = coreAlpha.shiftBank;
		coreModel.sigma2Bank = coreAlpha.sigma2Bank;
		coreModel.r = coreAlpha.r;
		coreModel.rg = coreAlpha.rg;

		return coreModel;
	}




	RegPixLstModel em2(double[] charX1, double[][] dataMatFIU, ArrayList<Integer[]> fiuPixList, boolean[][] seedAvailMap){
		int max_iter = 50;
		RegPixLstModel fiuModel = new RegPixLstModel();
		double[][] fiuCorrMap = new double[seedAvailMap.length][seedAvailMap[0].length];
		double maxCorr = 0;
		int[] seedPixel = new int[2];
		int numTimePoints = imageStack[0][0].length;
		double[] charX = new double[numTimePoints];
		double[] oldCharX = new double[numTimePoints];
		RegPixLstModel fiuAlpha;
		myBasicMath basicMathObj = new myBasicMath();


		//// Initiation
		// Seed for fitting process
		for(int i=0; i<seedAvailMap.length; i++){
			for(int j=0; j<seedAvailMap[0].length; j++){
				if(seedAvailMap[i][j]){
					fiuCorrMap[i][j] = corrZmapAll[i][j];
					if(fiuCorrMap[i][j] > maxCorr){
						maxCorr = fiuCorrMap[i][j];
						seedPixel[0] = i;
						seedPixel[1] = j;
					}
				}else{
					fiuCorrMap[i][j] = 0;
				}
			}
		}
		// charX
		for(int k=0; k<numTimePoints; k++){
			charX[k] = charX1[k];
			oldCharX[k] = charX[k];
		}

		// Stop criteria
		double delta = 1;
		int loop = 0;

		// Initiate updating order
		Update_Info updt_info = get_updt_info_simple(fiuPixList, seedPixel);		

		// Iterative process
		double delta_thr = 1e-3;
		loop = 0;
		while(!Double.isNaN(delta) && delta>delta_thr && loop<max_iter){
			loop ++;

			// 'Adaptive' threshold
			if(loop>20){
				delta_thr = 1e-2;
			}else if(loop>40){
				delta_thr = 1e-1;
			}else{
				delta_thr = 1e-3;
			}

			fiuAlpha = calAlpha_local_simple(dataMatFIU, charX, updt_info);

			charX = calCharX(dataMatFIU,fiuAlpha.betaBank,fiuAlpha.shiftBank,fiuAlpha.sigma2Bank);
			for(int k=0; k<numTimePoints; k++) oldCharX[k] = charX[k] - oldCharX[k];
			delta = Math.sqrt(basicMathObj.sampleVar(oldCharX)) / Math.sqrt(basicMathObj.sampleVar(charX));
			for(int k=0; k<numTimePoints; k++) oldCharX[k] = charX[k];
		}

		fiuAlpha = calAlpha_local_simple(dataMatFIU, charX, updt_info);

		fiuModel.charX = charX;
		fiuModel.seed = seedPixel;
		fiuModel.betaBank = fiuAlpha.betaBank;
		fiuModel.shiftBank = fiuAlpha.shiftBank;
		fiuModel.sigma2Bank = fiuAlpha.sigma2Bank;
		fiuModel.r = fiuAlpha.r;
		fiuModel.rg = fiuAlpha.rg;

		return fiuModel;
	}



	private double[] getInitCharX(int[] seedPixel){
		//	 For the first pixel in seedPixel, find the best correlation/propagation direction,
		//	 then get the average waveform in this direction (this pixel and two neighbors),
		//	 return the standardized average waveform.

		int x = seedPixel[0];
		int y = seedPixel[1];

		if(x <= 0) x = 1;
		if(x >= (imageStack.length-1)) x = imageStack.length-2;
		if(y <= 0) y = 1;
		if(y >= (imageStack[0].length-1)) y = imageStack[0].length-2;

		int numTimePoints = imageStack[0][0].length;
		double[] initCharX = new double[numTimePoints];
		double[] V0 = new double[numTimePoints];
		double[][] Vdir = new double[4][numTimePoints];
		double[] corrDir = new double[4];
		Random rNorm = new Random();
		PearsonCorr calCorObj = new PearsonCorr();
		double maxCorr;
		int dirID_maxCorr;
		double meanVal;
		double stdVal;
		int k;
		myBasicMath basicMathObj = new myBasicMath();

		for(k=0; k<numTimePoints; k++){
			V0[k] = (double)imageStack[x][y][k] + rNorm.nextGaussian()*0.0001;
		}

		// Direction 0
		for(k=0; k<numTimePoints; k++){
			Vdir[0][k] = ((double)imageStack[x-1][y][k] + (double)imageStack[x+1][y][k]) / 2;
		}
		// Direction 1
		for(k=0; k<numTimePoints; k++){
			Vdir[1][k] = ((double)imageStack[x][y-1][k] + (double)imageStack[x][y+1][k]) / 2;
		}
		// Direction 2
		for(k=0; k<numTimePoints; k++){
			Vdir[2][k] = ((double)imageStack[x-1][y+1][k] + (double)imageStack[x+1][y-1][k]) / 2;
		}
		// Direction 3
		for(k=0; k<numTimePoints; k++){
			Vdir[3][k] = ((double)imageStack[x-1][y-1][k] + (double)imageStack[x+1][y+1][k]) / 2;
		}

		maxCorr = 0;
		dirID_maxCorr = -1;
		// Find the direction with the max correlation with V0
		for(int i=0; i<4; i++){
			corrDir[i] = calCorObj.calculateCorr(V0, Vdir[i]);
			if(corrDir[i]>maxCorr){
				maxCorr = corrDir[i];
				dirID_maxCorr = i;
			}
		}

		for(k=0; k<numTimePoints; k++){
			initCharX[k] = (Vdir[dirID_maxCorr][k]*2 + V0[k]) / 3;
		}
		meanVal = basicMathObj.sampleMean(initCharX);
		stdVal = Math.sqrt(basicMathObj.sampleVar(initCharX, meanVal));

		for(k=0; k<numTimePoints; k++){
			initCharX[k] = (initCharX[k] - meanVal)/stdVal;
		}

		return initCharX;
	}



	//	Find order of the coordinate of pixels to update -----------
	//	give the already updated neighbor as well
	//	!! boundary check not done, not a issue since these points are not included
	public Update_Info get_updt_info_simple(ArrayList<Integer[]> curRegPixList, int[] seedPixel){
		int Nx = imageStack.length;
		int Ny = imageStack[0].length;
		int totalPix = Nx*Ny;
		int[] offset = new int[]{-1, 1, Nx, -Nx, Nx+1, -Nx-1, Nx-1, -Nx+1};
		int numRegPix = curRegPixList.size();
		int[][] mat = new int[Nx][Ny];
		// values: 0 - not interesting; 
		// 	       1 - unconsidered pixel in curRegPixList; 
		// 		   2 - pixels added to res[][] but its neighbors NOT checked;
		//	       3 - pixels added to res[][] AND its neighbors checked.
		for(int i=0; i<numRegPix; i++){
			mat[curRegPixList.get(i)[0]][curRegPixList.get(i)[1]] = 1;
		}

		int[][] res = new int[numRegPix][9];     // [:][0] for idx, [:][1:8] for already updated neighbors
		for(int i=0; i<numRegPix; i++){
			for(int j=0; j<9; j++){
				res[i][j] = -1;
			}
		}

		int x;   // 1-dim ID of the pixel   (starting from 0)
		int[] xnb = new int[offset.length];  // 1-dim IDs of x's 8 neighbors  (x1 in R code)
		int[] matValNb =new int[offset.length];   // Values of the 8 neighbors of x in mat  (xn in R code)
		int temp2Dx, temp2Dy;
		int n_updt;
		int n_next;
		ArrayList<Integer> nbr_updt = new ArrayList<Integer>();
		ArrayList<Integer> nbr_next = new ArrayList<Integer>();
		boolean flagBreakBreak = false;
		int pt2 = 1;
		res[0][0] = seedPixel[1]*Nx + seedPixel[0];      // res[0] is the seedPixel

		for(int ii=0; ii<numRegPix; ii++){
			x = res[ii][0];
			// When one region is clear at all, go to next disjoint region
			if(x==-1){
				flagBreakBreak = false;
				for(int i=0; i<mat.length; i++) {
					for(int j=0; j<mat[0].length; j++) {
						if(mat[i][j]==1){
							res[ii][0] = j*Nx + i;
							flagBreakBreak = true;
							break;
						}
					}
					if(flagBreakBreak) break;
				}
				x = res[ii][0];
				pt2 = pt2 + 1;

				temp2Dy = (int)Math.floor(x/Nx);
				temp2Dx = x - Nx*temp2Dy;
				mat[temp2Dx][temp2Dy] = 2;
			}

			nbr_updt.clear();
			nbr_next.clear();
			for(int j=0; j<offset.length; j++) {
				xnb[j] = x + offset[j];              // A neighbor of x
				if(xnb[j]<0 || xnb[j]>=totalPix){
					matValNb[j] = 0;
				}else{
					temp2Dy = (int)Math.floor(xnb[j]/Nx);
					temp2Dx = xnb[j] - Nx*temp2Dy;
					matValNb[j] = mat[temp2Dx][temp2Dy];       // This neighbor's mat value
					if(matValNb[j]==3) nbr_updt.add(xnb[j]);
					if(matValNb[j]==1) nbr_next.add(xnb[j]);
				}
			}
			n_updt = nbr_updt.size();
			n_next = nbr_next.size();
			if(n_updt>0){
				for(int j=1; j<9; j++){
					if(matValNb[j-1]==3){
						res[ii][j] = 1;
					}else{
						res[ii][j] = 0;
					}
				}
			}
			if(n_next>0){
				int k=0;
				for(int j=pt2; j<(pt2+n_next); j++){
					res[j][0] = nbr_next.get(k); 
					k++;
				}
				pt2 = pt2 + n_next;
			}
			for(int j=0; j<n_next; j++){
				temp2Dy = (int)Math.floor(nbr_next.get(j)/Nx);
				temp2Dx = nbr_next.get(j) - Nx*temp2Dy;
				mat[temp2Dx][temp2Dy] = 2;
			}
			temp2Dy = (int)Math.floor(x/Nx);
			temp2Dx = x - Nx*temp2Dy;
			mat[temp2Dx][temp2Dy] = 3;
		}

		//	Find the relationship between coordinate and curRegPixList order
		int[][] mat1 = new int[Nx][Ny];
		for(int i=0; i<Nx; i++){
			for(int j=0; j<Ny; j++){
				mat1[i][j] = -1;   // Background: -1
			}
		}
		for(int i=0; i<numRegPix; i++){
			mat1[curRegPixList.get(i)[0]][curRegPixList.get(i)[1]] = i;   // Note: starting from 0
		}


		Update_Info updt_info = new Update_Info();
		updt_info.updt_order = new int[numRegPix];
		for(int i=0; i<numRegPix; i++){
			updt_info.updt_order[i] = res[i][0];
		}
		updt_info.valid_pred = new int[numRegPix][8];
		for(int i=0; i<numRegPix; i++){
			for(int j=1; j<9; j++){
				updt_info.valid_pred[i][j-1] = res[i][j];
			}
		}
		updt_info.map_tbl = mat1;

		return updt_info;
	}



	//	core.alpha = calAlpha_local_simple(dataMatCore, charX, updt_info);
	private RegPixLstModel calAlpha_local_simple(double[][] dataMat, double[] charX, Update_Info updt_info){
		int Nx = imageStack.length;
		int[] offset = new int[]{-1, 1, Nx, -Nx, Nx+1, -Nx-1, Nx-1, -Nx+1};
		int numRegPix = dataMat.length;

		int[] u0 = updt_info.updt_order;  // order of update, coordinate index (1D-coded)
		int[][] p0 = updt_info.valid_pred;  // direction for prediction, coordinate index
		int[][] m0 = updt_info.map_tbl;  // matrix with element as its lstPixel index (starting from 0; background==-1)
		double[][] res0 = new double[numRegPix][4];    // results, ordered by lstPixel index
		//												  colnames(res0) : c('shift','coeff','sigma2','r')
		int idxi;    // 1D-coded coordinate index of pixel
		int temp2Dx, temp2Dy;
		int[] p0i;
		ArrayList<Integer> rg0 = new ArrayList<Integer>();
		boolean[] idx0 = new boolean[p0[0].length];
		int trueCount = 0;
		int cd0i; 
		int tempInt;
		int rg0Left, rg0Right;
		PixModel res00_pixMd;
		ArrayList<Integer[]> res_rg00 = new ArrayList<Integer[]>();
		Integer[] tempPt;


		for(int ii=0; ii<numRegPix; ii++){
			temp2Dy = (int)Math.floor(u0[ii]/Nx);
			temp2Dx = u0[ii] - Nx*temp2Dy;
			idxi = m0[temp2Dx][temp2Dy];    // index in curPixList or dataMat
			p0i = p0[ii];   // direction to look for prediction
			if(ii>8){
				trueCount = 0;
				for(int i=0; i<p0i.length; i++){
					if(p0i[i]>0){
						idx0[i] = true;
						trueCount ++;
					}else{
						idx0[i] = false;
					}
				}
				if(trueCount==0){
					rg0.clear();
					for(int i=-2; i<3; i++) rg0.add(i);
				}else{
					rg0.clear();
					for(int i=0; i<idx0.length; i++){
						if(idx0[i]){
							cd0i = u0[ii] + offset[i];
							temp2Dy = (int)Math.floor(cd0i/Nx);
							temp2Dx = cd0i - Nx*temp2Dy;
							tempInt = (int)Math.round(res0[m0[temp2Dx][temp2Dy]][1]);
							rg0.add(tempInt);
							rg0.add(tempInt - 1);
							rg0.add(tempInt + 1);
							rg0.add(tempInt - 2);
							rg0.add(tempInt + 2);
						}
					}
					rg0Left =  Collections.min(rg0);
					rg0Right =  Collections.max(rg0);
					if(rg0Left<-30) rg0Left = -30;
					if(rg0Right>30) rg0Right = 30;
					rg0.clear();
					for(int i=rg0Left; i<rg0Right+1; i++) rg0.add(i);
				}
			}else{
				rg0.clear();
				for(int i=-2; i<3; i++) rg0.add(i);
			}

			res00_pixMd = findBestMatch_ls(charX, dataMat[idxi], rg0);
			res0[idxi][0] = (double)res00_pixMd.shift;
			res0[idxi][1] = res00_pixMd.coeff;
			res0[idxi][2] = res00_pixMd.sigma2;
			res0[idxi][3] = res00_pixMd.r;

			tempPt = new Integer[rg0.size()];
			rg0.toArray(tempPt);
			res_rg00.add(tempPt);
		}

		double[] shiftBank = new double[numRegPix];
		for(int i=0; i<numRegPix; i++) shiftBank[i] = res0[i][0];
		double[] betaBank = new double[numRegPix];
		for(int i=0; i<numRegPix; i++) betaBank[i] = res0[i][1];
		double[] sigma2Bank = new double[numRegPix];
		for(int i=0; i<numRegPix; i++) sigma2Bank[i] = res0[i][2];
		double[] r = new double[numRegPix];
		for(int i=0; i<numRegPix; i++) r[i] = res0[i][3];

		RegPixLstModel outRes = new RegPixLstModel();
		outRes.shiftBank = shiftBank;
		outRes.betaBank = betaBank;
		outRes.sigma2Bank = sigma2Bank;
		outRes.r = r;
		outRes.rg = res_rg00;

		return outRes;
	}




	//	Fit the time shift model ------------
	//	Search the given shift range only
	private PixModel findBestMatch_ls(double[] charCurve, double[] yCurve, ArrayList<Integer> rg){
		PixModel resPixMd = new PixModel();

		//	shift: positive value means Y starts before X does
		ArrayList<Double> X = new ArrayList<Double>();
		ArrayList<Double> Y = new ArrayList<Double>();
		PearsonCorr corrCalObj = new PearsonCorr();
		double r = -2;
		GaussianLUT calGaussObj = new GaussianLUT();
		myBasicMath calBsMathObj = new myBasicMath();

		//  p-value of correlation r between current X and Y
		double maxZscore = -100;
		int shift = 0;
		double bestR = r;
		double coeff = -100;
		double sigma2 = -100;
		double[] tempVec = new double[charCurve.length];
		int vecLength = charCurve.length;
		double tempZscore;
		double sdxy = Math.sqrt(calBsMathObj.sampleVar(yCurve)) / Math.sqrt(calBsMathObj.sampleVar(charCurve));
		int iRg;
		for(int j=0; j<rg.size(); j++){
			iRg = rg.get(j);
			X.clear();
			Y.clear();
			if(iRg > 0){
				for(int k=0; k<(vecLength-iRg); k++) X.add(charCurve[k]);
				for(int k=iRg; k<vecLength; k++) Y.add(yCurve[k]);
			}else{
				for(int k=-iRg; k<vecLength; k++) X.add(charCurve[k]);
				for(int k=0; k<(vecLength+iRg); k++) Y.add(yCurve[k]);
			}
			r = corrCalObj.calculateCorr(X, Y);
			tempZscore = Math.sqrt(X.size()-3.0)/2.0 * Math.log((1.0+r)/(1.0-r));
			tempZscore = tempZscore + calGaussObj.myQnorm(1.0/((Math.abs((double)iRg)+1.0)*2.0));
			if(tempZscore > maxZscore){
				maxZscore = tempZscore;
				shift = iRg;
				bestR = r;
				coeff = bestR*sdxy;
				for(int i=0; i<tempVec.length; i++){
					if(i < X.size()){
						tempVec[i] = (Y.get(i) - coeff*X.get(i));
						tempVec[i] = tempVec[i]*tempVec[i];
					}else{
						tempVec[i] = 0;
					}
				}
				sigma2 = calBsMathObj.vectorSum(tempVec) / (X.size()-1);
			}
		}

		resPixMd.shift = shift;
		resPixMd.coeff = coeff;
		resPixMd.sigma2 = Math.max(sigma2, 1e-16);
		resPixMd.r = bestR;
		return resPixMd;
	}



	double[] calCharX(double[][] dataMat, double[] betaBank, double[] shiftBank, double[] sigma2Bank){
		int numPix = dataMat.length;
		int numTimePoints = dataMat[0].length;
		double[] relatSigma2Bank = new double[sigma2Bank.length];
		myBasicMath calBsMathObj = new myBasicMath();
		double tempVal = calBsMathObj.vectorMax(sigma2Bank);
		for(int i=0; i<relatSigma2Bank.length; i++) relatSigma2Bank[i] = sigma2Bank[i] / tempVal;

		double[] fittedX = new double[numTimePoints];		
		double[] weight = new double[numTimePoints];
		double beta;
		double sigma2;
		int shift;
		double meanVal;
		double stdVal;

		for(int ii=0; ii<numPix; ii++){
			beta = betaBank[ii];
			sigma2 = relatSigma2Bank[ii];
			shift = (int)Math.round(shiftBank[ii]);

			if(shift==0){
				// The characteristic waveform is the weighted average of estimated X(t)s from all waveforms.
				// And the weights are signal-to-noise rates
				for(int j=0; j<numTimePoints; j++){
					fittedX[j] = fittedX[j] + dataMat[ii][j]*beta/sigma2;
					weight[j] = weight[j] + beta*beta/sigma2;
				}

			}else if(shift > 0){
				for(int j=0; j<(numTimePoints-shift); j++){
					fittedX[j] = fittedX[j] + dataMat[ii][j+shift]*beta/sigma2;
					weight[j] = weight[j] + beta*beta/sigma2;
				}

			}else{   			// shift<0
				for(int j=-shift; j<numTimePoints; j++){
					fittedX[j] = fittedX[j] + dataMat[ii][j+shift]*beta/sigma2;
					weight[j] = weight[j] + beta*beta/sigma2;
				}
			}
		}

		for(int i=0; i<numTimePoints; i++){
			if(weight[i]==0) weight[i] = 1;
			fittedX[i] = fittedX[i]/weight[i];
		}
		meanVal = calBsMathObj.sampleMean(fittedX);
		stdVal = Math.sqrt(calBsMathObj.sampleVar(fittedX, meanVal));
		for(int i=0; i<numTimePoints; i++){
			fittedX[i] = (fittedX[i]-meanVal)/stdVal;
		}

		return fittedX;
	}



	private double[][] calResidu(double[][] dataMat, double[] charX, double[] betaBank, double[] shiftBank, double[] sigma2Bank){
		int numPix = dataMat.length;
		int numTimePoints = dataMat[0].length;
		double[][] resiMat = new double[numPix][numTimePoints];

		myBasicMath calBsMathObj = new myBasicMath();
		double beta;
		int shift;
		double sigma2;
		double stdVal;
		double[] yCurve;
		int k;

		for(int ii=0; ii<numPix; ii++){
			beta = betaBank[ii];
			sigma2 = sigma2Bank[ii];
			shift = (int)Math.round(shiftBank[ii]);
			yCurve = dataMat[ii];
			stdVal = Math.sqrt(calBsMathObj.sampleVar(yCurve));
			for(k=0; k<numTimePoints; k++){
				resiMat[ii][k] = yCurve[k]/stdVal * Math.sqrt(sigma2);
			}

			if(shift==0){
				for(k=0; k<numTimePoints; k++){
					resiMat[ii][k] = yCurve[k] - beta*charX[k];
				}

			}else if(shift > 0){
				for(k=shift; k<numTimePoints; k++){
					resiMat[ii][k] = yCurve[k] - beta*charX[k-shift];
				}

			}else{   			// shift<0
				for(k=0; k<(numTimePoints+shift); k++){
					resiMat[ii][k] = yCurve[k] - beta*charX[k-shift];
				}
			}
		}
		return resiMat;
	}



	private double[][] get_z_map_charX_residual(double[] charX, ArrayList<Integer[]> curRegPixList, ArrayList<Integer[]> curNeighborList,
			double[][] resiMatCore, double[][] resiMatNeib, double[] regShiftBank){

		double[][] diffZ = new double[imageStack.length][imageStack[0].length];
		double[][] charCorrMap;
		double[][] residCorrMap;
		int nTimePt = charX.length;
		double tempVal1, tempVal2;

		charCorrMap = calCorrMap_charX(charX, curRegPixList);
		residCorrMap = calCorrMap_res(curRegPixList, curNeighborList, resiMatCore, resiMatNeib);
		
		double constant1 = Math.sqrt(nTimePt-3);
		double constant2 = Math.sqrt(0.9088657);
		double constant3 = Math.sqrt(0.516);
		double constant4 = Math.sqrt(2);

		for(int x=0; x<diffZ.length; x++){
			for(int y=0; y<diffZ[0].length; y++){
				if(charCorrMap[x][y]==-100){
					diffZ[x][y] = 0;
				}else{
					tempVal1 = Math.log((1+charCorrMap[x][y])/(1-charCorrMap[x][y])) * constant1 / 2;
					tempVal1 = (tempVal1 - 0.385904) / constant2;
					tempVal2 = Math.log((1+residCorrMap[x][y])/(1-residCorrMap[x][y])) * constant1 / 2;
					tempVal2 = (tempVal2 - 1.02937) / constant3;
					diffZ[x][y] = (tempVal1 - tempVal2) / constant4;
				}
			}
		}

		return diffZ;
	}



	//  Correlation between pixels and charX ------------------
	private double[][] calCorrMap_charX(double[] charX, ArrayList<Integer[]> curRegPixList){
		double[][] charCorrMap = new double[imageStack.length][imageStack[0].length];
		int x,y,k;
		for(x=0; x<charCorrMap.length; x++){
			for(y=0; y<charCorrMap[0].length; y++){
				charCorrMap[x][y] = -100;
			}
		}
		int numRegPix = curRegPixList.size();
		int vecLength = imageStack[0][0].length; 
		double[] yCurve = new double[vecLength];
		PearsonCorr calCorObj = new PearsonCorr();

		for(int ii=0; ii<numRegPix; ii++){
			x = curRegPixList.get(ii)[0];
			y = curRegPixList.get(ii)[1];
			for(k=0; k<vecLength; k++) yCurve[k] = (double)imageStack[x][y][k];
			charCorrMap[x][y] = calCorObj.calculateCorr(charX, yCurve);
		}
		return charCorrMap;
	}


	//	Correlation between neighboring residuals ------------------
	private double[][] calCorrMap_res(ArrayList<Integer[]> curRegPixList, ArrayList<Integer[]> curNeighborList, 
			double[][] resiMatCore, double[][] resiMatNeib){

		double[][] residCorrMap = new double[imageStack.length][imageStack[0].length];
		// Initialize residCorrMap
		int x,y;
		for(x=0; x<residCorrMap.length; x++){
			for(y=0; y<residCorrMap[0].length; y++){
				residCorrMap[x][y] = -100;
			}
		}
		int numTimePoints = imageStack[0][0].length;
		double[] V0 = new double[numTimePoints];
		double[] Vnbavg = new double[numTimePoints];   // average of 8 neighbors
		PearsonCorr calCorObj = new PearsonCorr();
		int ii;
		int k;
		int idPixInLst;
		int deltx, delty;

		////	idInLstMap -----	 positive: (id in curRegPixList + 1) ;   negative: (-id in curNeighborList - 1)
		int[][] idInLstMap = new int[imageStack.length][imageStack[0].length]; 
		for(x=0; x<idInLstMap.length; x++){
			for(y=0; y<idInLstMap[0].length; y++){
				idInLstMap[x][y] = 0;
			}
		}
		for(ii=0; ii<curRegPixList.size(); ii++){
			x = curRegPixList.get(ii)[0];
			y = curRegPixList.get(ii)[1];
			idInLstMap[x][y] = ii + 1;
		}
		for(ii=0; ii<curNeighborList.size(); ii++){
			x = curNeighborList.get(ii)[0];
			y = curNeighborList.get(ii)[1];
			idInLstMap[x][y] = -ii - 1;
		}

		// Calculate residCorrMap: only consider pixels in the region
		for(ii=0; ii<curRegPixList.size(); ii++){
			x = curRegPixList.get(ii)[0];
			y = curRegPixList.get(ii)[1];

			for(k=0; k<numTimePoints; k++){
				// Reset the residual curve of current pixel (x,y)
				V0[k] = resiMatCore[ii][k];
				// Reset the average curve of 8 neighbors
				Vnbavg[k] = 0;
			}
			for(deltx=-1; deltx<=1; deltx++){
				for(delty=-1; delty<=1; delty++){
					if(deltx==0 && delty==0) continue;
					idPixInLst = idInLstMap[x+deltx][y+delty];
					if(idPixInLst > 0){   // if this neighbor is in the region (in curRegPixList)
						idPixInLst --;
						for(k=0; k<numTimePoints; k++){
							Vnbavg[k] += resiMatCore[idPixInLst][k];
						}
					}else{        // if this neighbor is out of the region (in curNeighborList)
						idPixInLst = - (idPixInLst + 1);
						for(k=0; k<numTimePoints; k++){
							Vnbavg[k] += resiMatNeib[idPixInLst][k];
						}
					}
				}
			}
			residCorrMap[x][y] = calCorObj.calculateCorr(V0, Vnbavg);
		}
		return residCorrMap;
	}
}
