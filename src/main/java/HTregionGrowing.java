import java.util.ArrayList;
import java.util.Collections;

/**
 * Region growing based on order statistics guided hypothesis testing (z-score map binarization method)
 *  
 * @author Yinxue
 */


class HTregionGrowing {
	double[][] zMap;
	double threshPixZ;
	double threshRegZ;
	int minSize = 10;
	int[][] connDmIDmap;
	int nConnDomain = 0;

	private boolean[][] availPixMap;   // This available map only used in this class
	private double ZcurReg = -100;
	private boolean[][] newRegPosMap;

	HTregionGrowing(){};



	HTregionGrowing(double[][] zMap1,boolean[][] availMap1,double threshPixZ1,double threshRegZ1,
			int minSize1,boolean isUnitSearch){
		zMap = zMap1;
		availPixMap = new boolean[zMap.length][zMap[0].length];
		for(int x=0; x<zMap.length; x++){
			for(int y=0; y<zMap[0].length; y++){
				availPixMap[x][y] = availMap1[x][y];
			}
		}
		threshPixZ = threshPixZ1;
		threshRegZ = threshRegZ1;
		minSize = minSize1;
		nConnDomain = 0;
		connDmIDmap = new int[zMap.length][zMap[0].length];
		newRegPosMap = new boolean[availPixMap.length][availPixMap[0].length];

		this.doBinarization(isUnitSearch);
	}



	void reDoBinarization(double[][] zMap1,boolean[][] availMap1,double threshPixZ1,double threshRegZ1,
			int minSize1,boolean isUnitSearch){
		zMap = zMap1;
		for(int x=0; x<zMap.length; x++){
			for(int y=0; y<zMap[0].length; y++){
				availPixMap[x][y] = availMap1[x][y];
			}
		}
		threshPixZ = threshPixZ1;
		threshRegZ = threshRegZ1;
		minSize = minSize1;
		nConnDomain = 0;
		for(int x=0; x<connDmIDmap.length; x++){
			for(int y=0; y<connDmIDmap[0].length; y++){
				connDmIDmap[x][y] = 0;
			}
		}

		this.doBinarization(isUnitSearch);
	}




	boolean[][] getMaxSizedConnDmMap(){
		boolean[][] maxConnDmMap = new boolean[connDmIDmap.length][connDmIDmap[0].length];
		int[] sizeConnDm = new int[nConnDomain + 1];
		int maxSize = 0;
		int idMaxConnDm = 0;

		if(nConnDomain==0) return null;

		// Count the size of each connected domain
		for(int k=0; k<=nConnDomain; k++){
			sizeConnDm[k] = 0;
		}
		for(int i=0; i<connDmIDmap.length; i++){
			for(int j=0; j<connDmIDmap[0].length; j++){
				sizeConnDm[connDmIDmap[i][j]]++;
			}
		}

		// Get the largest size
		for(int k=1; k<sizeConnDm.length; k++){
			if(sizeConnDm[k] > maxSize){
				maxSize = sizeConnDm[k];
				idMaxConnDm = k;
			}
		}

		// The boolean map of the largest domain
		for(int i=0; i<connDmIDmap.length; i++){
			for(int j=0; j<connDmIDmap[0].length; j++){
				if(connDmIDmap[i][j]==idMaxConnDm){
					maxConnDmMap[i][j] = true;	
				}else{
					maxConnDmMap[i][j] = false;
				}
			}
		}
		return maxConnDmMap;
	}




	private void doBinarization(boolean isUnitSearch){
		boolean flagContinue;
		double maxPixZ = -100;
		int remainTrueCount;
		int[] seedPix = new int[2];

		// Initialize the indicator variables and the seed pixel
		remainTrueCount = 0;
		for(int x=0; x<zMap.length; x++){
			for(int y=0; y<zMap[0].length; y++){
				if(availPixMap[x][y]){
					remainTrueCount += 1;
					if(zMap[x][y] > maxPixZ){
						maxPixZ = zMap[x][y];
						seedPix[0] = x;
						seedPix[1] = y;
					}
				}
			}
		}
		if(isUnitSearch){
			flagContinue = remainTrueCount > minSize;
		}else{
			flagContinue = (remainTrueCount > minSize) && (maxPixZ > threshPixZ);
		}


		while(flagContinue){
			// Find a new region starting from the seed, get its position map and region-wide z-score
			// Update the position map newRegPosMap and ZcurReg
			
			htrgSingleRun(seedPix,availPixMap);

			// Record the new founded "active region" (possibly not isolated)
			if(ZcurReg > threshRegZ){
				nConnDomain += 1;
				for(int x=0; x<availPixMap.length; x++){
					for(int y=0; y<availPixMap[0].length; y++){
						if(newRegPosMap[x][y]){
							availPixMap[x][y] = false;
							connDmIDmap[x][y] = nConnDomain;
							newRegPosMap[x][y] = false;
						}
					}
				}
			}else{
				for(int x=0; x<availPixMap.length; x++){
					for(int y=0; y<availPixMap[0].length; y++){
						if(newRegPosMap[x][y]){
							availPixMap[x][y] = false;
							newRegPosMap[x][y] = false;
						}
					}
				}
			}

			// Update the indicator variables and the seed pixel
			remainTrueCount = 0;
			maxPixZ = -100;
			for(int x=0; x<zMap.length; x++){
				for(int y=0; y<zMap[0].length; y++){
					if(availPixMap[x][y]){
						remainTrueCount += 1;
						if(zMap[x][y] > maxPixZ){
							maxPixZ = zMap[x][y];
							seedPix[0] = x;
							seedPix[1] = y;
						}
					}
				}
			}	

			if(isUnitSearch){
				flagContinue = remainTrueCount > minSize;
			}else{
				flagContinue = (remainTrueCount > minSize) && (maxPixZ > threshPixZ);
			}
		}

		// Merge regions that are connected to each other, updating connDmIDmap and nConnDomain
		this.updateConnDomains();
	}





	private void htrgSingleRun(int[] seedPix,boolean[][] availMap){
		// Find a new region starting from the seed, get its position map newRegPosMap and region-wide z-score ZcurReg
		// Return the position map and update ZcurReg

		ArrayList<Double[]> curRegPixList = new ArrayList<Double[]>();
		ArrayList<Double[]> curNeighborList = new ArrayList<Double[]>();
		boolean[][] neighbAvailMap = new boolean[availMap.length][availMap[0].length];
		ArrayList<Double> Z_reg_Nb = new ArrayList<Double>();   // z-score list of pixels in curRegPixList and curNeighborList 
		ArrayList<Double> Z_nb = new ArrayList<Double>();
		ArrayList<Double> Z_nb_sorted = new ArrayList<Double>();  // ranked z-scores of neighbors (descending order)
		ArrayList<Integer> indInOri_Zsorted = new ArrayList<Integer>();  // index of ordered z-scores in neighbors (descending order)
		ArrayList<Double[]> curNeighborList_sorted = new ArrayList<Double[]>(); // list of current neighbors being ranked by z-score
		int tempInd = 0;
		ArrayList<Double> cumSumZ_reg_Nb = new ArrayList<Double>(); // cumulative sum of z-score in set C
		double tempSumZ = 0;
		ArrayList<Boolean> flagNbConsidered = new ArrayList<Boolean>();
		ArrayList<Double[]> curNeighborList_enclose = new ArrayList<Double[]>();
		Double[] tempPix;
		ArrayList<Double> Z_reg_Nb_sorted = new ArrayList<Double>();
		ArrayList<Integer> orderOfZ_pixInZregNb = new ArrayList<Integer>();  // order of pixels in Z_reg_Nb (ascending order)
		ArrayList<Integer> subOrderOfZ_pixInZregNb = new ArrayList<Integer>();  // a sublist of orderOfZ_pixInZregNb
		int totN;
		ArrayList<Double> candidateRegZ = new ArrayList<Double>();
		int numEquivPix;
		int numPixToAdd;
		int selN;
		double biasBestZ;
		double varBestZ;
		double tempBestZ = -100;
		double bestZ = -100;
		
		////   Initiation
		ZcurReg = 0;   // clear the current region's z-score  (the output of this function)
		for(int x=1; x<(neighbAvailMap.length-1); x++){   // (Exclude the borders)
			for(int y=1; y<(neighbAvailMap[0].length-1); y++){
				neighbAvailMap[x][y] = availMap[x][y];
			}
		}

		// The initial region contains only the seed pixel
		curRegPixList.add(new Double[]{(double)seedPix[0],(double)seedPix[1],zMap[seedPix[0]][seedPix[1]]});
		neighbAvailMap[seedPix[0]][seedPix[1]] = false;      // Exclude pixels already in the region pixel list
		numPixToAdd = 1;

		while(numPixToAdd!=0){
			// Find neighbors that can be considered in this iteration by re-deciding curNeighborList
			findNeighbors(curNeighborList,curRegPixList,neighbAvailMap);  

			// The first part of Z_reg_Nb is consist of z of pixels in curRegPixList
			Z_reg_Nb.clear();
			for(int i=0; i<curRegPixList.size(); i++){
				Z_reg_Nb.add(curRegPixList.get(i)[2]);
			}

			//  Sort curNeighborList by z's, get the ordered z's as Z_nb_sorted (descending order)
			Z_nb.clear();
			Z_nb_sorted.clear();
			indInOri_Zsorted.clear();
			curNeighborList_sorted.clear();
			for(int i=0; i<curNeighborList.size(); i++){
				Z_nb.add(curNeighborList.get(i)[2]);
				Z_nb_sorted.add(-curNeighborList.get(i)[2]);
			} 
			Collections.sort(Z_nb_sorted); // negative z in ascending order (z in descending order)
			for(int i=0; i<Z_nb_sorted.size(); i++){
				Z_nb_sorted.set(i, -Z_nb_sorted.get(i));  // the i-th max z in Z_nb
				tempInd = Z_nb.indexOf(Z_nb_sorted.get(i));  // index of the i-th max z in Z_nb (also in curNeighborList)
				indInOri_Zsorted.add(tempInd);
				curNeighborList_sorted.add(curNeighborList.get(tempInd));  // sort the pixels only by redirecting pointers
			}

			// Let cumSumZ_reg_Nb[0] be the sum of z of all pixels in curRegPixList
			cumSumZ_reg_Nb.clear();
			tempSumZ = 0;
			for(int i=0; i<curRegPixList.size(); i++){
				tempSumZ += curRegPixList.get(i)[2];
			}
			cumSumZ_reg_Nb.add(tempSumZ);

			flagNbConsidered.clear();
			for(int i=0; i<curNeighborList_sorted.size(); i++){
				flagNbConsidered.add(false);
			}


			// The following part mainly considers the "enclosing" problem
			curNeighborList_enclose.clear();
			for(int i=0; i<curNeighborList_sorted.size(); i++){  // for each currently potential neighbor (in lstNeighbor_sorted) i
				// 		cumSumZ_reg_Nb is an accumulating series: cumSumZ_reg_Nb[k]=cumSumZ_reg_Nb[k-1]+z[k] 
				// 		(k is an index in curNeighborList_enclose)

				if(flagNbConsidered.get(i)) continue;

				curNeighborList_enclose.add(new Double[]{curNeighborList_sorted.get(i)[0],curNeighborList_sorted.get(i)[1],
						curNeighborList_sorted.get(i)[2]});
				tempSumZ = cumSumZ_reg_Nb.get(cumSumZ_reg_Nb.size()-1) + curNeighborList_sorted.get(i)[2];
				flagNbConsidered.set(i, true);

				// 		The second part of Z_reg_Nb are z of current neighbors (curNeighborList_enclose)
				Z_reg_Nb.add(curNeighborList_sorted.get(i)[2]);

				// 		Find the enclosed pixels after we take into account the above pixel
				if(i!=(curNeighborList_sorted.size()-1)){
					for(int j=i+1; j<curNeighborList_sorted.size(); j++){
						if(!flagNbConsidered.get(j)){
							tempPix = curNeighborList_sorted.get(j);
							if(checkEnclose(tempPix,curRegPixList,curNeighborList_enclose,flagNbConsidered)){ 
								// If pixel j IS enclosed by the currently already-considered pixels
								curNeighborList_enclose.add(new Double[]{tempPix[0],tempPix[1],tempPix[2]});
								tempSumZ += tempPix[2];
								cumSumZ_reg_Nb.add(0.0);
								flagNbConsidered.set(j, true);
								Z_reg_Nb.add(tempPix[2]);
							}
						}
					}
				}
				cumSumZ_reg_Nb.add(tempSumZ);
			}

			// The score (without correction) mu*sqrt(n) 
			//    candidateRegZ is the linear function L of ordered variable 
			//    Now numEquivPix is the real number of corresponding pixels
			candidateRegZ.clear();
			for(int i=0; i<cumSumZ_reg_Nb.size(); i++){
				numEquivPix = curRegPixList.size() + i;
				candidateRegZ.add(cumSumZ_reg_Nb.get(i) / Math.sqrt(numEquivPix));
			}

			// Find the value and index of the max mu*sqrt(n) in candidateRegZ
			tempBestZ = Collections.max(candidateRegZ);
			numPixToAdd = candidateRegZ.indexOf(tempBestZ);

			// The number of pixels in the new region (to be used to replace the old)
			selN = curRegPixList.size() + numPixToAdd;   //n_C


			// Correct for the order statistics
			//// We need to find the orders of the selected pixels in the whole C (including A and B')

			////  Order all pixels in curNeighborList_enclose (including both curRegPixList and curNeighborPixList)
			// curNeighborList_enclose and Z_reg_Nb are corresponding in the indices of pixels
			Z_reg_Nb_sorted.clear();
			for(int i=0; i<Z_reg_Nb.size(); i++){
				Z_reg_Nb_sorted.add(-Z_reg_Nb.get(i));
			}
			Collections.sort(Z_reg_Nb_sorted);
			//  Find the orders of pixels in curNeighborList_enclose in terms of their z-scores
			// (orderOfZ_pixInZregNb indices indicate the indices in curNeighborList_enclose (i.e. in Z_reg_Nb),
			//   while values indicate the ranks)
			for(int i=0; i<Z_reg_Nb_sorted.size(); i++){
				Z_reg_Nb_sorted.set(i, -Z_reg_Nb_sorted.get(i));  // the i-th max z in Z_reg_Nb
			}
			orderOfZ_pixInZregNb.clear();
			for(int i=0; i<Z_reg_Nb.size(); i++){
				tempInd = Z_reg_Nb_sorted.indexOf(Z_reg_Nb.get(i));    // (descending order rank - 1)
				tempInd = Z_reg_Nb.size()-tempInd;
				orderOfZ_pixInZregNb.add(tempInd);      // ascending order of this pixel's z
			}

			// The total number of pixels in the whole current region + neighborhood
			totN = curRegPixList.size() + curNeighborList_sorted.size();   // (n_A + n_B) 

			// Mean of linear function L of order statistics
			subOrderOfZ_pixInZregNb.clear();
			for(int i=0; i<selN; i++){
				subOrderOfZ_pixInZregNb.add(orderOfZ_pixInZregNb.get(i));  // sublist of orderOfZ_pixInZregNb
			}
			biasBestZ = calBias_old(totN, subOrderOfZ_pixInZregNb);

			// Variance of linear function L of order statistics
			varBestZ = calVar_old(totN,subOrderOfZ_pixInZregNb);

			subOrderOfZ_pixInZregNb.clear();
			for(int i=0; i<totN; i++) subOrderOfZ_pixInZregNb.add(i+1);
			varBestZ = varBestZ * totN / calVar_old(totN,subOrderOfZ_pixInZregNb);

			tempBestZ = (Collections.max(cumSumZ_reg_Nb) - biasBestZ) / Math.sqrt(varBestZ);

			if(tempBestZ < bestZ) break;     // Stop searching when regZ reaches the max

			bestZ = tempBestZ;
			////   Now, update the region!
			if(numPixToAdd!=0){    // Already got to the maximum Zreg
				for(int i=0; i<numPixToAdd; i++){
					tempPix = curNeighborList_enclose.get(i);
					curRegPixList.add(tempPix);
					neighbAvailMap[(int)Math.round(tempPix[0])][(int)Math.round(tempPix[1])] = false;
				}
			}
		}

		for(int i=0; i<curRegPixList.size(); i++){
			tempPix = curRegPixList.get(i);
			newRegPosMap[(int)Math.round(tempPix[0])][(int)Math.round(tempPix[1])] = true;
		}
		ZcurReg = bestZ;
	}



	private void findNeighbors(ArrayList<Double[]> curNeighborList, ArrayList<Double[]> curRegPixList, boolean[][] neighbAvailMap){
		// Update curNeighborList1 and neighbAvailMap1 (by using the pointers)
		Double[] curPix;
		int neighbX, neighbY;

		curNeighborList.clear();

		for(int i=0; i<curRegPixList.size(); i++){
			curPix = curRegPixList.get(i);

			for(int deltx=-1; deltx<=1; deltx++){
				for(int delty=-1; delty<=1; delty++){
					if(deltx==0 && delty==0) continue;
					neighbX = (int)Math.round(curPix[0]+deltx);
					neighbY = (int)Math.round(curPix[1]+delty);
					if(neighbX<0 || neighbX>=neighbAvailMap.length || neighbY<0 || neighbY>=neighbAvailMap[0].length) continue;
					if(neighbAvailMap[neighbX][neighbY]){
						curNeighborList.add(new Double[]{(double)neighbX, (double)neighbY, zMap[neighbX][neighbY]});
						neighbAvailMap[neighbX][neighbY] = false;
					}
				}
			}
		}
	}
	
	
	
	private void findNeighborsInt(ArrayList<Integer[]> curNeighborList, ArrayList<Integer[]> curRegPixList, boolean[][] neighbAvailMap){
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
	
	


	private boolean checkEnclose(Double[] thePix, ArrayList<Double[]> curRegPixList, ArrayList<Double[]> curNeighborList_enclose,
			ArrayList<Boolean> flagNbConsidered){
		// In two cases we think of thePix as enclosed by the currently already-considered pixels:
		//  a) its all 4-connect neighbors have been considered
		//  b) its 5 or more 8-connect neighbors have been considered

		int x = (int)Math.round(thePix[0]);
		int y = (int)Math.round(thePix[1]);
		int nEnclosed4cnnNb = 0;
		int nEnclosed8cnnNb = 0;
		Double[] tempPix;
		int xTemp, yTemp;

		for(int i=0; i<curRegPixList.size(); i++){
			tempPix = curRegPixList.get(i);
			xTemp = (int)Math.round(tempPix[0]);
			yTemp = (int)Math.round(tempPix[1]);
			if((xTemp>=x-1) & (xTemp<=x+1) & (yTemp>=y-1) & (yTemp<=y+1)){
				nEnclosed8cnnNb ++;
				if((xTemp==x) | (yTemp==y)) nEnclosed4cnnNb ++;
			}
		}
		for(int i=0; i<curNeighborList_enclose.size(); i++){
			tempPix = curNeighborList_enclose.get(i);
			xTemp = (int)Math.round(tempPix[0]);
			yTemp = (int)Math.round(tempPix[1]);
			if((xTemp>=x-1) & (xTemp<=x+1) & (yTemp>=y-1) & (yTemp<=y+1)){
				nEnclosed8cnnNb ++;
				if((xTemp==x) | (yTemp==y)) nEnclosed4cnnNb ++;
			}
		}

		if(nEnclosed4cnnNb>=3 | nEnclosed8cnnNb>4){
			return true;
		}else{
			return false;
		}
	}


	private double calVar(int n_AB, ArrayList<Integer> order_pixInRegion){
		double var = 0;
		double tempIJterm = 0;
		double tempJterm = 0;
		double v_i;
		double v_j;
		double[] coef_of_eachOrder = new double[order_pixInRegion.size()];  // coef = 1/phi(Phi^-1(v_i))
		GaussianLUT calGaussObj = new GaussianLUT();

		Collections.sort(order_pixInRegion);

		for(int i=0; i<order_pixInRegion.size(); i++){
			v_i = ((double)order_pixInRegion.get(i) + 0.5) / ((double)n_AB + 1.0);
			coef_of_eachOrder[i] = 1 / calGaussObj.myDnorm( calGaussObj.myQnorm(v_i) );
		}

		tempIJterm = 0;
		for(int i=0; i<order_pixInRegion.size(); i++){
			tempJterm = 0;
			for(int j=i+1; j<order_pixInRegion.size(); j++){
				v_j = ((double)order_pixInRegion.get(j) + 0.5) / ((double)n_AB + 1.0);
				tempJterm += (1-v_j) * coef_of_eachOrder[j];  
			}

			v_i = ((double)order_pixInRegion.get(i) + 0.5) / ((double)n_AB + 1.0);
			tempIJterm += v_i*coef_of_eachOrder[i] * tempJterm;
		}
		tempIJterm *= 2.0;

		for(int i=0; i<order_pixInRegion.size(); i++){
			v_i = ((double)order_pixInRegion.get(i) + 0.5) / ((double)n_AB + 1.0);
			tempIJterm += v_i*coef_of_eachOrder[i] * (1-v_i)*coef_of_eachOrder[i];
		}

		var = tempIJterm / ((double)order_pixInRegion.size()* (double)n_AB);
		return var;

	}


	private double calVar_old(int n_AB, ArrayList<Integer> order_pixInRegion){
		ArrayList<Integer> order_pixInRegion_sorted = new ArrayList<Integer>();
		for(int i=0; i<order_pixInRegion.size(); i++)  order_pixInRegion_sorted.add(n_AB + 1 -order_pixInRegion.get(i));
		Collections.sort(order_pixInRegion_sorted); 
		GaussianLUT calGaussObj = new GaussianLUT();

		double tempVar = 0;
		double tempTemp = 0;
		double[] Finv = new double[n_AB+1];
		for(int i=0; i<=n_AB; i++){
			Finv[i] = - calGaussObj.myQnorm(((double)i+1.0)/((double)n_AB+2.0));
		}
		double[] difFJ = new double[n_AB];
		for(int i=0; i<n_AB; i++){
			difFJ[i] = Finv[i] - Finv[i+1];
		}

		double difFI;
		for(int i=0; i<order_pixInRegion_sorted.size(); i++){
			difFI = difFJ[order_pixInRegion_sorted.get(i)-1];
			tempTemp = 0;
			for(int j=i; j<order_pixInRegion_sorted.size(); j++){
				tempTemp +=  (1 - (double)order_pixInRegion_sorted.get(j)/((double)n_AB+1.0) ) * difFJ[order_pixInRegion_sorted.get(j)-1];
			}
			tempVar += 2 * ( (double)order_pixInRegion_sorted.get(i)/((double)n_AB+1.0) ) * difFI * tempTemp ;
		}
		return tempVar*(double)n_AB;
	}



	private double calBias(int n_AB, ArrayList<Integer> order_pixInRegion){
		double bias = 0;
		GaussianLUT calGaussObj = new GaussianLUT();

		for(int i=0; i<order_pixInRegion.size(); i++){
			bias += calGaussObj.myQnorm((double)order_pixInRegion.get(i)/((double)n_AB+1.0));
		}
		bias *= 1.0/Math.sqrt((double)order_pixInRegion.size());
		return bias;
	}



	private double calBias_old(int n_AB, ArrayList<Integer> order_pixInRegion){
		ArrayList<Integer> order_pixInRegion_sorted = new ArrayList<Integer>();
		for(int i=0; i<order_pixInRegion.size(); i++)  order_pixInRegion_sorted.add(n_AB + 1 -order_pixInRegion.get(i));
		double bias = 0;
		GaussianLUT calGaussObj = new GaussianLUT();

		for(int i=0; i<order_pixInRegion_sorted.size(); i++){
			bias -= calGaussObj.myQnorm((double)order_pixInRegion_sorted.get(i)/((double)n_AB+1.0));
		}
		return bias;
	}


	
	// Merge regions that are connected to each other, updating connDmIDmap and nConnDomain
	private void updateConnDomains(){
		// Turn the original connected domain map into binary
		boolean[][] binConnMap = new boolean[connDmIDmap.length][connDmIDmap[0].length];
		int x,y;
		for(x=0; x<binConnMap.length; x++){
			for(y=0; y<binConnMap[0].length; y++){
				binConnMap[x][y] = (connDmIDmap[x][y]>0);
			}
		}
		int[][] newConnDmIDmap = new int[connDmIDmap.length][connDmIDmap[0].length];
		ArrayList<Integer[]> badSeedList = new ArrayList<Integer[]>();
		ArrayList<Integer[]> badSeedHoleList = new ArrayList<Integer[]>();
		int tempPixCount;
		int newRegID = 0;
		for(x=0; x<binConnMap.length; x++){
			for(y=0; y<binConnMap[0].length; y++){
				if(binConnMap[x][y]){	  // once a unconsidered available pixel is found, search from it
					badSeedList.clear();
					badSeedHoleList.clear();
					badSeedList.add(new Integer[]{x,y});
					binConnMap[x][y] = false;
					tempPixCount = 1;
					findNeighborsInt(badSeedHoleList, badSeedList, binConnMap);
					while(!badSeedHoleList.isEmpty()){   // search until it cannot be enlarged any more
						for(int j=0; j<badSeedHoleList.size(); j++){
							badSeedList.add(badSeedHoleList.get(j));
							binConnMap[badSeedHoleList.get(j)[0]][badSeedHoleList.get(j)[1]] = false;
							tempPixCount ++;
						}
						findNeighborsInt(badSeedHoleList, badSeedList, binConnMap);
					}
					if(tempPixCount >= minSize){
						newRegID ++;
						for(int j=0; j<badSeedList.size(); j++){
							newConnDmIDmap[badSeedList.get(j)[0]][badSeedList.get(j)[1]] = newRegID;
						}
					}
				}
			}
		}
		for(x=0; x<connDmIDmap.length; x++){
			for(y=0; y<connDmIDmap[0].length; y++){
				connDmIDmap[x][y] = newConnDmIDmap[x][y];
			}
		}
		nConnDomain = newRegID;
	}
}



