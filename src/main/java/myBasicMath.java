
class myBasicMath {
	
	double vectorSum(double[] curveA){
		double sumVal = 0;
		for(int k=0; k<curveA.length; k++){
			sumVal += curveA[k];
		}
		return sumVal;
	}
	
	double vectorMax(double[] curveA){
		double maxVal = curveA[0];
		for(int k=1; k<curveA.length; k++){
			if(curveA[k] > maxVal) maxVal = curveA[k];
		}
		return maxVal;
	}
	
	double matrix2DMax(double[][] matA){
		double maxVal = matA[0][0];
		for(int i=0; i<matA.length; i++){
			for(int j=0; j<matA[0].length; j++){
				if(matA[i][j] > maxVal) maxVal = matA[i][j];
			}
		}
		return maxVal;
	}
	
	int matrix2DMax(int[][] matA){
		int maxVal = matA[0][0];
		for(int i=0; i<matA.length; i++){
			for(int j=0; j<matA[0].length; j++){
				if(matA[i][j] > maxVal) maxVal = matA[i][j];
			}
		}
		return maxVal;
	}
	
	double vectorMin(double[] curveA){
		double minVal = curveA[0];
		for(int k=1; k<curveA.length; k++){
			if(curveA[k] < minVal) minVal = curveA[k];
		}
		return minVal;
	}
	
	double matrix2DMin(double[][] matA){
		double minVal = matA[0][0];
		for(int i=0; i<matA.length; i++){
			for(int j=0; j<matA[0].length; j++){
				if(matA[i][j] < minVal) minVal = matA[i][j];
			}
		}
		return minVal;
	}
	
	int matrix2DMin(int[][] matA){
		int minVal = matA[0][0];
		for(int i=0; i<matA.length; i++){
			for(int j=0; j<matA[0].length; j++){
				if(matA[i][j] < minVal) minVal = matA[i][j];
			}
		}
		return minVal;
	}
	
	
	
	double sampleMean(double[] curveA){
		double meanVal = 0;
		for(int k=0; k<curveA.length; k++){
			meanVal += curveA[k];
		}
		meanVal /= curveA.length;
		return meanVal;
	}
	
	double sampleVar(double[] curveA){
		double meanVal = 0;
		for(int k=0; k<curveA.length; k++){
			meanVal += curveA[k];
		}
		meanVal /= curveA.length;
		
		double smpVarVal = 0;
		for(int k=0; k<curveA.length; k++){
			smpVarVal += (curveA[k] - meanVal)*(curveA[k] - meanVal);
		}
		smpVarVal /= (curveA.length-1);
		
		return smpVarVal;
	}
	
	double sampleVar(double[] curveA, double meanVal){
		double smpVarVal = 0;
		for(int k=0; k<curveA.length; k++){
			smpVarVal += (curveA[k] - meanVal)*(curveA[k] - meanVal);
		}
		smpVarVal /= (curveA.length-1);
		
		return smpVarVal;
	}
	
	
}
