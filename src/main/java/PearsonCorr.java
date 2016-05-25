import java.util.ArrayList;

class PearsonCorr {
	double correlationAB = -2;


	double calculateCorr(double[] curveA1, double[] curveB1){

		if(curveA1.length != curveB1.length){
			//			System.out.println("Error: Curves for correlation calculation do not match in length");
			return -2;
		}

		double[] curveA = new double[curveA1.length];
		double[] curveB = new double[curveB1.length];
		for(int k=0; k<curveA1.length;k++){
			curveA[k] = curveA1[k];
			curveB[k] = curveB1[k];
		}

		myBasicMath matCalObj = new myBasicMath();

		double meanA = matCalObj.sampleMean(curveA);
		double meanB = matCalObj.sampleMean(curveB);
		double stdA = Math.sqrt(matCalObj.sampleVar(curveA, meanA));
		double stdB = Math.sqrt(matCalObj.sampleVar(curveB, meanB));
		correlationAB = 0;

		for(int k=0; k<curveA.length; k++){
			curveA[k] -= meanA;
			curveA[k] /= stdA;
			curveB[k] -= meanB;
			curveB[k] /= stdB;
			correlationAB += curveA[k]*curveB[k];
		}
		correlationAB /= (curveA.length-1);

		return correlationAB;
	}


	double calculateCorr(ArrayList<Double> curveA1, ArrayList<Double> curveB1){

		if(curveA1.size() != curveB1.size()){
			//			System.out.println("Error: Curves for correlation calculation do not match in length");
			return -2;
		}

		double[] curveA = new double[curveA1.size()];
		double[] curveB = new double[curveB1.size()];
		for(int k=0; k<curveA1.size();k++){
			curveA[k] = curveA1.get(k);
			curveB[k] = curveB1.get(k);
		}

		myBasicMath matCalObj = new myBasicMath();

		double meanA = matCalObj.sampleMean(curveA);
		double meanB = matCalObj.sampleMean(curveB);
		double stdA = Math.sqrt(matCalObj.sampleVar(curveA, meanA));
		double stdB = Math.sqrt(matCalObj.sampleVar(curveB, meanB));
		correlationAB = 0;

		for(int k=0; k<curveA.length; k++){
			curveA[k] -= meanA;
			curveA[k] /= stdA;
			curveB[k] -= meanB;
			curveB[k] /= stdB;
			correlationAB += curveA[k]*curveB[k];
		}
		correlationAB /= (curveA.length-1);

		return correlationAB;
	}

}
