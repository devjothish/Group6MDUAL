package MDUAL;
import java.io.IOException;

import Coreprocessing.QueryGenerator;
import Coreprocessing.Simulator;

public class Final{
	public static void main(String[] args) throws IOException {
		String dataset = "STK";
		int numQueries = 10;
		double changedQRatio = 0.2;
		String[] varyingParams = new String[] {"R","K","S","W"};
		int defaultW = 1000;
		int gcdS = 50;
		int defaultK = 5;
		int nW = 10;
		int repeatNum = 5;
		int variationTimes = 10;
		QueryGenerator queryGen = new QueryGenerator(dataset, defaultW, gcdS, defaultK, variationTimes);
		
		System.out.println(String.format("%-10s %10s %10s %10s %10s %10s %10s %10s", "Dataset", "Queryset", "ChgQRatio", "Time", "AvgMem", "PeakMem", "#Out", "#OutQ"));
		for (int i = 0; i < repeatNum; i++) {
				String queryset = queryGen.generate(numQueries, nW, varyingParams);
				Simulator sim = new Simulator(dataset, queryset);
				sim.run(nW, numQueries, changedQRatio);
			}
		System.out.println();
	}
}
