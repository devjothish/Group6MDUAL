package MDUAL;
import java.io.*;
import java.util.*;

public class MergedClassFunctional {
	private double[] minValues;
	private double[] maxValues;
	private ArrayList<Integer> priorityList;
	public int dim;
	public int subDim;
	public int id;
	public double R;
	public int K;
	public int W;
	public int S;
	public HashSet<Tuple> outliers;
	public int nnToFind;
	private BufferedReader br; 
	private String filePath;
	public int maxW;
	public int gcdS;
	public double minR;
	public int slideID;		
	public double[] value;
	public ArrayList<Short> subDimCellIdx, fullDimCellIdx;
	public int nn;
	public HashSet<Integer> outlierQueryIDs;
	/** Variables for SOP **/
	public LinkedHashMap<Tuple,Double> LSky; // "old" <points, normalized distance> are placed as "last"
	public int[] layerCount;
	public HashMap<Integer,Boolean> safeDue; // for query id, until when this tuple is safe inlier with respect to the slide id of newest slide
	/** Variables for pMCSKY **/
	public int mc;

	// Contents from DataLoader.java
	public void DataLoader(String dataset) throws IOException {
		filePath = "datasets/"+dataset+".csv";
		br = new BufferedReader(new FileReader(filePath));
		String line = br.readLine();
		String[] rawValues = line.split(",");
		dim = rawValues.length;
		subDim = (dim > 15 ? 3 : rawValues.length); //default sub-dimensionality of high-dim(>15) data set is 3 (i.e., cells are created by the three dimensionalities).
		minValues = new double[dim];
		maxValues = new double[dim];
		priorityList = new ArrayList<Integer>();
		
		for(int i = 0; i < dim; i++) {
			minValues[i] = Double.MAX_VALUE;
			maxValues[i] = Double.MIN_VALUE;
			priorityList.add(i);
		}
		while(line!=null) {
			rawValues = line.split(",");
			for (int i = 0; i < dim; i++) {
				double value = Double.parseDouble(rawValues[i]);
				if(minValues[i]>value) minValues[i] = value;
				if(maxValues[i]<value) maxValues[i] = value;
			}
			line = br.readLine();
		}
	}

	public ArrayList<Tuple> getNewSlideTuples(int itr, int S) throws IOException {
		ArrayList<Tuple> newSlide = new ArrayList<Tuple>();
		BufferedReader br = new BufferedReader(new FileReader(filePath));
		String line = br.readLine();
		int tid = 0;
		
		while(line!=null) {
			if(tid>=itr*S) {
				String[] rawValues = line.split(",");
				double[] value = new double[dim];
				for(int i = 0; i<dim; i++) value[i] = Double.parseDouble(rawValues[priorityList.get(i)]);

				Tuple tuple = new Tuple(tid, itr, value);
				newSlide.add(tuple);
			}
			tid++;
			if(tid==(itr+1)*S) break;
			line = br.readLine();
		}
		return newSlide;
	}
	
	public double[] getMinValues() {
		
		return minValues;
	}
	
	public double[] getMaxValues() {
		
		return maxValues;
	}

	// Contents from Query.java
	public void Query(int id, double R, int K, int W, int S) {
		this.id = id;
		this.R = R;
		this.K = K;
		this.W = W;
		this.S = S;
	}
	
	@Override
	public boolean equals(Object obj){
		Query _obj = (Query) obj;
		return _obj.id == this.id;
	}
	// Contents from QueryLoader.java
	public void QueryLoader(String queryset) throws IOException {
		filePath = "Generatedqueryset/"+queryset+".csv";
		br = new BufferedReader(new FileReader(filePath));
		maxW = Integer.MIN_VALUE;
		gcdS = Integer.MAX_VALUE;
		minR = Integer.MAX_VALUE;
		
		String line = br.readLine();		
		while(line!=null) {
			String[] rawValues = line.split(",");
			double R = Double.parseDouble(rawValues[3]);
			int W = Integer.parseInt(rawValues[5]);
			int S = Integer.parseInt(rawValues[6]);
			if(maxW<W) maxW = W;
			if(gcdS>S) gcdS = S;
			if(minR>R) minR = R;
			line = br.readLine();
		}
	}
	
	public HashMap<Integer,Query> getQuerySet(int curr_itr) throws IOException {
		HashMap<Integer,Query> querySet = new HashMap<Integer,Query>();
		br = new BufferedReader(new FileReader(filePath));
		String line = br.readLine();		
		
		while(line!=null) {
			String[] rawValues = line.split(",");
			int id = Integer.parseInt(rawValues[0]); 
			int s_time = Integer.parseInt(rawValues[1]);
			int e_time = Integer.parseInt(rawValues[2]);
			if(s_time <= curr_itr && curr_itr < e_time) {
				double R = Double.parseDouble(rawValues[3]);
				int K = Integer.parseInt(rawValues[4]);
				int W = Integer.parseInt(rawValues[5]);
				int S = Integer.parseInt(rawValues[6]);
				
				Query query = new Query(id, R,K,W,S);
				querySet.put(id, query);
			}
			line = br.readLine();
		}
		return querySet;
	}
	
	public HashMap<Integer,Query> getQuerySetByQID(int fromQID, int numQueries) throws IOException {
		HashMap<Integer,Query> querySet = new HashMap<Integer,Query>();
		br = new BufferedReader(new FileReader(filePath));
		String line = br.readLine();		
		
		while(line!=null) {
			String[] rawValues = line.split(",");
			int id = Integer.parseInt(rawValues[0]); 
			if(id >= fromQID && id < fromQID+numQueries) {
				double R = Double.parseDouble(rawValues[3]);
				int K = Integer.parseInt(rawValues[4]);
				int W = Integer.parseInt(rawValues[5]);
				int S = Integer.parseInt(rawValues[6]);
				
				Query query = new Query(id, R,K,W,S);
				querySet.put(id, query);
			}
			line = br.readLine();
		}
		return querySet;
	}

}