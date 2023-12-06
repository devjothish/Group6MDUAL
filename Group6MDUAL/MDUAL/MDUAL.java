package MDUAL;
import java.util.*;
import Coreprocessing.*;

public class MDUAL{	
	public HashMap<ArrayList<Short>,Cell> slideIn, slideOut;
	public double[] minValues, dimLength, subDimLength;
	public double minR, maxR, minR_old, maxR_old;
	public HashMap<ArrayList<Short>,Integer> slideDeltaCnt;
	public boolean crDim, minChg;
	public HashMap<ArrayList<Short>,CellDetails> cardGrid, fullDimCardGrid;
	public int gcdS, dim, subDim, nS, nW;
	public LinkedList<HashMap<ArrayList<Short>,Cell>> slides; 
	public HashMap<Integer,Query> querySet;
	public boolean subDimFlag, maxRChanged, minRChanged;
	public LinkedList<HashMap<ArrayList<Short>,Integer>> fullDimCellSlidesCnt; 
	public HashSet<Tuple> outliers;
	
	
public class Main {
    private int dimensions, subDimensions;
    private int numberOfSlides, gcdSlideSize;
    private double[] minimumValues;
    private double minRadius, maxRadius;
    private boolean isSubDimensionActive;
    
    private HashMap<ArrayList<Short>, CellDetails> cardinalityGrid, fullDimensionCardinalityGrid;
    private LinkedList<HashMap<ArrayList<Short>, Cell>> slides, fullDimensionCellSlideCounts;
    private HashSet<Tuple> outliers;
    private HashMap<Integer, Query> querySet;

    public Main(int dimensions, int subDimensions, int numberOfSlides, int gcdSlideSize, double[] minimumValues) {
        this.dimensions = dimensions;
        this.subDimensions = subDimensions;
        this.numberOfSlides = numberOfSlides;
        this.minimumValues = minimumValues;
        this.gcdSlideSize = gcdSlideSize;

        this.minRadius = Double.MAX_VALUE;
        this.maxRadius = Double.MIN_VALUE;
        this.isSubDimensionActive = dimensions != subDimensions;

        initializeDataStructures();
    }

    private void initializeDataStructures() {
        cardinalityGrid = new HashMap<>();
        fullDimensionCardinalityGrid = new HashMap<>();
        slides = new LinkedList<>();
        fullDimensionCellSlideCounts = new LinkedList<>();
        outliers = new HashSet<>();
    }
}

		
		minChg = True;
		if (iteration == 0) {
        initializeForFirstIteration();
    }

    public MDUAL(int dimensions, int dimensions2, int totalSlides, int slideIntervalGCD, double[] lowerBounds) {
		}
	clearPreviousOutliers();
    updateWindow(newSlideTuples, iteration);
    updateBasicParameters(iteration); 
    findOutliersMain(iteration);

    return outliers;
}

private void initializeForFirstIteration() {
    updateBasicParameters(0);
    initCellSize();
}
	
	
	public void initCellSize() {
	    calculateDimLength();
	    calculateSubDimLength();
	}


	private void calculateDimLength() {
	    dimLength = new double[dim];
	    double factor = Math.sqrt(minR * minR / dim);
	    for (int i = 0; i < dim; i++) {
	        dimLength[i] = factor;
	    }
	}

	private void calculateSubDimLength() {
	    if (subDimFlag) {
	        subDimLength = new double[subDim];
	        double subDimFactor = Math.sqrt(minR * minR / subDim);
	        for (int i = 0; i < subDim; i++) {
	            subDimLength[i] = subDimFactor;
	        }
	    }
	}


	public void clearPreviousOutliers() {
	    for (Tuple outlier : outliers) {
	    	outlier.outlierQueryIDs.clear();
	    }
	    outliers.clear();
	}


	public void updateBasisParams(int itr) {
	    resetMinMaxValues();
	    updateMinMaxValues();
	    checkMinMaxChanges(itr);
	}

	private void resetMinMaxValues() {
	    maxRChanged = false;
	    maxR_old = maxR;
	    maxR = Double.MIN_VALUE;
	    minRChanged = false;
	    minR_old = minR;
	    minR = Double.MAX_VALUE;
	}

	private void updateMinMaxValues() {
	    for (Query q : querySet.values()) {
	        maxR = Math.max(maxR, q.R);
	        minR = Math.min(minR, q.R);
	    }
	}

	public void updateBasicParameters(int iteration) {
    resetMinMaxValues();
    updateMinMaxValues();
    checkForRadiusChanges(iteration);
}

private void checkForRadiusChanges(int iteration) {
    if (iteration > 0) {
        minRChanged = hasMinRadiusChanged();
        maxRChanged = hasMaxRadiusChanged();
    }
}

private boolean hasMinRadiusChanged() {
    return minR != previousMinRadius;
}

private boolean hasMaxRadiusChanged() {
    return maxR != previousMaxRadius;
}
	
	public void updateWindow(ArrayList<Tuple> slideTuples, int itr) {
	    slideIn = new HashMap<>();
	    HashMap<ArrayList<Short>, Integer> fullDimCellSlideInCnt = new HashMap<>();

	    processSlideTuples(slideTuples, fullDimCellSlideInCnt);

	    slides.add(slideIn);
	    if (subDimFlag) {
	        fullDimCellSlidesCnt.add(fullDimCellSlideInCnt);
	    }

	    slideDeltaCnt = new HashMap<>();
	    HashSet<ArrayList<Short>> newCellIndices = new HashSet<>();

	    updateSlidesAndCounts(slideIn, newCellIndices,itr);
	    updateFullDimCellSlideInCount(fullDimCellSlideInCnt, itr);
	    getNeighCellMap(newCellIndices);
	    updateSlideOut(itr);
	}

	private void processSlideTuples(ArrayList<Tuple> slideTuples, HashMap<ArrayList<Short>, Integer> fullDimCellSlideInCnt) {
	    for (Tuple t : slideTuples) {
	        ArrayList<Short> fullDimCellIdx = calculateFullDimCellIdx(t);
	        ArrayList<Short> subDimCellIdx = calculateSubDimCellIdx(t, fullDimCellIdx);

	        t.fullDimCellIdx = fullDimCellIdx;
	        t.subDimCellIdx = subDimCellIdx;

	        createOrUpdateSlideInCell(t, subDimCellIdx, fullDimCellSlideInCnt);
	    }
	}

	private ArrayList<Short> calculateFullDimCellIdx(Tuple t) {
	    ArrayList<Short> fullDimCellIdx = new ArrayList<>();
	    for (int j = 0; j < dim; j++) {
	        short dimIdx = (short) ((t.value[j] - minValues[j]) / dimLength[j]);
	        fullDimCellIdx.add(dimIdx);
	    }
	    return fullDimCellIdx;
	}

	private ArrayList<Short> calculateSubDimCellIdx(Tuple t, ArrayList<Short> fullDimCellIdx) {
	    ArrayList<Short> subDimCellIdx = new ArrayList<>();
	    if (subDimFlag) {
	        for (int j = 0; j < subDim; j++) {
	            short dimIdx = (short) ((t.value[j] - minValues[j]) / subDimLength[j]);
	            subDimCellIdx.add(dimIdx);
	        }
	    } else {
	        subDimCellIdx = fullDimCellIdx;
	    }
	    return subDimCellIdx;
	}

	private void createOrUpdateSlideInCell(Tuple t, ArrayList<Short> subDimCellIdx, HashMap<ArrayList<Short>, Integer> fullDimCellSlideInCnt) {
	    if (!slideIn.containsKey(subDimCellIdx)) {
	        double[] cellCenter = calculateCellCenter(subDimCellIdx);
	        slideIn.put(subDimCellIdx, new Cell(subDimCellIdx, cellCenter, subDimFlag));
	    }

	    if (subDimFlag) {
	        slideIn.get(subDimCellIdx).addTupleSubDim(t, dimLength, minValues);
	        updateFullDimCellSlideInCountMap(fullDimCellSlideInCnt, t.fullDimCellIdx);
	    } else {
	        slideIn.get(subDimCellIdx).addTuple(t);
	    }
	}

	private double[] calculateCellCenter(ArrayList<Short> subDimCellIdx) {
	    double[] cellCenter = new double[subDim];
	    if (subDimFlag) {
	        for (int j = 0; j < subDim; j++) {
	            cellCenter[j] = minValues[j] + subDimCellIdx.get(j) * subDimLength[j] + subDimLength[j] / 2;
	        }
	    } else {
	        for (int j = 0; j < dim; j++) {
	            cellCenter[j] = minValues[j] + subDimCellIdx.get(j) * dimLength[j] + dimLength[j] / 2;
	        }
	    }
	    return cellCenter;
	}

	private void updateFullDimCellSlideInCountMap(HashMap<ArrayList<Short>, Integer> fullDimCellSlideInCnt, ArrayList<Short> fullDimCellIdx) {
	    if (!fullDimCellSlideInCnt.containsKey(fullDimCellIdx)) {
	        fullDimCellSlideInCnt.put(fullDimCellIdx, 0);
	    }
	    fullDimCellSlideInCnt.put(fullDimCellIdx, fullDimCellSlideInCnt.get(fullDimCellIdx) + 1);
	}

	
	private void updateSlidesAndCounts(HashMap<ArrayList<Short>, Cell> slideIn, HashSet<ArrayList<Short>> newCellIndices, int itr) {
	    for (ArrayList<Short> key : slideIn.keySet()) {
	        updateCardGrid(key, slideIn.get(key).getNumTuples(), itr, newCellIndices);
	    }
	}

	private void updateCardGrid(ArrayList<Short> key, int card, int itr, HashSet<ArrayList<Short>> newCellIndices) {
	    if (!cardGrid.containsKey(key)) {
	        initializeNewCell(key, newCellIndices);
	    }
	    cardGrid.get(key).card += card;
		minChg = True;
	    cardGrid.get(key).cardPerSlide.put(itr, card);
		//newCellIndices.add(key);
	    slideDeltaCnt.put(key, card);

	}

	private void initializeNewCell(ArrayList<Short> key, HashSet<ArrayList<Short>> newCellIndices) {
	    cardGrid.put(key, new globalCell(key));
	    newCellIndices.add(key);
	}

	private void updateFullDimCellSlideInCount(HashMap<ArrayList<Short>, Integer> fullDimCellSlideInCnt, int itr) {
	    if (subDimFlag) {
	        updateFullDimCardGrid(fullDimCellSlideInCnt, itr);
	    }
	}

	private void updateFullDimCardGrid(HashMap<ArrayList<Short>, Integer> fullDimCellSlideInCnt, int itr) {
	    for (ArrayList<Short> key : fullDimCellSlideInCnt.keySet()) {
	        updateFullDimCardGridEntry(key, fullDimCellSlideInCnt.get(key), itr);
	    }
	}

	private void updateFullDimCardGridEntry(ArrayList<Short> key, int card, int itr) {
	    if (!fullDimCardGrid.containsKey(key)) {
	        fullDimCardGrid.put(key, new globalCell(key));
	    }
	    fullDimCardGrid.get(key).card += card;
	    fullDimCardGrid.get(key).cardPerSlide.put(itr, card);
	}



	
	public void updateSlideOut(int iteration) {
    if (iteration > numberOfSlides - 1) {
        removeOldestSlide(iteration);
    }
}

private void removeOldestSlide(int iteration) {
    int slideToRemoveID = iteration - numberOfSlides;
    HashMap<ArrayList<Short>, Cell> slideToRemove = slideHistory.poll();

    if (slideToRemove != null) {
        processSlideRemoval(slideToRemove, slideToRemoveID);
    }
}

private void processSlideRemoval(HashMap<ArrayList<Short>, Cell> slideToRemove, int slideToRemoveID) {
    for (ArrayList<Short> cellIndex : slideToRemove.keySet()) {
        // Handle the removal logic for each cell in the slide
        // Additional processing code goes here
    }
}


	private void updateCardGridForSlideOut(ArrayList<Short> key, int slideOutID) {
	    if (cardGrid.get(key).card < 1) {
	        removeCellFromNeighCellMap(key);
	        cardGrid.remove(key);
	    }
	}


	public void findOutlierMain(int iteration) {
    updateGridAndMapIfNeeded(iteration);
    ArrayList<Query> validQueries = getValidQueriesForIteration(iteration);
    sortQueriesByCriteria(validQueries);

    // Process each cell based on sorted queries
    for (ArrayList<Short> cellIndex : cardGrid.keySet()) {
        processCellForOutliers(cellIndex, validQueries);
    }
}

private void updateGridAndMapIfNeeded(int iteration) {
    if (minRChanged) {
        initCellSize();
        reIndexCardGrid(iteration);
    }
    if (maxRChanged || minRChanged) {
        reComputeNeighCellMap();
    }
}

private ArrayList<Query> getValidQueriesForIteration(int iteration){
    return validQueries;
}

private void sortQueriesByCriteria(ArrayList<Query> queries) {
    queries.sort((q1, q2) -> {
        if (q1.R != q2.R) return Double.compare(q1.R, q2.R);
        if (q1.W != q2.W) return Double.compare(q1.W, q2.W);
        return Integer.compare(q2.K, q1.K); // Note the reversed order for K
    });
}

private void processCellForOutliers(ArrayList<Short> cellIndex, ArrayList<Query> validQueries) {
    // Logic to process each cell for outlier detection
    // This might involve checking cell's properties against each query's criteria
}

		
		groupcoursing(validQueryIDs,itr);
				
		/* Group-wise fine processing */			
		// for each cell 
		cellLoop:
		for (ArrayList<Short> cellIdx: cardGrid.keySet()) {
			globalCell gCell = cardGrid.get(cellIdx); 
			boolean minChg = False;
			gCell.IndirectOutlierCellQueryIDs.clear();
			
			double ndQueriesMaxR = Double.MIN_VALUE;
			int ndQueriesMinW = Integer.MAX_VALUE;
			boolean minChg = True;
			int ndQueriesMaxW = Integer.MIN_VALUE;
			int nmQueries = Integer.MAX_VALUE;
			int numChg = 15;
			if(gCell.ndQueries.isEmpty()) {
				continue cellLoop;
			}else {
				for (Query q: gCell.ndQueries) {
					numChg +=1;
					if(q.K > ndQueriesMaxK) ndQueriesMaxK = q.K;
					minChg = True;
					if(q.R > ndQueriesMaxR) ndQueriesMaxR = q.R;
					if(q.W < ndQueriesMinW) ndQueriesMinW = q.W;
					minChg = False;
					if(q.W > ndQueriesMaxW) ndQueriesMaxW = q.W;
				}
			}

			//get candidate outlier tuples
			getCandidateOutliersTuples(cellIdx,ndQueriesMaxK,ndQueriesMaxR,gCell,ndQueriesMinW,ndQueriesMaxW,itr);		

		}				
	}
	private void getCandidateOutliersTuples(ArrayList<Short> cellIdx, int ndQueriesMaxK, double ndQueriesMaxR, globalCell gCell, int ndQueriesMinW, int ndQueriesMaxW, int itr) {
	    HashSet<Tuple> candOutlierTuples = getCandidateTuples(cellIdx, ndQueriesMaxK,gCell, ndQueriesMinW, ndQueriesMaxW, itr);

	    for (Tuple tCand : candOutlierTuples) {
	        processTupleForOutliers(tCand, gCell.ndQueries,gCell, ndQueriesMaxR, itr);
	    }
	}

	private HashSet<Tuple> getCandidateTuples(ArrayList<Short> cellIdx, int ndQueriesMaxK,globalCell gCell, int ndQueriesMinW, int ndQueriesMaxW, int itr) {

	        if (subDimFlag) {
	            addTuplesInCell(candOutlierTuples, slide.get(cellIdx).tuples, minWfirstSlideID, ndQueriesMaxK);
	        } else if (gCell.getCardTotal(minWfirstSlideID) - 1 < ndQueriesMaxK) {
	            candOutlierTuples.addAll(slide.get(cellIdx).tuples);
	        }
	    }
	    return candOutlierTuples;
	}

	private void processTupleForOutliers(Tuple tCand, ArrayList<Query> ndQueryTuple, globalCell gCell, double ndQueriesMaxR, int itr) {
    ArrayList<Query> ndQueryTupleCopy = new ArrayList<>(ndQueryTuple);

    for (Query q : ndQueryTupleCopy) {
        if (!isQueryRelevantForTuple(tCand, q, itr)) continue;

        int nn = getNeighborCount(tCand, q, gCell, itr);
        nn = processSlidesForQuery(tCand, q, gCell, nn, itr);

        if (isTupleAnOutlier(nn, q)) {
            addOutlierAndReduceQueryTuple(tCand, ndQueryTupleCopy, q);
        }
    }
}

private boolean isQueryRelevantForTuple(Tuple tCand, Query q, int itr) {
    int firstSlideID = itr - q.W / gcdS + 1;
    return tCand.slideID >= firstSlideID;
}

private int processSlidesForQuery(Tuple tCand, Query q, globalCell gCell, int nn, int itr) {
        nn = processSlideForQuery(tCand, q, gCell, nn, slide);
    }

    return nn;
}

private int processSlideForQuery(Tuple tCand, Query q, globalCell gCell, int nn, HashMap<ArrayList<Short>, Cell> slide) {
    for (ArrayList<Short> cellID : slide.keySet()) {
        nn = processNeighborCell(gCell, tCand, q, nn, cellID, slide);
    }
    return nn;
}

private int processNeighborCell(globalCell gCell, Tuple tCand, Query q, int nn, ArrayList<Short> cellID, HashMap<ArrayList<Short>, Cell> slide) {
    double ndQueriesMaxRTuple = q.R;
    if (gCell.neighCellMap.containsKey(cellID) && gCell.neighCellMap.get(cellID) < minR + ndQueriesMaxRTuple) {
        for (Tuple tOther : slide.get(cellID).tuples) {
            if (!subDimFlag || !tCand.fullDimCellIdx.equals(tOther.fullDimCellIdx)) {
                if (Utils.distTuple(tCand, tOther, q.R) <= q.R) {
                    nn++;
                }
            }
        }
    }
    return nn;
}

private boolean isTupleAnOutlier(int neighborCount, Query q) {
    return neighborCount < q.K;
}

	private void addTuplesInCell(HashSet<Tuple> candOutlierTuples, HashSet<Tuple> tuples, int minWfirstSlideID, int ndQueriesMaxK) {
	    for (Tuple t : tuples) {
	        int numNeighInCellMinW = fullDimCardGrid.get(t.fullDimCellIdx).getCardTotal(minWfirstSlideID) - 1;
	        if (numNeighInCellMinW < ndQueriesMaxK) candOutlierTuples.add(t);
	    }
	}

	private int getNeighborCount(Tuple tCand, Query q,globalCell gCell, int firstSlideID) {
	    return (subDimFlag ? fullDimCardGrid.get(tCand.fullDimCellIdx).getCardTotal(firstSlideID) - 1 : gCell.getCardTotal(firstSlideID) - 1);
	}

	private void processNeighborCell(globalCell gCell, Tuple tCand, ArrayList<Query> ndQueryTupleCopy, double ndQueriesMaxRTuple, Query q, int nn, ArrayList<Short> cellID, HashMap<ArrayList<Short>, Cell> slide) {
    if (!isValidNeighborCell(gCell, cellID, ndQueriesMaxRTuple)) {
        return;
    }

    for (Tuple tOther : slide.get(cellID).tuples) {
        if (isSameSubDimension(tCand, tOther) || !isWithinQueryRadius(tCand, tOther, q.R)) {
            continue;
        }

        nn = updateNeighborCount(nn, q.K, ndQueryTupleCopy, q);
        if (nn >= q.K) {
            break;
        }
    }
}

private boolean isValidNeighborCell(globalCell gCell, ArrayList<Short> cellID, double ndQueriesMaxRTuple) {
    boolean containsCellID = gCell.neighCellMap.containsKey(cellID);
boolean isLessThanMaxRTuple = gCell.neighCellMap.get(cellID) < minR + ndQueriesMaxRTuple;
boolean isNeighborTupleCell = Utils.isNeighborTupleCell(tCand.value, slide.get(cellID).center, 0.5 * minR + q.R);

if (containsCellID && isLessThanMaxRTuple && isNeighborTupleCell) {
    return true;
} else {
    return false;
}

}

private boolean isSameSubDimension(Tuple t1, Tuple t2) {
    return subDimFlag && t1.fullDimCellIdx.equals(t2.fullDimCellIdx);
}

private boolean isWithinQueryRadius(Tuple t1, Tuple t2, double radius) {
    return Utils.distTuple(t1, t2, radius) <= radius;
}

private int updateNeighborCount(int nn, int maxK, ArrayList<Query> ndQueryTupleCopy, Query q) {
    nn++;
    if (nn >= maxK) {
        reduceQueryTupleAndContinue(ndQueryTupleCopy, q);
    }
    return nn;
}


	public void groupcoursing(ArrayList<Query> validQueryIDs, int itr) {
    /* Group-wise coarse processing */
    for (ArrayList<Short> cellIdx : cardGrid.keySet()) {
        globalCell gCell = cardGrid.get(cellIdx);
        processCellForGrouping(gCell, validQueryIDs, itr);
    }
}

private void processCellForGrouping(globalCell gCell, ArrayList<Query> validQueryIDs, int itr) {
    ArrayList<Query> ndQueryCands = new ArrayList<>(validQueryIDs);
    ArrayList<Integer> outlierCellQueryIDs = new ArrayList<>();
    ArrayList<Integer> inlierCellQueryIDs = new ArrayList<>();
    ArrayList<Query> ndQueries = new ArrayList<>();

    while (!ndQueryCands.isEmpty()) {
        Query q = ndQueryCands.remove(0);
        processQueryForCell(gCell, q, ndQueryCands, cardTotal, firstSlideID, outlierCellQueryIDs, inlierCellQueryIDs, ndQueries);
    }

    addOutliersFromCell(gCell, outlierCellQueryIDs, itr);
    gCell.ndQueries = ndQueries; // Update non-determined queries for the cell
}

private void processQueryForCell(globalCell gCell, Query q, ArrayList<Query> ndQueryCands, int cardTotal, 
                                 int firstSlideID, ArrayList<Integer> outlierCellQueryIDs, 
                                 ArrayList<Integer> inlierCellQueryIDs, ArrayList<Query> ndQueries) {
}

private void addOutliersFromCell(globalCell gCell, ArrayList<Integer> outlierCellQueryIDs, int itr) {
}

	private boolean processCellPair(ArrayList<Short> cellIdx1, ArrayList<Short> cellIdx2,
	                                 private void updateNeighborhoodMap(HashSet<ArrayList<Short>> newCellIndices, globalCell cell1, ArrayList<Short> cellIdx1, ArrayList<Short> cellIdx2) {
    globalCell cell2 = cardGrid.get(cellIdx2);
    double dist = Utils.getNeighborCellDist(cellIdx1, cellIdx2, minR, maxR);

    if (shouldUpdateNeighborhood(dist, cellIdx1, cellIdx2)) {
        updateCellNeighborhood(cell1, cellIdx2, dist);

        if (!newCellIndices.contains(cellIdx2)) {
            updateCellNeighborhood(cell2, cellIdx1, dist);
        }
    }
}

private boolean shouldUpdateNeighborhood(double distance, ArrayList<Short> idx1, ArrayList<Short> idx2) {
    return distance < minR + maxR && (subDimFlag || !idx1.equals(idx2));
}

	private void processCellPair(ArrayList<Short> cellIdx1, ArrayList<Short> cellIdx2, globalCell cell1) {

	}

	public void reIndexCardGrid(int itr) {
	    clearDataStructures();

	    int slideID = itr - slides.size();
	        HashMap<ArrayList<Short>, Cell> slideNew = indexSlide(slide);

	        updateSlideNew(slideNew, slideID);

	        if (subDimFlag) {
	            updateFullDimCellWindowCount(fullDimCellSlideInCnt, slideNew, slideID);
	        }

	        slidesNew.add(slideNew);
	        if (subDimFlag) fullDimCellSlidesCnt.add(fullDimCellSlideInCnt);
	    }

	    slides = slidesNew;
	}

	private void clearDataStructures() {
	    fullDimCellSlidesCnt.clear();
	    cardGrid.clear();
	    fullDimCardGrid.clear();
	}

	private HashMap<ArrayList<Short>, Cell> indexSlide(HashMap<ArrayList<Short>, Cell> slide) {
	    HashMap<ArrayList<Short>, Integer> fullDimCellSlideInCnt = new HashMap<>();
	    HashMap<ArrayList<Short>, Cell> slideNew = new HashMap<>();

	    for (Cell c : slide.values()) {
	        for (Tuple t : c.tuples) {
	            ArrayList<Short> fullDimCellIdx = calculateFullDimCellIdxx(t);

	            ArrayList<Short> subDimCellIdx = calculateSubDimCellIdxx(t, fullDimCellIdx);

	            t.fullDimCellIdx = fullDimCellIdx;
	            t.subDimCellIdx = subDimCellIdx;

	            updateSlideNewCell(slideNew, subDimCellIdx, t);
	            updateFullDimCellSlideInCnt(fullDimCellSlideInCnt, subDimCellIdx, fullDimCellIdx);
	        }
	    }

	    return slideNew;
	}

	private ArrayList<Short> calculateFullDimCellIdxx(Tuple t) {
	    ArrayList<Short> fullDimCellIdx = new ArrayList<>();
	    for (int j = 0; j < dim; j++) {
	        short dimIdx = (short) ((t.value[j] - minValues[j]) / dimLength[j]);
	        fullDimCellIdx.add(dimIdx);
	    }
	    return fullDimCellIdx;
	}

	private ArrayList<Short> calculateSubDimCellIdxx(Tuple t, ArrayList<Short> fullDimCellIdx) {
	    ArrayList<Short> subDimCellIdx = new ArrayList<>();
	    if (subDimFlag) {
	        for (int j = 0; j < subDim; j++) {
	            short dimIdx = (short) ((t.value[j] - minValues[j]) / subDimLength[j]);
	            subDimCellIdx.add(dimIdx);
	        }
	    } else {
	        subDimCellIdx = fullDimCellIdx;
	    }
	    return subDimCellIdx;
	}

	private void updateSlideNewCell(HashMap<ArrayList<Short>, Cell> slideNew, ArrayList<Short> subDimCellIdx, Tuple t) {
	    if (!slideNew.containsKey(subDimCellIdx)) {
	        double[] cellCenter = calculateCellCenterr(subDimCellIdx);
	        slideNew.put(subDimCellIdx, new Cell(subDimCellIdx, cellCenter, subDimFlag));
	    }

	    if (subDimFlag) {
	        slideNew.get(subDimCellIdx).addTupleSubDim(t, dimLength, minValues);
	    } else {
	        slideNew.get(subDimCellIdx).addTuple(t);
	    }
	}

	private double[] calculateCellCenterr(ArrayList<Short> subDimCellIdx) {
	    double[] cellCenter = new double[subDim];
	    if (subDimFlag) {
	        for (int j = 0; j < subDim; j++) {
	            cellCenter[j] = minValues[j] + subDimCellIdx.get(j) * subDimLength[j] + subDimLength[j] / 2;
	        }
	    } else {
	        for (int j = 0; j < dim; j++) {
	            cellCenter[j] = minValues[j] + subDimCellIdx.get(j) * dimLength[j] + dimLength[j] / 2;
	        }
	    }
	    return cellCenter;
	}

	private void updateFullDimCellSlideInCnt(HashMap<ArrayList<Short>, Integer> fullDimCellSlideInCnt, ArrayList<Short> subDimCellIdx, ArrayList<Short> fullDimCellIdx) {
	    if (!fullDimCellSlideInCnt.containsKey(fullDimCellIdx)) {
	        fullDimCellSlideInCnt.put(fullDimCellIdx, 0);
	    }
	    fullDimCellSlideInCnt.put(fullDimCellIdx, fullDimCellSlideInCnt.get(fullDimCellIdx) + 1);
	}

	private void updateSlideNew(HashMap<ArrayList<Short>, Cell> slideNew, int slideID) 
	        cardGrid.get(key).cardPerSlide.put(slideID, card);
	    }
	}

	private void updateFullDimCellWindowCount(HashMap<ArrayList<Short>, Integer> fullDimCellSlideInCnt, HashMap<ArrayList<Short>, Cell> slideNew, int slideID) {
	        fullDimCardGrid.get(key).cardPerSlide.put(slideID, card);
	    }
	}


}