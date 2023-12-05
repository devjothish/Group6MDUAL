package MDUAL;
import java.util.*;

public class MDUAL {
	public double minR, maxR, minR_old, maxR_old;
	public int gcdS;
	public int dim, subDim;
	public boolean subDimFlag;
	public int nS, nW;
	public double[] minValues;
	public double[] dimLength, subDimLength;
	
	public HashMap<ArrayList<Short>,Cell> slideIn, slideOut;
	public HashMap<ArrayList<Short>,Integer> slideDeltaCnt;
	public HashMap<ArrayList<Short>,globalCell> cardGrid, fullDimCardGrid;
	public LinkedList<HashMap<ArrayList<Short>,Cell>> slides; 
	public LinkedList<HashMap<ArrayList<Short>,Integer>> fullDimCellSlidesCnt; 
	public HashSet<Tuple> outliers;

	public HashMap<Integer,Query> querySet;
	
	public boolean maxRChanged, minRChanged;
	
	public MDUAL(int dim, int subDim, int nS, int gcdS, double[] minValues) {
		this.dim = dim;
		this.subDim = subDim;
		this.subDimFlag = dim != subDim;
		this.minR = Double.MAX_VALUE;
		this.maxR = Double.MIN_VALUE;
		this.gcdS = gcdS;
		this.nS = nS;
		this.minValues = minValues;
		
		this.cardGrid = new HashMap<ArrayList<Short>,globalCell>();
		this.fullDimCardGrid = new HashMap<ArrayList<Short>,globalCell>();
		this.fullDimCellSlidesCnt = new LinkedList<HashMap<ArrayList<Short>,Integer>>(); 
		this.slides = new LinkedList<HashMap<ArrayList<Short>,Cell>>();
		this.slideOut = new HashMap<ArrayList<Short>,Cell>();
		this.outliers = new HashSet<Tuple>();
	}

	public HashSet<Tuple> findOutlier(ArrayList<Tuple> newSlideTuples, HashMap<Integer,Query> newQuerySet, int itr) {
		this.querySet = newQuerySet;
		
		if(itr==0) {
			this.updateBasisParams(itr); 
			this.initCellSize();
		}
		this.clearPreviousOutliers();
		this.updateWindow(newSlideTuples,itr);
		this.updateBasisParams(itr); 
		this.findOutlierMain(itr);

		return this.outliers;
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
	    minRChanged = false;
	    minR_old = minR;
	    maxR_old = maxR;
	    minR = Double.MAX_VALUE;
	    maxR = Double.MIN_VALUE;
	}

	private void updateMinMaxValues() {
	    for (Query q : querySet.values()) {
	        maxR = Math.max(maxR, q.R);
	        minR = Math.min(minR, q.R);
	    }
	}

	private void checkMinMaxChanges(int itr) {
	    if (itr > 0 && minR != minR_old) {
	        minRChanged = true;
	    }
	    if (itr > 0 && maxR != maxR_old) {
	        maxRChanged = true;
	    }
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
	    cardGrid.get(key).cardPerSlide.put(itr, card);
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



	
	private void updateSlideOut(int itr) {
	    if (itr > nS - 1) {
	        int slideOutID = itr - nS;
	        slideOut = slides.poll();

	        for (ArrayList<Short> key : slideOut.keySet()) {
	            updateCardGridForSlideOut(key, slideOutID);
	            updateSlideDeltaCount(key, slideOut.get(key).getNumTuples());
	        }

	        if (subDimFlag) {
	            updateFullDimCellSlideOutCount(slideOutID);
	        }
	    }
	}

	private void updateCardGridForSlideOut(ArrayList<Short> key, int slideOutID) {
	    int card = slideOut.get(key).getNumTuples();
	    cardGrid.get(key).card -= card;
	    cardGrid.get(key).cardPerSlide.remove(slideOutID);

	    if (cardGrid.get(key).card < 1) {
	        removeCellFromNeighCellMap(key);
	        cardGrid.remove(key);
	    }
	}

	private void removeCellFromNeighCellMap(ArrayList<Short> key) {
	    for (ArrayList<Short> neighCellIdx : cardGrid.get(key).neighCellMap.keySet()) {
	        if (key.equals(neighCellIdx)) continue;
	        cardGrid.get(neighCellIdx).neighCellMap.remove(key);
	    }
	}

	private void updateSlideDeltaCount(ArrayList<Short> key, int card) {
	    if (slideDeltaCnt.containsKey(key)) {
	        slideDeltaCnt.put(key, slideDeltaCnt.get(key) - card);
	    } else {
	        slideDeltaCnt.put(key, card * -1);
	    }
	}

	private void updateFullDimCellSlideOutCount(int slideOutID) {
	    HashMap<ArrayList<Short>, Integer> fullDimCellSlideOutCnt = fullDimCellSlidesCnt.poll();
	    for (ArrayList<Short> key : fullDimCellSlideOutCnt.keySet()) {
	        updateFullDimCardGridForSlideOut(key, fullDimCellSlideOutCnt.get(key), slideOutID);
	    }
	}

	private void updateFullDimCardGridForSlideOut(ArrayList<Short> key, int card, int slideOutID) {
	    fullDimCardGrid.get(key).card -= card;
	    fullDimCardGrid.get(key).cardPerSlide.remove(slideOutID);

	    if (fullDimCardGrid.get(key).card < 1) {
	        fullDimCardGrid.remove(key);
	    }
	}


	public void findOutlierMain(int itr) {
		if(minRChanged) {
			initCellSize();
			reIndexCardGrid(itr);
		}
		if(maxRChanged || minRChanged) this.reComputeNeighCellMap();
		
		ArrayList<Query> validQueryIDs = new ArrayList<Query>();
		for(Query q: querySet.values()) {
			if((itr+1) % (q.S/gcdS) == 0) validQueryIDs.add(q); //check if slide condition is met
		}

		//inlier-first
		Collections.sort(validQueryIDs, new Comparator<Query>(){ //Sort by order of smaller W -> smaller R -> larger K 
			@Override
			public int compare(Query q1, Query q2) {
				if(q1.R>q2.R) {
					return 1;
				}else if(q1.R==q2.R){
					if(q1.W>q2.W) return 1;
					else if(q1.W==q2.W){
						if(q1.K<=q2.K) return 1;
						else return -1;
					}else return -1;
				}else return -1;
			}
		});
		
		groupcoursing(validQueryIDs,itr);
				
		/* Group-wise fine processing */			
		// for each cell 
		cellLoop:
		for (ArrayList<Short> cellIdx: cardGrid.keySet()) {
			globalCell gCell = cardGrid.get(cellIdx); 
			gCell.IndirectOutlierCellQueryIDs.clear();
			
			int ndQueriesMinW = Integer.MAX_VALUE;
			int ndQueriesMaxW = Integer.MIN_VALUE;
			int ndQueriesMaxK = Integer.MIN_VALUE;
			double ndQueriesMaxR = Double.MIN_VALUE;
			if(gCell.ndQueries.isEmpty()) {
				continue cellLoop;
			}else {
				for (Query q: gCell.ndQueries) {
					if(q.K > ndQueriesMaxK) ndQueriesMaxK = q.K;
					if(q.R > ndQueriesMaxR) ndQueriesMaxR = q.R;
					if(q.W < ndQueriesMinW) ndQueriesMinW = q.W;
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
	    HashSet<Tuple> candOutlierTuples = new HashSet<>();
	    int minWfirstSlideID = itr - ndQueriesMinW / gcdS + 1;
	    int maxWfirstSlideID = itr - ndQueriesMaxW / gcdS + 1;
	    int slideID = itr - slides.size();

	    for (HashMap<ArrayList<Short>, Cell> slide : slides) {
	        slideID++;
	        if (slideID < maxWfirstSlideID || !slide.containsKey(cellIdx)) continue;

	        if (subDimFlag) {
	            addTuplesInCell(candOutlierTuples, slide.get(cellIdx).tuples, minWfirstSlideID, ndQueriesMaxK);
	        } else if (gCell.getCardTotal(minWfirstSlideID) - 1 < ndQueriesMaxK) {
	            candOutlierTuples.addAll(slide.get(cellIdx).tuples);
	        }
	    }
	    return candOutlierTuples;
	}

	private void processTupleForOutliers(Tuple tCand, ArrayList<Query> ndQueryTuple,globalCell gCell, double ndQueriesMaxR, int itr) {
	    ArrayList<Query> ndQueryTupleCopy = new ArrayList<>(ndQueryTuple);
	    double ndQueriesMaxRTuple = ndQueriesMaxR;

	    queryLoop:
	    while (!ndQueryTupleCopy.isEmpty()) {
	        Query q = ndQueryTupleCopy.iterator().next();
	        ndQueryTupleCopy.remove(q);

	        int firstSlideID = itr - q.W / gcdS + 1;
	        if (tCand.slideID < firstSlideID) continue queryLoop;

	        int nn = getNeighborCount(tCand, q, gCell, firstSlideID);
	        int slideID = itr - slides.size();

	        for (HashMap<ArrayList<Short>, Cell> slide : slides) {
	            slideID++;
	            if (slideID < firstSlideID) continue;
	            for (ArrayList<Short> cellID : slide.keySet()) {
	                processNeighborCell(gCell, tCand, ndQueryTupleCopy, ndQueriesMaxRTuple, q, nn, cellID,slide);
	            }
	        }

	        if (nn < q.K) {
	            addOutlierAndReduceQueryTuple(tCand, ndQueryTupleCopy, q);
	        }
	    }
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

	private void processNeighborCell(globalCell gCell, Tuple tCand, ArrayList<Query> ndQueryTupleCopy, double ndQueriesMaxRTuple, Query q, int nn, ArrayList<Short> cellID,HashMap<ArrayList<Short>, Cell> slide) {
	    if (gCell.neighCellMap.containsKey(cellID)
	            && gCell.neighCellMap.get(cellID) < minR + ndQueriesMaxRTuple
	            && Utils.isNeighborTupleCell(tCand.value, slide.get(cellID).center, 0.5 * minR + q.R)) {
	        for (Tuple tOther : slide.get(cellID).tuples) {
	            if (subDimFlag && tCand.fullDimCellIdx.equals(tOther.fullDimCellIdx)) continue;
	            if (Utils.distTuple(tCand, tOther, q.R) <= q.R) {
	                nn++;
	                if (nn >= q.K) {
	                    reduceQueryTupleAndContinue(ndQueryTupleCopy, q);
	                }
	            }
	        }
	    }
	}

	private void addOutlierAndReduceQueryTuple(Tuple tCand, ArrayList<Query> ndQueryTupleCopy, Query q) {
	    tCand.outlierQueryIDs.add(q.id);
	    outliers.add(tCand);

	    Iterator<Query> ndQueryItr = ndQueryTupleCopy.iterator();
	    double ndQueriesMaxRTuple = Double.MIN_VALUE;
	    while (ndQueryItr.hasNext()) {
	        Query q2 = ndQueryItr.next();
	        if (q2.W <= q.W && q2.R <= q.R && q2.K >= q.K) {
	            ndQueryItr.remove();
	            tCand.outlierQueryIDs.add(q2.id);
	        } else if (q2.R > ndQueriesMaxRTuple) {
	            ndQueriesMaxRTuple = q2.R;
	        }
	    }
	}

	private void reduceQueryTupleAndContinue(ArrayList<Query> ndQueryTupleCopy, Query q) {
	    Iterator<Query> ndQueryTupleItr = ndQueryTupleCopy.iterator();
	    double ndQueriesMaxRTuple = Double.MIN_VALUE;
	    while (ndQueryTupleItr.hasNext()) {
	        Query q2 = ndQueryTupleItr.next();
	        if (q2.W >= q.W && q2.R >= q.R && q2.K <= q.K) {
	            ndQueryTupleItr.remove();
	        } else if (q2.R > ndQueriesMaxRTuple) {
	            ndQueriesMaxRTuple = q2.R;
	        }
	    }
	}

	public void groupcoursing(ArrayList<Query> validQueryIDs,int itr) {
		/* Group-wise coarse processing */ 
		for (ArrayList<Short> cellIdx: cardGrid.keySet()) {
			globalCell gCell = cardGrid.get(cellIdx); 
			ArrayList<Query> ndQueryCands = new ArrayList<Query>(validQueryIDs);
			
			ArrayList<Integer> outlierCellQueryIDs = new ArrayList<Integer>();
			ArrayList<Integer> inlierCellQueryIDs = new ArrayList<Integer>();
			ArrayList<Query> ndQueries = new ArrayList<Query>();
			//Verify outlier/intlier cell query IDs
			while(!ndQueryCands.isEmpty()) {
				Query q = ndQueryCands.iterator().next();
				ndQueryCands.remove(q);
				if(gCell.IndirectOutlierCellQueryIDs.contains(q.id)) {
					outlierCellQueryIDs.add(q.id);
					continue;
				}
				int firstSlideID = itr - q.W/gcdS + 1;
				int cardTotal = gCell.getCardTotal(firstSlideID);
				
				if(!subDimFlag && cardTotal > q.K){ //inlier cell for q
					inlierCellQueryIDs.add(q.id);
					Iterator<Query> ndQueryItr = ndQueryCands.iterator();
					while(ndQueryItr.hasNext()) { //propagation
						Query q2 = ndQueryItr.next();
						if(q2.W >= q.W && q2.R >= q.R && q2.K <= q.K) {
							ndQueryItr.remove();
						}
					}
				}else{ 
					int thredNeighCellCardTotal = 0;
					for(ArrayList<Short> neighCellIdx: gCell.getThredNeighCellsIn(q.R-minR)) thredNeighCellCardTotal += cardGrid.get(neighCellIdx).getCardTotal(firstSlideID);					
					if(!subDimFlag && thredNeighCellCardTotal + cardTotal > q.K){ //inlier cell for q
						inlierCellQueryIDs.add(q.id); 
						Iterator<Query> ndQueryItr = ndQueryCands.iterator();
						while(ndQueryItr.hasNext()) { //propagation
							Query q2 = ndQueryItr.next();
							if(q2.W >= q.W && q2.R >= q.R && q2.K <= q.K) {
								ndQueryItr.remove();
							}
						}
					}else {
						//Get total cards of neighbor cells. 
						//if not sub dim, add cell card since it is not contained in neighbor cells
						int neighCellCardTotal = (subDimFlag? 0: cardTotal);
						for(ArrayList<Short> neighCellIdx: gCell.getThredNeighCellsOut(q.R+minR)) neighCellCardTotal += cardGrid.get(neighCellIdx).getCardTotal(firstSlideID);

						if(neighCellCardTotal <= q.K) { //outlier cell for q
							outlierCellQueryIDs.add(q.id);
							Iterator<Query> ndQueryItr = ndQueryCands.iterator();
							while(ndQueryItr.hasNext()) { //propagation
								Query q2 = ndQueryItr.next();
								if(q2.W <= q.W && q2.R <= q.R && q2.K >= q.K) {
									ndQueryItr.remove();
									outlierCellQueryIDs.add(q2.id);
								}
							}
						}else{
							ndQueries.add(q);
						}
											
						
					}
				}
			}
			
			//Add outlier tuples by outlierCellQueryIDs
			for (int qid: outlierCellQueryIDs) {
				int firstSlideID = itr - querySet.get(qid).W/gcdS + 1;
				int slideID = itr - slides.size();
				for(HashMap<ArrayList<Short>, Cell> slide: slides) {
					slideID++;
					if(slideID < firstSlideID || !slide.containsKey(cellIdx)) continue; //check if slide is inside the query window condition OR contains the cell
					for(Tuple t: slide.get(cellIdx).tuples) {
						t.outlierQueryIDs.add(qid);
						outliers.add(t);
					}
				}
			}
			gCell.ndQueries = ndQueries;
		}
	}
	
	public void getNeighCellMap(HashSet<ArrayList<Short>> newCellIndices) {
	    for (ArrayList<Short> newCellIdx : newCellIndices) {
	        globalCell newCell = cardGrid.get(newCellIdx);

	        for (ArrayList<Short> candCellIdx : cardGrid.keySet()) {
	            if (processCellPair(newCellIdx, candCellIdx, newCellIndices, newCell)) {
	                processCellPair(candCellIdx, newCellIdx, newCellIndices, cardGrid.get(candCellIdx));
	            }
	        }
	    }
	}

	private boolean processCellPair(ArrayList<Short> cellIdx1, ArrayList<Short> cellIdx2,
	                                 HashSet<ArrayList<Short>> newCellIndices, globalCell cell1) {
	    globalCell cell2 = cardGrid.get(cellIdx2);
	    double dist = Utils.getNeighborCellDist(cellIdx1, cellIdx2, minR, minR + maxR);

	    if (dist < minR + maxR && (subDimFlag || !cellIdx1.equals(cellIdx2))) {
	        cell1.neighCellMap.put(cellIdx2, dist);

	        if (!newCellIndices.contains(cellIdx2)) {
	            cell2.neighCellMap.put(cellIdx1, dist);
	        }
	        return true;
	    }

	    return false;
	}

	
	public void reComputeNeighCellMap() {
	    for (ArrayList<Short> cellIdx : cardGrid.keySet()) {
	        globalCell cell = cardGrid.get(cellIdx);
	        cell.neighCellMap.clear();

	        for (ArrayList<Short> otherCellIdx : cardGrid.keySet()) {
	            if (!cell.neighCellMap.containsKey(otherCellIdx)) {
	                processCellPair(cellIdx, otherCellIdx, cell);
	            }
	        }
	    }
	}

	private void processCellPair(ArrayList<Short> cellIdx1, ArrayList<Short> cellIdx2, globalCell cell1) {
	    double dist = Utils.getNeighborCellDist(cellIdx1, cellIdx2, minR, minR + maxR);

	    if (dist < minR + maxR && (subDimFlag || !cellIdx1.equals(cellIdx2))) {
	        cell1.neighCellMap.put(cellIdx2, dist);
	        cardGrid.get(cellIdx2).neighCellMap.put(cellIdx1, dist);
	    }
	}

	public void reIndexCardGrid(int itr) {
	    clearDataStructures();

	    int slideID = itr - slides.size();
	    LinkedList<HashMap<ArrayList<Short>,Cell>> slidesNew = new LinkedList<HashMap<ArrayList<Short>,Cell>>();
	    for (HashMap<ArrayList<Short>, Cell> slide : slides) {
	        slideID++;
	        HashMap<ArrayList<Short>, Integer> fullDimCellSlideInCnt = new HashMap<>();
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

	private void updateSlideNew(HashMap<ArrayList<Short>, Cell> slideNew, int slideID) {
	    for (ArrayList<Short> key : slideNew.keySet()) {
	        int card = slideNew.get(key).getNumTuples();
	        if (!cardGrid.containsKey(key)) cardGrid.put(key, new globalCell(key));
	        cardGrid.get(key).card += card;
	        cardGrid.get(key).cardPerSlide.put(slideID, card);
	    }
	}

	private void updateFullDimCellWindowCount(HashMap<ArrayList<Short>, Integer> fullDimCellSlideInCnt, HashMap<ArrayList<Short>, Cell> slideNew, int slideID) {
	    for (ArrayList<Short> key : fullDimCellSlideInCnt.keySet()) {
	        int card = fullDimCellSlideInCnt.get(key);
	        if (!fullDimCardGrid.containsKey(key)) {
	            fullDimCardGrid.put(key, new globalCell(key));
	        }
	        fullDimCardGrid.get(key).card += card;
	        fullDimCardGrid.get(key).cardPerSlide.put(slideID, card);
	    }
	}


}