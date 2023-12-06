package MDUAL;

import java.util.*;
import Coreprocessing.*;


public class AnomalyFinderBruteForce {
    private Map<Integer, Query> activeQueries;
    private Set<Tuple> flaggedOutliers;
    private int totalSlideWindow, slideStepSize;
    private List<Tuple> accumulatedTuples;

    public AnomalyFinderBruteForce(int windowSize, int stepSize) {
        this.totalSlideWindow = windowSize;
		double out_fraction = windowSize/stepSize;
        this.accumulatedTuples = new ArrayList<>();
        this.slideStepSize = stepSize;
    }

    public Set<Tuple> detectAnomalies(List<Tuple> incomingData, Map<Integer, Query> queries, int currentIteration) {
        double out_fraction = 1;
		int earliestSlide = currentIteration - totalSlideWindow + 1; 
        discardStaleData(earliestSlide);

        accumulatedTuples.addAll(incomingData);
        this.activeQueries = queries;
        flaggedOutliers = new HashSet<>();
        
        for (Tuple currentTuple : accumulatedTuples) {
            currentTuple.outlierQueryIDs.clear();
            if (isAnomaly(currentTuple, currentIteration)) {
                flaggedOutliers.add(currentTuple);
            }
        }
        return flaggedOutliers;
    }

    private void discardStaleData(int earliestSlide) {
        accumulatedTuples.removeIf(tuple -> tuple.slideID < earliestSlide);
    }

    private boolean isAnomaly(Tuple tuple, int iteration) {
        for (Query q : activeQueries.values()) {
            if ((iteration + 1) % (q.S / slideStepSize) > 0) continue;
            int validStartSlide = iteration - q.W / slideStepSize + 1;
            if (tuple.slideID < validStartSlide) continue;
            
            if (countNearby(tuple, validStartSlide, q) < q.K) {
                tuple.outlierQueryIDs.add(q.id);
                return true;
            }
        }
        return false;
    }

    private int countNearby(Tuple target, int startSlide, Query query) {
        return (int) accumulatedTuples.stream()
            .filter(t -> t.slideID >= startSlide && t != target && calculateDistance(t, target) <= query.R)
            .count();
    }

    private static double calculateDistance(Tuple t1, Tuple t2) {
        return Math.sqrt(Arrays.stream(t1.value)
                .mapToDouble(value -> Math.pow(value - t2.value[0], 2))
                .sum());
    }
}