package MDUAL;

import java.util.*;
import java.util.stream.*;

public class CellDetails {
    private List<Short> cellCoordinates, neighboringCells;
    private Map<List<Short>, Double> neighboringCellsDistance;
    private boolean isRecentlyUpdated;
    private int cellCardinality;
    
    private List<Integer> outlierQueryIDs = new ArrayList<>();
    private Map<Integer, Integer> cardinalityBySlide;
    private List<Query> nonDirectQueries = new ArrayList<>();

    public CellDetails(List<Short> coordinates) {
        this.cellCoordinates = new ArrayList<>(coordinates);
        this.cellCardinality = 0;
        this.neighboringCellsDistance = new HashMap<>();
        this.isRecentlyUpdated = false;
        this.cardinalityBySlide = new HashMap<>();
    }

    public int calculateTotalCardinality(int startingSlideId) {
        return cardinalityBySlide.entrySet().stream()
                .filter(entry -> entry.getKey() >= startingSlideId)
                .mapToInt(Map.Entry::getValue)
                .sum();
    }
    
    public List<List<Short>> fetchNeighborCellsWithinThreshold(double threshold, boolean inclusive) {
        return neighboringCellsDistance.entrySet().stream()
                .filter(entry -> inclusive ? entry.getValue() <= threshold : entry.getValue() < threshold)
                .map(Map.Entry::getKey)
                .collect(Collectors.toList());
    }

    public List<Short> getCellCoordinates() {
        return new ArrayList<>(cellCoordinates);
    }
    public void setCellCoordinates(List<Short> coordinates) {
        this.cellCoordinates = new ArrayList<>(coordinates);
    }
}
