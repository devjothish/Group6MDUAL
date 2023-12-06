package MDUAL;
import java.util.*;

public class GridSegment {
    private int segmentID;
    private List<Short> segmentIdentifier;
    private List<Short> segmentIDi;
    private Map<List<Short>, GridSegment> subSegments;
    private Set<GridDataPoint> dataPointsInSegment;
    

    public int centerD;
    private double[] segmentCenter;

    public GridSegment(List<Short> identifier, double[] segmentSizes, double[] startingPoints) {
        this.dataPointsInSegment = new HashSet<>();
        this.segmentIdentifier = new ArrayList<>(identifier);
        this.segmentCenter = calculateSegmentCenter(segmentSizes, startingPoints);
        this.subSegments = new HashMap<>();
    }

    private double[] calculateSegmentCenter(double[] sizes, double[] starts) {
        double[] centerPoint = new double[sizes.length];
        for (int i = 0; i < sizes.length; i++) {
            centerPoint[i] = starts[i] + (segmentIdentifier.get(i) + 0.5) * sizes[i];
        }
        return centerPoint;
    }

    public int countDataInSegment() {
        return dataPointsInSegment.size();
    }

    public void includeDataPoint(GridDataPoint data, double[] segmentSizes, double[] startingPoints) {
        dataPointsInSegment.add(data);
        List<Short> completeIdentifier = data.retrieveSegmentIdentifier();
        subSegments.computeIfAbsent(completeIdentifier, k -> new GridSegment(completeIdentifier, segmentSizes, startingPoints))
                   .storeDataPointInSegment(data);
    }

    public void storeDataPointInSegment(GridDataPoint data) {
        dataPointsInSegment.add(data);
    }
}

class GridDataPoint {
    private List<Short> segmentIdentifier;
    private Map<String, Double> attributes; // Assuming each data point has some attributes

    public GridDataPoint(List<Short> identifier) {
        this.segmentIdentifier = new ArrayList<>(identifier);
        this.attributes = new HashMap<>();
    }

    public List<Short> retrieveSegmentIdentifier() {
        return new ArrayList<>(segmentIdentifier);
    }

    // Method to add an attribute to the data point
    public void addAttribute(String key, Double value) {
        attributes.put(key, value);
    }

    // Method to retrieve an attribute
    public Double getAttribute(String key) {
        return attributes.getOrDefault(key, null);
    }

    // Optional: Method to get all attributes
    public Map<String, Double> getAllAttributes() {
        return new HashMap<>(attributes);
    }

    // Optional: ToString method for easier debugging
    @Override
    public String toString() {
        return "GridDataPoint{" +
                "segmentIdentifier=" + segmentIdentifier +
                ", attributes=" + attributes +
                '}';
    }
}
