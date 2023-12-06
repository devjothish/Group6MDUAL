package Coreprocessing;

import MDUAL.*;
import java.io.IOException;
import java.util.*;

public class Processor {
    private DataQueryManager datasetManager = new DataQueryManager();
    private DataQueryManager queryManager = new DataQueryManager();
    private double cumulativeExecutionTime;
    private double cumulativeMemoryUsage;
    private MemoryCalculationStrategy performanceMonitor;
    private String sourceDataset;
    private String querySetIdentifier;
    private int maxWindowSize;
    private int slideIntervalGCD;
    private int totalSlides;
    public int outlierQueryCount;
    public int detectedOutliers;

    public Processor(String sourceData, String querySet) throws IOException {
        initializeDataset(sourceData);
        initializeQuerySet(querySet);
        this.sourceDataset = sourceData;
        this.querySetIdentifier = querySet;
        this.slideIntervalGCD = queryManager.getSlideIntervalGCD();
        this.totalSlides = queryManager.getMaxWindowSize() / queryManager.getSlideIntervalGCD();
        performanceMonitor = new MemoryCalculationStrategy();         
    }

    public void executeAnalysis(int windowSpan, int queryVolume, double queryChangeFactor) throws IOException {
        MDUAL analysisAlgorithm = new MDUAL(datasetManager.getDimensions(), datasetManager.getDimensions(), totalSlides, slideIntervalGCD, datasetManager.getLowerBounds());
        performanceMonitor.beginMonitoring();

        int processedWindows = 0;
        int adjustedQueryVolume = (int) (queryVolume * queryChangeFactor);
        int detectedOutliers = 0;
        int outlierQueryCount = 0;
        
        for (int cycle = 0; cycle < windowSpan + totalSlides - 1; cycle++) {
            Map<Integer, Query> activeQueries = queryManager.retrieveQueriesByRange(cycle * adjustedQueryVolume, queryVolume);
            
            if (activeQueries.isEmpty()) break;
            List<Tuple> slideData = datasetManager.fetchSlideData(cycle, slideIntervalGCD);
            if (slideData.isEmpty()) break;
            Set<Tuple> identifiedOutliers = new HashSet<>();

            long startMoment = PerformanceMetric.getCPUTime();
            identifiedOutliers = analysisAlgorithm.identifyOutliers(slideData, activeQueries, cycle);

            long endMoment = Measure.getCPUTime();
            long memoryConsumed = Measure.getMemory();
            
            if (cycle >= totalSlides - 1) {
                detectedOutliers += identifiedOutliers.size();
                for (Tuple t : identifiedOutliers) outlierQueryCount += t.getOutlierQueryIDs().size();
                processedWindows++;
                cumulativeExecutionTime += (endMoment - startMoment) / 1000000.0; // Convert to ms
                cumulativeMemoryUsage += memoryConsumed;
            }
        }

        outputResults(windowSpan, queryVolume, queryChangeFactor, processedWindows);
        performanceMonitor.beginMonitoring();
    }

    private void initializeDataset(String data) throws IOException {
        datasetManager.prepareData(data);
    }

    private void initializeQuerySet(String querySet) throws IOException {
        queryManager.setupQueries(querySet);
    }

    private void outputResults(int windowSpan, int queryVolume, double queryChangeFactor, int processedWindows) {
        int detectedOutliers = 0;
        int outlierQueryCount = 0;
        System.out.println(String.format("%-10s %10s %10.1f %10.2f %10.1f %10.1f %10d %10d", sourceDataset,
                querySetIdentifier, queryChangeFactor, cumulativeExecutionTime / processedWindows,
                cumulativeMemoryUsage / processedWindows, performanceMonitor.getPeakMemoryUsage(),
                detectedOutliers / processedWindows, outlierQueryCount / processedWindows));
    }

}
