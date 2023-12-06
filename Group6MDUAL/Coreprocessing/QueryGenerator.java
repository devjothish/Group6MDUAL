package Coreprocessing;

import java.io.*;
import java.util.Random;
import java.util.logging.*;
import MDUAL.Query;
public class QueryGenerator {

    private final QueryConfig config;
    private final Random random;
    private final Logger logger;

    public QueryGenerator(QueryConfig config) {
        this.config = config;
        this.random = new Random();
        this.logger = Logger.getLogger(QueryGenerator.class.getName());
    }

    public void generateAndSaveQueries(int numQueries, String[] varyingParams, String filePath) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
            for (int i = 0; i < numQueries; i++) {
                Query query = generateOne(i, varyingParams);
                writer.write(query.toString());
                writer.newLine();
            }
        } catch (IOException e) {
            logger.log(Level.SEVERE, "Error writing to file", e);
        }
    }

    private Query generateOne(int qID, String[] varyingParams) {
        QueryParams queryParams = new QueryParams(config.defaultR, config.defaultK, config.defaultW, config.gcdS);
        for (String param : varyingParams) {
            adjustQueryParam(param, queryParams);
        }
        return new Query(qID, queryParams.R, queryParams.K, queryParams.W, queryParams.S);
    }

    private void adjustQueryParam(String param, QueryParams queryParams) {
        switch (param) {
            case "R":
                queryParams.R = roundToTwoDecimals(1 + random.nextDouble() * (config.variationTimes - 1) * config.defaultR);
                break;
            case "K":
                queryParams.K = (int) (config.defaultK * (1 + random.nextDouble() * config.variationTimes));
                break;
            case "S":
                queryParams.S = (int) (config.gcdS * (1 + random.nextDouble() * config.variationTimes));
                break;
            case "W":
                queryParams.W = queryParams.S + queryParams.S * (int) (random.nextDouble() * (config.defaultW * config.variationTimes / queryParams.S));
                break;
        }
    }

    private double roundToTwoDecimals(double value) {
        return Math.round(value * 100) / 100.0;
    }
}

class QueryConfig {
    final String dataset;
    final int defaultW, gcdS, defaultK, variationTimes;
    final double defaultR;

    QueryConfig(String dataset, int defaultW, int gcdS, int defaultK, int variationTimes) {
        this.dataset = dataset;
        this.defaultW = defaultW;
        this.gcdS = gcdS;
        this.defaultK = defaultK;
        this.variationTimes = variationTimes;
        this.defaultR = determineDefaultR(dataset);
    }

    private double determineDefaultR(String dataset) {
        // Logic to determine defaultR based on dataset
        return 0.5; // Default return for this example
    }
}

class QueryParams {
    double R;
    int K, W, S;

    QueryParams(double R, int K, int W, int S) {
        this.R = R;
        this.K = K;
        this.W = W;
        this.S = S;
    }
}
