package MDUAL;

import java.util.Scanner;
import Coreprocessing.*;
import MDUAL.*;

public class Generator {
    public static void main(String[] args) {
        UserInput input = gatherInput();
        QueryBuilder queryBuilder = new QueryBuilder(input.dataSource, input.windowSize, input.slideStep, input.kThreshold, input.variationScale);

        System.out.println(String.format("%10s %12s %10s %11s %8s %10s %11s", "QuerySet", "ChangeRate", "MemAvg", "MemPeak", "CPULoad", "Detected", "QueriesOut"));

        runSimulations(queryBuilder, input);
    }

    private static UserInput gatherInput() {
        Scanner scanner = new Scanner(System.in);

        System.out.print("Enter data source: ");
        String dataSource = scanner.nextLine();

        System.out.print("Enter window size (W): ");
        int windowSize = Integer.parseInt(scanner.nextLine());

        System.out.print("Enter slide step (S): ");
        int slideStep = Integer.parseInt(scanner.nextLine());

        System.out.print("Enter K threshold: ");
        int kThreshold = Integer.parseInt(scanner.nextLine());

        System.out.print("Number of windows (N): ");
        int numWindows = Integer.parseInt(scanner.nextLine());

        System.out.println("Simulation will start. Press CTRL+C to stop.");
        scanner.close();

        return new UserInput(dataSource, windowSize, slideStep, kThreshold, numWindows, 10, 0.2, 5);
    }

    private static void runSimulations(QueryBuilder queryBuilder, UserInput input) {
        int runIndex = 0;
        while (runIndex < input.totalRuns) {
            String querySet = queryBuilder.create(input.numWindows, input.queryCount, input.paramsToVary);
            SimulationCoordinator simulation = new SimulationCoordinator();
            simulation.run(input.numWindows, input.queryCount, input.queryChangeRate);
            runIndex++;
        }
    }

    private static class UserInput {
        String dataSource;
        int windowSize;
        int slideStep;
        int kThreshold;
        int numWindows;
        int queryCount;
        double queryChangeRate;
        int totalRuns;
        String[] paramsToVary;

        UserInput(String dataSource, int windowSize, int slideStep, int kThreshold, int numWindows, int queryCount, double queryChangeRate, int totalRuns) {
            this.dataSource = dataSource;
            this.windowSize = windowSize;
            this.slideStep = slideStep;
            this.kThreshold = kThreshold;
            this.numWindows = numWindows;
            this.queryCount = queryCount;
            this.queryChangeRate = queryChangeRate;
            this.totalRuns = totalRuns;
            this.paramsToVary = new String[] {"R", "K", "S", "W"};
        }
    }
}
