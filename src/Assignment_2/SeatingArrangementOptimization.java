package Assignment_2;

import java.util.*;


public class SeatingArrangementOptimization {
    static double[][] dislikes = {
        {0.0, 0.68, 0.55, 0.30, 0.82, 0.48, 0.33, 0.10, 0.76, 0.43},
        {0.68, 0.0, 0.90, 0.11, 0.76, 0.20, 0.55, 0.17, 0.62, 0.99},
        {0.55, 0.90, 0.0, 0.70, 0.63, 0.96, 0.51, 0.90, 0.88, 0.64},
        {0.30, 0.11, 0.70, 0.0, 0.91, 0.86, 0.78, 0.99, 0.53, 0.92},
        {0.82, 0.76, 0.63, 0.91, 0.0, 0.43, 0.88, 0.53, 0.42, 0.75},
        {0.48, 0.20, 0.96, 0.86, 0.43, 0.0, 0.63, 0.97, 0.37, 0.26},
        {0.33, 0.55, 0.51, 0.78, 0.88, 0.63, 0.0, 0.92, 0.87, 0.81},
        {0.10, 0.17, 0.90, 0.99, 0.53, 0.97, 0.92, 0.0, 0.81, 0.78},
        {0.76, 0.62, 0.88, 0.53, 0.42, 0.37, 0.87, 0.81, 0.0, 0.45},
        {0.43, 0.99, 0.64, 0.92, 0.75, 0.26, 0.81, 0.78, 0.45, 0.0}
    };

    static String[] guests = {
        "Ahmad", "Salem", "Ayman", "Hani", "Kamal", "Samir", "Hakem", "Fuad", "Ibrahim", "Khalid"
    };

    static int numGuests = guests.length;
    static Random random = new Random();

    public static void main(String[] args) {
        Graph graph = new Graph(numGuests);
        // Construct graph with dislike percentages as weights
        for (int i = 0; i < numGuests; i++) {
            for (int j = i + 1; j < numGuests; j++) {
                double dislikePercentage = dislikes[i][j];
                graph.addEdge(i, j, dislikePercentage);
            }
        }

        // Genetic Algorithm
        List<String> arrangementGenetic = geneticAlgorithm(graph, 100, 1000, 0.1);
        double costGenetic = calculateCost(arrangementGenetic, graph);
        System.out.println("Genetic Algorithm \nBest seating arrangement: " + arrangementGenetic);
        System.out.println("Total Cost: " + costGenetic);

        // Simulated Annealing
        List<String> arrangementSimulatedAnnealing = simulatedAnnealing(graph, 1000, 0.99, 10000);
        double costSimulatedAnnealing = calculateCost(arrangementSimulatedAnnealing, graph);
        System.out.println("\nSimulated Annealing \nBest seating arrangement: " + arrangementSimulatedAnnealing);
        System.out.println("Total Cost: " + costSimulatedAnnealing);

        // Hill Climbing
        List<String> arrangementHillClimbing = hillClimbing(graph, 100);
        double costHillClimbing = calculateCost(arrangementHillClimbing, graph);
        System.out.println("\nHill Climbing \nBest seating arrangement: " + arrangementHillClimbing);
        System.out.println("Total Cost: " + costHillClimbing);
    }

    static double calculateCost(List<String> arrangement, Graph graph) {
        double totalCost = 0;
        for (int i = 0; i < arrangement.size() - 1; i++) {
            String guest1 = arrangement.get(i);
            String guest2 = arrangement.get(i + 1);
            int index1 = Arrays.asList(guests).indexOf(guest1);
            int index2 = Arrays.asList(guests).indexOf(guest2);
            totalCost += graph.adjMatrix[index1][index2];
        }
        // Add the cost between the last and the first guest to make it a round table
        int lastIndex = Arrays.asList(guests).indexOf(arrangement.get(arrangement.size() - 1));
        int firstIndex = Arrays.asList(guests).indexOf(arrangement.get(0));
        totalCost += graph.adjMatrix[lastIndex][firstIndex];
        return totalCost;
    }

    static List<String> geneticAlgorithm(Graph graph, int populationSize, int numGenerations, double mutationRate) {
        List<List<String>> population = new ArrayList<>();

        // Initialize population with random arrangements
        for (int i = 0; i < populationSize; i++) {
            List<String> arrangement = new ArrayList<>(Arrays.asList(guests));
            Collections.shuffle(arrangement);
            population.add(arrangement);
        }

        for (int gen = 0; gen < numGenerations; gen++) {
            // Evaluate fitness
            population.sort(Comparator.comparingDouble(arr -> calculateCost(arr, graph)));

            // Selection: keep the top 60%
            int cutoff = (int) (0.6 * populationSize);
            List<List<String>> newPopulation = new ArrayList<>(population.subList(0, cutoff));

            // Crossover 
            while (newPopulation.size() < populationSize) {
                List<String> parent1 = population.get(random.nextInt(cutoff));
                List<String> parent2 = population.get(random.nextInt(cutoff));
                List<String> child = crossover(parent1, parent2);
                mutate(child, mutationRate);
                newPopulation.add(child);
            }
            population = newPopulation;
        }

        // Return the best arrangement found
        return population.stream().min(Comparator.comparingDouble(arr -> calculateCost(arr, graph))).orElse(new ArrayList<>());
    }

    static List<String> crossover(List<String> parent1, List<String> parent2) {
        int crossoverPoint = random.nextInt(parent1.size());
        Set<String> childSet = new HashSet<>(parent1.subList(0, crossoverPoint));
        List<String> child = new ArrayList<>(parent1.subList(0, crossoverPoint));
        for (String guest : parent2) {
            if (!childSet.contains(guest)) {
                child.add(guest);
                childSet.add(guest);
            }
        }
        return child;
    }

    static void mutate(List<String> arrangement, double mutationRate) {
        if (random.nextDouble() < mutationRate) {
            int index1 = random.nextInt(arrangement.size());
            int index2 = random.nextInt(arrangement.size());
            Collections.swap(arrangement, index1, index2);
        }
    }

    static List<String> simulatedAnnealing(Graph graph, double initialTemperature, double coolingRate, int numIterations) {
        List<String> currentArrangement = new ArrayList<>(Arrays.asList(guests));
        Collections.shuffle(currentArrangement);
        List<String> bestArrangement = new ArrayList<>(currentArrangement);

        double currentCost = calculateCost(currentArrangement, graph);
        double bestCost = currentCost;

        double temperature = initialTemperature;

        for (int i = 0; i < numIterations; i++) {
            List<String> newArrangement = new ArrayList<>(currentArrangement);
            int index1 = random.nextInt(newArrangement.size());
            int index2 = random.nextInt(newArrangement.size());
            Collections.swap(newArrangement, index1, index2);

            double newCost = calculateCost(newArrangement, graph);

            if (acceptanceProbability(currentCost, newCost, temperature) > random.nextDouble()) {
                currentArrangement = newArrangement;
                currentCost = newCost;
            }

            if (newCost < bestCost) {
                bestArrangement = newArrangement;
                bestCost = newCost;
            }

            temperature *= coolingRate;
        }

        return bestArrangement;
    }

    static double acceptanceProbability(double currentCost, double newCost, double temperature) {
        if (newCost < currentCost) {
            return 1.0;
        }
        return Math.exp((currentCost - newCost) / temperature);
    }

    static List<String> hillClimbing(Graph graph, int numRestarts) {
        List<String> bestArrangement = new ArrayList<>(Arrays.asList(guests));
        Collections.shuffle(bestArrangement);
        double bestCost = calculateCost(bestArrangement, graph);

        for (int restart = 0; restart < numRestarts; restart++) {
            List<String> currentArrangement = new ArrayList<>(Arrays.asList(guests));
            Collections.shuffle(currentArrangement);
            double currentCost = calculateCost(currentArrangement, graph);

            boolean improved;
            do {
                improved = false;
                for (int i = 0; i < currentArrangement.size() - 1; i++) {
                    for (int j = i + 1; j < currentArrangement.size(); j++) {
                        List<String> newArrangement = new ArrayList<>(currentArrangement);
                        Collections.swap(newArrangement, i, j);
                        double newCost = calculateCost(newArrangement, graph);
                        if (newCost < currentCost) {
                            currentArrangement = newArrangement;
                            currentCost = newCost;
                            improved = true;
                        }
                    }
                }
            } while (improved);

            if (currentCost < bestCost) {
                bestArrangement = currentArrangement;
                bestCost = currentCost;
            }
        }

        return bestArrangement;
    }

    static class Graph {
        int V;
        double[][] adjMatrix;

        public Graph(int V) {
            this.V = V;
            this.adjMatrix = new double[V][V];
        }

        public void addEdge(int u, int v, double weight) {
            adjMatrix[u][v] = weight;
            adjMatrix[v][u] = weight;
        }
    }
}
