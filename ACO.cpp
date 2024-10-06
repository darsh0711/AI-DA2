// ALPHA and BETA are the parameters that control the influence of pheromone and distance respectively.
// RHO is the evaporation rate of pheromone and Q is the total pheromone deposited by each ant.
// These parameters can be tuned to get better results for different problems.
// Larger the value of ALPHA, more importance is given to pheromone.
// Larger the value of BETA, more importance is given to distance.
// Larger the value of RHO, faster the pheromone evaporates.
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <ctime>

using namespace std;

const double ALPHA = 1.0;   // Influence of pheromone
const double BETA = 5.0;    // Influence of distance (visibility)
const double RHO = 0.5;     // Evaporation rate
const double Q = 100.0;     // Total pheromone deposited

// City class to store coordinates of each city
class City {
public:
    int id;
    double x, y;

    City(int id, double x, double y) : id(id), x(x), y(y) {}

    // Function to calculate Euclidean distance between two cities
    double distanceTo(City other) {
        return sqrt(pow(x - other.x, 2) + pow(y - other.y, 2));
    }
};

// Ant class to represent each ant and its tour
class Ant {
public:
    vector<int> tour;           // Stores the cities visited
    vector<bool> visited;       // Keeps track of visited cities
    double tourLength;          // Total distance of the tour
    int currentCity;            // Current city where the ant is

    Ant(int numCities) {
        visited.resize(numCities, false);
        tourLength = 0;
    }

    // Choose the next city based on pheromone and visibility
    // Choose the next city based on pheromone and visibility
int chooseNextCity(vector<vector<double>>& pheromone, vector<vector<double>>& visibility) {
    int nextCity = -1;
    double sumProbabilities = 0.0;
    vector<double> probabilities(visited.size(), 0.0);  // Fix: use visited.size() instead of tour.size()

    // Calculate probability for each unvisited city
    for (int i = 0; i < visited.size(); ++i) {
        if (!visited[i]) {
            double pheromoneLevel = pow(pheromone[currentCity][i], ALPHA);
            double visibilityLevel = pow(visibility[currentCity][i], BETA);
            probabilities[i] = pheromoneLevel * visibilityLevel;
            sumProbabilities += probabilities[i];
        }
    }

    // Check if sumProbabilities is greater than 0 (prevents division by zero)
    if (sumProbabilities == 0) {
        cerr << "Error: All cities seem visited or no valid probabilities" << endl;
        return nextCity;
    }

    // Roulette wheel selection for the next city
    double randomValue = ((double) rand() / RAND_MAX) * sumProbabilities;
        for (int i = 0; i < visited.size(); ++i) {
            if (!visited[i]) {
                randomValue -= probabilities[i];
                if (randomValue <= 0.0) {
                    nextCity = i;
                    break;
                }
            }
        }

        return nextCity;
    }


    // Reset ant for the next iteration
    void reset(int startCity) {
        fill(visited.begin(), visited.end(), false);
        tour.clear();
        tourLength = 0;
        currentCity = startCity;
        visited[startCity] = true;
        tour.push_back(startCity);
    }

    // Add city to the tour and mark it as visited
    void visitCity(int city, double distance) {
        currentCity = city;
        tour.push_back(city);
        visited[city] = true;
        tourLength += distance;
    }
};

// ACO class to manage the algorithm
class ACO {
public:
    vector<City> cities;                   // List of cities
    vector<vector<double>> pheromone;      // Pheromone levels between cities
    vector<vector<double>> visibility;     // Inverse of distances between cities
    vector<Ant> ants;                      // List of ants
    int numCities, numAnts;                // Number of cities and ants

    ACO(int numCities, int numAnts) : numCities(numCities), numAnts(numAnts) {
        pheromone.resize(numCities, vector<double>(numCities, 1.0));  // Initialize pheromone to 1
        visibility.resize(numCities, vector<double>(numCities, 0.0)); // Visibility (1/distance)
        ants.resize(numAnts, Ant(numCities));
    }

    // Add city to the list of cities
    void addCity(int id, double x, double y) {
        cities.push_back(City(id, x, y));
    }

    // Calculate the distance matrix and visibility (1/distance)
    void calculateVisibility() {
        for (int i = 0; i < numCities; ++i) {
            for (int j = 0; j < numCities; ++j) {
                if (i != j) {
                    double dist = cities[i].distanceTo(cities[j]);
                    visibility[i][j] = 1.0 / dist;
                }
            }
        }
    }

    // Run the ACO algorithm
    void runACO(int iterations) {
        srand(time(0));  // Initialize random seed

        for (int iter = 0; iter < iterations; ++iter) {
            // Place ants randomly in the cities
            for (int i = 0; i < numAnts; ++i) {
                int startCity = rand() % numCities;
                ants[i].reset(startCity);
            }

            // Each ant constructs a tour
            for (int step = 0; step < numCities - 1; ++step) {
                for (int i = 0; i < numAnts; ++i) {
                    Ant& ant = ants[i];
                    int nextCity = ant.chooseNextCity(pheromone, visibility);
                    double distance = cities[ant.currentCity].distanceTo(cities[nextCity]);
                    ant.visitCity(nextCity, distance);
                }
            }

            // Update pheromone levels
            updatePheromones();
        }

        // Find the best solution
        double bestLength = numeric_limits<double>::max();
        vector<int> bestTour;
        for (int i = 0; i < numAnts; ++i) {
            if (ants[i].tourLength < bestLength) {
                bestLength = ants[i].tourLength;
                bestTour = ants[i].tour;
            }
        }

        // Print the best solution
        cout << "Best Tour: ";
        for (int city : bestTour) {
            cout << city << " -> ";
        }
        cout << bestTour[0] << endl;
        cout << "Tour Length: " << bestLength << endl;
    }

    // Update pheromones based on the tours taken by ants
    void updatePheromones() {
        // Evaporate pheromones
        for (int i = 0; i < numCities; ++i) {
            for (int j = 0; j < numCities; ++j) {
                pheromone[i][j] *= (1.0 - RHO);
            }
        }

        // Add pheromones for each ant's tour
        for (int i = 0; i < numAnts; ++i) {
            Ant& ant = ants[i];
            double contribution = Q / ant.tourLength;
            for (int j = 0; j < ant.tour.size() - 1; ++j) {
                int from = ant.tour[j];
                int to = ant.tour[j + 1];
                pheromone[from][to] += contribution;
                pheromone[to][from] += contribution;
            }
        }
    }
};

int main() {
    int numCities = 5;
    int numAnts = 5;
    int iterations = 50;

    ACO aco(numCities, numAnts);

    // Add cities (id, x, y)
    aco.addCity(0, 1.0, 1.0);
    aco.addCity(1, 4.0, 1.0);
    aco.addCity(2, 4.0, 5.0);
    aco.addCity(3, 1.0, 5.0);
    aco.addCity(4, 2.5, 3.0);

    // Calculate distances and visibility
    aco.calculateVisibility();

    // Run the ACO algorithm
    aco.runACO(iterations);

    return 0;
}
