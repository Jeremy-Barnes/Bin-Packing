
package genetic_algo;

/**
 * @author Jeremy Barnes
 * Bin packing Simple Genetic Algorithm
 * 
 * Explanation of chromosome: A chromosome is an integer string saved in an array. 
 * The int string is order-based, so that chromosome[0]
 * refers to the same package each time (Package 0). 
 * The value of the array at that index is the bin that package has been placed in. 
 * So chromosome[0] = 3 indicates that the first package in the list has been placed in bin 3. 
 * The weights of those packages are stored in a
 * different array, packages[]. This means that the packages are in a constant order, 
 * and the packages array never needs modification.
 * Instead of moving packages, we're moving bins.
 * 
 * To deal with infeasible chromosomes, we honor the initial packing as much as possible, 
 * chromosome[0] will go in Bin 3 (see above), unless we find
 * that Bin 3 has already been filled, in which case we will default to 
 * Best-First packing strategy to place that package and repair the
 * chromosome. The repaired chromosome will be the one used to get a fitness value.
 * 
 * This class has functions for Simulated Annealing (SA), Genetic Algorithm (GA), and Foolish Hill Climbing (FHC/HC)
 */
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

import javax.swing.JOptionPane;


public class Core {

	private static Random rand = new Random();

	/* ************ General Bin Packing Parameters - most set by interface ************* */
	protected int binCount;
	protected int packageCount;
	protected int binSize;
	protected int maxPackageSize;
	protected int[] packageWeights;
	String crossover;

	/* ************ GA Parameters - most set by interface ************* */
	protected String mutationMethod;
	protected float selectionRate = .6f; //not set by interface
	float eliteFitness = 0;
	int[] eliteChromosome;

	/* ************ SA Parameters - most set by interface ************* */
	int timeElapsed;
	int timeMax;
	float temperature;
	int[] currentSolution;



	public static void main(String[] args) {
		new Core();
	}

	public Core() {
		WindowInterface win = new WindowInterface(this);
	}

	public void runGeneticAlgorithm(String sel, String mut, float mutRate, String cross, int gens, int[] pkgs, int pkgCt, int binCt, int binSize) {
		eliteFitness = binCt * 10;
		this.binSize = binSize;
		maxPackageSize = binSize / 2;
		packageCount = pkgCt;
		packageWeights = pkgs;
		binCount = binCt;
		crossover = cross;

		if (pkgs.length != pkgCt) { //data sanity check
			JOptionPane.showMessageDialog(null, "You have invalid data - the package count specified does not match the number of packages given! \nExiting.");
			return;
		}
		for (int i = 0; i < pkgs.length; i++) {
			if (pkgs[i] > binSize) {
				JOptionPane.showMessageDialog(null, "You have invalid data - the package at position " + i + " is larger than all the bins! \n(Package size is : " + pkgs[i]
					+ " and bin size is " + binSize + ") \nExiting.");
				return;
			}
			if (pkgs[i] > binSize / 2) {
				JOptionPane.showMessageDialog(null, "WARNING: The problem specification says \"each package has positive size ranging from 1 unit to half the size of each bin\". \n\n "
					+ "The package at position " + i + " is larger than the problem's max size! (Package size is : " + pkgs[i] + " and bin size is " + binSize + ")\n"
					+ "The algorithm will run anyways, but may not function correctly. Good luck! ");
				break;
			}
		}

		int[][] population = generateInitialPopulation(binCt, binSize, pkgCt, pkgCt * 20); //initial population
		eliteChromosome = population[0];

		for (int generation = 0; generation < gens; generation++) {
			System.out.println("Generation " + generation);

			eliteChromosome = getElite(population);
			if (canReduceBins()) { //try to halve the problem space with a binary chopper
				binCount = binCount / 2;
				population = generateInitialPopulation(binCount, binSize, pkgCt, pkgCt * 20);
				System.out.println("Binary chop - " + binCount);
			} else { //remove KNOWN unnecessary bins from consideration
				int minBins = getFitness(eliteChromosome);
				if (totalInfeasibility(eliteChromosome) == 0 && minBins < binCount) {
					System.out.println("Reduced problem space to " + (minBins));
					binCount = minBins;
				}
			}

			//select breeders
			ArrayList<int[]> breeders = null;
			if (sel.equalsIgnoreCase("Roulette"))
				breeders = selectionRoulette(population, selectionRate);
			else if (sel.equalsIgnoreCase("Rank"))
				breeders = selectionRank(population, selectionRate);

			ArrayList<int[]> replacements = new ArrayList<int[]>();
			while (breeders.size() > 0) {//reproduction loop

				int rand1 = rand.nextInt(breeders.size()); //randomly choose breeding pair
				int rand2 = rand.nextInt(breeders.size());
				while (rand2 == rand1 && breeders.size() > 1) { //try to sure no parents breeding with themselves.
					rand2 = rand.nextInt(breeders.size());
				}
				if (breeders.size() == 1) { //but single parents have to mate with themselves.
					replacements.add(breeders.get(0));
					breeders.remove(0);
					continue;
				}

				int[] child1 = null; //crossover
				int[] child2 = null;
				int[][] children = null;
				if (crossover.equalsIgnoreCase("One Point Crossover")) {
					children = crossoverOnePoint(breeders.get(rand1), breeders.get(rand2));
					child1 = children[0];
					child2 = children[1];
				} else if (crossover.equalsIgnoreCase("Two Point Crossover")) {
					children = crossoverTwoPoint(breeders.get(rand1), breeders.get(rand2));
					child1 = children[0];
					child2 = children[1];
				}

				if (rand.nextDouble() <= mutRate) { //mutate child1 
					if (mut.equalsIgnoreCase("Pairwise Exchange")) {
						child1 = mutationPairwiseExchange(child1);
					} else {
						child1 = mutationSingleMove(child1);
					}
				}
				if (rand.nextDouble() <= mutRate) { //mutate child2
					if (mut.equalsIgnoreCase("Pairwise Exchange")) {
						child1 = mutationPairwiseExchange(child1);
					} else {
						child1 = mutationSingleMove(child1);
					}
				}

				replacements.add(child1);
				replacements.add(child2);
				breeders.remove(rand1);
				if (rand2 > rand1) //deal with shifted arraylist
					breeders.remove(rand2 - 1);
				else
					breeders.remove(rand2);

			} //end reproduction

			population = replaceLowestFitness(population, replacements); //put in the children

		} //END GA



		System.out.println();
		System.out.println();

		System.out.println("Bins Used:");
		eliteChromosome = getElite(population);
		System.out.println(getFitness(eliteChromosome));



		for (int i = 0; i < binCount; i++) {
			int binWeight = 0;
			for (int j = 0; j < eliteChromosome.length; j++) {
				if (eliteChromosome[j] == i) {
					binWeight += packageWeights[j];
				}
			}
			System.out.println("There are " + binWeight + " weight units in Bin " + i);
		}

		
		System.out.println("Chromosome: ");
		for (int j = 0; j < eliteChromosome.length; j++) {
				System.out.print(eliteChromosome[j] + ",");
		}
	}

	public void runSimulatedAnnealing(String perturb, int timeToRun, int[] pkgs, int pkgCt, int binCt, int binSize, int innerIterations, boolean foolish) {
		if (pkgs.length != pkgCt) { //data sanity check
			JOptionPane.showMessageDialog(null, "You have invalid data - the package count specified does not match the number of packages given! \nExiting.");
			return;
		}
		for (int i = 0; i < pkgs.length; i++) {
			if (pkgs[i] > binSize) {
				JOptionPane.showMessageDialog(null, "You have invalid data - the package at position " + i + " is larger than all the bins! \n(Package size is : " + pkgs[i]
					+ " and bin size is " + binSize + ") \nExiting.");
				return;
			}
		}

		this.binSize = binSize;
		packageCount = pkgCt;
		packageWeights = pkgs;
		binCount = binCt;

		temperature = rand.nextFloat(); //initial temp
		currentSolution = generateInitialPopulation(binCount, binSize, pkgCt, 1)[0];

		for (int generation = 0; generation < timeToRun; generation++) {
			System.out.println(generation + " out of " + timeToRun + " inner loop iterations - " + innerIterations);
			for (int iterations = 0; iterations < innerIterations; iterations++) {
				int[] newSolution = currentSolution.clone(); //mutation uses side-effects, clone to avoid
				if (perturb.equalsIgnoreCase("Pairwise Exchange")) {
					newSolution = mutationPairwiseExchange(newSolution);
				} else {
					newSolution = mutationSingleMove(newSolution);
				}
				int hNew = h(newSolution);
				int hCurrent = h(currentSolution);

				if (!foolish)
					if (hNew < hCurrent || (rand.nextDouble() < Math.pow(Math.E, (hCurrent - hNew) / temperature))) {
						currentSolution = newSolution;
					} else if (hNew < hCurrent) {
						currentSolution = newSolution;
					}

			}
			temperature = (temperature * 0.8f); //alpha
			innerIterations *= 1.2; //beta

		}
		System.out.println("Bins used: " + h(currentSolution));
		
		for (int i = 0; i < binCount; i++) {
			int binWeight = 0;
			for (int j = 0; j < currentSolution.length; j++) {
				if (currentSolution[j] == i) {
					binWeight += packageWeights[j];
				}
			}
			System.out.println("There are " + binWeight + " weight units in Bin " + i);
		}

		
		System.out.println("Solution: ");
		for (int j = 0; j < currentSolution.length; j++) {
				System.out.print(currentSolution[j] + ",");
		}
		
	}

	/**
	 * Generates a list of chromosomes entirely randomly by placing the packages into random bins. No concern for feasibility here.
	 */
	private int[][] generateInitialPopulation(int binCt, int binSize, int pkgCt, int populationSize) {
		int[][] population = new int[populationSize][pkgCt];
		for (int i = 0; i < populationSize; i++) {
			for (int j = 0; j < pkgCt; j++) {
				population[i][j] = rand.nextInt(binCt);
			}
		}
		return population;
	}

	/**
	 * Evaluates fitness of chromosomes, chooses some for breeding (minimization roulette) returns chromosomes that won breeding lottery.
	 */
	private ArrayList<int[]> selectionRoulette(int[][] chromosomes, float percentToSelect) {

		ArrayList<int[]> selectedChromosomes = new ArrayList<int[]>();

		//sum up fitnesses for population
		int fitnessSum = 0;
		for (int[] chromosome : chromosomes) {
			fitnessSum += getFitness(chromosome);
		}

		//minimization problem - invert fitness
		int invertFitnessSum = 0;
		for (int[] chromosome : chromosomes) {
			invertFitnessSum += Math.round((fitnessSum * 1.0) / getFitness(chromosome)); //small fitness values (few bins) are favored
		}

		//get fitnesses of chromosomes, save in array to avoid recalculating
		int[] fitnesses = new int[chromosomes.length];
		for (int i = 0; i < chromosomes.length; i++) {
			fitnesses[i] = (int) Math.round((fitnessSum * 1.0) / getFitness(chromosomes[i]));
		}

		//spin the roulette wheel
		for (int i = 0; i < chromosomes.length * percentToSelect; i++) {
			int currentFitnessTotal = 0;
			int fitnessCutoff = rand.nextInt(invertFitnessSum); //where the wheel stops

			for (int j = 0; j < chromosomes.length; j++) { //go through all chromosomes to see who is selected this spin
				if (currentFitnessTotal < fitnessCutoff)
					currentFitnessTotal += fitnesses[j];

				if (currentFitnessTotal >= fitnessCutoff) { //if this chromosome pushed us over the edge, save it
					selectedChromosomes.add(chromosomes[j]);
					break;
				}
			}
		}//end roulette spinning
		return selectedChromosomes;
	}

	/**
	 * Evaluates the fitness of chromosomes, chooses some for breeding (min rank) returns chromosomes that won breeding lottery.
	 */
	private ArrayList<int[]> selectionRank(int[][] chromosomes, float percentToSelect) {

		ArrayList<int[]> selectedChromosomes = new ArrayList<int[]>();

		//get fitnesses and ranks of chromosomes, save in array to avoid recalculating
		int[] fitnesses = new int[chromosomes.length];
		ArrayList<Integer> ranks = new ArrayList<Integer>(); //holds indices of chromosomes, position in arraylist indicates rank

		for (int i = 0; i < chromosomes.length; i++) { //set up list of fitnesses and ranks
			fitnesses[i] = getFitness(chromosomes[i]);


			boolean added = false;
			for (int keyIndex = 0; keyIndex < ranks.size(); keyIndex++) { //only search the second half
				if (fitnesses[ranks.get(keyIndex)] > fitnesses[i]) { //if this one is worse
					ranks.add(keyIndex, i); //inverted rank
					added = true;
					break;
				}
			}
			if (!added) {
				ranks.add(i);
			}

		}// end fitness ranking (low fitness is best, last index is worst)


		//spin the roulette wheel
		for (int i = 0; i < chromosomes.length * percentToSelect; i++) {
			int currentFitnessTotal = 0;
			int fitSelector = rand.nextInt(chromosomes.length); //where the wheel stops
			int fitnessCutoff = (fitSelector * (fitSelector + 1)) / 2;

			for (int j = 0; j < ranks.size(); j++) { //go through all ranks to see who is selected this spin
				if (currentFitnessTotal < fitnessCutoff)
					currentFitnessTotal += chromosomes.length - ranks.get(j); //reorder fitness properly

				if (currentFitnessTotal >= fitnessCutoff) { //if this rank pushed us over the edge, save it
					selectedChromosomes.add(chromosomes[ranks.get(j)]);
					break;
				}
			}
		}//end roulette spinning

		return selectedChromosomes;
	}

	/**
	 * Crossover of two parents using two random crossover points, crossover takes form - P2 | P1 | P2
	 */
	private int[][] crossoverTwoPoint(int[] parent1, int[] parent2) {
		int crossPoint1 = rand.nextInt(parent1.length);
		int crossLength = rand.nextInt(parent1.length + 1 - crossPoint1);
		int crossPoint2 = crossPoint1 + crossLength;

		int[] child = new int[parent1.length];
		int[] child2 = new int[parent1.length];
		for (int i = crossPoint1; i < crossPoint2; i++) {
			child[i] = parent1[i]; //items between cross points are direct copied from p1
			child2[i] = parent2[i]; //opposite for child 2
		}
		for (int i = 0; i < crossPoint1; i++) {
			child[i] = parent2[i]; //items on either side copied direct from p2
			child2[i] = parent1[i]; //opposite for child 2
		}
		for (int i = crossPoint2; i < parent1.length; i++) {
			child[i] = parent2[i]; //items on either side copied direct from p2
			child2[i] = parent1[i]; //opposite for child 2
		}
		return new int[][] { child, child2 };
	}

	/**
	 * Crossover of two parents using two random crossover points - crossover takes form P1|P2
	 */
	private int[][] crossoverOnePoint(int[] parent1, int[] parent2) {
		int crossPoint1 = rand.nextInt(parent1.length);

		int[] child = new int[parent1.length];
		int[] child2 = new int[parent1.length];
		for (int i = 0; i < crossPoint1; i++) {
			child[i] = parent1[i]; //before cross point = p1
			child2[i] = parent2[i]; //opposite for child2
		}

		for (int i = crossPoint1; i < parent1.length; i++) {
			child[i] = parent2[i]; //after and including crosspoint = p2
			child2[i] = parent1[i]; //opposite for child2
		}
		return new int[][] { child, child2 };

	}

	/**
	 * Mutation - Chooses two elements in different bins and swaps them (acts like partition problem)
	 */
	private int[] mutationPairwiseExchange(int[] mutant) {
		int pairPoint1 = rand.nextInt(mutant.length);
		int pairPoint2 = pairPoint1;
		while (mutant[pairPoint2] == mutant[pairPoint1]) {
			pairPoint2 = rand.nextInt(mutant.length); //find element in different bucket
		}
		int temp = mutant[pairPoint1];
		mutant[pairPoint1] = mutant[pairPoint2];
		mutant[pairPoint2] = temp;
		return mutant;
	}

	/**
	 * Mutation - chooses one element, places it in a different bin. Behaves like a partition problem.
	 */
	private int[] mutationSingleMove(int[] mutant) {
		int pairPoint1 = rand.nextInt(mutant.length);
		int newBin = mutant[pairPoint1];
		while (newBin == mutant[pairPoint1]) {
			newBin = rand.nextInt(binCount);
		}
		mutant[pairPoint1] = newBin; //swap to a different bucket
		return mutant;
	}

	/* ******************************* GA HELPER FUNCTIONS ******************************* */

	/**
	 * See top comment on chromosome structure, this tries to honor packing order if possible. If bin is full, best-first packing gets used.
	 * This method also contains a mutation operation - bin dumping.
	 */
	private int[] bestFitModified(int[] chromosome) {
		int[] bins = new int[binCount];
		ArrayList<Integer> cannotPackList = new ArrayList<Integer>();

		//honor chromosome listing when possible
		for (int i = 0; i < chromosome.length; i++) {
			//if this bin isnt full & won't be after placement & this chrom hasnt gone off the end 
			//(VERY RARE, only happens when binCount is reduced w/ elitism)
			if (chromosome[i] < bins.length && bins[chromosome[i]] + packageWeights[i] <= binSize) {
				bins[chromosome[i]] += packageWeights[i]; //put package into bin
			} else { //store those we can't honor
				chromosome[i] = -1;
				cannotPackList.add(new Integer(i)); //save all those packages whos original bins were full
			}
		}

		//best fit on all the rest whos bins were full up
		int tries = 0;
		for (int j = 0; j < cannotPackList.size(); j++) {
			int i = cannotPackList.get(j);
			int bestBinIndex = -1;

			while (bestBinIndex < 0 && tries < bins.length * 2) { //while we still havent found a bin to put this in and time isn't up

				int bestBinLeftOverSpace = binSize * 100; //how much space is left in the bin when the packages is placed in it (ideally 0 = full bin)

				for (int bin = 0; bin < bins.length; bin++) {
					//if this bin is closer to full than the best bin and not overfilled save it
					if ((binSize - (bins[bin] + packageWeights[i])) >= 0 && (binSize - (bins[bin] + packageWeights[i])) < bestBinLeftOverSpace) {
						bestBinIndex = bin;
						bestBinLeftOverSpace = binSize - (bins[bin] + packageWeights[i]);
					}
				}

				if (bestBinIndex < 0) { //if could not pack, dump out a bin (mutation)
					j = 0;
					int dumpBin = rand.nextInt(binCount);
					for (int binIndex = dumpBin; binIndex < binCount; binIndex++) {
						bins[binIndex] = 0;
						for (int chromIndex = 0; chromIndex < chromosome.length; chromIndex++) {
							if (chromosome[chromIndex] == binIndex) {
								chromosome[chromIndex] = -1;
								cannotPackList.add(new Integer(chromIndex));
							}
						}
					}

				}
				tries++;
			}

			if (tries >= bins.length * 2) {
				return null; //this chromosome is unsalvageable because it is impossible to fit the packages in the bins
			}


			bins[bestBinIndex] += packageWeights[i];
			chromosome[i] = bestBinIndex;
			cannotPackList.remove(j);
			j--;
		}

		return chromosome;
	}

	/**
	 * Inserts a set of chromosomes into a given population by removing the worst members and replacing them with children.
	 */
	private int[][] replaceLowestFitness(int[][] population, ArrayList<int[]> replacers) {
		float[] fitnesses = new float[population.length];
		for (int i = 0; i < population.length; i++) {
			fitnesses[i] = getFitness(population[i]);
		}

		//find worst chromosomes around
		ArrayList<Integer> worstFitnessesIndexes = new ArrayList<Integer>();
		float[] fitnessesSorted = fitnesses.clone();
		Arrays.sort(fitnessesSorted);
		for (int i = 0; i < population.length; i++) {
			for (int j = 0; j < replacers.size() + 1; j++) {
				if (fitnesses[i] <= fitnessesSorted[j]) {
					worstFitnessesIndexes.add(new Integer(i));
					break;
				}
			}
			if (worstFitnessesIndexes.size() > replacers.size() + 1) //only find as many bad chromosomes as we have replacements
				break;
		}

		for (int i = 0; i < replacers.size(); i++) {
			int[] replacer = replacers.get(i);
			population[worstFitnessesIndexes.get(i)] = replacer; //replace worst fitnesses with replacements
		}

		population[worstFitnessesIndexes.get(replacers.size())] = eliteChromosome; //make sure elite stays
		return population;
	}

	/**
	 * Checks population for elite chromosome, or returns old elite if there are no challengers.
	 */
	private int[] getElite(int[][] population) {
		float[] fitnesses = new float[population.length];

		for (int i = 0; i < population.length; i++) {
			fitnesses[i] = getFitness(population[i]);
			if (fitnesses[i] < eliteFitness) {
				System.out.println(getFitness(eliteChromosome) + " old" + getFitness(population[i]) + " new");
				eliteFitness = fitnesses[i];
				eliteChromosome = population[i];
			}
		}
		return eliteChromosome;
	}

	/**
	 * Finds fitness of a chromosome - the number of bins (minimization). All chromosomes are fixed if infeasible using side-effects.
	 */
	private Integer getFitness(int[] chromosome) {
		int[] bins = new int[binCount]; //list of found bins
		for (int binIndex = 0; binIndex < bins.length; binIndex++) { //initialize list of bins to unfound
			bins[binIndex] = -1;
		}
		int[] newChrom = bestFitModified(chromosome.clone()); //fix if infeasible
		if (newChrom == null) { //this is an unfixable chromosome, inform caller
			return binCount * 2 + totalInfeasibility(chromosome);
		}

		for (int i = 0; i < newChrom.length; i++) {
			chromosome[i] = newChrom[i]; //copy over to take advantage of side-effects
		}
		int uniqueBins = 0;

		for (int i = 0; i < chromosome.length; i++) {
			if (bins[chromosome[i]] == -1) { //if we haven't seen this bin before track it and increment bin counter
				bins[chromosome[i]] = chromosome[i];
				uniqueBins++;
			} else
				continue;
		}

		return uniqueBins;
	}

	/**
	 * If a chromosome cannot be fixed in a reasonable time we mark it as unfit by finding how much it overfills each bin.
	 */
	private int totalInfeasibility(int[] chromosome) {
		int infeasibility = 0;
		int[] bins = new int[binCount];
		for (int i = 0; i < chromosome.length; i++) {
			if (chromosome[i] >= bins.length) {
				infeasibility += 10;
				continue;
			}
			try {
				bins[chromosome[i]] += packageWeights[i]; //put package into bin
			} catch (Exception e) {
				System.out.println("noooo");
			}

		}

		for (int i = 0; i < bins.length; i++) {
			if (bins[i] > binSize)
				infeasibility += bins[i] - binSize;
		}
		return infeasibility;
	}

	/**
	 * Tests to see if problem can be solved with half as many bins
	 */
	private boolean canReduceBins() {
		int oldBinCount = binCount;
		boolean canReduce = false;
		binCount = binCount / 2;
		int[] chromosome = eliteChromosome.clone();

		if (getFitness(chromosome) <= binCount)
			canReduce = true;
		else
			canReduce = false;
		binCount = oldBinCount; //no side effects
		return canReduce;
	}

	/* ******************************* SA HELPER FUNCTIONS ******************************* */

	/**
	 * SA Fitness. SA does not repair infeasibles, so we try to take
	 * infeasibility into account here.
	 */
	public int h(int[] solution) {
		int[] bins = new int[binCount];
		int nonEmptyBins = 0;
		int overFilledBins = 0;
		int overFillWeight = 0;

		for (int i = 0; i < solution.length; i++) {
			int bin = solution[i];
			int weight = packageWeights[i];
			bins[bin] += weight;
		}

		for (int i = 0; i < bins.length; i++) {
			if (bins[i] != 0) {
				nonEmptyBins++;
			}
			if (bins[i] > binSize) {
				overFilledBins++;
				overFillWeight += (bins[i] - binSize);
			}
		}
		if (overFilledBins > 0)
			return (nonEmptyBins + binSize * (overFillWeight));
		else
			return nonEmptyBins;
	}
}
