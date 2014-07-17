
package genetic_algo;

/**
 * @author Jeremy Barnes
 * 
 * This class is solely for the Swing Interface because I don't like command-line applications.
 * No significant logic happens in this class, only some button pressing and file reading.
 * The actual Bin-Packing algorithm begins in Core.java at the runAlgorithm method.
 * 
 * Note - the Interface is ugly as sin because I don't care, this is about the GA.
 */

import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.JSpinner;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;


public class WindowInterface extends JFrame implements ActionListener, ChangeListener {

	Core computationEngine;

	/* *********** Control Panel Items *************** */
	public JComboBox<String> algoChooser;
	public JComboBox<String> selectionChooser;
	public JComboBox<String> crossoverChooser;
	public JComboBox<String> mutationChooser;
	public JSlider mutationRateSlider;
	public JLabel mutationValueLabel;
	public JSpinner generationCounter;
	public JTextField filePath;
	public JButton fileButton;
	public JLabel selLabel = new JLabel("Choose selection method:");
	public JLabel mutLabel = new JLabel("Choose mutation method and rate:");
	public JLabel crossLabel = new JLabel("Choose crossover method:");
	public JLabel genLabel = new JLabel("Choose number of generations:");

	/* *********** Output Panel Items **************** */
	public JButton startButton;
	public JTextArea outputTextbox;

	/* *********** Non-GUI components ************** */
	public File dataSet;
	public int binCount;
	public int packageCount;
	public int binSize;
	public int[] packages;

	//JFrame set up, ignore unless you're really into boring Swing coding
	public WindowInterface(Core computationEngine) {
		super("Bin Packer - Jeremy Barnes");
		this.setDefaultCloseOperation(EXIT_ON_CLOSE);
		super.setSize(430, 675);
		this.computationEngine = computationEngine;
		this.setLayout(new GridLayout(2, 1, 0, 10));

		/* *********** Set up control panel for GA (controlPanel) *************** */
		JPanel controlPanel = new JPanel();

		controlPanel.setLayout(new GridLayout(15, 2, 20, 2));

		algoChooser = new JComboBox<String>(new String[] { "Genetic Algorithm", "Simulated Annealing", "Foolish Hill Climbing" });
		algoChooser.addActionListener(this);

		selectionChooser = new JComboBox<String>(new String[] { "Roulette", "Rank" });
		selectionChooser.addActionListener(this);

		crossoverChooser = new JComboBox<String>(new String[] { "One Point Crossover", "Two Point Crossover" });
		crossoverChooser.addActionListener(this);

		mutationChooser = new JComboBox<String>(new String[] { "Pairwise Exchange", "Single Move" });
		mutationChooser.addActionListener(this);

		mutationRateSlider = new JSlider(0, 10000, 500);
		mutationRateSlider.addChangeListener(this);
		mutationValueLabel = new JLabel("5.0");

		generationCounter = new JSpinner(new SpinnerNumberModel(1, 0, Integer.MAX_VALUE, 1));
		generationCounter.addChangeListener(this);
		generationCounter.setVisible(true);

		controlPanel.add(new JLabel("Choose Algorithm:"));
		controlPanel.add(new JLabel(" "));

		controlPanel.add(algoChooser);
		controlPanel.add(new JLabel(" "));

		controlPanel.add(new JLabel(" "));
		controlPanel.add(new JLabel(" "));

		controlPanel.add(new JLabel(" "));
		controlPanel.add(new JLabel(" "));

		controlPanel.add(selLabel);
		controlPanel.add(crossLabel);

		controlPanel.add(selectionChooser);
		controlPanel.add(crossoverChooser);

		controlPanel.add(new JLabel(" ")); //spacers
		controlPanel.add(new JLabel(" "));

		controlPanel.add(mutLabel);
		controlPanel.add(mutationValueLabel);
		controlPanel.add(mutationChooser);
		controlPanel.add(mutationRateSlider);

		controlPanel.add(new JLabel(" ")); //spacers
		controlPanel.add(new JLabel(" "));

		controlPanel.add(genLabel);
		controlPanel.add(new JLabel(" "));
		controlPanel.add(generationCounter);
		generationCounter.setValue(new Integer(200));
		controlPanel.add(new JLabel(" "));

		controlPanel.add(new JLabel(" ")); //spacers
		controlPanel.add(new JLabel(" "));

		controlPanel.add(new JLabel("Load Data File"));
		controlPanel.add(new JLabel(" "));

		filePath = new JTextField();
		filePath.setEditable(false);
		fileButton = new JButton("Select Data File");
		fileButton.addActionListener(this);
		controlPanel.add(fileButton);
		controlPanel.add(filePath);

		this.add(controlPanel);

		/* *********** End Set up control panel for GA (contentPane) *************** */

		/* *********** Set up output panel (outputPanel) *************** */
		JPanel outputPanel = new JPanel(new GridLayout(2, 3, 0, 40));

		outputPanel.add(new JLabel("  ")); //startbutton and spacers
		startButton = new JButton("Start!");
		startButton.addActionListener(this);
		outputPanel.add(startButton);
		outputPanel.add(new JLabel("  "));

//		outputPanel.add(new JLabel("          GA Output:"));
//		outputTextbox = new JTextArea();
//		outputTextbox.setLineWrap(true);
//		outputPanel.add(outputTextbox);
//		outputPanel.add(new JLabel("  "));


		this.add(outputPanel);

		/* *********** End Set up output panel (outputPanel) *************** */

		this.setVisible(true);
	}

	@Override
	public void stateChanged(ChangeEvent arg0) {
		if (arg0.getSource().equals(mutationRateSlider)) {
			mutationValueLabel.setText(mutationRateSlider.getValue() / 100. + "");
		}
	}

	@Override
	public void actionPerformed(ActionEvent arg0) {
		if (arg0.getSource().equals(algoChooser)) { //choose algorithm and adjust interface
			if (algoChooser.getSelectedItem().equals("Genetic Algorithm")) {
				selectionChooser.setVisible(true);
				crossoverChooser.setVisible(true);
				mutationRateSlider.setVisible(true);
				mutationValueLabel.setVisible(true);
				mutLabel.setText("Choose mutation method and rate:");
				selLabel.setVisible(true);
				crossLabel.setVisible(true);
				genLabel.setText("Choose number of generations:");
				mutationChooser.setSelectedIndex(0);
				generationCounter.setValue(new Integer(200));
			} else if (algoChooser.getSelectedItem().equals("Simulated Annealing") || algoChooser.getSelectedItem().equals("Foolish Hill Climbing")) {
				selectionChooser.setVisible(false);
				crossoverChooser.setVisible(false);
				mutationRateSlider.setVisible(false);
				mutationValueLabel.setVisible(false);
				mutLabel.setText("Choose perturbation function:");
				selLabel.setVisible(false);
				crossLabel.setVisible(false);
				genLabel.setText("Choose number of outer iterations:");
				mutationChooser.setSelectedIndex(1);
				generationCounter.setValue(new Integer(70));
			}
		}


		if (arg0.getSource().equals(fileButton)) { //choosing data file
			JFileChooser jfc = new JFileChooser();
			if (jfc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION) {
				dataSet = jfc.getSelectedFile();
				filePath.setText(dataSet.getAbsolutePath());
				setUpFileParams(dataSet);
				return;
			}
		}

		if (arg0.getSource().equals(startButton)) { //START
			String selection = (String) selectionChooser.getSelectedItem();
			String mutation = (String) mutationChooser.getSelectedItem();
			String crossover = (String) crossoverChooser.getSelectedItem();
			float mutationRate = mutationRateSlider.getValue() / 10000.0f;
			int generations = Integer.parseInt(generationCounter.getValue().toString());
			if (packages == null || packages.length != packageCount || packageCount == 0 || binCount == 0 || binSize == 0 || generations < 1) {
				JOptionPane.showMessageDialog(null, "Invalid Parameters Chosen - check control panel and data file!");
				return;
			}
			if (algoChooser.getSelectedItem().equals("Simulated Annealing")) {
				computationEngine.runSimulatedAnnealing(mutation, generations, packages, packageCount, binCount, binSize, 8, false);
			} else if (algoChooser.getSelectedItem().equals("Genetic Algorithm")) {
				computationEngine.runGeneticAlgorithm(selection, mutation, mutationRate, crossover, generations, packages, packageCount, binCount, binSize);
			} else {
				computationEngine.runSimulatedAnnealing(mutation, generations, packages, packageCount, binCount, binSize, 8, true);
			}
			
			return;
		}
	}

	/**
	 * Side effect method loads data set, bin count, package count, and bin size from data set TXT file. Lots of string parsing.
	 * @param dataSet - TXT file with data
	 */
	public void setUpFileParams(File dataSet) {
		try {
			BufferedReader fileIn = new BufferedReader(new FileReader(dataSet));
			String currentLine = null;
			while ((currentLine = fileIn.readLine()) != null) { //parse text file
				if (currentLine.startsWith("//") || currentLine.length() < 2)
					continue; //allow whole-line comments

				if (currentLine.startsWith("Bin Count:")) {
					String bins = currentLine.substring(10).replace(" ", "");
					this.binCount = Integer.parseInt(bins);
				}
				if (currentLine.startsWith("Package Count:")) {
					String packages = currentLine.substring(14).replace(" ", "");
					this.packageCount = Integer.parseInt(packages);
				}
				if (currentLine.startsWith("Bin Size:")) {
					String binsize = currentLine.substring(9).replace(" ", "");
					this.binSize = Integer.parseInt(binsize);

				}
				if (currentLine.startsWith("Packages:")) {
					String packageList = currentLine.substring(9).replace(" ", "");
					String[] packageStringList = packageList.split(",");
					packages = new int[packageStringList.length];
					for (int i = 0; i < packageStringList.length; i++) {
						packages[i] = Integer.parseInt(packageStringList[i]);
					}

				}
			}//end while
			fileIn.close();
			//check for invalid parameters
			if (packages == null || packages.length != packageCount || packageCount == 0 || binCount == 0 || binSize == 0) {
				throw new Exception();
			}

		} catch (Exception e) {
			JOptionPane.showMessageDialog(null, "Invalid Data File Chosen!");
		}
	}

}
