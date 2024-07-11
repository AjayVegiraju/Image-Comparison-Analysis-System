import java.awt.*;
import javax.swing.*;
import java.awt.event.*;
import javax.imageio.ImageIO;
import java.io.*;
import java.util.*;
import java.awt.image.BufferedImage;
import java.nio.charset.StandardCharsets;
import java.util.List;

public class ImageBrowser extends JFrame {

    private final JPanel thumbnailPanel;
    private final JLabel displayLabel;

    public int weight = 1 / 89;
    Comparator<File> comparator = (file1, file2) -> {
        int number1 = extractImageNumber(file1.getName());
        int number2 = extractImageNumber(file2.getName());
        return Integer.compare(number1, number2);
    };

    private final ArrayList<File> imageFiles = new ArrayList<>();

    // change the folderPath variable to the absolute path of the GUI Demo directory
    // in your system
    static String folderPath = "/Users/ajayvegiraju/Documents/CSS484/GUI Demo 3/images";

    private static Integer[][] intensityMatrix;

    File currentFile;
    // Stores the distance between image files (Intensity)
    Double[][] manhattanDistances;
    // Stores the distance between image files (Color Code)
    Double[][] colorManhattanDistances;
    Double[][] weightedManhattanDistances;
    File[] curSortedFiles;
    File[] curColorSortedFiles;
    File[] WeightSortedFiles;
    private final List<File> relevantFiles = new ArrayList<>();
    private static Double[][] rfMatrix;
    Map<Integer, Double[]> imageFeatures = new HashMap<>();
    private static double[] weights;

    Double[][] rows;
    Double[] row;

    private static Integer[][] colorCodeMatrix;
    Double[] columnAverages;
    Double[] columnStdDevs;

    public ImageBrowser() {
        setTitle("Image Browser");
        setSize(1000, 600);
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        loadImageFiles(folderPath);
        int rows = imageFiles.size() / 3;

        GridLayout gridLayout = new GridLayout(rows, 3); // creating GridLayout object
        gridLayout.setHgap(10); // setting horizontal gap to 10
        gridLayout.setVgap(10);// setting vertical gap to 10

        thumbnailPanel = new JPanel(gridLayout);

        JScrollPane scrollPane = new JScrollPane(thumbnailPanel);
        scrollPane.getVerticalScrollBar().setUnitIncrement(10);
        scrollPane.getHorizontalScrollBar().setUnitIncrement(10);

        displayLabel = new JLabel();
        displayLabel.setHorizontalAlignment(JLabel.CENTER);

        JButton colorCodeBasedRetrieval = new JButton("Retrieve by Color Code");
        colorCodeBasedRetrieval.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                colorCodeCompareImages(); // This should populate curColorSortedFiles
                displayThumbnails(Arrays.asList(curColorSortedFiles)); // This should display the thumbnails of
                                                                       // curSortedFiles
            }
        });

        JButton intensityBasedRetrieval = new JButton("Retrieve by Intensity");
        intensityBasedRetrieval.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                compareImages(); // This should populate curSortedFiles
                displayThumbnails(Arrays.asList(curSortedFiles)); // This should display the thumbnails of
                                                                  // curSortedFiles
            }
        });

        JButton rfBasedRetrieval = new JButton("Retrieve by Relevance Feedback");
        rfBasedRetrieval.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {

                getRelevantFiles();
                rfMatrix();
                loadRfComparisons();
                rfCompareImages();
                displayThumbnails(Arrays.asList(WeightSortedFiles));
            }
        });

        JPanel buttonsPanel = new JPanel();
        buttonsPanel.setLayout(new BoxLayout(buttonsPanel, BoxLayout.Y_AXIS));
        buttonsPanel.add(colorCodeBasedRetrieval);
        buttonsPanel.add(rfBasedRetrieval);
        buttonsPanel.add(intensityBasedRetrieval);

        JPanel rightPanel = new JPanel(new BorderLayout());
        rightPanel.add(displayLabel, BorderLayout.CENTER);
        rightPanel.add(buttonsPanel, BorderLayout.SOUTH);

        JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, scrollPane, rightPanel);
        splitPane.setDividerLocation(500);
        add(splitPane);
        displayThumbnails(imageFiles);

    }

    void populateWeights() {
        if (weights == null) {
            weights = new double[89];

            for (int i = 0; i < weights.length; i++) {
                weights[i] = 1 / 89;
            }
        }
    }

    public void getRelevantFiles() {
        System.out.println("got here to the relevant files");
        System.out.println("this is the size of the relevantFiles list: " + relevantFiles.size());
        for (File file : relevantFiles) {
            System.out.println("This is a relevant file: " + file.getName());
        }
    }

    private void displayThumbnails(List<File> filesToDisplay) {

        thumbnailPanel.removeAll(); // Clear existing thumbnails

        if (filesToDisplay == null || filesToDisplay.isEmpty()) {
            filesToDisplay = imageFiles;
        }

        for (File file : filesToDisplay) {
            ImageIcon imageIcon = new ImageIcon(new ImageIcon(file.getAbsolutePath()).getImage()
                    .getScaledInstance(100, 100, Image.SCALE_DEFAULT));
            JLabel label = new JLabel(imageIcon);
            label.addMouseListener(new MouseAdapter() {
                @Override
                public void mouseClicked(MouseEvent e) {
                    relevantFiles.clear();
                    currentFile = file;
                    displayImage(file);
                    relevantFiles.add(currentFile);
                    System.out.println(currentFile.getName());
                }
            });

            JPanel container = new JPanel();
            container.setLayout(new BoxLayout(container, BoxLayout.PAGE_AXIS));
            container.setPreferredSize(new Dimension(110, 130));
            container.add(label);

            JCheckBox checkBox = new JCheckBox("Relevant", relevantFiles.contains(file));
            checkBox.addItemListener(e -> {
                if (e.getStateChange() == ItemEvent.SELECTED) {
                    // Add to relevant files if not already present
                    if (!relevantFiles.contains(file)) {
                        relevantFiles.add(file);
                    }
                } else if (e.getStateChange() == ItemEvent.DESELECTED) {
                    // Remove from relevant files
                    relevantFiles.remove(file);
                }
            });

            container.add(checkBox);
            thumbnailPanel.add(container);
        }

        thumbnailPanel.revalidate();
        thumbnailPanel.repaint();
    }

    /*
     * to get the int value of the image file name for example 1.jpg will return 1
     * so i can use it to
     * sort the image files arraylist using a custom comparator
     */
    private static int extractImageNumber(String filename) {
        try {
            return Integer.parseInt(filename.substring(0, filename.lastIndexOf('.')));
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Invalid image filename format: " + filename);
        }
    }

    /* Loads all od the image files and populates an arrayList called imageFiles */
    private void loadImageFiles(String folderPath) {
        File folder = new File(folderPath);
        if (!folder.exists() || !folder.isDirectory()) {
            System.out.println("Invalid folder path: " + folderPath);
            return;
        }
        File[] files = folder.listFiles((dir, name) -> name.toLowerCase().endsWith("jpg"));
        if (files != null) {
            for (File file : files) {
                if (file.isFile()) {
                    imageFiles.add(file);
                }
            }
            Collections.sort(imageFiles, comparator);
            int numOfImages = imageFiles.size();
            intensityMatrix = new Integer[numOfImages][25];
            colorCodeMatrix = new Integer[numOfImages][64];

            for (int i = 0; i < intensityMatrix.length; i++) {
                for (int j = 0; j < intensityMatrix[i].length; j++) {
                    intensityMatrix[i][j] = 0;
                }
            }

            for (int i = 0; i < colorCodeMatrix.length; i++) {
                for (int j = 0; j < colorCodeMatrix[i].length; j++) {
                    colorCodeMatrix[i][j] = 0;
                }
            }
        } else {
            System.out.println("No JPG files found in: " + folderPath);
        }
    }

    /* This populates the color code feature matrix */
    public void loadColorCode(String folderPath) {

        if (imageFiles != null) {

            try (FileOutputStream output = new FileOutputStream("ColorCodeOutput.txt");
                    Writer colorOutputStreamWriter = new OutputStreamWriter(output, StandardCharsets.UTF_8)) {
                int imageCount = 0;
                for (File file : imageFiles) {

                    try {
                        BufferedImage image = ImageIO.read(file);
                        int height = image.getHeight();
                        int width = image.getWidth();

                        colorOutputStreamWriter.write("Image " + file.getName() + "\n\n");

                        Integer[] bin = new Integer[64];
                        int pixelCount = 0;

                        for (int i = 0; i < bin.length; i++) {
                            bin[i] = 0;
                        }

                        for (int i = 0; i < height; i++) {

                            for (int j = 0; j < width; j++) {

                                int rgb = image.getRGB(j, i);
                                int red = (rgb >> 16) & 0xFF;
                                int green = (rgb >> 8) & 0xFF;
                                int blue = rgb & 0xFF;

                                bin[rgbToDecimal(red, green, blue)]++;
                                pixelCount++;
                            }
                        }

                        System.arraycopy(bin, 0, colorCodeMatrix[imageCount], 0, 64);

                        for (int i = 0; i < bin.length; i++) {
                            colorOutputStreamWriter.write("bin " + i + " = " + bin[i] + "\n");
                        }

                        colorOutputStreamWriter.write("\n" +
                                "\n-------------------------------------------------------------\n");
                        imageCount++;

                    } catch (IOException e) {
                        e.printStackTrace();
                    }

                }

            } catch (IOException e) {
                System.out.println(e);
            }
        }
    }

    /* a helper method for loadColorCode to return an int value for and rgb value */
    public int rgbToDecimal(int red, int green, int blue) {
        int rr = red >> 6;
        int gg = green >> 6;
        int bb = blue >> 6;
        return (rr << 4) | (gg << 2) | bb;
    }

    public void rfMatrix() {
        if(colorCodeMatrix == null ) loadColorCode(folderPath);
        if(intensityMatrix == null) loadIntensity(folderPath);
        int size = colorCodeMatrix[0].length + intensityMatrix[0].length;
        System.out.println("This is the size of the intensity length: " + intensityMatrix[0].length);
        System.out.println("This is the size of the Color code length: " + colorCodeMatrix[0].length);

        System.out.println("This is the size of the rf matrix: " + size);
        rfMatrix = new Double[100][size];
        int rfIndex = 0;
        for (int i = 0; i < 100; i++) {
            for (int j = 0; j < size; j++) {
                if (j < intensityMatrix[i].length) {
                    double value = (double) (intensityMatrix[i][j]) / 98304;

                    rfMatrix[i][j] = Math.round(value * 100000.0) / 100000.0;
                } else {

                    double value = (double) (colorCodeMatrix[i][rfIndex++]) / 98304;
                    rfMatrix[i][j] = Math.round(value * 100000.0) / 100000.0;
                }
            }

            rfIndex = 0;

        }
        normalizeRfMatrix(rfMatrix);
        printRfMatrix();
    }

    public void normalizeRfMatrix(Double[][] matrix) {
        int numRows = matrix.length;
        int numCols = matrix[0].length;

        Double[] columnAverages = new Double[numCols];
        Double[] columnStdDevs = new Double[numCols];
        Double minNonZeroStdDev = Double.MAX_VALUE;

        // Calculate averages and standard deviations for each column
        for (int j = 0; j < numCols; j++) {
            double sum = 0;
            for (int i = 0; i < numRows; i++) {
                sum += matrix[i][j];
            }
            columnAverages[j] = sum / numRows;
        }

        // Calculate standard deviations and find the minimum non-zero standard
        // deviation
        for (int j = 0; j < numCols; j++) {
            double sumOfSquares = 0;
            for (int i = 0; i < numRows; i++) {
                double deviation = matrix[i][j] - columnAverages[j];
                sumOfSquares += deviation * deviation;
            }
            columnStdDevs[j] = Math.sqrt(sumOfSquares / (numRows - 1));
            if (columnStdDevs[j] > 0 && columnStdDevs[j] < minNonZeroStdDev) {
                minNonZeroStdDev = columnStdDevs[j];
            }
        }

        // Normalize each value in the matrix
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                if (columnStdDevs[j] != 0) {
                    matrix[i][j] = (matrix[i][j] - columnAverages[j]) / columnStdDevs[j];
                } else {
                    // Handle the case where the standard deviation is 0.

                    matrix[i][j] = 0.0;
                }
            }
        }
    }

    public void getRows() {
        rows = new Double[relevantFiles.size()][89];
        for (int i = 0; i < rows.length; i++) {
            row = rfMatrix[extractImageNumber(relevantFiles.get(i).getName()) - 1];
            for (int j = 0; j < rows[i].length; j++) {
                rows[i][j] = row[j];
            }

        }
    }

    public void printRfMatrix() {
        try (FileOutputStream output = new FileOutputStream("rfMatrix.txt");
                Writer outputStreamWriter = new OutputStreamWriter(output, StandardCharsets.UTF_8)) {
            for (int i = 0; i < rfMatrix.length; i++) {
                for (int j = 0; j < rfMatrix[i].length; j++) {
                    outputStreamWriter.write(String.format("%7.4f | ", rfMatrix[i][j]));
                }

                outputStreamWriter.write("\n");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    /* this method populates the intensity matrix */
    public void loadIntensity(String folderPath) {

        if (imageFiles != null) {

            try (FileOutputStream output = new FileOutputStream("output.txt");
                    Writer outputStreamWriter = new OutputStreamWriter(output, StandardCharsets.UTF_8)) {
                int imageCount = 0;
                for (File file : imageFiles) {

                    try {
                        BufferedImage image = ImageIO.read(file);
                        int height = image.getHeight();
                        int width = image.getWidth();

                        outputStreamWriter.write("Image " + file.getName() + "\n\n");

                        Integer[] bin = new Integer[25];
                        int pixelCount = 0;

                        for (int i = 0; i < bin.length; i++) {
                            bin[i] = 0;
                        }

                        for (int i = 0; i < height; i++) {

                            for (int j = 0; j < width; j++) {

                                int rgb = image.getRGB(j, i);
                                int red = (rgb >> 16) & 0xFF;
                                int green = (rgb >> 8) & 0xFF;
                                int blue = rgb & 0xFF;

                                double intensity = (0.299 * red) + (0.587 * green) + (0.114 * blue);
                                // gets the index to increment
                                int histbin = histBin(intensity);
                                bin[histbin]++;

                                pixelCount++;
                            }
                        }

                        System.arraycopy(bin, 0, intensityMatrix[imageCount], 0, 25);

                        for (int i = 0; i < bin.length; i++) {
                            outputStreamWriter.write("bin " + i + " = " + bin[i] + "\n");
                        }

                        outputStreamWriter.write("\n" +
                                "\n-------------------------------------------------------------\n");
                        imageCount++;

                    } catch (IOException e) {
                        e.printStackTrace();
                    }

                }
            } catch (IOException e) {
                System.out.println(e);
            }
        }
    }



    public void loadRfComparisons() {
        getRows();
        populateWeights();

        if (relevantFiles.isEmpty()) {
            // If no checkboxes are checked, set all weights to 1/number of features.
            Arrays.fill(weights, 1.0 / 89);
            for (int i = 0; i < weights.length; i++) {
                System.out.println("Default weight at index " + i + ": " + weights[i]);
            }
        } else {
            double[] avgs = calculateColumnMeans(rows);
            double[] std = calculateColumnStandardDeviations(rows, avgs);

            double stdMin = Double.MAX_VALUE;

            // Find the minimum non-zero standard deviation
            for (double d : std) {
                if (d != 0 && d < stdMin) {
                    stdMin = d;
                }
            }

            // Print the minimum non-zero std deviation
            System.out.println("Minimum non-zero std deviation: " + stdMin);

            // Adjust weights based on standard deviations
            for (int i = 0; i < weights.length; i++) {
                if (std[i] != 0) {
                    weights[i] = 1 / std[i];
                } else if (avgs[i] != 0) { // Check if average is non-zero for this column
                    if (stdMin != 0) { // Check to avoid division by zero
                        weights[i] = 1 / (0.5 * stdMin); // Set weight to half of the minimum non-zero standard deviation
                    } else {
                        weights[i] = 0; // If all std are zero, then set weight to zero
                    }
                } else {
                    weights[i] = 0; // If std[i] is zero and avgs[i] is zero
                }
                // Print the weight after adjustment
                System.out.println("Weight at index " + i + ": " + weights[i]);
            }

            // Normalize the weights so that they sum to 1
            double sumOfWeights = 0.0;
            for (double weight : weights) {
                sumOfWeights += weight;
            }
            System.out.println("Sum of weights before normalization: " + sumOfWeights);

            for (int i = 0; i < weights.length; i++) {
                weights[i] /= sumOfWeights; // Divide each weight by the sum of all weights to normalize
                // Print the normalized weight
                System.out.println("Normalized weight at index " + i + ": " + weights[i]);
            }
        }

        // Perform weighted comparisons
        weightedComparisons(weights);

        // Write the weighted Manhattan distances to a file
        try (FileOutputStream ot = new FileOutputStream("weightManDist.txt");
             OutputStreamWriter writer = new OutputStreamWriter(ot)) {
            for (int i = 0; i < weightedManhattanDistances.length; i++) {
                for (int j = 0; j < weightedManhattanDistances[i].length; j++) {
                    writer.write(weightedManhattanDistances[i][j] + ", ");
                }
                writer.write("\n");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    // 1. Method to calculate the mean of each column in the matrix
    private double[] calculateColumnMeans(Double[][] matrix) {
        double[] means = new double[matrix[0].length];
        for (int j = 0; j < matrix[0].length; j++) {
            double sum = 0;
            for (Double[] doubles : matrix) {
                sum += doubles[j];
            }
            means[j] = sum / matrix.length;
        }
        return means;
    }

    // 2. Method to calculate the standard deviation of each column in the matrix
    private double[] calculateColumnStandardDeviations(Double[][] matrix, double[] means) {
        double[] stdDevs = new double[matrix[0].length];
        for (int j = 0; j < matrix[0].length; j++) {
            double sumOfSquaredDiffs = 0;
            for (Double[] doubles : matrix) {
                sumOfSquaredDiffs += Math.pow(doubles[j] - means[j], 2);
            }
            stdDevs[j] = Math.sqrt(sumOfSquaredDiffs / matrix.length);
        }
        return stdDevs;
    }

    private double calculateWeightedManhattanDistance(Double[] vectorA, Double[] vectorB, double[] weights) {
        if (vectorA.length != vectorB.length || vectorA.length != weights.length) {
            throw new IllegalArgumentException("All arrays must be of the same length.");
        }

        double sum = 0.0;

        for (int i = 0; i < vectorA.length; i++) {
            double normalizedDifference = Math.abs(vectorA[i] - vectorB[i]);
            sum += normalizedDifference * weights[i];

            // Debugging output
            // System.out.println("Index: " + i + ", vectorA: " + vectorA[i] + ", vectorB: "
            // + vectorB[i] + ", weight: " + weights[i] + ", normalizedDifference: " +
            // normalizedDifference + ", partial sum: " + sum);
        }

        // double roundedSum = Math.round(1000.0 * sum) / 1000.0;

        return sum;
    }

    public void rfCompareImages() {
        weightedComparisons(weights);

        // gets the row number of the specific image in the colorCode manhattan distance
        // matrix
        int row = extractImageNumber(currentFile.getName()) - 1;

        // copy of the specific row retrieved by the previous line
        Double[] imageCount = new Double[imageFiles.size()];

        // copies the values of the original row into the duplicate array imageCount
        System.arraycopy(weightedManhattanDistances[row], 0, imageCount, 0, imageFiles.size());

        // sorts the duplicate array image count
        Arrays.sort(imageCount);

        // defines the file array that stores the files in sorted order
        WeightSortedFiles = new File[imageCount.length];
        Double[] original = weightedManhattanDistances[row];

        for (int i = 0; i < WeightSortedFiles.length; i++) {
            // gets the original index of the file from the sorted list of distances
            int originalIndex = findInOriginal(imageCount[i], original);
            WeightSortedFiles[i] = imageFiles.get(originalIndex);
        }
    }

    private void weightedComparisons(double[] weights) {
        int numOfImages = imageFiles.size();
        weightedManhattanDistances = new Double[numOfImages][numOfImages];

        loadIntensity(folderPath); // Assuming this method loads your intensityMatrix

        for (int i = 0; i < numOfImages; i++) {
            for (int j = i; j < numOfImages; j++) { // Adjusting the start index of j
                if (i == j) {
                    weightedManhattanDistances[i][j] = 0.0;
                } else {
                    if (rfMatrix[i].length == rfMatrix[j].length) {

                        double distance = calculateWeightedManhattanDistance(rfMatrix[i], rfMatrix[j], weights);

                        // Since distance[i][j] == distance[j][i]
                        weightedManhattanDistances[i][j] = distance;
                        weightedManhattanDistances[j][i] = distance;
                    }

                }
            }
        }

        // Assuming this method prints your matrix to a file
    }

    private void comparisons() {
        int numOfImages = imageFiles.size();
        manhattanDistances = new Double[numOfImages][numOfImages];

        loadIntensity(folderPath);
        for (int i = 0; i < numOfImages; i++) {
            for (int j = i; j < numOfImages; j++) { // Adjusting the start index of j

                if (i == j) {
                    manhattanDistances[i][j] = 0.0;
                } else {
                    double distance = calculateManhattanDistance(intensityMatrix[i], intensityMatrix[j]);
                    // Since distance[i][j] == distance[j][i]
                    manhattanDistances[i][j] = distance;
                    manhattanDistances[j][i] = distance;
                }
            }
        }

        printMatrixToFile(manhattanDistances);

    }

    /* helper method that find the index of the original file */
    public int findInOriginal(Double distance, Double[] distanceValues) {
        int index = -1;

        for (int i = 0; i < imageFiles.size(); i++) {
            if (distanceValues != null) {
                if (distance == distanceValues[i]) {
                    index = i;
                }
            }

        }

        return index;
    }

    /* Prints the Intensity matrix Manhattan distances to a file */
    public static void printMatrixToFile(Double[][] matrix) {
        try (FileOutputStream out = new FileOutputStream("manhattanDistanceOutput.txt");
                Writer outputStreamWriter = new OutputStreamWriter(out, StandardCharsets.UTF_8)) {

            for (int i = 0; i < matrix.length; i++) {
                for (int j = 0; j < matrix[i].length; j++) {
                    // Format number with 3 decimal places and at least 7 characters wide
                    outputStreamWriter.write(String.format("%7.3f | ", matrix[i][j]));
                }
                outputStreamWriter.write("\n");
            }
            System.out.println(matrix.length);

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /* method to determine which index (bin) to put the intensity value in */
    private static int histBin(Double intensity) {
        if (intensity >= 0 && intensity <= 10) {
            return 0;
        } else if (intensity > 10 && intensity <= 20) {
            return 1;
        } else if (intensity > 20 && intensity <= 30) {
            return 2;
        } else if (intensity > 30 && intensity <= 40) {
            return 3;
        } else if (intensity > 40 && intensity <= 50) {
            return 4;
        } else if (intensity > 50 && intensity <= 60) {
            return 5;
        } else if (intensity > 60 && intensity <= 70) {
            return 6;
        } else if (intensity > 70 && intensity <= 80) {
            return 7;
        } else if (intensity > 80 && intensity <= 90) {
            return 8;
        } else if (intensity > 90 && intensity <= 100) {
            return 9;
        } else if (intensity > 100 && intensity <= 110) {
            return 10;
        } else if (intensity > 110 && intensity <= 120) {
            return 11;
        } else if (intensity > 120 && intensity <= 130) {
            return 12;
        } else if (intensity > 130 && intensity <= 140) {
            return 13;
        } else if (intensity > 140 && intensity <= 150) {
            return 14;
        } else if (intensity > 150 && intensity <= 160) {
            return 15;
        } else if (intensity > 160 && intensity <= 170) {
            return 16;
        } else if (intensity > 170 && intensity <= 180) {
            return 17;
        } else if (intensity > 180 && intensity <= 190) {
            return 18;
        } else if (intensity > 190 && intensity <= 200) {
            return 19;
        } else if (intensity > 200 && intensity <= 210) {
            return 20;
        } else if (intensity > 210 && intensity <= 220) {
            return 21;
        } else if (intensity > 220 && intensity <= 230) {
            return 22;
        } else if (intensity > 230 && intensity <= 240) {
            return 23;
        } else {
            return 24;
        }
    }

    /* populates the intensity manhattan distance matrix with values */

    /* populates the color code manhattan distance matrix with values */
    private void colorComparisons() {
        int numOfImages = imageFiles.size();
        colorManhattanDistances = new Double[numOfImages][numOfImages];

        loadColorCode(folderPath);

        for (int i = 0; i < numOfImages; i++) {
            for (int j = i; j < numOfImages; j++) { // Adjusting the start index of j

                if (i == j) {
                    colorManhattanDistances[i][j] = 0.0;
                } else {
                    double distance = calculateManhattanDistance(colorCodeMatrix[i], colorCodeMatrix[j]);
                    // Since distance[i][j] == distance[j][i]
                    colorManhattanDistances[i][j] = distance;
                    colorManhattanDistances[j][i] = distance;
                }
            }
        }

        printColorMatrixToFile(colorManhattanDistances);
    }

    /*
     * populates the curSorted files array with a ranked list of files that are
     * similar to the current selected file
     * (Ranked by Intensity)
     */
    public void compareImages() {
        comparisons();
        // gets the row number of the specific image in the intensity manhattan distance
        // matrix
        int row = extractImageNumber(currentFile.getName()) - 1;
        System.out.println(row);

        // copy of the specific row retrieved by the previous line
        Double[] imageCount = new Double[imageFiles.size()];

        // copies the values of the original row into the duplicate array imageCount
        for (int i = 0; i < imageFiles.size(); i++) {
            if (manhattanDistances[row][i] == null) {
                imageCount[i] = 0.0;
            } else {
                imageCount[i] = manhattanDistances[row][i];
            }

        }

        // sorts the duplicate array image count
        Arrays.sort(imageCount);

        // defines the file array that stores the files in sorted order
        curSortedFiles = new File[imageCount.length];
        Double[] original = manhattanDistances[row];

        for (int i = 0; i < curSortedFiles.length; i++) {
            // gets the original index of the file from the sorted list of distances
            int originalIndex = findInOriginal(imageCount[i], original);
            curSortedFiles[i] = imageFiles.get(originalIndex);
        }

    }

    public static void printMatrix(Double[][] matrix) {

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                System.out.print(matrix[i][j] + " | ");
            }
            System.out.println();
        }

        System.out.println(matrix.length);
    }

    /*
     * Populates the curColorSortedFiles array with a sorted list of files that are
     * similar to the selected file
     * (Ranked by Color Code)
     */
    public void colorCodeCompareImages() {
        colorComparisons();

        // gets the row number of the specific image in the colorCode manhattan distance
        // matrix
        int row = extractImageNumber(currentFile.getName()) - 1;

        // copy of the specific row retrieved by the previous line
        Double[] imageCount = new Double[imageFiles.size()];

        // copies the values of the original row into the duplicate array imageCount
        System.arraycopy(colorManhattanDistances[row], 0, imageCount, 0, imageFiles.size());

        // sorts the duplicate array image count
        Arrays.sort(imageCount);

        // defines the file array that stores the files in sorted order
        curColorSortedFiles = new File[imageCount.length];
        Double[] original = colorManhattanDistances[row];

        for (int i = 0; i < curColorSortedFiles.length; i++) {
            // gets the original index of the file from the sorted list of distances
            int originalIndex = findInOriginal(imageCount[i], original);
            curColorSortedFiles[i] = imageFiles.get(originalIndex);
        }
    }

    /* prints the color code matrix to an external file */
    private void printColorMatrixToFile(Double[][] matrix) {
        try (FileOutputStream out = new FileOutputStream("colorMan.txt", true); // Appending, optional
                Writer outputStreamWriter = new OutputStreamWriter(out, StandardCharsets.UTF_8)) {

            if (matrix == null || matrix.length == 0) {
                System.out.println("Matrix is empty or null.");
                return;
            }

            for (int i = 0; i < matrix.length; i++) {
                if (matrix[i] == null) {
                    System.out.println("Row " + i + " is null.");
                    continue;
                }
                for (int j = 0; j < matrix[i].length; j++) {
                    if (matrix[i][j] == null) {
                        outputStreamWriter.write(String.format("%7s | ", "NULL"));
                    } else {
                        outputStreamWriter.write(String.format("%7.3f | ", matrix[i][j]));
                    }
                }
                outputStreamWriter.write("\n");
            }
            outputStreamWriter.flush(); // Explicitly flushing
            System.out.println("Matrix printed. Size: " + matrix.length);

        } catch (IOException e) {
            e.printStackTrace();
            System.out.println("An error occurred while writing to the file.");
        }
    }

    /*
     * helper method to calculate manhattan distance between two images' intensity
     * values
     */
    private double calculateManhattanDistance(Integer[] vectorA, Integer[] vectorB) {
        double sum = 0.0;
        int pixels = 256 * 384; // Assuming constant image size, adjust as per actual size

        for (int i = 0; i < vectorA.length; i++) {
            sum += Math.abs(((double) vectorA[i] / pixels) - ((double) vectorB[i] / pixels));
        }

        return Math.round(1000.0 * sum) / 1000.0;
    }



    private void displayImage(File file) {
        currentFile = file;
        if (file != null && file.exists()) {
            ImageIcon imageIcon = new ImageIcon(new ImageIcon(file.getAbsolutePath()).getImage()
                    .getScaledInstance(480, 360, Image.SCALE_DEFAULT));
            displayLabel.setIcon(imageIcon);

            // Set the label text to the file name
            displayLabel.setText(file.getName());
            // If you want the text to appear below the image, you might need to use HTML tags
            displayLabel.setVerticalTextPosition(SwingConstants.BOTTOM);
            displayLabel.setHorizontalTextPosition(SwingConstants.CENTER);
        } else {
            displayLabel.setIcon(null);
            displayLabel.setText("File does not exist.");
        }
    }


    public static void main(String[] args) {

        try {
            ImageBrowser a = new ImageBrowser();

            SwingUtilities.invokeLater(() -> new ImageBrowser().setVisible(true));
            a.colorCodeCompareImages();
            a.rfMatrix();

            System.out.println("rfMatrix Length" + rfMatrix.length + "rfMatrix Row Length: " + rfMatrix[0].length);
            a.compareImages();
            a.loadRfComparisons();
        } catch (NullPointerException e) {
            System.out.println("Please select an image from the left side of the GUI window (click on a thumbnail)");
        } catch (Exception e) {
            e.printStackTrace();
        }

    }
}
