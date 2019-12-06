package views;

import filtres.Filter_3D;
import filtres.FilterType_3D;
import filtres.Gabor3DFilter_3D;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.ImageWindow;
import ij.gui.StackWindow;
import ij.io.OpenDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

import java.util.ArrayList;
import java.util.Scanner;

import process.Double2DArray_3D;
import process.Double3DArray_3D;
import process.Vsnr_denoiser_3D;

import java.awt.Font;
import java.io.File;

/**
 * 
 * @author benjamin.font
 *
 * Plugin developer : Benjamin Font, Léo Mouly
 * Algorithm : Jérôme Fehrenbach, Pierre Weiss, Corinne Lorenzo 
 *
 *         This plugin is used for denoise 3D stacks. To denoise 2d images, use
 *         the VSNR_2D plugin.
 */
public class Vsnr_Plugin_3D implements PlugInFilter {

	// Used for the GUI and the filters creation
	private String inputMethod;
	private String filterType;
	private ImagePlus image;
	private ImageStack stack;
	private ImagePlus noiseImage;
	private int size;
	private int width;
	private int height;
	private String type;
	private ArrayList<Filter_3D> filterList = new ArrayList<Filter_3D>();

	// Filters parameters
	private double noiseLevel = 1;
	private double sigmax = 10;
	private double sigmay = 10;
	private double sigmaz = 10;
	private double thetax = 0;
	private double thetay = 0;
	private double thetaz = 0;

	// Denoising parameters
	private String algorithm;
	private double nit = 50;

	// error for the text file
	private boolean errorInTxtFile = false;

	/**
	 * Run function
	 */
	@Override
	public void run(ImageProcessor ip) {
		if (ip.equals(null)) {
			GenericDialog g = new GenericDialog("Open an image please !");
			g.setLocationRelativeTo(null);
			g.pack();
			g.setVisible(true);
		}
		if (this.image.getImageStack().getSize() == 1) {
			GenericDialog g = new GenericDialog("Error");
			g.addMessage("Use VSNR_Macro_2D for images !");
			g.setLocationRelativeTo(null);
			g.pack();
			g.setVisible(true);
		} else {
			if (this.askMethod()) {
				if (this.inputMethod.equals("Graphic User Interface (GUI)")) {
					if (this.askFilter()) {
						if (this.askParameters()) {
							// the user want to add more than 1 filter
							while (this.askOtherFilter()) {
								if (this.askFilter()) {
									this.askParameters();
								}
							}
							if (this.askIterationsAndAlgorithm()) {
								printParameters();
								// the algorithm is TV
								if (this.algorithm.equals("TV")) {
									StackWindow wimg = new StackWindow(
											this.image);
									wimg.setTitle("Image");
									wimg.show();

									StackWindow wtv3d = new StackWindow(
											denoiseTV_3D());
									wtv3d.setTitle("Denoised Image - TV");
									wtv3d.show();

									StackWindow wflt = new StackWindow(
											this.noiseImage);
									wflt.setTitle("Filtre");
									wflt.show();
								}
								// the algorithm is H1
								else {
									StackWindow wimg = new StackWindow(
											this.image);
									wimg.setTitle("Image");
									wimg.show();

									StackWindow wh13d = new StackWindow(
											denoiseH1_3D());
									wh13d.setTitle("Denoised Image - H1");
									wh13d.show();

									StackWindow wflt = new StackWindow(
											this.noiseImage);
									wflt.setTitle("Filtre");
									wflt.show();
								}
							}
						}
					}
				}
				// the user want to use texts files (.txt)
				else {
					// we read the file
					if (this.readFile()) {

						// the algorithm is TV
						if (this.algorithm.equals("TV")) {
							ImageWindow wimg = new StackWindow(this.image);
							wimg.setTitle("Image");
							wimg.show();

							ImageWindow wh13d = new StackWindow(denoiseTV_3D());
							wh13d.setTitle("Denoised Image - TV");
							wh13d.show();

							ImageWindow wflt = new StackWindow(this.noiseImage);
							wflt.setTitle("Filtre");
							wflt.show();
						} else {
							// the algorithm is H1
							ImageWindow wimg = new StackWindow(this.image);
							wimg.setTitle("Image");
							wimg.show();

							ImageWindow wh13d = new StackWindow(denoiseH1_3D());
							wh13d.setTitle("Denoised Image - H1");
							wh13d.show();

							ImageWindow wflt = new StackWindow(this.noiseImage);
							wflt.setTitle("Filtre");
							wflt.show();
						}
					}
				}
			}
		}
	}

	/**
	 * setup function
	 */
	@Override
	public int setup(String arg0, ImagePlus img) {
		if (img == null) {
			IJ.error("Open an image please !");
		}
		this.image = img;
		this.size = img.getStackSize();
		this.width = img.getWidth();
		this.height = img.getHeight();
		this.type = this.image.getBitDepth()+"-bit";
		return DOES_ALL;
	}

	/**
	 * this method print the text file parameters in a Fiji log
	 */
	public void printParameters() {
		IJ.log("Denoising_Algorithm: " + this.algorithm);
		IJ.log("Iteration_Number: " + this.nit);
		IJ.log("***");
		for (Filter_3D f : this.filterList) {
			if (f instanceof Gabor3DFilter_3D) {
				Gabor3DFilter_3D f2 = (Gabor3DFilter_3D) f;
				IJ.log("Filter_Type: Gabor");
				IJ.log("Noise_Level: " + f2.getAlpha());
				IJ.log("Sigma_X: " + f2.getSigmax());
				IJ.log("Sigma_Y: " + f2.getSigmay());
				IJ.log("Sigma_Z: " + f2.getSigmaz());
				IJ.log("Angle_X: " + f2.getThetax());
				IJ.log("Angle_Y: " + f2.getThetay());
				IJ.log("Angle_Z: " + f2.getThetaz());
				IJ.log("***");
			} else {
				IJ.log("Filter_Type: Dirac");
				IJ.log("Noise_Level: " + f.getAlpha());
				IJ.log("***");
			}
		}
	}

	/**
	 * this function read a text file and get the plugin's parameters from it
	 * 
	 * @return false if an error is detected in the text file, else return true
	 */
	public boolean readFile() {
		OpenDialog od = new OpenDialog("Choose the file to read", "");

		String path = od.getDirectory() + od.getFileName();
		try {
			// this scanner read the lines
			Scanner scanFile = new Scanner(new File(path));

			// while the file is not end
			while (scanFile.hasNextLine() && this.errorInTxtFile == false) {
				// this scanner read the words
				Scanner scanLine = new Scanner(scanFile.nextLine());
				switch (scanLine.next()) {

				case "Denoising_Algorithm:":
					this.algorithm = scanLine.next().toString();
					// if error
					if (!(this.algorithm.equals("TV"))
							&& !(this.algorithm.equals("H1"))) {

						this.errorInTxtFile = true;
					}
					break;

				case "Iteration_Number:":
					this.nit = Double.parseDouble(scanLine.next());
					break;

				case "Filter_Type:":
					this.filterType = scanLine.next().toString();
					// if error
					if (!(this.filterType.equals("Dirac"))
							&& !(this.filterType.equals("Gabor"))) {
						this.errorInTxtFile = true;
					}
					break;

				case "Noise_Level:":
					this.noiseLevel = Double.parseDouble(scanLine.next());
					if (this.filterType.equals("Dirac")) {
						this.createDirac();
					}
					break;

				case "Sigma_X:":
					this.sigmax = Double.parseDouble(scanLine.next());
					break;

				case "Sigma_Y:":
					this.sigmay = Double.parseDouble(scanLine.next());
					break;

				case "Sigma_Z:":
					this.sigmaz = Double.parseDouble(scanLine.next());
					break;

				case "Angle_X:":
					this.thetax = Double.parseDouble(scanLine.next());
					break;

				case "Angle_Y:":
					this.thetay = Double.parseDouble(scanLine.next());
					break;

				case "Angle_Z:":
					this.thetaz = Double.parseDouble(scanLine.next());
					if (this.filterType.equals("Gabor")) {
						this.createGabor();
					}
					break;

				case "***":
					break;
				}
				scanLine.close();
			}
			scanFile.close();

		} catch (Exception e) {
			e.printStackTrace();
			IJ.log("Error : the text file is not conform !");
			return false;
		}
		if (this.errorInTxtFile == true) {
			IJ.log("Error : the text file is not conform !");
			return false;
		}
		return true;
	}

	/**
	 * this function ask the user the way to use the plugin, using the GUI or
	 * the Text File (.txt) mode
	 * 
	 * @return false if the window is closed by the user, else return true
	 */
	public boolean askMethod() {
		GenericDialog g = new GenericDialog("Choose the way to use the plugin");

		// signature
		String authors = "Welcome to VSNR plugin !\n \nIn case you use this algorithm, please cite:\nJerome Fehrenbach, Pierre Weiss, Corinne Lorenzo - ITAV\nPlugin Developer : Benjamin Font, Leo Mouly\n \n";
		g.addMessage(authors, new Font(authors,Font.CENTER_BASELINE,13));

		
		String s = "Where are the parameters from ?";
		String[] tabChoice = { "Graphic User Interface (GUI)",
				"Text File (.txt)" };
		String defaultItem = "Choose the input type";
		g.addChoice(s, tabChoice, defaultItem);
		g.pack();
		g.showDialog();

		this.inputMethod = g.getNextChoice();

		if (g.wasCanceled()) {
			return false;
		} else {
			return true;
		}
	}

	/**
	 * thus function ask the user the filter type he want to use
	 * 
	 * @return false if the window is closed by the user, else return true
	 */
	public boolean askFilter() {
		GenericDialog g = new GenericDialog("Filter to use !");

		String[] filterChoice = { "Dirac", "Gabor" };
		g.addChoice("Filter type ?", filterChoice, "Filter...");
		g.pack();
		g.showDialog();

		this.filterType = g.getNextChoice();

		if (g.wasCanceled()) {
			return false;
		} else {
			return true;
		}
	}

	/**
	 * this function ask the user the parameters of the filter
	 * 
	 * @return false if the window is closed by the user, else return true
	 */
	public boolean askParameters() {
		// the filter is dirac
		if (this.filterType.equals("Dirac")) {
			GenericDialog g = new GenericDialog(
					"Setting the Dirac's parameters !");

			g.addNumericField("Noise level :", this.noiseLevel, 0);
			g.pack();
			g.showDialog();

			this.noiseLevel = g.getNextNumber();

			if (g.wasCanceled()) {
				return false;
			} else {
				createDirac();
				return true;
			}
			// the filter is a gabor
		} else {
			GenericDialog g = new GenericDialog(
					"Setting the Gabor's parameters !");

			g.addNumericField("Noise level :", this.noiseLevel, 0);
			g.addNumericField("Sigma X :", this.sigmax, 1);
			g.addNumericField("Sigma Y :", this.sigmay, 1);
			g.addNumericField("Sigma Z :", this.sigmaz, 1);
			g.addNumericField("Theta X :", this.thetax, 1);
			g.addNumericField("Theta Y :", this.thetay, 1);
			g.addNumericField("Theta Z :", this.thetaz, 1);

			g.pack();
			g.showDialog();

			this.noiseLevel = g.getNextNumber();
			this.sigmax = g.getNextNumber();
			this.sigmay = g.getNextNumber();
			this.sigmaz = g.getNextNumber();
			this.thetax = g.getNextNumber();
			this.thetay = g.getNextNumber();
			this.thetaz = g.getNextNumber();

			if (g.wasCanceled()) {
				return false;
			} else {
				createGabor();
				return true;
			}

		}
	}

	/**
	 * this function ask the user if he want to set another filter
	 * 
	 * @return false if the window is closed by the user, else return true
	 */
	public boolean askOtherFilter() {
		GenericDialog g = new GenericDialog("Add another filter ?");
		g.addMessage("Do you want to add an other filter ?");

		g.pack();
		g.showDialog();

		if (g.wasCanceled()) {
			return false;
		} else {
			return true;
		}
	}

	/**
	 * this function ask the user the algorithm to use and the number of
	 * iterations
	 * 
	 * @return false if the window is closed by the user, else return true
	 */
	public boolean askIterationsAndAlgorithm() {
		GenericDialog g = new GenericDialog(
				"Number of iterations and Algorithm ?");

		String[] denoisingChoice = { "TV", "H1" };
		g.addChoice("Choose the denoising algorithm", denoisingChoice, "TV");
		g.addNumericField("Iterations (TV only) :", this.nit, 0);
		IterationListener_3D il = new IterationListener_3D();
		g.addDialogListener(il);
		g.pack();
		g.showDialog();

		this.algorithm = g.getNextChoice();
		this.nit = g.getNextNumber();

		if (g.wasCanceled()) {
			return false;
		} else {
			return true;
		}
	}

	/**
	 * this function create a dirac filter using the parameters entered by the
	 * user
	 */
	public void createDirac() {

		Filter_3D dirac = new Filter_3D(FilterType_3D.DIRAC, this.width, this.height,
				this.size);
		dirac.setAlpha(this.noiseLevel);
		// we add the filter to the array
		this.filterList.add(dirac);

	}

	/**
	 * this function create a gabor filter using the parameters entered by the
	 * user
	 */
	public void createGabor() {
		Gabor3DFilter_3D gabor3d = new Gabor3DFilter_3D(this.image,
				FilterType_3D.GABOR3D, this.width, this.height, this.size);
		gabor3d.setFlt(new Double2DArray_3D(this.height, this.width));
		gabor3d.setAlpha(this.noiseLevel);
		gabor3d.setSigmax(sigmax);
		gabor3d.setSigmay(sigmay);
		gabor3d.setSigmaz(sigmaz);
		gabor3d.setThetax(thetax);
		gabor3d.setThetay(thetay);
		gabor3d.setThetaz(thetaz);
		// we add the filter to the array
		this.filterList.add(gabor3d);
	}

	/**
	 * this function denoise the stack using the Vsnr_denoiser.denoiseTV_3D()
	 * function
	 * 
	 * @return the denoised stack
	 */
	public ImagePlus denoiseTV_3D() {

		this.stack = image.getStack();

		// security
		if (filterList.isEmpty()) {

			IJ.log("Unable to process (no filters set)!\nTry to add some filters first \n");
		}
		if (this.image == null) {
			IJ.log("Something bad happened. You probably closed the image.\nPlease reload VSNR! \n");
		}

		// preparation of the parameters of the denoising algorithm
		Double3DArray_3D[] Psis = new Double3DArray_3D[filterList.size()];
		double[] alphas = new double[filterList.size()];

		// we prepare the differents images (Double3DArray)
		Double3DArray_3D img3DToDenoise = new Double3DArray_3D(this.size,
				this.height, this.width);
		Double3DArray_3D denoised = new Double3DArray_3D(this.size, this.height,
				this.width);
		Double3DArray_3D noise = new Double3DArray_3D(this.size, this.height,
				this.width);

		// we get the pixels of the image inside a Double3DArray
		ImageProcessor sliceProcessor;

		for (int k = 0; k < size; k++) {
			sliceProcessor = stack.getProcessor(k + 1);
			for (int i = 0; i < height; i++) {
				for (int j = 0; j < stack.getWidth(); j++) {
					img3DToDenoise.setValue(sliceProcessor.getPixel(j, i), k,
							i, j, false);
				}
			}
		}

		// this int will help for filling alphas and psis
		int position = 0;

		// for each filter
		for (Filter_3D filter_3D : this.filterList) {
			// if the filter is a gabor, we use the method computeNewGabor()
			if (filter_3D instanceof Gabor3DFilter_3D) {
				Gabor3DFilter_3D gabor = (Gabor3DFilter_3D) filter_3D;
				gabor.computeNewGabor3D();
				alphas[position] = gabor.getAlpha();
				Psis[position] = gabor.getStackFilter();
			} else {
				alphas[position] = filter_3D.getAlpha();
				Psis[position] = filter_3D.getStackFilter();
			}
			position++;
		}

		// we use the denoising algorithm
		denoised = Vsnr_denoiser_3D.denoiseTV_3D(img3DToDenoise, Psis, alphas,
				(int) this.nit, noise);

		// Pierre's code
		Double3DArray_3D denoised3D = new Double3DArray_3D(size, height, width);
		Double3DArray_3D noise3D = new Double3DArray_3D(size, height, width);

		for (int k = 0; k < size; k++) {
			for (int i = 0; i < denoised3D.getRows(); i++) {
				for (int j = 0; j < denoised3D.getColumns(); j++) {
					denoised3D.setValue(denoised.getValue(k, i, j, false), k,
							i, j, false);
					noise3D.setValue(noise.getValue(k, i, j, false), k, i, j,
							false);
				}
			}
		}

		ImagePlus resultImg = IJ.createImage("vsnr_result_"
				+ Vsnr_Plugin_3D.this.image.getTitle(), this.type, width,
				height, 0);
		ImageStack resultStack = IJ.createImage(
				"vsnr_result_" + Vsnr_Plugin_3D.this.image.getTitle(),
				this.type, width, height, 0).getImageStack();
		resultStack.deleteLastSlice();

		ImagePlus resultNoise = IJ.createImage("vsnr_noise_"
				+ Vsnr_Plugin_3D.this.image.getTitle(), this.type, width,
				height, 0);
		ImageStack resultStackNoise = IJ.createImage(
				"vsnr_noise_" + Vsnr_Plugin_3D.this.image.getTitle(), this.type,
				width, height, 0).getImageStack();
		resultStackNoise.deleteLastSlice();


		for (int k = 0; k < size; k++) {
			resultImg = IJ.createImage("vsnr_result_"
					+ Vsnr_Plugin_3D.this.image.getTitle(), this.type, width,
					height, size);
			resultNoise = IJ.createImage("vsnr_noise_"
					+ Vsnr_Plugin_3D.this.image.getTitle(), this.type, width,
					height, size);
			/*
			double[] minnmax1 = new double[2];
			
			double[] minnmax2 = new double[2];
			minnmax1 = denoised3D.minNmax();
			minnmax2 = noise3D.minNmax();
			for (int i = 0; i < height; i++) {
				for (int j = 0; j < width; j++) {
					double val = (denoised3D.getValue(k, i, j, false) - minnmax1[0])
							* 65535.0D / (minnmax1[1] - minnmax1[0]);
					resultImg.getProcessor().putPixel(j, i,
							(int) Math.floor(val));
					val = (noise3D.getValue(k, i, j, false) - minnmax2[0])
							* 65535.0D / (minnmax2[1] - minnmax2[0]);
					resultNoise.getProcessor().putPixel(j, i,
							(int) Math.floor(val));
				}
			}
			*/
			for (int i = 0; i < height; i++) {
				for (int j = 0; j < width; j++) {
					double val = denoised3D.getValue(k, i, j, false);
					resultImg.getProcessor().putPixel(j, i,
							(int) Math.floor(val));
					val = noise3D.getValue(k, i, j, false);
					resultNoise.getProcessor().putPixel(j, i,
							(int) Math.floor(val));
				}
			}
			
			
			resultStack.addSlice("Slice #" + k, resultImg.getProcessor());
			resultStackNoise.addSlice("Slice noise #" + k,
					resultNoise.getProcessor());
		}
		ImagePlus imresult = new ImagePlus("Denoised Image", resultStack);
		//imresult.setDisplayRange(0.0D, 65535.0D);
		ImagePlus imnoise = new ImagePlus("Noise", resultStackNoise);
		//imnoise.setDisplayRange(0.0D, 65535.0D);

		this.noiseImage = imnoise;

		this.noiseImage.setDisplayRange(this.image.getDisplayRangeMin(),this.image.getDisplayRangeMax());
		imresult.setDisplayRange(this.image.getDisplayRangeMin(),this.image.getDisplayRangeMax());
		
		return imresult;

	}

	/**
	 * this function denoise the stack using the Vsnr_denoiser.denoiseH1_3D()
	 * function
	 * 
	 * @return the denoised stack
	 */
	public ImagePlus denoiseH1_3D() {

		this.stack = image.getStack();

		// security
		if (filterList.isEmpty()) {

			IJ.log("Unable to process (no filters set)!\nTry to add some filters first \n");
		}
		if (this.image == null) {
			IJ.log("Something bad happened. You probably closed the image.\nPlease reload VSNR! \n");
		}

		// preparation of the parameters of the denoising algorithm
		Double3DArray_3D[] Psis = new Double3DArray_3D[filterList.size()];
		double[] alphas = new double[filterList.size()];

		// we prepare the differents images (Double3DArray)
		Double3DArray_3D img3DToDenoise = new Double3DArray_3D(this.size,
				this.height, this.width);
		Double3DArray_3D denoised = new Double3DArray_3D(this.size, this.height,
				this.width);
		Double3DArray_3D noise = new Double3DArray_3D(this.size, this.height,
				this.width);

		// we get the pixels of the image inside a Double3DArray
		ImageProcessor sliceProcessor;

		for (int k = 0; k < size; k++) {
			sliceProcessor = stack.getProcessor(k + 1);
			for (int i = 0; i < height; i++) {
				for (int j = 0; j < stack.getWidth(); j++) {
					img3DToDenoise.setValue(sliceProcessor.getPixel(j, i), k,
							i, j, false);
				}
			}
		}

		// this int will help for filling alphas and psis
		int position = 0;

		// for each filter
		for (Filter_3D filter_3D : this.filterList) {
			// if the filter is a gabor, we use the method computeNewGabor()
			if (filter_3D instanceof Gabor3DFilter_3D) {
				Gabor3DFilter_3D gabor = (Gabor3DFilter_3D) filter_3D;
				gabor.computeNewGabor3D();
				alphas[position] = gabor.getAlpha();
				Psis[position] = gabor.getStackFilter();
			} else {
				alphas[position] = filter_3D.getAlpha();
				Psis[position] = filter_3D.getStackFilter();
			}
			position++;
		}

		// we use the denoising algorithm
		denoised = Vsnr_denoiser_3D.denoiseH1_3D(img3DToDenoise, Psis, alphas,
				noise);

		// Pierre's code
		Double3DArray_3D denoised3D = new Double3DArray_3D(size, height, width);
		Double3DArray_3D noise3D = new Double3DArray_3D(size, height, width);

		for (int k = 0; k < size; k++) {
			for (int i = 0; i < denoised3D.getRows(); i++) {
				for (int j = 0; j < denoised3D.getColumns(); j++) {
					denoised3D.setValue(denoised.getValue(k, i, j, false), k,
							i, j, false);
					noise3D.setValue(noise.getValue(k, i, j, false), k, i, j,
							false);
				}
			}
		}

		ImagePlus resultImg = IJ.createImage("vsnr_result_"
				+ Vsnr_Plugin_3D.this.image.getTitle(), this.type, width,
				height, 0);
		ImageStack resultStack = IJ.createImage(
				"vsnr_result_" + Vsnr_Plugin_3D.this.image.getTitle(),
				this.type, width, height, 0).getImageStack();
		resultStack.deleteLastSlice();

		ImagePlus resultNoise = IJ.createImage("vsnr_noise_"
				+ Vsnr_Plugin_3D.this.image.getTitle(), this.type, width,
				height, 0);
		ImageStack resultStackNoise = IJ.createImage(
				"vsnr_noise_" + Vsnr_Plugin_3D.this.image.getTitle(), this.type,
				width, height, 0).getImageStack();
		resultStackNoise.deleteLastSlice();

		for (int k = 0; k < size; k++) {
			resultImg = IJ.createImage("vsnr_result_"
					+ Vsnr_Plugin_3D.this.image.getTitle(), this.type, width,
					height, size);
			resultNoise = IJ.createImage("vsnr_noise_"
					+ Vsnr_Plugin_3D.this.image.getTitle(), this.type, width,
					height, size);

			/*
			double[] minnmax1 = new double[2];
			
			double[] minnmax2 = new double[2];
			minnmax1 = denoised3D.minNmax();
			minnmax2 = noise3D.minNmax();
			for (int i = 0; i < height; i++) {
				for (int j = 0; j < width; j++) {
					double val = (denoised3D.getValue(k, i, j, false) - minnmax1[0])
							* 65535.0D / (minnmax1[1] - minnmax1[0]);
					resultImg.getProcessor().putPixel(j, i,
							(int) Math.floor(val));
					val = (noise3D.getValue(k, i, j, false) - minnmax2[0])
							* 65535.0D / (minnmax2[1] - minnmax2[0]);
					resultNoise.getProcessor().putPixel(j, i,
							(int) Math.floor(val));
				}
			}
			*/
			
			for (int i = 0; i < height; i++) {
				for (int j = 0; j < width; j++) {
					double val = denoised3D.getValue(k, i, j, false);
					resultImg.getProcessor().putPixel(j, i,
							(int) Math.floor(val));
					val = noise3D.getValue(k, i, j, false);
					resultNoise.getProcessor().putPixel(j, i,
							(int) Math.floor(val));
				}
			}
			
			resultStack.addSlice("Slice #" + k, resultImg.getProcessor());
			resultStackNoise.addSlice("Slice noise #" + k,
					resultNoise.getProcessor());
		}
		ImagePlus imresult = new ImagePlus("Denoised Image", resultStack);
		//imresult.setDisplayRange(0.0D, 65535.0D);
		ImagePlus imnoise = new ImagePlus("Noise", resultStackNoise);
		//imnoise.setDisplayRange(0.0D, 65535.0D);

		this.noiseImage = imnoise;
		
		this.noiseImage.setDisplayRange(this.image.getDisplayRangeMin(),this.image.getDisplayRangeMax());
		imresult.setDisplayRange(this.image.getDisplayRangeMin(),this.image.getDisplayRangeMax());

		
		return imresult;

	}
}
