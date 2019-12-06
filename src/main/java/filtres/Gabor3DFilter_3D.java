package filtres;

import java.io.Serializable;
import process.Double3DArray_3D;

import ij.ImagePlus;

public class Gabor3DFilter_3D extends Filter_3D implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = -7402740491860517689L;
	private double sigmax, sigmay, sigmaz, thetax, thetay, thetaz;
	private Double3DArray_3D gabor3D;

	public Gabor3DFilter_3D(ImagePlus imgF, FilterType_3D tF, int widthF, int heightF) {
		super(imgF, tF, widthF, heightF);
		this.sigmax = 10;
		this.sigmay = 10;
		this.sigmaz = 10;
		this.thetax = 0;
		this.thetay = 0;
		this.thetaz = 0;
		gabor3D = new Double3DArray_3D(imgF.getNSlices(), imgF.getHeight(),
				imgF.getWidth());
	}

	public Gabor3DFilter_3D(ImagePlus imgF, FilterType_3D typeF, int widthF,
			int heightF, int s) {
		super(imgF, typeF, 0, 0.1, widthF, heightF, s);
		this.sigmax = 10;
		this.sigmay = 10;
		this.sigmaz = 10;
		this.thetax = 0;
		this.thetay = 0;
		this.thetaz = 0;
		gabor3D = new Double3DArray_3D(imgF.getNSlices(), imgF.getHeight(),
				imgF.getWidth());
	}

	public double getSigmax() {
		return sigmax;
	}

	public void setSigmax(double sigmax) {
		this.sigmax = sigmax;
	}

	public double getSigmay() {
		return sigmay;
	}

	public void setSigmay(double sigmay) {
		this.sigmay = sigmay;
	}

	public double getSigmaz() {
		return sigmaz;
	}

	public void setSigmaz(double sigmaz) {
		this.sigmaz = sigmaz;
	}

	public double getThetax() {
		return thetax;
	}

	public void setThetax(double thetax) {
		this.thetax = thetax;
	}

	public double getThetay() {
		return thetay;
	}

	public void setThetay(double thetay) {
		this.thetay = thetay;
	}

	public double getThetaz() {
		return thetaz;
	}

	public void setThetaz(double thetaz) {
		this.thetaz = thetaz;
	}

	public Double3DArray_3D getStackFilter() {
		return gabor3D;
	}

	public void computeNewGabor3D() {
		int width = this.getWidth();
		int height = this.getHeight();
		int slices = this.getSlices();

		// thetas must be in radians !
		double thetax = (double) (this.getThetax() / 180) * Math.PI;
		double thetay = (double) (this.getThetay() / 180) * Math.PI;
		double thetaz = (double) (this.getThetaz() / 180) * Math.PI;

		double sigmax = this.getSigmax();
		double sigmay = this.getSigmay();
		double sigmaz = this.getSigmaz();

		double c1 = Math.cos(thetax);
		double s1 = Math.sin(thetax);

		double c2 = Math.cos(thetay);
		double s2 = Math.sin(thetay);

		double c3 = Math.cos(thetaz);
		double s3 = Math.sin(thetaz);

		double Rx[][] = { { 1, 0, 0 }, { 0, c1, -s1 }, { 0, s1, c1 } };

		double Ry[][] = { { c2, 0, s2 }, { 0, 1, 0 }, { -s2, 0, c2 } };

		double Rz[][] = { { c3, -s3, 0 }, { s3, c3, 0 }, { 0, 0, 1 } };

		double RR[][] = new double[3][3];

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					RR[i][j] = RR[i][j] + Rx[i][k] * Ry[k][j];
				}
			}
		}

		double RRR[][] = new double[3][3];

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					RRR[i][j] = RRR[i][j] + RR[i][k] * Rz[k][j];

				}
			}
		}

		double s = 0, x, y, z, value;
		double res[] = new double[3];

		// computing gabor pixels
		for (int i = 0; i < width; i++) {
			x = (i + 1) - Math.floor(width / 2);
			for (int j = 0; j < height; j++) {
				y = (j + 1) - Math.floor(height / 2);
				for (int k = 0; k < slices; k++) {
					z = (k + 1) - Math.floor(slices / 2);

					for (int l = 0; l < 3; l++) {
						res[l] = RRR[l][0] * x + RRR[l][1] * y + RRR[l][2] * z;
						// IJ.log("res ="+res[l]);
					}

					value = Math
							.exp(-((res[0] / sigmax) * (res[0] / sigmax)
									+ (res[1] / sigmay) * (res[1] / sigmay) + (res[2] / sigmaz)
									* (res[2] / sigmaz)));
					gabor3D.setValue(value, k, j, i, false);
					s = s + value;

				}
			}
		}

		double pixel;
		double min = 999999999;
		double max = -999999999;
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				for (int k = 0; k < slices; k++) {
					pixel = (gabor3D.getValue(k, i, j, false) / s);
					gabor3D.setValue(pixel, k, i, j, false);
					if (pixel > max) {
						max = pixel;
					} else if (pixel < min) // && pixel != 0 ?
					{
						min = pixel;
					}

				}
			}
		}
		// rescale to 0->255
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				for (int k = 0; k < slices; k++) {
					pixel = Math
							.floor(((gabor3D.getValue(k, i, j, false) - min) * 255)
									/ (max - min));
					gabor3D.setValue(pixel, k, i, j, false);

				}
			}
		}
	}

	public Double3DArray_3D getGabor3D() {
		return gabor3D;
	}

	public void setGabor3D(Double3DArray_3D gabor3d) {
		gabor3D = gabor3d;
	}

}
