package process;
import ij.*;

public class Vsnr_denoiser_3D {

	/**
	 * Warning : this is not entireley my code so the parameters naming is maybe
	 * wrong
	 * 
	 * This function is denoising a stack using the TV algorithm
	 * 
	 * @param u0
	 *            stack to denoise
	 * @param psis
	 *            array of filters
	 * @param etas
	 *            array of noise levels
	 * @param nit
	 *            number of iterations
	 * @param noise
	 *            the final noise array
	 * @return the denoised stack
	 */
	public static Double3DArray_3D denoiseTV_3D(Double3DArray_3D u0,
			Double3DArray_3D[] psis, double[] etas, int nit, Double3DArray_3D noise) {
		try {
			IJ.showStatus("Starting denoising ...");

			int alpha = 1;

			int nbSlices = u0.getSlices();
			int nbRows = u0.getRows();
			int nbCols = u0.getColumns();

			Double3DArray_3D d1 = new Double3DArray_3D(nbSlices, nbRows, nbCols);
			d1.setValue(1, 0, 0, 0, false);
			d1.setValue(-1, 0, d1.getRows() - 1, 0, false);

			Double3DArray_3D d2 = new Double3DArray_3D(nbSlices, nbRows, nbCols);
			d2.setValue(1, 0, 0, 0, false);
			d2.setValue(-1, 0, 0, d2.getColumns() - 1, false);

			Double3DArray_3D d3 = new Double3DArray_3D(nbSlices, nbRows, nbCols);
			d3.setValue(1, 0, 0, 0, false);
			d3.setValue(-1, d1.getSlices() - 1, 0, 0, false);

			d1 = d1.getFFTn();

			d2 = d2.getFFTn();

			d3 = d3.getFFTn();

			Double3DArray_3D fu0 = u0.getFFTn();

			Double3DArray_3D[] fPsis = new Double3DArray_3D[psis.length];
			for (int i = 0; i < psis.length; i++) {
				fPsis[i] = psis[i].getFFTn();
			}
			double[] normesInf = new double[psis.length];
			double h = 0;
			Double3DArray_3D hh = new Double3DArray_3D(nbSlices, nbRows, nbCols);
			for (int i = 0; i < psis.length; i++) {
				normesInf[i] = 0;
				hh.compute_hh(fPsis[i]);
				h = hh.compute_h(d1);
				if (normesInf[i] < h) {
					normesInf[i] = h;
				}

				h = hh.compute_h(d2);
				if (normesInf[i] < h) {
					normesInf[i] = h;
				}
				h = hh.compute_h(d3);
				if (normesInf[i] < h) {
					normesInf[i] = h;
				}
			}

			double[] alphas = u0.compute_alphas(etas, normesInf);

			Double3DArray_3D fPsi = new Double3DArray_3D(nbSlices, nbRows, nbCols);
			for (int i = 0; i < psis.length; i++) {
				fPsi.compute_fPsi3(fPsis[i], alphas[i]);
			}

			fPsi.compute_fPsi4(alpha);
			double[] CF = new double[nit];

			Double3DArray_3D y1 = new Double3DArray_3D(nbSlices, nbRows, nbCols);
			Double3DArray_3D y2 = new Double3DArray_3D(nbSlices, nbRows, nbCols);
			Double3DArray_3D y3 = new Double3DArray_3D(nbSlices, nbRows, nbCols);

			Double3DArray_3D q1 = new Double3DArray_3D(nbSlices, nbRows, nbCols);
			Double3DArray_3D q2 = new Double3DArray_3D(nbSlices, nbRows, nbCols);
			Double3DArray_3D q3 = new Double3DArray_3D(nbSlices, nbRows, nbCols);

			Double3DArray_3D d1u0 = new Double3DArray_3D(nbSlices, nbRows, nbCols);
			d1u0 = d1u0.compute_dXu0(d1, fu0);

			Double3DArray_3D d2u0 = new Double3DArray_3D(nbSlices, nbRows, nbCols);
			d2u0 = d2u0.compute_dXu0(d2, fu0);
			Double3DArray_3D d3u0 = new Double3DArray_3D(nbSlices, nbRows, nbCols);
			d3u0 = d3u0.compute_dXu0(d3, fu0);
			double L = -1;

			L = fPsi.compute_L(alpha, d1, d2, d3);
			Double3DArray_3D fy1, fy2, fy3, fAy = new Double3DArray_3D(nbSlices,
					nbRows, nbCols), aY = new Double3DArray_3D(nbSlices, nbRows,
					nbCols), nablaF1, nablaF2, nablaF3, q1p, q2p, q3p, b, nq = new Double3DArray_3D(
					nbSlices, nbRows, nbCols);

			for (int i = 1; i <= nit; i++) {

				fy1 = y1.getFFTn();
				fy2 = y2.getFFTn();
				fy3 = y3.getFFTn();
				fAy = fAy.compute_fAy(fPsi, fy1, fy2, fy3, d1, d2, d3);
				aY = fAy.getIFFTn(true);

				CF[i - 1] = aY.compute_CF(d1u0, y1, d2u0, y2, d3u0, y3, alpha);

				nablaF1 = d1u0.compute_nablaFX(alpha, d1, fPsi, fAy);
				nablaF2 = d2u0.compute_nablaFX(alpha, d2, fPsi, fAy);
				nablaF3 = d3u0.compute_nablaFX(alpha, d3, fPsi, fAy);
				q1p = new Double3DArray_3D(q1.getArray());
				q2p = new Double3DArray_3D(q2.getArray());
				q3p = new Double3DArray_3D(q3.getArray());

				q1 = q1.compute_qX(y1, L, nablaF1);
				q2 = q2.compute_qX(y2, L, nablaF2);
				q3 = q3.compute_qX(y3, L, nablaF3);

				nq = q1.compute_nq(q2, q3);

				q1.compute_qX2(nq);
				q2.compute_qX2(nq);
				q3.compute_qX2(nq);

				y1 = y1.compute_yX(q1, i, q1p);
				y2 = y2.compute_yX(q2, i, q2p);
				y3 = y3.compute_yX(q3, i, q3p);

				System.gc();

				// we warn the user that we are denoising
				IJ.showStatus("Denoising... Iteration n°" + i);
				IJ.showProgress(i, (nit + 1));
			}

			Double3DArray_3D u = new Double3DArray_3D(nbSlices, nbRows, nbCols);

			u.compute_u(u0, fPsi, fAy, alpha, noise);

			// We should clear the memory
			y1 = null;
			y2 = null;
			y3 = null;
			d1 = null;
			d2 = null;
			d3 = null;
			nq = null;
			q1 = null;
			q2 = null;
			q3 = null;
			q1p = null;
			q2p = null;
			q3p = null;
			nablaF1 = null;
			nablaF2 = null;
			nablaF3 = null;
			d1u0 = null;
			d2u0 = null;
			d3u0 = null;
			fy1 = null;
			fy2 = null;
			fy3 = null;
			fAy = null;
			aY = null;
			hh = null;
			CF = null;
			fu0 = null;
			fPsi = null;
			fPsis = null;
			System.gc();

			// everything is fine, so we return the denoised image !
			IJ.showStatus("Denoising ended !");
			return u;
		} catch (NullPointerException e) {

			IJ.handleException(e);
			IJ.log("Denoising aborded.");
			return null;
		}
	}

	/**
	 * Warning : this is not entireley my code so the parameters naming is maybe
	 * wrong
	 * 
	 * This function is denoising a stack using the TV algorithm
	 * 
	 * @param u0
	 *            stack to denoise
	 * @param psis
	 *            array of filters
	 * @param etas
	 *            array of noise levels
	 * @param noise
	 *            the final noise array
	 * @return the denoised stack
	 */
	public static Double3DArray_3D denoiseH1_3D(Double3DArray_3D u0,
			Double3DArray_3D[] psis, double[] etas, Double3DArray_3D noise) {
		try {

			// we warn th user that we are denoising
			IJ.showStatus("Starting denoising ...");

			int nbSlices = u0.getSlices();
			int nbRows = u0.getRows();
			int nbCols = u0.getColumns();

			Double3DArray_3D d1 = new Double3DArray_3D(nbSlices, nbRows, nbCols);
			d1.setValue(1, 0, 0, 0, false);
			d1.setValue(-1, 0, d1.getRows() - 1, 0, false);

			Double3DArray_3D d2 = new Double3DArray_3D(nbSlices, nbRows, nbCols);
			d2.setValue(1, 0, 0, 0, false);
			d2.setValue(-1, 0, 0, d2.getColumns() - 1, false);

			Double3DArray_3D d3 = new Double3DArray_3D(nbSlices, nbRows, nbCols);
			d3.setValue(1, 0, 0, 0, false);
			d3.setValue(-1, d1.getSlices() - 1, 0, 0, false);

			d1 = d1.getFFTn();
			d2 = d2.getFFTn();
			d3 = d3.getFFTn();

			Double3DArray_3D fu0 = new Double3DArray_3D(nbSlices, nbRows, nbCols);
			fu0 = u0.getFFTn(); // fu0=fftn(u0);
			double norm_fu0 = fu0.getNorm();

			Double3DArray_3D[] fPsis = new Double3DArray_3D[psis.length];
			for (int i = 0; i < psis.length; i++) {
				fPsis[i] = psis[i].getFFTn();
			}

			int m = psis.length;
			double[] nbruits;
			nbruits = new double[psis.length]; // morozov's estimated noises

			for (int i = 0; i < psis.length; i++) {
				nbruits[i] = etas[i] * norm_fu0;// each component of the Morozov
												// estimated noises
			}

			Double3DArray_3D d = new Double3DArray_3D(nbSlices, nbRows, nbCols);
			d.computeD(d1, d2, d3); // d=abs(fd1).^2+abs(fd2).^2+abs(fd3).^2;
			d1 = null;
			d2 = null;
			d3 = null;

			Double3DArray_3D temp = new Double3DArray_3D(nbSlices, nbRows, nbCols);
			double[] alphas = new double[m];

			for (int i = 0; i < m; i++) {
				temp = (fPsis[i].getAbs()).getArrayPowerOf2();
				temp.multiplyByComplexArray(d);
				temp.multiplyByComplexArray(fu0);
				alphas[i] = temp.getNorm();
				alphas[i] = alphas[i] / nbruits[i];
			}

			Double3DArray_3D sumpsi = new Double3DArray_3D(nbSlices, nbRows, nbCols);
			Double3DArray_3D[] bruiths = new Double3DArray_3D[m];
			double[] nbruiths = new double[m];

			for (int i = 0; i < m; i++) {
				sumpsi.compute_sumpsi(fPsis[i], d, alphas[i]);
			}

			Double3DArray_3D d_mutliply_fu0 = d.multiply(fu0);
			temp.add1(sumpsi);
			for (int i = 0; i < m; i++) {
				bruiths[i] = new Double3DArray_3D(nbSlices, nbRows, nbCols);
				bruiths[i] = (fPsis[i].getAbs()).getArrayPowerOf2();
				bruiths[i].multiplyByComplexArray(d_mutliply_fu0);
				bruiths[i].divideByRealArray(temp);
				nbruiths[i] = (bruiths[i].getNorm()) / alphas[i];
			}

			boolean b = true;
			int nit = 0;
			int maxiter = 5;
			while (b && nit <= maxiter) {
				alphas[0] = alphas[0] * nbruiths[0] / nbruits[0];
				for (int i = 0; i < m - 1; i++) {
					sumpsi.ReInitializeToZero();
					for (int j = 0; j < m; j++) {
						sumpsi.compute_sumpsi(fPsis[j], d, alphas[j]);
					}
					bruiths[i + 1] = (fPsis[i + 1].getAbs()).getArrayPowerOf2();
					bruiths[i + 1].multiplyByComplexArray(d_mutliply_fu0);
					temp.add1(sumpsi);
					bruiths[i + 1].divideByRealArray(temp);
					nbruiths[i + 1] = (bruiths[i + 1].getNorm())
							/ alphas[i + 1];
					alphas[i + 1] = (nbruiths[i + 1] / nbruits[i + 1])
							* alphas[i + 1];
				}
				sumpsi.ReInitializeToZero();
				for (int j = 0; j < m; j++) {
					sumpsi.compute_sumpsi(fPsis[j], d, alphas[j]);
				}

				temp.add1(sumpsi);
				for (int i = 0; i < m; i++) {
					bruiths[i] = (fPsis[i].getAbs()).getArrayPowerOf2();
					bruiths[i].multiplyByComplexArray(d_mutliply_fu0);
					bruiths[i].divideByReal(alphas[i]);
					bruiths[i].divideByRealArray(temp);
					nbruiths[i] = (bruiths[i].getNorm());
				}
				int i = 0;
				int j = 0;
				while (i < m && b) {
					if ((nbruits[i] - nbruiths[i]) / nbruits[i] < (1 / 10)) {
						j++;
					}
					i++;
				}
				if (j == m)
					b = false;

				// we warn the user that we are denoising
				IJ.showStatus("Denoising ...");
				IJ.showProgress(nit, (maxiter + 1));
				nit++;
			}

			Double3DArray_3D fnoise = new Double3DArray_3D(nbSlices, nbRows, nbCols);
			for (int i = 0; i < m; i++) {
				fnoise.compute_fnoise(bruiths[i]);
			}
			fnoise = fnoise.getIFFTn(true);
			noise.setArray(fnoise.getArray());

			Double3DArray_3D u = new Double3DArray_3D(nbSlices, nbRows, nbCols);
			u.denoiseImage(u0, noise);

			IJ.showStatus("Denoising done !");
			return u;
		} catch (Exception e) {
			IJ.handleException(e);
		}
		return null;
	}

}
