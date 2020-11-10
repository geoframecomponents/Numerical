package it.geoframe.blogspot.numerical.newtonalgorithm;

import java.util.List;

import it.geoframe.blogspot.numerical.linearsystemsolver.Thomas;
import it.geoframe.blogspot.closureequation.equationstate.EquationState;

/*
 * GNU GPL v3 License
 *
 * Copyright 2020  Niccolo` Tubini
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/**
 * This class carries out the Nested-Newton algorithm
 * (A NESTED NEWTON-TYPE ALGORITHM FOR FINITE VOLUME METHODS SOLVING RICHARDS' EQUATION IN MIXED FORM, Casulli V., Zanolli P., Journal Scientific Computing, 2010)
 *  @author Niccolo' Tubini
 */

public class NestedNewtonThomas {
	private double outerResidual;
	private double innerResidual;

	private int nestedNewton;
	private int MAXITER_NEWT;
	private int KMAX;

	private double newtonTolerance;

	private double[] x;
	private double[] y;
	private double[] dx;
	private double[] mainDiagonal;
	private double[] upperDiagonal;
	private double[] lowerDiagonal;
	private double[] rhss;


	private double[] fs;
	private double[] fks;
	private double[] bb;
	private double[] cc;
	private double[] dis;
	private int[] rheologyID;
	private int[] parameterID;
	private double[] x_outer;

	private double delta;

	private List<EquationState> equationState;
	private Thomas thomasAlg = new Thomas();



	public NestedNewtonThomas(int nestedNewton, double newtonTolerance, int MAXITER_NEWT, int VECTOR_LENGTH, List<EquationState> equationState, double delta){

		this.nestedNewton = nestedNewton;
		this.newtonTolerance = newtonTolerance;
		this.MAXITER_NEWT = MAXITER_NEWT;
		//		this.NUM_CONTROL_VOLUMES = VECTOR_LENGTH;
		this.equationState = equationState;
		this.delta = delta;


		x = new double[VECTOR_LENGTH];
		dx = new double[VECTOR_LENGTH];
		fs = new double[VECTOR_LENGTH];
		fks = new double[VECTOR_LENGTH];
		bb = new double[VECTOR_LENGTH]; 
		cc = new double[VECTOR_LENGTH];
		dis = new double[VECTOR_LENGTH];
		x_outer = new double[VECTOR_LENGTH];
	}




	public void set(double[] x, double[] y, double[] mainDiagonal, double[] upperDiagonal, double[] lowerDiagonal, double[] rhss, int KMAX, int[] parameterID, int[] rheologyID){

		this.x = x;
		this.y = y;
		this.mainDiagonal = mainDiagonal;
		this.upperDiagonal = upperDiagonal;
		this.lowerDiagonal = lowerDiagonal;
		this.rhss = rhss;

		this.KMAX = KMAX;
		
		this.parameterID = parameterID;
		this.rheologyID = rheologyID;

	}



	public double[] solver(){

		/*
		 *  Initial guess
		 */
		if(nestedNewton == 0) {

		} else {
			for(int element = 0; element < KMAX; element++) {
				//x[element] = Math.min(x[element], xStar[element]-1 );
				x[element] = equationState.get(rheologyID[element]).initialGuess(x[element],parameterID[element],element);
			}
		}

		/* 
		 * OUTER CYCLE
		 */
		for(int i = 0; i < MAXITER_NEWT; i++) {
			// I have to assign 0 to outerResidual otherwise I will take into account of the previous error
			outerResidual = 0.0;
			for(int element = 0; element < KMAX; element++) {
				if(element==0) {
					fs[element] = equationState.get(rheologyID[element]).equationState(x[element],y[element],parameterID[element],element) - rhss[element]  + mainDiagonal[element]*x[element] + upperDiagonal[element]*x[element+1];
					dis[element] = equationState.get(rheologyID[element]).dEquationState(x[element],y[element],parameterID[element],element);
//										System.out.println(element+" "+fs[element]);
				} else if(element==KMAX-1) {
					fs[element] = equationState.get(rheologyID[element]).equationState(x[element],y[element],parameterID[element],element) - rhss[element] + lowerDiagonal[element]*x[element-1] + mainDiagonal[element]*x[element];
					dis[element] = equationState.get(rheologyID[element]).dEquationState(x[element],y[element],parameterID[element],element);
//										System.out.println(element+" "+fs[element]);
				} else {
					fs[element] = equationState.get(rheologyID[element]).equationState(x[element],y[element],parameterID[element],element) - rhss[element] + lowerDiagonal[element]*x[element-1] + mainDiagonal[element]*x[element] + upperDiagonal[element]*x[element+1];
					dis[element] = equationState.get(rheologyID[element]).dEquationState(x[element],y[element],parameterID[element],element);
//										System.out.println(element+" "+fs[element]);
				}

				outerResidual += fs[element]*fs[element];
			}
			outerResidual = Math.pow(outerResidual,0.5);  
//									System.out.println("\tOuter iteration " + i + " with residual " +  outerResidual);
			if(outerResidual < newtonTolerance) {
				break;
			}

			if(nestedNewton == 0){
				
				bb = mainDiagonal.clone();
				cc = upperDiagonal.clone();
				for(int y = 0; y < KMAX; y++) {
					bb[y] += dis[y];
				}
				thomasAlg.set(cc,bb,lowerDiagonal,fs,KMAX);
				dx = thomasAlg.solver();

				//// UPDATE SOLUTION////
				for(int element = 0; element < KMAX; element++) {
					x[element] = x[element] - dx[element]*delta; //if multiply by delta you get the globally convergent Newton aka Newton damped.
				}
				
			} else {
				x_outer = x.clone();
				//				for(int ii = 0; ii < KMAX; ii++) {
				//					x[ii] = Math.max(x[ii], xStar[ii] );
				//				}


				//// INNER CYCLE ////
				for(int j = 0; j < MAXITER_NEWT; j++) {
					// I have to assign 0 to innerResidual otherwise I will take into account of the previous error
					innerResidual = 0.0; 
					for(int element=0; element < KMAX; element++) {
						if(element==0) {
							fks[element] = equationState.get(rheologyID[element]).pIntegral(x[element],y[element],parameterID[element],element) - ( equationState.get(rheologyID[element]).qIntegral(x_outer[element],y[element],parameterID[element],element) + equationState.get(rheologyID[element]).q(x_outer[element],y[element],parameterID[element],element)*(x[element]-x_outer[element]) ) - this.rhss[element] + mainDiagonal[element]*x[element] + upperDiagonal[element]*x[element+1];
							dis[element] = ( equationState.get(rheologyID[element]).p(x[element],y[element],parameterID[element],element) - equationState.get(rheologyID[element]).q(x_outer[element],y[element],parameterID[element],element) );
//														System.out.println(element+" "+fks[element]);
//														System.out.println(element+" "+dis[element]);
						} else if(element==KMAX-1) {
							fks[element] = equationState.get(rheologyID[element]).pIntegral(x[element],y[element],parameterID[element],element) - ( equationState.get(rheologyID[element]).qIntegral(x_outer[element],y[element],parameterID[element],element) + equationState.get(rheologyID[element]).q(x_outer[element],y[element],parameterID[element],element)*(x[element]-x_outer[element]) ) - this.rhss[element] + lowerDiagonal[element]*x[element-1] + mainDiagonal[element]*x[element];
							dis[element] = ( equationState.get(rheologyID[element]).p(x[element],y[element],parameterID[element],element) - equationState.get(rheologyID[element]).q(x_outer[element],y[element],parameterID[element],element) );
//														System.out.println(element+" "+fks[element]);
//														System.out.println(element+" "+dis[element]);
						} else {
							fks[element] = equationState.get(rheologyID[element]).pIntegral(x[element],y[element],parameterID[element],element) - ( equationState.get(rheologyID[element]).qIntegral(x_outer[element],y[element],parameterID[element],element) + equationState.get(rheologyID[element]).q(x_outer[element],y[element],parameterID[element],element)*(x[element]-x_outer[element]) ) - this.rhss[element]  + lowerDiagonal[element]*x[element-1] + mainDiagonal[element]*x[element] + upperDiagonal[element]*x[element+1];
							dis[element] = ( equationState.get(rheologyID[element]).p(x[element],y[element],parameterID[element],element) - equationState.get(rheologyID[element]).q(x_outer[element],y[element],parameterID[element],element) );
//														System.out.println(element+" "+fks[element]);
//														System.out.println(element+" "+dis[element]);
						}

						innerResidual += fks[element]*fks[element];
					}
					innerResidual = Math.pow(innerResidual,0.5);

//															System.out.println("\t\t-Inner iteration " + j + " with residual " +  innerResidual);    

					if(innerResidual < newtonTolerance) {
						break;
					}

					//// THOMAS ALGORITHM////
					// Attention: the main diagonal of the coefficient matrix must not change!! The same for the upper diagonal

					bb = mainDiagonal.clone();
					cc = upperDiagonal.clone();
					for(int element = 0; element < KMAX; element++) {
						bb[element] += dis[element];
					}
					thomasAlg.set(cc,bb,lowerDiagonal,fks,KMAX);
					dx = thomasAlg.solver();

					//// UPDATE solution ////
					for(int element = 0; element < KMAX; element++) {
						x[element] = x[element] - dx[element];
						//						System.out.println(dx[s]);
					}
				}
			} //// INNER CYCLE END ////
		}
		return x;
	}

}
