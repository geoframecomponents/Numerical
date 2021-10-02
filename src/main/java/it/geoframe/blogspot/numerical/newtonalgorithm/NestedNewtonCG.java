/*
 * GNU GPL v3 License
 *
 * Copyright 2021 Niccolo` Tubini
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

package it.geoframe.blogspot.numerical.newtonalgorithm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import it.geoframe.blogspot.closureequation.equationstate.EquationState;
import it.geoframe.blogspot.numerical.linearsystemsolver.*;
import it.geoframe.blogspot.numerical.matop.*;

public class NestedNewtonCG {

	private double outerResidual;
	private double innerResidual;

	private int MAXITER_NEWT;

	private double newtonTolerance;
	private double tmp;


	private List<Double> rhss;
	private List<Double> mainDiagonal;
	private List<Double> x;
	private List<Double> y;
	private List<Integer> elementEquationStateID;
	private List<Integer> elementParameterID;
	private List<EquationState> equationState;
	private List<Double> Apsi;
	private List<Double> fs;
	private List<Double> fks;
	private List<Double> dis;
	private List<Double> dx;
	private List<Double> x_outer;


	private Matop matop;
	private ConjugateGradientMethod cg;



	/**
	 * @param nestedNewton control parameter to choose between simple Newton method (0), or the nested Newton one (1)
	 * @param newtonTolerance prefixed tolerance representing the maximum mass balance error allowed  
	 * @param MAXITER_NEWT prefixed maximum number of iteration
	 * @param NUM_CONTROL_VOLUMES number of control volumes
	 * @param soilPar is the class to compute the soil hydraulic properties
	 * @param totalDepth is the class to compute the total water depth
	 * @param par1SWRC vector containing the first parameter of the SWRC, it is a vector of length NUM_CONTROL_VOLUMES-1
	 * @param par2SWRC vector containing the second parameter of the SWRC, it is a vector of length NUM_CONTROL_VOLUMES-1
	 * @param thetaR vector containing the adimensional residual water contentfor each control volume, it is a vector of length NUM_CONTROL_VOLUMES-1
	 * @param thetaS vector containing the adimensional water content at saturation for each control volume, it is a vector of length NUM_CONTROL_VOLUMES-1
	 */
	public NestedNewtonCG(double newtonTolerance, int MAXITER_NEWT, List<EquationState> equationState, Matop matop, double cgTolerance,
			List<Integer> elementParameterID, List<Integer> elementEquationStateID) {

		this.newtonTolerance = newtonTolerance;
		this.MAXITER_NEWT = MAXITER_NEWT;
		this.equationState = equationState;
		this.elementParameterID = elementParameterID;
		this.elementEquationStateID = elementEquationStateID;
		this.cg = new ConjugateGradientMethod(matop, cgTolerance, elementParameterID);

		this.matop = matop;


		fs		  = new ArrayList<Double>(Arrays.asList(new Double[elementParameterID.size()]));
		fks		  = new ArrayList<Double>(Arrays.asList(new Double[elementParameterID.size()]));
		dis		  = new ArrayList<Double>(Arrays.asList(new Double[elementParameterID.size()]));
		dx		  = new ArrayList<Double>(Arrays.asList(new Double[elementParameterID.size()]));
		x_outer	  = new ArrayList<Double>(Arrays.asList(new Double[elementParameterID.size()]));

	}



	/**
	 * @param psis vector contains the suction values, it is a vector of length NUM_CONTROL_VOLUMES
	 * @param upperDiagonal upper diagonal of the coefficient matrix A of the linear system, it is a vector of length NUM_CONTROL_VOLUMES
	 * @param mainDiagonal main diagonal of the coefficient matrix A of the linear system, it is a vector of length NUM_CONTROL_VOLUMES
	 * @param lowerDiagonal lower diagonal of the coefficient matrix A of the linear system, it is a vector of length NUM_CONTROL_VOLUMES
	 * @param rhss right hand side term of the linear system, it is a vector of length NUM_CONTROL_VOLUMES
	 */
	public void set(List<Double> x, List<Double> y, List<Double> rhss, List<Double> mainDiagonal) {

		this.x = new ArrayList<Double>(x);
		this.y = new ArrayList<Double>(y);
		this.rhss = new ArrayList<Double>(rhss);
		this.mainDiagonal = new ArrayList<Double>(mainDiagonal);

	}



	public List<Double> solver() {



		// Initial guess for the outer iteration
		for(int element=1; element<elementParameterID.size(); element++) {
			tmp = equationState.get(elementEquationStateID.get(element)).initialGuess(x.get(element), elementParameterID.get(element), element);
			x.set(element, tmp);
		}

		//// OUTER CYCLE ////
		for(int i = 0; i < MAXITER_NEWT; i++) {
			// I have to assign 0 to outerResidual otherwise I will take into account of the previous error
			outerResidual = 0.0;
			for(int element=1; element<elementParameterID.size(); element++) {

				dis.set(element, 0.0);

			}

			Apsi = matop.solve(dis, x);
//						System.out.println("Apsi outer:");
//						for(int element=1; element<elementParameterID.size(); element++) {
//							System.out.println("\t"+ element + "\t" + Apsi.get(element));
//						}
//
//						System.out.println("fs :");
			for(int element=1; element<elementParameterID.size(); element++) {

				tmp = equationState.get(elementEquationStateID.get(element)).dEquationState(x.get(element), y.get(element), elementParameterID.get(element), element);
				dis.set(element, tmp);

				tmp = equationState.get(elementEquationStateID.get(element)).equationState(x.get(element), y.get(element), elementParameterID.get(element), element) - rhss.get(element) + Apsi.get(element);
				fs.set(element, tmp);
//								System.out.println("\t"+ element + "\t" + fs.get(element));

				outerResidual += tmp*tmp;
			}

			outerResidual = Math.pow(outerResidual,0.5);  
//			System.out.println("\t\t-Outer iteration " + i + " with residual " +  outerResidual);
			if(outerResidual < newtonTolerance) {

				break;

			}


			for(int element=1; element<elementParameterID.size(); element++) {

				x_outer.set(element, x.get(element));

			}

			//// INNER CYCLE ////
			for(int j = 0; j < MAXITER_NEWT; j++) {
				// I have to assign 0 to innerResidual otherwise I will take into account of the previous error
				innerResidual = 0.0; 
				for(int element=1; element<elementParameterID.size(); element++) {
					dis.set(element, 0.0);
				}
				Apsi = matop.solve(dis, x);
//									System.out.println("Apsi inner:");
//									for(int element=1; element<elementParameterID.size(); element++) {
//										System.out.println("\t"+ element + "\t" + Apsi.get(element));
//									}
//
//									System.out.println("fks :");
				for(int element=1; element<elementParameterID.size(); element++) {

					tmp = equationState.get(elementEquationStateID.get(element)).p(x.get(element), y.get(element), elementParameterID.get(element), element) - equationState.get(elementEquationStateID.get(element)).q(x_outer.get(element), y.get(element), elementParameterID.get(element), element); 
//					System.out.println("\t\telement "+element+"\t"+x.get(element)+"\t"+equationState.get(elementEquationStateID.get(element)).p(x.get(element), y.get(element), elementParameterID.get(element), element)+"\t"+equationState.get(elementEquationStateID.get(element)).q(x_outer.get(element), y.get(element), elementParameterID.get(element), element));
					dis.set(element, tmp);

					tmp = equationState.get(elementEquationStateID.get(element)).pIntegral(x.get(element), y.get(element), elementParameterID.get(element), element)
							- ( equationState.get(elementEquationStateID.get(element)).qIntegral(x_outer.get(element), y.get(element), elementParameterID.get(element), element) + equationState.get(elementEquationStateID.get(element)).q(x_outer.get(element), y.get(element), elementParameterID.get(element), element)*(x.get(element) - x_outer.get(element)) )
							- rhss.get(element) + Apsi.get(element);
					fks.set(element, tmp);
//											System.out.println("\t"+ element + "\t" + fks.get(element));

					innerResidual += tmp*tmp;
				}

				innerResidual = Math.pow(innerResidual,0.5);
//				System.out.println("\t\t\t-Inner iteration " + j + " with residual " +  innerResidual);    
				if(innerResidual < newtonTolerance) {

					break;

				}


				/* 
				 * CONJUGATE GRADIENT METHOD
				 */
				dx = cg.solve(dis, fks, mainDiagonal);
//					System.out.println("CG");
//					System.out.println("dx psi");
				for(int element=1; element<elementParameterID.size(); element++) {

					tmp  = x.get(element)-dx.get(element);
					x.set(element, tmp);
//					System.out.println(dx.get(element) +" "+ tmp);

				}

			} //// OUTER CYCLE END ////

		}
		
		return x;


	}

}

