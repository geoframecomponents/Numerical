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

package it.geoframe.blogspot.numerical.linearsystemsolver;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

//import bidimensionalDomain.Geometry;
//import bidimensionalDomain.Topology;
import it.geoframe.blogspot.numerical.matop.*;

/**
 * <h1>Conjugate gradient method </h1>
 * Matrix-free conjugate method for the solution of A*x=b
 * A is symmetric and positive definite.
 * <p>
 * 
 * We do not need to store the matrix, just to provide a class which computes
 * the matrix-vector product A*x
 * This implementation works for the 1D, 2D as well as 3D case. What changes, is the 
 * class that computes the matrix-vector.
 * 
 * With preconditioner
 * 
 * @author Niccolo' Tubini, Riccardo Rigon, Michael Dumbser
 * @version 0.1
 * @since 2018-11-23
 * @see <a href="https://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf">J. R. Shewchuk, An Introduction to the
 *	Conjugate Gradient Method Without the Agonizing Pain</a>
 */

public class ConjugateGradientMethod {
	
	private double alpha;
	private double alphak;
	private double lambda;
	private double tmp;
	private double cgTolerance;
	
	private List<Double> residual;
	private List<Double> x;
	private List<Double> p;
	private List<Double> Apsi;
	
	private List<Integer> elements;

	private Matop matop;
	
	public ConjugateGradientMethod(Matop matop, double cgTolerance, List<Integer> elements) {
		
		residual = new ArrayList<Double>(Arrays.asList(new Double[elements.size()]));
		x = new ArrayList<Double>(Arrays.asList(new Double[elements.size()]));
		p = new ArrayList<Double>(Arrays.asList(new Double[elements.size()]));
		Apsi = new ArrayList<Double>(Arrays.asList(new Double[elements.size()]));
		this.elements = new ArrayList<Integer>(elements);
		this.matop = matop;
		this.cgTolerance = cgTolerance;
	}

	public List<Double> solve(List<Double> dis, List<Double> rhs, List<Double> mainDiagonal){

		x = new ArrayList<Double>(rhs);
//				System.out.println("x:");
//				for(int element=1; element<elements.size(); element++) {
//					System.out.println("\t"+ element + "\t" + x.get(element));
//				}
		Apsi = matop.solve(dis, x);
//				System.out.println("Apsi primo:");
//				for(int element=1; element<elements.size(); element++) {
//					System.out.println("\t"+ element + "\t" + Apsi.get(element)+ "\t" + dis.get(element)+ "\t" + x.get(element) );
//				}
		alpha = 0.0;
		for(int element=1; element<elements.size(); element++) {
			residual.set(element, rhs.get(element)-Apsi.get(element)); 
			// 			no preconditioner
			//			p.put(element, residual.get(element));
			// With preconditioner
			p.set(element, residual.get(element)/(mainDiagonal.get(element)+dis.get(element)));
			alpha += p.get(element)*residual.get(element);
		}
		int iter =1;
		for(int k=1; k<4*x.size();k++) {

			if(Math.sqrt(alpha)<=cgTolerance) {
//				System.out.println("\t\t\t\t\tk: " + k +"\tsqrt(alpha): " +Math.sqrt(alpha));
				break;
			}

			tmp = 0.0;
//					System.out.println("p:");
//					for(int element=1; element<elements.size(); element++) {
//						System.out.println("\t"+ element + "\t" + p.get(element));
//					}
			Apsi = matop.solve(dis, p);
//					System.out.println("Apsi :");
//					for(int element=1; element<elements.size(); element++) {
//						System.out.println("\t"+ element + "\t" + Apsi.get(element));
//					}
			tmp = 0.0;
			for(int element=1; element<elements.size(); element++) {
				tmp += p.get(element)*Apsi.get(element);
			}
			lambda = alpha/tmp;
//			System.out.println("\t\t\t\ttmp: " + tmp + "\tlambda: " + lambda);

			alphak = alpha;
			alpha = 0.0;
			for(int element=1; element<elements.size(); element++) {
				tmp = x.get(element) + lambda*p.get(element);
				x.set(element, tmp);
				tmp = residual.get(element) - lambda*Apsi.get(element);
				residual.set(element, tmp);
				// 				no preconditioner
				//				alpha += tmp*tmp;
				// with preconditioner
				alpha += tmp/(mainDiagonal.get(element)+dis.get(element))*tmp;
			}
			
//			System.out.println("\t\t\t\t\talpha: " + alpha + "\talphak: " + alphak);
			for(int element=1; element<elements.size(); element++) {
				// 				no preconditioner
				//				tmp = residual.get(element)+alpha/alphak*p.get(element);
				// with preconditioner
				tmp = residual.get(element)/(mainDiagonal.get(element)+dis.get(element))+alpha/alphak*p.get(element);
				p.set(element, tmp );
			}

			iter++;
		}
//		System.out.println("\t\t\t\t\tk: " + iter +"\tsqrt(alpha): " +Math.sqrt(alpha));

		return x;
	}

}
