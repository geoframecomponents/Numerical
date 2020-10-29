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

package linearsystemsolver;

import java.util.HashMap;
import java.util.Map;

import bidimensionalDomain.Geometry;
import bidimensionalDomain.Topology;
import matop.*;

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
	
	double alpha;
	double alphak;
	double lambda;
	double tmp;
	double cgTolerance;
	Map<Integer, Double> residual;
	Map<Integer, Double> x;
	Map<Integer, Double> p;
	Map<Integer, Double> Apsi;
	Topology topology;
	Geometry geometry;
	Matop matop;
	
	public ConjugateGradientMethod(Matop matop, double cgTolerance) {
		
		residual = new HashMap<Integer, Double>();
		x = new HashMap<Integer, Double>();
		p = new HashMap<Integer, Double>();
		Apsi = new HashMap<Integer, Double>();
		this.matop = matop;
		this.cgTolerance = cgTolerance;
	}

	public Map<Integer, Double> solve(Map<Integer, Double> dis, Map<Integer, Double> b, Map<Integer, Double> mainDiagonal){

		x = b;

		Apsi = matop.solve(dis, x);
		alpha = 0.0;
		for(Integer element : topology.s_i.keySet()) {
			residual.put(element, b.get(element)-Apsi.get(element)); 
			// 			no preconditioner
			//			p.put(element, residual.get(element));
			// With preconditioner
			p.put(element, residual.get(element)/(mainDiagonal.get(element)+dis.get(element)));
			alpha += p.get(element)*residual.get(element);
		}
		int iter =1;
		for(int k=1; k<4*x.size();k++) {

			if(Math.sqrt(alpha)<=cgTolerance) {
//				System.out.println("\t\t\t\t\tk: " + k +"\tsqrt(alpha): " +Math.sqrt(alpha));
				break;
			}

			tmp = 0.0;

			Apsi = matop.solve(dis, p);

			tmp = 0.0;
			for(Integer element : topology.s_i.keySet()) {
				tmp += p.get(element)*Apsi.get(element);
			}
			lambda = alpha/tmp;

			alphak = alpha;
			alpha = 0.0;
			for(Integer element : topology.s_i.keySet()) {
				tmp = x.get(element) + lambda*p.get(element);
				x.put(element, tmp);
				tmp = residual.get(element) - lambda*Apsi.get(element);
				residual.put(element, tmp);
				// 				no preconditioner
				//				alpha += tmp*tmp;
				// with preconditioner
				alpha += tmp/(mainDiagonal.get(element)+dis.get(element))*tmp;
			}
			
			for(Integer element : topology.s_i.keySet()) {
				// 				no preconditioner
				//				tmp = residual.get(element)+alpha/alphak*p.get(element);
				// with preconditioner
				tmp = residual.get(element)/(mainDiagonal.get(element)+dis.get(element))+alpha/alphak*p.get(element);
				p.put(element, tmp );
			}

			iter++;
		}
		System.out.println("\t\t\t\t\tk: " + iter +"\tsqrt(alpha): " +Math.sqrt(alpha));

		return x;
	}

}
