/*
 * GNU GPL v3 License
 *
 * Copyright 2021 Niccolò Tubini, Giuseppe Formetta, Riccardo Rigon
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

package it.geoframe.blogspot.numerical.ode;


/**
 * @author Niccolò Tubini, Giuseppe Formetta
 * 
 */
public class NestedNewton {

	private double x0;
	
	private double x_k;
	
	private double tol = 1e-10;
	
	private double error;
	
	private int iteration;
	
	private int maxIteration = 100;
	
	private double f;
	
	private double ff;
	
	public double solve( double x, OrdinaryDifferentialEquation ode ) {
		
		this.x0 = x;
		
		error = 1.0;
		
		iteration = 1;
				
		
		for(int outerIteration=0; outerIteration<maxIteration; outerIteration++) {
			
			f = ode.compute(x0) - ode.computeRHS();
			
			if(Math.abs(f)<tol) {
				
				return x0;
				
			}
			
			x_k = x0;
			
			for(int innerrIteration=0; innerrIteration<maxIteration; innerrIteration++) {
				
				ff = ode.computePIntegral(x0) - (ode.computeQIntegral(x_k)+ode.computeQ(x_k)*(x0-x_k)) - ode.computeRHS();
				
				if(Math.abs(ff)<tol) {
					
					return x0;
					
				}
				
				x0 = x0 - 1/(ode.computeP(x0)-ode.computeQ(x_k))*ff;
			}
		}

		return x0;
	}
	
}
