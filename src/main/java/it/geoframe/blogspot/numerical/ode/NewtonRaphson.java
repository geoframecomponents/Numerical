/*
 * GNU GPL v3 License
 *
 * Copyright 2021 Niccolo Tubini, Giuseppe Formetta, Riccardo Rigon
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
 * @author Niccolo Tubini, Giuseppe Formetta
 * 
 */
public class NewtonRaphson implements NonLinearEquationSolver {

	private double x0;
	
	private double xNew;
	
	private double tol = 1e-10;
	
	private double error;
	
	private int iteration;
	
	private int maxIteration = 100;
	
	
	
	public double solve( double x, OrdinaryDifferentialEquation ode ) {
		
		this.x0 = x;
		
		error = 1.0;
		
		iteration = 1;
		
		
		while(error>tol) {
		
			if(iteration>maxIteration) {
				
				System.out.println("\t NewtonRaphson reach maximum numeber of iteration " + ode.toString());
				break;
				
			}
			
			xNew = -ode.compute(this.x0)/ode.computeDerivative(this.x0) + this.x0;
			
			error = Math.abs(xNew-this.x0)/this.x0;
			
			this.x0 = xNew;
			
			iteration++;
				
		}
		
		return xNew;
	}
	
}
