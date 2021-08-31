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

public interface OrdinaryDifferentialEquation {
	
		
	public abstract double compute(double x);
	
	public abstract double computeDerivative(double x);
	
	public abstract double computeP(double x);
		
	public abstract double computePIntegral(double x);
	
	public abstract double computeRHS();
	
	public default  double computeQ(double x) {
		
		return computeP(x) - computeDerivative(x);
				
	}

	public default double computeQIntegral(double x) {
		
		return computePIntegral(x) - compute(x);
		
	}

}
