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
public class ODESolidMass implements OrdinaryDifferentialEquation {

	private double initialCondition;
	private double snowfall;
	private double freezing;
	private double melting;
	
	
	public ODESolidMass(double initialCondition, double snowfall, double freezing, double melting) {
		
		this.initialCondition = initialCondition;
		this.snowfall = snowfall;
		this.freezing = freezing;
		this.melting = melting;
		
	}
	
	@Override
	public double compute(double solidWater) {
		// TODO Auto-generated method stub
		
		return solidWater - initialCondition - snowfall - freezing + melting;
	}

	@Override
	public double computeDerivative(double x) {
		// TODO Auto-generated method stub
		return 1;
	}

	@Override
	public double computeP(double x) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double computePIntegral(double x) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double computeRHS() {
		// TODO Auto-generated method stub
		return 0;
	}

}
