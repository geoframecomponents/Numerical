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
public class ODELiquidMass implements OrdinaryDifferentialEquation {

	private double initialCondition;
	private double rainfall;
	private double freezing;
	private double melting;
	private double snowPorosity;
	
	
	public ODELiquidMass(double initialCondition, double rainfall, double freezing, double melting, double snowPorosity) {
		
		this.initialCondition = initialCondition;
		this.rainfall = rainfall;
		this.freezing = freezing;
		this.melting = melting;
		this.snowPorosity = snowPorosity;
		
	}
	
	@Override
	public double compute(double liquidWater) {
		// TODO Auto-generated method stub

		return Math.min(liquidWater,snowPorosity) + Math.max(liquidWater-snowPorosity, 0) - initialCondition - rainfall + freezing - melting;
		
	}

	@Override
	public double computeDerivative(double liquidWater) {
		// TODO Auto-generated method stub
		
//		if(liquidWater>snowPorosity) {
//			
//			return 2;
//			
//		} else {
			
			return 1;
			
//		}
		
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
