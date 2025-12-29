/*
 * GNU GPL v3 License
 *
 * Copyright 2021 Niccolo Tubini, Giuseppe Formetta
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
public class SolveODE {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
	
		double snowPorosity = -9999.0;
		double solidWater = -9999.0;
		double liquidWater = -9999.0;
		double q = -9999.0;
		
		double solidWaterInitialCondition = 0.1;
		
		double liquidWaterInitialCondition = 0.0;
		
		double rainfall = 0.05;
		
		double snowfall = 0.95;
				
		double freezing = 0;
		
		double melting = 0.498843205767106;
		
		

		NewtonRaphson newton = new NewtonRaphson();
		
		
		OrdinaryDifferentialEquation odeSolidWater = new ODESolidMass(solidWaterInitialCondition, snowfall, freezing, melting);
		
		solidWater = newton.solve(0, odeSolidWater);
		System.out.println("error ODE solid water " + (solidWater  - solidWaterInitialCondition - snowfall - freezing + melting));

		snowPorosity = solidWater * 0.3504315;
				
		OrdinaryDifferentialEquation odeLiquidWater = new ODELiquidMass(liquidWaterInitialCondition, rainfall, freezing, melting, snowPorosity);

		liquidWater = newton.solve(0, odeLiquidWater);
		System.out.println("error ODE liquid water " + (Math.min(liquidWater,snowPorosity) + Math.max(liquidWater-snowPorosity, 0) - liquidWaterInitialCondition - rainfall + freezing - melting));

		if(liquidWater>snowPorosity) {
			q = liquidWater-snowPorosity;
			liquidWater = snowPorosity;
		}
		
		
		System.out.println("solidWater " + solidWater + "\t" + " liquidWater " + "\t" + liquidWater + "\t" + " melting discharge" + "\t" + q);
		System.out.println("error " + (rainfall+snowfall - (solidWater-solidWaterInitialCondition) - (liquidWater-liquidWaterInitialCondition) - q));

		
	}

}
