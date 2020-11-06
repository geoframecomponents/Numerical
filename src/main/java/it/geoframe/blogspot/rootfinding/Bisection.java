/*
 * GNU GPL v3 License
 *
 * Copyright 2019 Niccolo` Tubini
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
 * 
 */
package it.geoframe.blogspot.rootfinding;


import stateequation.StateEquation;

/**
 * @author Niccolo` Tubini
 *
 */

public class Bisection {


	private StateEquation stateEquation;
	private double fa;
	private double fb;
	private double fc;
	private double c;
	private double tolerance = 1e-11;
	private int counter;

	public Bisection(StateEquation stateEquation) {

		this.stateEquation = stateEquation;		

	}



	public double findZero(double a, double b, double y, int id, int element) {

		counter = 1;

		fa = stateEquation.ddStateEquation(a, y, id, element);
		fb = stateEquation.ddStateEquation(b, y, id, element);
		c = (a+b)/2;
		fc = stateEquation.ddStateEquation(c, y, id, element);

		if(fc == 0.0) {
			return c;
		} else {
			while(Math.abs(a-b) > tolerance) {


				if(fa*fb > 0.0) {
					System.out.println("\tBISECTION: error finding zero.");
				}
//				System.out.println("a "+a+" b "+b+" fa*fb "+fa*fb);
				c = (a+b)/2;
				fc = stateEquation.ddStateEquation(c, y, id, element);

				if(fc == 0.0) {
					return c;
				} else {
					if(fa*fc < 0.0) {
						b = c;
						fb = fc;
					} else {
						a = c;
						fa = fc;
					}
				}

				if(counter>150) {
					System.out.println("\tBISECTION: reached 250 iteration |a-b| = "+Math.abs(a-b));
					return c;
				}
				counter ++;
			}
		}
//		System.out.println("element: " +element+ "counter: "+counter+", intervall: "+Math.abs(a-b));

		return c;

	}
}
