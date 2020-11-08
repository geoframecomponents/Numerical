package it.geoframe.blogspot.numerical.newtonalgorithm;

import java.util.HashMap;
import java.util.Map;

import bidimensionalDomain.Geometry;
import bidimensionalDomain.Topology;
import it.geoframe.blogspot.numerical.linearsystemsolver.*;
import it.geoframe.blogspot.numerical.matop.*;

public class NestedNewtonCG {

	private double outerResidual;
	private double innerResidual;

	int nestedNewton;
	int MAXITER_NEWT;

	double newtonTolerance;
	double tmp;
	double timeDelta;


	Map<Integer, Double> rhss;
	Map<Integer, Double> mainDiagonal;
	Map<Integer, Double> variable;
	
	private Map<Integer, Double> Apsi;
	private Map<Integer, Double> fs;
	private Map<Integer, Double> fks;
	private Map<Integer, Double> dis;
	private Map<Integer, Double> dpsis;
	private Map<Integer, Double> psis_outer;
	private Map<Integer, Double> psism;

	Topology topology;
	Geometry geometry;
//	SoilWaterRetentionCurve swrc;
//	//TotalDepth totalDepth;
//	//Thomas thomasAlg = new Thomas();
//
	Matop matop;
	ConjugateGradientMethod cg;



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
	public NestedNewtonCG(int nestedNewton, double newtonTolerance, int MAXITER_NEWT, SoilWaterRetentionCurve swrc, Matop matop, double cgTolerance){

		this.nestedNewton = nestedNewton;
		this.newtonTolerance = newtonTolerance;
		this.MAXITER_NEWT = MAXITER_NEWT;
		this.topology = Topology.getInstance();
		this.geometry = Geometry.getInstance();
//		this.swrc = swrc;
//		this.matop = matop;
		this.cg = new ConjugateGradientMethod(matop, cgTolerance);


		fs			  = new HashMap<Integer, Double>();
		fks			  = new HashMap<Integer, Double>();
		dis			  = new HashMap<Integer, Double>();
		dpsis		  = new HashMap<Integer, Double>();
		psis_outer	  = new HashMap<Integer, Double>();
		psism          = new HashMap<Integer, Double>();
	}



	/**
	 * @param psis vector contains the suction values, it is a vector of length NUM_CONTROL_VOLUMES
	 * @param upperDiagonal upper diagonal of the coefficient matrix A of the linear system, it is a vector of length NUM_CONTROL_VOLUMES
	 * @param mainDiagonal main diagonal of the coefficient matrix A of the linear system, it is a vector of length NUM_CONTROL_VOLUMES
	 * @param lowerDiagonal lower diagonal of the coefficient matrix A of the linear system, it is a vector of length NUM_CONTROL_VOLUMES
	 * @param rhss right hand side term of the linear system, it is a vector of length NUM_CONTROL_VOLUMES
	 */
	public void set(Map<Integer, Double> variable, Map<Integer, Double> rhss, Map<Integer, Double> mainDiagonal){

		this.variable = variable;
		this.rhss = rhss;
		this.mainDiagonal = mainDiagonal;

	}



	public Map<Integer, Double> solver(){



		// Initial guess for the outer iteration
		for(Integer element : topology.s_i.keySet()) {
			tmp = Math.min(variable.get(element), SoilParameters.psiStar1[SoilParameters.elementsLabel.get(element)]);
//			Variables.waterSuctions.put(element, tmp);
		}

		//// OUTER CYCLE ////
		for(int i = 0; i < MAXITER_NEWT; i++) {
			// I have to assign 0 to outerResidual otherwise I will take into account of the previous error
			outerResidual = 0.0;
			for(Integer element : topology.s_i.keySet()) {
				dis.put(element, 0.0);
			}
			Apsi = matop.solve(dis, variable);
//			System.out.println("Apsi outer:");
//			for(Integer element : Topology.s_i.keySet()) {
//				System.out.println("\t"+ element + "\t" + Apsi.get(element));
//			}
			
//			System.out.println("fs :");
			for(Integer element : topology.s_i.keySet()) {
				tmp =  swrc.dWaterContent(variable.get(element),element)*geometry.elementsArea.get(element);
				dis.put(element, tmp );
				tmp = swrc.waterContent(variable.get(element),element)*geometry.elementsArea.get(element) - rhss.get(element) + Apsi.get(element);
				fs.put(element, tmp);
//				System.out.println(element +": " +swrc.waterContent(Variables.waterSuctions.get(element),element) +" "+Geometry.elementsArea.get(element) +" "+rhss.get(element)+" "+ Apsi.get(element)+" " +fs.get(element) + " " + dis.get(element));
//				System.out.println("\t"+ element + "\t" + fs.get(element));
		
				outerResidual += tmp*tmp;
			}
			outerResidual = Math.pow(outerResidual,0.5);  
			System.out.println("\t\t-Outer iteration " + i + " with residual " +  outerResidual);
			if(outerResidual < newtonTolerance) {
				break;
			}
			if(nestedNewton == 0){

			}else{

				// Initial guess for the inner iteration (optional)
				for(Integer element : topology.s_i.keySet()) {
//					psis_outer.put(element, Variables.waterSuctions.get(element));
//					Variables.waterSuctions.put(element, Math.max(Variables.waterSuctions.get(element), SoilParameters.psiStar1[SoilParameters.elementsLabel.get(element)]));
				}

				//// INNER CYCLE ////
				for(int j = 0; j < MAXITER_NEWT; j++) {
					// I have to assign 0 to innerResidual otherwise I will take into account of the previous error
					innerResidual = 0.0; 
					for(Integer element : topology.s_i.keySet()) {
						dis.put(element, 0.0);
					}
					Apsi = matop.solve(dis, variable);
//					System.out.println("Apsi inner:");
//					for(Integer element : Topology.s_i.keySet()) {
//						System.out.println("\t"+ element + "\t" + Apsi.get(element));
//					}
//					System.out.println("inner psi-psi_outer :");
					for(Integer element : topology.s_i.keySet()) {

						tmp = (swrc.p(variable.get(element),element) - swrc.q(psis_outer.get(element),element))*geometry.elementsArea.get(element);
						dis.put(element, tmp);
						tmp = swrc.pIntegral(variable.get(element),element)*geometry.elementsArea.get(element)
								- ( swrc.qIntegral(psis_outer.get(element),element) + swrc.q(psis_outer.get(element),element)*(variable.get(element) - psis_outer.get(element)) )*Geometry.elementsArea.get(element)
								- rhss.get(element) + Apsi.get(element);
						fks.put(element, tmp);

//						System.out.println(element +" "+Variables.waterSuctions.get(element) + " " + psis_outer.get(element));
//						System.out.println("\t"+ element + "\t" + fks.get(element));

						innerResidual += tmp*tmp;
					}

					innerResidual = Math.pow(innerResidual,0.5);
					System.out.println("\t\t\t-Inner iteration " + j + " with residual " +  innerResidual);    
					if(innerResidual < newtonTolerance) {
//						System.out.println("Psi:");
//						for(Integer element : Topology.s_i.keySet()) {
//							System.out.println("\t"+ element + "\t" + Variables.waterSuctions.get(element));
//
//						}
//
//						System.out.println("\n\n");
//						System.out.println("dpsi:");
//						for(Integer element : Topology.s_i.keySet()) {
//							System.out.println("\t"+ element + "\t" + dpsis.get(element));
//
//						}
//						
						break;
					}
					//// CONJUGATE GRADIENT METHOD////
					dpsis = cg.solve(dis, fks, mainDiagonal);
//					dpsis = cg.solve(dis, fks);

					for(Integer element : topology.s_i.keySet()) {
//						psism.put(element, Variables.waterSuctions.get(element));
						tmp  = variable.get(element)-dpsis.get(element);
						variable.put(element, tmp );
					}

//					if(j>1) {
//						for(Integer element : Topology.s_i.keySet()) {
//							Variables.waterSuctions.put(element, Math.min( Variables.waterSuctions.get(element), psism.get(element)) );
//							if(i>1) {
//									Variables.waterSuctions.put(element, Math.max( Variables.waterSuctions.get(element), psis_outer.get(element)) );
//							}
//						}
//					}
				} //// INNER CYCLE END ////

			}
			//		return variable;
//			if(i>1) {
//				for(Integer element : Topology.s_i.keySet()) {
//					Variables.waterSuctions.put(element, Math.max( Variables.waterSuctions.get(element), psis_outer.get(element)) );
//				}
//			}
			return variable;
		} //// OUTER CYCLE END ////

	}


}

