// Author: 			Jennifer Brown
// Revision History: 	Created Oct. 2012
//
// Purpose:			Takes matrices A, B, and calculates:
// 					exp(-t(A+B))= e^(-tA) e^(-tB) e^(t^2 [B,A]/2)
// 					where [B,A] = BA-AB
// 					This is a Trotter Approximation, see Raedt 1983 
// 					
