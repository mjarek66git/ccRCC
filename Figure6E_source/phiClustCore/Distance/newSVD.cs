using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Distance
{
 

    class newSVD
    {
        const double thr = 1e-10;
        static double[] y = new double[3];
        static double[,] AA = new double[3, 3];
        static double[,] LDU = new double[3, 3];
        static int[] P = new int[3];
        static private double [,] U=new double [3,3];
        static double [] S=new double [3];
        static private double [,] V=new double [3,3];

        public double [,] GetV(){ return V;}
        public double [,] GetU(){return U;}
        void cross(double [,] o,int k, double[,] x,int i, double[,] y,int j)
        {

            o[k,0]=x[i,1]*y[j,2]-x[i,2]*y[j,1];
            o[k,1] = -(x[i,0] * y[j,2] - x[i,2] * y[j,0]);
            o[k,2] = x[i,0] * y[j,1] - x[i,1] * y[j,0];

        }
        void cross(double[,] o, int k, double[] x, double[,] y, int j)
        {

            o[k, 0] = x[1] * y[j, 2] - x[2] * y[j, 1];
            o[k, 1] = -(x[0] * y[j, 2] - x[2] * y[j, 0]);
            o[k, 2] = x[0] * y[j, 1] - x[1] * y[j, 0];

        }
        
        void sort3(double [] x) 
        {
            double tmp;

            if (x[0] < x[1]) 
            {
                tmp = x[0];
                x[0] = x[1];
                x[1] = tmp;
            }
            if (x[1] < x[2]) 
            {
                if (x[0] < x[2]) 
                {
                        tmp = x[2];
                        x[2] = x[1];
                        x[1] = x[0];
                        x[0] = tmp;
                }
                else 
                {
                        tmp = x[1];
                        x[1] = x[2];
                        x[2] = tmp;
                }
            }
        }
        void NormRows(double[,] x)
        {
            for (int i = 0; i < 3; i++)
            {
                double tmp = Math.Sqrt(x[i, 0] * x[i, 0] + x[i, 1] * x[i, 1] + x[i, 2] * x[i, 2]);
                x[i, 0] /= tmp;
                x[i, 1] /= tmp;
                x[i, 2] /= tmp;
            }
        }
        double [] ldubsolve3(double []y,double [,] LDU, int [] P) 
        {
  
            double [] x=new double [3];
            x[P[2]] = y[2];
            x[P[1]] = y[1] - LDU[P[2],1]*x[P[2]];
            x[P[0]] = y[0] - LDU[P[2],0]*x[P[2]] - LDU[P[1],0]*x[P[1]];

            return x;
        }
        double [,] matmul3(double [,]A, double [,] B) 
        {
            double [,]C=new double [3,3];
            C[0,0] = A[0,0]*B[0,0] + A[1,0]*B[0,1] + A[2,0]*B[0,2];
            C[1,0] = A[0,0]*B[1,0] + A[1,0]*B[1,1] + A[2,0]*B[1,2];
            C[2,0] = A[0,0]*B[2,0] + A[1,0]*B[2,1] + A[2,0]*B[2,2];
            C[0,1] = A[0,1]*B[0,0] + A[1,1]*B[0,1] + A[2,1]*B[0,2];
            C[1,1] = A[0,1]*B[1,0] + A[1,1]*B[1,1] + A[2,1]*B[1,2];
            C[2,1] = A[0,1]*B[2,0] + A[1,1]*B[2,1] + A[2,1]*B[2,2];
            C[0,2] = A[0,2]*B[0,0] + A[1,2]*B[0,1] + A[2,2]*B[0,2];
            C[1,2] = A[0,2]*B[1,0] + A[1,2]*B[1,1] + A[2,2]*B[1,2];
            C[2,2] = A[0,2]*B[2,0] + A[1,2]*B[2,1] + A[2,2]*B[2,2];

            return C;
        }

        double [] matvec3(double [,] A,double [,] x,int i) 
        {
            double [] y=new double [3];
            y[0] = A[0,0]*x[i,0] + A[1,0]*x[i,1] + A[2,0]*x[i,2];
            y[1] = A[0,1]*x[i,0] + A[1,1]*x[i,1] + A[2,1]*x[i,2];
            y[2] = A[0,2]*x[i,0] + A[1,2]*x[i,1] + A[2,2]*x[i,2];

            return y;
        }
        double [,] ata3(double [,] A) 
        {
            double [,] AA=new double [3,3];
            AA[0,0] = A[0,0]*A[0,0] + A[0,1]*A[0,1] + A[0,2]*A[0,2];
            AA[1,0] = A[0,0]*A[1,0] + A[0,1]*A[1,1] + A[0,2]*A[1,2];
            AA[2,0] = A[0,0]*A[2,0] + A[0,1]*A[2,1] + A[0,2]*A[2,2];
            AA[0,1] = AA[1,0];
            AA[1,1] = A[1,0]*A[1,0] + A[1,1]*A[1,1] + A[1,2]*A[1,2];
            AA[2,1] = A[1,0]*A[2,0] + A[1,1]*A[2,1] + A[1,2]*A[2,2];
            AA[0,2] = AA[2,0];
            AA[1,2] = AA[2,1];
            AA[2,2] = A[2,0]*A[2,0] + A[2,1]*A[2,1] + A[2,2]*A[2,2];
    
            return AA;
        }
        double [,] aat3(double [,] A) 
        {
            double [,] AA=new double [3,3];
            AA[0,0] = A[0,0]*A[0,0] + A[1,0]*A[1,0] + A[2,0]*A[2,0];
            AA[1,0] = A[0,0]*A[0,1] + A[1,0]*A[1,1] + A[2,0]*A[2,1];
            AA[2,0] = A[0,0]*A[0,2] + A[1,0]*A[1,2] + A[2,0]*A[2,2];
            AA[0,1] = AA[1,0];
            AA[1,1] = A[0,1]*A[0,1] + A[1,1]*A[1,1] + A[2,1]*A[2,1];
            AA[2,1] = A[0,1]*A[0,2] + A[1,1]*A[1,2] + A[2,1]*A[2,2];
            AA[0,2] = AA[2,0];
            AA[1,2] = AA[2,1];
            AA[2,2] = A[0,2]*A[0,2] + A[1,2]*A[1,2] + A[2,2]*A[2,2];

            return AA;
        }   
        void trans3(double [,] A) 
        {
            double tmp;
            tmp = A[1,0];
            A[1,0] = A[0,1];
            A[0,1] = tmp;
            tmp = A[2,0];
            A[2,0] = A[0,2];
            A[0,2] = tmp;
            tmp = A[2,1];
            A[2,1] = A[1,2];
            A[1,2] = tmp;
        }
        void solvecubic(double [] c) 
        {
            double sq3d2 = 0.86602540378443864676;
            double c2d3 = c[2]/3;
            double c2sq = c[2]*c[2];
            double Q = (3*c[1]-c2sq)/9;
            double R = (c[2]*(9*c[1]-2*c2sq)-27*c[0])/54;
            double tmp, t, sint, cost;
            if (Q < 0) 
            {
                tmp = 2*Math.Sqrt(-Q);
                t = Math.Acos(R/Math.Sqrt(-Q*Q*Q))/3;
                cost = tmp*Math.Cos(t);
                sint = tmp*Math.Sin(t);
                c[0] = cost - c2d3;
                cost = -0.5*cost - c2d3;
                sint = sq3d2*sint;
                c[1] = cost - sint;
                c[2] = cost + sint;
            }
            else 
            {
                tmp = Math.Pow(R, 1/3.0); 
                c[0] = -c2d3 + 2*tmp;
                c[1] = c[2] = -c2d3 - tmp;
            }
    }
    int [] ldu3(double [,] A) 
    {
	    int tmp;
        int []P=new int [3];
	    P[1] = 1;
	    P[2] = 2;

	    P[0] = Math.Abs(A[1,0]) > Math.Abs(A[0,0]) ? 
		(Math.Abs(A[2,0]) > Math.Abs(A[1,0]) ? 2 : 1) : 
		(Math.Abs(A[2,0]) > Math.Abs(A[0,0]) ? 2 : 0);
	    P[P[0]] = 0;

	    if (Math.Abs(A[P[2],1]) > Math.Abs(A[P[1],1])) 
        {
		    tmp = P[1];
		    P[1] = P[2];
		    P[2] = tmp;
    	}

	    if (A[P[0],0] != 0) 
        {
		    A[P[1],0] = A[P[1],0]/A[P[0],0];
		    A[P[2],0] = A[P[2],0]/A[P[0],0];
		    A[P[0],1] = A[P[0],1]/A[P[0],0];
		    A[P[0],2] = A[P[0],2]/A[P[0],0];
	    }

	    A[P[1],1] = A[P[1],1] - A[P[0],1]*A[P[1],0]*A[P[0],0];

	    if (A[P[1],1] != 0) 
        {
		    A[P[2],1] = (A[P[2],1] - A[P[0],1]*A[P[2],0]*A[P[0],0])/A[P[1],1];
		    A[P[1],2] = (A[P[1],2] - A[P[0],2]*A[P[1],0]*A[P[0],0])/A[P[1],1];
	    }

	    A[P[2],2] = A[P[2],2] - A[P[0],2]*A[P[2],0]*A[P[0],0] - A[P[1],2]*A[P[2],1]*A[P[1],1];

        return P;
    }
    public newSVD(double [,] A) 
    {        
        int k;
        double[] sTmp;
	/*
	 * Steps:
	 * 1) Use eigendecomposition on A^T A to compute V.
	 * Since A = U S V^T then A^T A = V S^T S V^T with D = S^T S and V the 
	 * eigenvalues and eigenvectors respectively (V is orthogonal).
	 * 2) Compute U from A and V.
	 * 3) Normalize columns of U and V and root the eigenvalues to obtain 
	 * the singular values.
	 */

	/* Compute AA = A^T A */
	AA=ata3(A);

	/* Form the monic characteristic polynomial */
	S[2] = -AA[0,0] - AA[1,1] - AA[2,2];
	S[1] = AA[0,0]*AA[1,1] + AA[2,2]*AA[0,0] + AA[2,2]*AA[1,1] - 
		AA[2,1]*AA[1,2] - AA[2,0]*AA[0,2] - AA[1,0]*AA[0,1];
	S[0] = AA[2,1]*AA[1,2]*AA[0,0] + AA[2,0]*AA[0,2]*AA[1,1] + AA[1,0]*AA[0,1]*AA[2,2] -
		AA[0,0]*AA[1,1]*AA[2,2] - AA[1,0]*AA[2,1]*AA[0,2] - AA[2,0]*AA[0,1]*AA[1,2];

	/* Solve the cubic equation. */
	solvecubic(S);

	/* All roots should be positive */
	if (S[0] < 0)
		S[0] = 0;
	if (S[1] < 0)
		S[1] = 0;
	if (S[2] < 0)
		S[2] = 0;

	/* Sort from greatest to least */
	sort3(S);

	/* Form the eigenvector system for the first (largest) eigenvalue */
    LDU=(double [,])AA.Clone();
	LDU[0,0] -= S[0];
	LDU[1,1] -= S[0];
	LDU[2,2] -= S[0];

	/* Perform LDUP decomposition */
	P=ldu3(LDU);

	/* 
	 * Write LDU = AA-I*lambda.  Then an eigenvector can be
	 * found by solving LDU x = LD y = L z = 0
	 * L is invertible, so L z = 0 implies z = 0
	 * D is singular since det(AA-I*lambda) = 0 and so 
	 * D y = z = 0 has a non-unique solution.
	 * Pick k so that D_kk = 0 and set y = e_k, the k'th column
	 * of the identity matrix.
	 * U is invertible so U x = y has a unique solution for a given y.
	 * The solution for U x = y is an eigenvector.
	 */

	/* Pick the component of D nearest to 0 */
	y[0] = y[1] = y[2] = 0;
	k = Math.Abs(LDU[P[1],1]) < Math.Abs(LDU[P[0],0]) ?
		(Math.Abs(LDU[P[2],2]) < Math.Abs(LDU[P[1],1]) ? 2 : 1) :
		(Math.Abs(LDU[P[2],2]) < Math.Abs(LDU[P[0],0]) ? 2 : 0);
	y[k] = 1;

	/* Do a backward solve for the eigenvector */
	sTmp=ldubsolve3(y, LDU, P);

    for (int i = 0; i < 3; i++)
        V[0, i] = sTmp[i];
	/* Form the eigenvector system for the last (smallest) eigenvalue */
    LDU=(double [,])AA.Clone();
	LDU[0,0] -= S[2];
	LDU[1,1] -= S[2];
	LDU[2,2] -= S[2];

	/* Perform LDUP decomposition */
	P=ldu3(LDU);

	/* 
	 * NOTE: The arrangement of the ternary operator output is IMPORTANT!
	 * It ensures a different system is solved if there are 3 repeat eigenvalues.
	 */

	/* Pick the component of D nearest to 0 */
	y[0] = y[1] = y[2] = 0;
	k = Math.Abs(LDU[P[0],0]) < Math.Abs(LDU[P[2],2]) ?
		(Math.Abs(LDU[P[0],0]) < Math.Abs(LDU[P[1],1]) ? 0 : 1) :
		(Math.Abs(LDU[P[1],1]) < Math.Abs(LDU[P[2],2]) ? 1 : 2);
	y[k] = 1;

	/* Do a backward solve for the eigenvector */
	sTmp=ldubsolve3(y, LDU, P);

    for (int i = 0; i < 3; i++)
        V[2, i] = sTmp[i];
	 /* The remaining column must be orthogonal (AA is symmetric) */
	cross(V,1,V,2, V,0);
	/* Count the rank */
    k = 0;
    for (int i = 0; i < 3; i++)
        if (S[i] > thr)
            k++;
	
	switch (k) {
		case 0:
			/*
			 * Zero matrix. 
			 * Since V is already orthogonal, just copy it into U.
			 */
            U = (double [,])V.Clone();
			break;
		case 1:
			/* 
			 * The first singular value is non-zero.
			 * Since A = U S V^T, then A V = U S.
			 * A V_1 = S_11 U_1 is non-zero. Here V_1 and U_1 are
			 * column vectors. Since V_1 is known, we may compute
			 * U_1 = A V_1.  The S_11 factor is not important as
			 * U_1 will be normalized later.
			 */
			    sTmp=matvec3(A, V,0);
                for (int i = 0; i < 3; i++)
                    U[0, i] = sTmp[i];

			/* 
			 * The other columns of U do not contribute to the expansion
			 * and we may arbitrarily choose them (but they do need to be
			 * orthogonal). To ensure the first cross product does not fail,
			 * pick k so that U_k1 is nearest 0 and then cross with e_k to
			 * obtain an orthogonal vector to U_1.
			 */
			y[0] = y[1] = y[2] = 0;
			k = Math.Abs(U[0,0]) < Math.Abs(U[0,2]) ?
				(Math.Abs(U[0,0]) < Math.Abs(U[0,1]) ? 0 : 1) :
				(Math.Abs(U[0,1]) < Math.Abs(U[0,2]) ? 1 : 2);
			y[k] = 1;


			cross(U,1,y,U,0);

			/* Cross the first two to obtain the remaining column */
			cross(U,2, U,0, U,1);
			break;
		case 2:
			/*
			 * The first two singular values are non-zero.
			 * Compute U_1 = A V_1 and U_2 = A V_2. See case 1
			 * for more information.
			 */
			sTmp=matvec3(A, V,0);
            for (int i = 0; i < 3; i++)
                U[0, i] = sTmp[i];
            sTmp=matvec3(A, V, 1);
            for (int i = 0; i < 3; i++)
                U[1, i] = sTmp[i];

			/* Cross the first two to obtain the remaining column */
			cross(U,2, U,0, U,1);
			break;
		case 3:
			/*
			 * All singular values are non-zero.
			 * We may compute U = A V. See case 1 for more information.
			 */
			U=matmul3(A, V);
			break;

	}

	/* Normalize the columns of U and V */
    NormRows(V);
    NormRows(U);
	/* S was initially the eigenvalues of A^T A = V S^T S V^T which are squared. */
	S[0] = Math.Sqrt(S[0]);
	S[1] = Math.Sqrt(S[1]);
	S[2] = Math.Sqrt(S[2]);
}

  }
}
