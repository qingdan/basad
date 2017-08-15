//
//  tPrior
//  Created by Qingyan Xiang on 4/19/17.
//  Copyright © 2017 Qingyan Xiang. All rights reserved.
//

#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <math.h>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "utiliti.h"


#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define P1 10e-4
#define P2 10e-4
#define beta1 2
#define beta2 5

using namespace std;
using namespace Eigen;
using namespace Rcpp;




extern "C"{
  void basadT(double* X_, double* Y_, double* Z0, double* B0, double* sig_, double* pr_,
             int* n_, int* p_, double* nu_,double* s0_,double* s1_, int* nburn_, int* niter_,
             int* nsplit_, double* outZ,  double* outbeta, double* outSig, int *Fast){

      
    int n=*n_, p=*p_, nburn=*nburn_, niter=*niter_;
    double pr=*pr_, sig =*sig_, nu = *nu_;
    
    MatrixXd X(n,p+1);
    VectorXd Y(n);
    
    int i, j;
    
    for (i=0;i<n;i++){
      for (j=1;j<p+1;j++) {
        X(i,j)=X_[j*n+i];
      }
      Y(i)=Y_[i];
      X(i,0) = 1.0;
    }
    
    long idum = -time(0);
    
    int itr, nsplit =*nsplit_,vsize = (p+1)/nsplit;
    double s0 = *s0_, s1 = *s1_, a, b;
    
    
    
    MatrixXd B(niter+nburn, p+1), Z(niter+nburn, p+1), G(p+1, p+1), COV(vsize, vsize), COVsq(vsize, vsize), tempG(vsize, p+1-vsize), COVeg(vsize, vsize), COV2(p+1, p+1), COVeg2(p+1, p+1), COVsq2(p+1, p+1);
	
    VectorXd tmu(p+1), T1(p+1), tempgas2(p+1), tempgas2n(n), delta(n),  prob(p+1), sigma(niter + nburn), vec1(n), vec2(p+1), T1Inver(p+1),  tempDelta(p+1);
    
      
    
    //allocation for new faster algorithm
      
    MatrixXd Phi(n, p + 1), diagN(n, n), tempA(n,n), D(p+1, p+1);
    VectorXd D_diag(p+1), alpha(n), u(p+1), v(n), w(n), tempB(n);
	
	
	//allocation for tPrior
	
	MatrixXd t(niter + nburn, p+1);
      
    //MatrixXd
      
    //Initial Values

    for(j=1;j<p+1;j++){
		B(0,j)=B0[j-1];
		Z(0,j)=Z0[j-1];
    }
    B(0,0) = 0, Z(0,0) = 1;
	
	sigma(0) = sig;
    for( itr = 1; itr < (niter + nburn) ; itr++ )
        sigma(itr) = 1.0;
	for( itr = 0; itr < (niter + nburn) ; itr++ )
        for( i = 0; i < p+1; i++ )
			t(itr, i) = 1.0;
    
    
    //Gibss
      
      
    G = X.transpose() * X;
    tmu = X.transpose() * Y;
      
    
      
    for( itr = 1; itr <(nburn + niter); itr++ ){
        
		//cout << itr << endl;
		sig = sigma(itr - 1);
      
        //updating B
      
	    //original way
        if( *Fast == 0 ){
	   
			B.row(itr) = B.row(itr - 1);
			for( j = 0; j < p+1; j++)
				T1(j) = 1/( Z(itr-1, j) * s1 * t(itr-1, j) + ( 1 - Z(itr - 1, j) ) * s0 * t(itr-1, j));

			COV2 = G;
      
			for(i=0;i<p+1;i++)
				COV2(i,i) += T1(i);
      
			SelfAdjointEigenSolver<MatrixXd> eigensolver2(COV2);
      
			if (eigensolver2.info() != Success) abort();
				COVeg2 = eigensolver2.eigenvalues().asDiagonal();
			for( i = 0; i < p + 1; i++ )
				COVeg2(i,i)=1/sqrt(eigensolver2.eigenvalues()(i));
      
			COVsq2 = eigensolver2.eigenvectors() * COVeg2 * eigensolver2.eigenvectors().transpose();
			for( i = 0; i < p + 1; i++ )
				tempgas2(i) = gasdev(&idum);
      
			B.row(itr) = COVsq2 * COVsq2 * tmu + sqrt(sig) * COVsq2 *tempgas2;
		} 
        else{  // The Faster way
        
			B.row(itr) = B.row(itr - 1);
			for( j = 0; j < p+1; j++)
				T1(j) = ( Z(itr - 1, j) * s1 * t(itr-1, j) + (1 - Z(itr - 1, j)) * s0 * t(itr-1, j)  );
        
			D_diag =  T1 * sig;
			Phi = X/sqrt(sig);
			alpha = Y/sqrt(sig);
        
        //sample u
			for( i = 0; i < p+1; i++ )
				tempgas2(i) = gasdev(&idum);
        
			u = tempgas2.cwiseProduct( D_diag.cwiseSqrt() );
        
			for( i = 0; i < n; i++ )
				delta(i) = gasdev(&idum);
        
			v = (Phi * u) + delta;
			D = D_diag.asDiagonal();
        
			for( i = 0; i < n; i++ )
				tempgas2n(i) = 1.00;
        
			diagN = tempgas2n.asDiagonal();
			tempA = Phi * D * Phi.transpose() + diagN;
			tempB = alpha - v;
        
			w = tempA.colPivHouseholderQr().solve( tempB );
			B.row(itr) = u + ( D * Phi.transpose()) * w;
			//cout << "run" << endl;
        }

    

  
        //updating Z
        Z(itr,0)=1;
		for(i=1;i<p+1;i++){
            double s1sq = sig * s1 * t(itr-1, i);
			double s0sq = sig * s0 * t(itr-1, i);
            prob(i)=pr* mydnorm(B(itr,i), (double)0, s1sq)/( pr* mydnorm(B(itr,i), (double)0, s1sq) + (1-pr)* mydnorm(B(itr,i), (double)0, s0sq ));
            Z(itr,i) = (ran1(&idum) < prob(i))? 1:0;
        }
        
        //updating Sigma
        
        for( j = 0; j < p + 1; j++ )
            T1(j) =  1/ ( Z(itr, j) * s1 * t(itr-1, j) + ( 1 - Z(itr, j) ) * s0 * t(itr-1, j) );
      
		vec1 = Y - X * B.row(itr).transpose();
		a = P1 + n * 0.5 + p * 0.5;
		b = P2 + 0.5 * vec1.dot(vec1) + 0.5 * B.row(itr) * T1.asDiagonal() * B.row(itr).transpose();
        
		sig  = 1/gamdev(a, b, &idum);
		sigma(itr) = sig;
	  
	  //updating Tau
	 
        
		for( j = 0; j < p + 1; j++ ){
			T1(j) =  ( Z(itr, j) * s1  + ( 1 - Z(itr, j) ) * s0  );
      
			a = nu * 0.5 +  0.5;
			b = nu * 0.5  + 0.5 * B(itr, j) * B(itr, j) / (2 * sig * T1(j));
          
			t(itr, j) = 1/gamdev(a, b, &idum);
		}
	  
		if( itr % 400 == 0)
			cout << itr << endl;
        
    }
    
	//cout << *Fast << endl;
		
		//Return output value
	for( i=0;i<(nburn+niter);i++){
		for( j=0;j<(p+1);j++){
			outZ[(nburn+niter)*j+i]=(double)Z(i,j);
			outbeta[(nburn+niter)*j+i]=B(i,j);
		}
	}
	for( i=0;i<(nburn+niter);i++)
		outSig[ i ] = (double)sigma(i);
		
      
	}
}

extern "C"{
    void basadTPr(double* X_, double* Y_, double* Z0, double* B0, double* sig_, double* pr_,
                int* n_, int* p_, double* nu_,double* s0_,double* s1_, int* nburn_, int* niter_,
                int* nsplit_, double* outZ,  double* outbeta, double* outSig, double *outPr, int *Fast){
        
        
        int n=*n_, p=*p_, nburn=*nburn_, niter=*niter_;
        double pr=*pr_, sig =*sig_, nu = *nu_, prTemp;
        
        MatrixXd X(n,p+1);
        VectorXd Y(n);
        VectorXd prV(niter + nburn);
        
        int i, j;
        
        for (i=0;i<n;i++){
            for (j=1;j<p+1;j++) {
                X(i,j)=X_[j*n+i];
            }
            Y(i)=Y_[i];
            X(i,0) = 1.0;
        }
        
        long idum = -time(0);
        
        int itr, nsplit =*nsplit_,vsize = (p+1)/nsplit;
        double s0 = *s0_, s1 = *s1_, a, b;
        
        
        
        MatrixXd B(niter+nburn, p+1), Z(niter+nburn, p+1), G(p+1, p+1), COV(vsize, vsize), COVsq(vsize, vsize), tempG(vsize, p+1-vsize), COVeg(vsize, vsize), COV2(p+1, p+1), COVeg2(p+1, p+1), COVsq2(p+1, p+1);
        
        VectorXd tmu(p+1), T1(p+1), tempgas2(p+1), tempgas2n(n), delta(n),  prob(p+1), sigma(niter + nburn), vec1(n), vec2(p+1), T1Inver(p+1),  tempDelta(p+1);
        
        
        
        //allocation for new faster algorithm
        
        MatrixXd Phi(n, p + 1), diagN(n, n), tempA(n,n), D(p+1, p+1);
        VectorXd D_diag(p+1), alpha(n), u(p+1), v(n), w(n), tempB(n);
        
        
        //allocation for tPrior
        
        MatrixXd t(niter + nburn, p+1);
        
        //MatrixXd
        
        //Initial Values
        
        for(j=1;j<p+1;j++){
            B(0,j)=B0[j-1];
            Z(0,j)=Z0[j-1];
        }
        B(0,0) = 0;
        Z(0,0) = 1;
        sigma(0) = sig;
        prV(0) = pr;
        for( itr = 1; itr < (niter + nburn) ; itr++ )
            sigma(itr) = 1.0;
        for( itr = 0; itr < (niter + nburn) ; itr++ )
            for( i = 0; i < p+1; i++ )
                t(itr, i) = 1.0;
        
        
        //Gibss
        cout << "prior probability that a coefficient is nonzero is estimated by Gibbs sampling" << endl;
        
        G = X.transpose() * X;
        tmu = X.transpose() * Y;
        
        
        
        for( itr = 1; itr <(nburn + niter); itr++ ){
            
            sig = sigma(itr - 1);
            prTemp = prV(itr - 1);
            
            //updating B
            
            //original way
            if( *Fast == 0 ){
                
                B.row(itr) = B.row(itr - 1);
                for( j = 0; j < p+1; j++)
                    T1(j) = 1/( Z(itr-1, j) * s1 * t(itr-1, j) + ( 1 - Z(itr - 1, j) ) * s0 * t(itr-1, j));
                
                COV2 = G;
                
                for(i=0;i<p+1;i++)
                    COV2(i,i) += T1(i);
                
                SelfAdjointEigenSolver<MatrixXd> eigensolver2(COV2);
                
                if (eigensolver2.info() != Success) abort();
                COVeg2 = eigensolver2.eigenvalues().asDiagonal();
                for( i = 0; i < p + 1; i++ )
                    COVeg2(i,i)=1/sqrt(eigensolver2.eigenvalues()(i));
                
                COVsq2 = eigensolver2.eigenvectors() * COVeg2 * eigensolver2.eigenvectors().transpose();
                for( i = 0; i < p + 1; i++ )
                    tempgas2(i) = gasdev(&idum);
                
                B.row(itr) = COVsq2 * COVsq2 * tmu + sqrt(sig) * COVsq2 *tempgas2;
            }
            else{  // The Faster way
                
                B.row(itr) = B.row(itr - 1);
                for( j = 0; j < p+1; j++)
                    T1(j) = ( Z(itr - 1, j) * s1 * t(itr-1, j) + (1 - Z(itr - 1, j)) * s0 * t(itr-1, j)  );
                
                D_diag =  T1 * sig;
                Phi = X/sqrt(sig);
                alpha = Y/sqrt(sig);
                
                //sample u
                for( i = 0; i < p+1; i++ )
                    tempgas2(i) = gasdev(&idum);
                
                u = tempgas2.cwiseProduct( D_diag.cwiseSqrt() );
                
                for( i = 0; i < n; i++ )
                    delta(i) = gasdev(&idum);
                
                v = (Phi * u) + delta;
                D = D_diag.asDiagonal();
                
                for( i = 0; i < n; i++ )
                    tempgas2n(i) = 1.00;
                
                diagN = tempgas2n.asDiagonal();
                tempA = Phi * D * Phi.transpose() + diagN;
                tempB = alpha - v;
                
                w = tempA.colPivHouseholderQr().solve( tempB );
                B.row(itr) = u + ( D * Phi.transpose()) * w;

            }
            
            
            
            
            //updating Z
            Z(itr,0)=1;
            for(i=1;i<p+1;i++){
                double s1sq = sig * s1 * t(itr-1, i);
                double s0sq = sig * s0 * t(itr-1, i);
                prob(i)=prTemp* mydnorm(B(itr,i), (double)0, s1sq)/( prTemp* mydnorm(B(itr,i), (double)0, s1sq) + (1-prTemp)* mydnorm(B(itr,i), (double)0, s0sq ));
                Z(itr,i) = (ran1(&idum) < prob(i))? 1:0;
            }
            
            //updating Sigma
            
            for( j = 0; j < p + 1; j++ )
                T1(j) =  1/ ( Z(itr, j) * s1 * t(itr-1, j) + ( 1 - Z(itr, j) ) * s0 * t(itr-1, j) );
            
            vec1 = Y - X * B.row(itr).transpose();
            a = P1 + n * 0.5 + p * 0.5;
            b = P2 + 0.5 * vec1.dot(vec1) + 0.5 * B.row(itr) * T1.asDiagonal() * B.row(itr).transpose();
            
            sig  = 1/gamdev(a, b, &idum);
            sigma(itr) = sig;
            
            //updating Tau
            
            
            for( j = 0; j < p + 1; j++ ){
                T1(j) =  ( Z(itr, j) * s1  + ( 1 - Z(itr, j) ) * s0  );
                
                a = nu * 0.5 +  0.5;
                b = nu * 0.5  + 0.5 * B(itr, j) * B(itr, j) / (2 * sig * T1(j));
                
                t(itr, j) = 1/gamdev(a, b, &idum);
            }
            
   
            
            //updating Pr
            double a1 = beta1 + Z.row(itr).sum();
            double b1 = beta2 + p + 1 - Z.row(itr).sum();
            prV(itr) = betadev(a1, b1, &idum);
            
            if( itr % 400 == 0)
                cout << itr << endl;
            
        }
        
        //cout << *Fast << endl;
        
        //Return output value
        for( i=0;i<(nburn+niter);i++){
            for( j=0;j<(p+1);j++){
                outZ[(nburn+niter)*j+i]=(double)Z(i,j);
                outbeta[(nburn+niter)*j+i]=B(i,j);
            }
        }
        for( i=0;i<(nburn+niter);i++){
            outSig[ i ] = (double)sigma(i);
            outPr[ i ] = (double)prV(i);
        }
    }
}




RcppExport SEXP basadFunctionT(SEXP X, SEXP Y, SEXP Z0, SEXP B0, SEXP sig, SEXP pr, SEXP n, SEXP p, SEXP nu, SEXP s0, SEXP s1, SEXP nburn, SEXP niter, SEXP nsplit, SEXP Fast  ){
  
	Rcpp::NumericMatrix XX(X);
	Rcpp::NumericVector YY(Y);
	Rcpp::NumericVector ZZ(Z0);
	Rcpp::NumericVector BB(B0);
	double sighat = Rcpp::as<double>(sig);
	double ppr = Rcpp::as<double>(pr);
	int nn = Rcpp::as<int>(n);
	int pp = Rcpp::as<int>(p);
	double nunu = Rcpp::as<double>(nu);
	double ss0 = Rcpp::as<double>(s0);
	double ss1 = Rcpp::as<double>(s1);
	int nnburn = Rcpp::as<int>(nburn);
	int nniter = Rcpp::as<int>(niter);
	int nnsplit = Rcpp::as<int>(nsplit);
	int fast = Rcpp::as<int>(Fast);
	int i,j;

  
	double *XXX = new double[nn * (pp+1)];
	for(i =0; i < nn; i++ ){
		for(j = 0; j < (pp + 1); j++ ){
		XXX[j * nn + i  ] = XX(i, j);
		}
	} 
    
	double *YYY = new double[nn];
	for( i = 0; i < nn; i++ )
		YYY[i] = YY(i);
  
  
	double *ZZZ = new double[pp+1];
	for( i = 0; i < pp + 1; i++ )
		ZZZ[i] = ZZ(i);
    
	double *BBB = new double[pp+1];
	for( i = 0; i < pp + 1; i++ )
		BBB[i] = BB(i);

	double *outZZ = new double[  (nnburn + nniter) * (pp + 1) ];
	double *outBB = new double[  (nnburn + nniter) * (pp + 1) ];
	double *outSig = new double[ (nnburn + nniter) ];
    double *outPr = new double[ (nnburn + nniter) ];
  
   if( ppr >0 )
	basadT( XXX, YYY, ZZZ, BBB, &sighat, &ppr, &nn, &pp, &nunu, &ss0, &ss1, &nnburn, &nniter, &nnsplit, outZZ, outBB, outSig, &fast);
   else{
       ppr = 0.5;
       	basadTPr( XXX, YYY, ZZZ, BBB, &sighat, &ppr, &nn, &pp, &nunu, &ss0, &ss1, &nnburn, &nniter, &nnsplit, outZZ, outBB, outSig, outPr, &fast);
   }
  
  
	int nnrow;
	nnrow = nnburn + nniter;
	int nncol;
	nncol = pp + 1;
    
   // cout << nnrow << nncol << endl;
  
    Rcpp::NumericMatrix realOutB( nnrow , nncol );
    Rcpp::NumericMatrix realOutZ( nnrow , nncol );

  
  
	for( int k = 0; k < nnrow; k++ ){
		for( int h = 0; h < nncol; h++ ){
			realOutB(k, h) = outBB[ h *(nnburn + nniter) + k];
			realOutZ(k, h) = outZZ[ h *(nnburn + nniter) + k];
		}
	}
   
    Rcpp::NumericVector realOutSig( nnburn + nniter );
    Rcpp::NumericVector realOutPr(  nnburn + nniter );
   
    for( int k = 0; k < nnburn + nniter; k++ ){
        realOutSig(k) = outSig[k];
        realOutPr(k) = outPr[k];
    }

  
    delete [] XXX;
    delete [] YYY;
    delete [] ZZZ;
    delete [] BBB;
    delete [] outZZ;
    delete [] outBB;
    delete [] outSig;
  
	return List::create(
		Named("B") = realOutB,
		Named("Z") = realOutZ
		//Named("Sig") = realOutSig
	);
   
  
  //return wrap(10);
}

