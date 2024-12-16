#include <cstdio>
#include <cstddef>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <iterator>

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>

//R CMD SHLIB -lRlapack fn-surv7-C-mu-k1.cpp
extern "C"
{
	double fn_lgamma(double x)
	{
		double x0,x2,xp,gl,gl0;
		int n,k;
		static double a[] = {
			8.333333333333333e-02,
			-2.777777777777778e-03,
			7.936507936507937e-04,
			-5.952380952380952e-04,
			8.417508417508418e-04,
			-1.917526917526918e-03,
			6.410256410256410e-03,
			-2.955065359477124e-02,
			1.796443723688307e-01,
			-1.39243221690590};
		
		x0 = x;
		if (x <= 0.0) return 1e308;
		else if ((x == 1.0) || (x == 2.0)) return 0.0;
		else if (x <= 7.0) {
			n = (int)(7-x);
			x0 = x+n;
		}
		
		x2 = 1.0/(x0*x0);
		xp = 2.0*M_PI;
		gl0 = a[9];
		
		for (k=8;k>=0;k--) {
			gl0 = gl0*x2 + a[k];
		}
		
		gl = gl0/x0+0.5*log(xp)+(x0-0.5)*log(x0)-x0;
		if (x <= 7.0) {
			for (k=1;k<=n;k++) {
				gl -= log(x0-1.0);
				x0 -= 1.0;
			}
		}
		return gl;
	}
	
 
    void fn_compute_BeZ(double *ZZ, double *Beta_j, int Tot_N, int PP, double *BeZ_j)
    {
        int p, i_sam;
        double BeZ;
        
        // \sum_p beta_jp*Z_ip ==> Tot_N-dim vector
        for (i_sam=0; i_sam < Tot_N; i_sam++) {
            
            BeZ = 0.0;
            //#### Z: Tot_N*p covariates
            for (p=0; p < PP; p++) {
                BeZ = BeZ + ZZ[p*Tot_N + i_sam]*Beta_j[p];
            }
            
            BeZ_j[i_sam] = BeZ;
        }
    }
    
    
/*
################################################################################################################################
## Beta : beta_{pj}: effect of covariate p on OTU j
## ### \Beta: p*J matrix.
### MMu:  Tot_N*J matrix: J otus and Tot_N samples ** Note ** some rows are repitition due to replicate.
################################################################################################################################

 
 ### Y: Tot_N*J matrix: J otus and Tot_N samples

 
 #### K: Tot_N * M where  M: the number of basis points
 #### Th0: J-dim vector
 #### Th: M*J matrix
 ####  ==> KTh: Tot_N*J matrix
 
 ### r: Tot_N
 #### s: J-dim vector
 
 #### Z: Tot_N*p covariates   *** Note ***: Z is the same for all K_t replicates at time t
 #### BeZ: (p*J)*(Tot_N*p)  ==> Tot_N*J
 #### Change Beta  ===> Change BeZ =====> Chage Mu
 
 ####  \beta_jp \sim N(0, \sig^2_j*phi_jp)
 ####  \phi: p*J matrix.
 */
    
    void fn_update_Beta(int JJ, int PP, int Tot_N, double *ZZ, double *Th0, double *KTh, double *YY, double *SS, double *RR, double *sig2, double *phi, double *Beta_cur, double *BeZ_cur, double *MMu_cur)
    {
        int j, p, ind1, ind2, i_sam;
        double A_cur, A_pro, u;
        double *mu_cur, *mu_pro, *m_pro, *m_cur, *Beta_j_cur, *Beta_j_pro, *BeZ_j_cur, *BeZ_j_pro;
        
        //Allocate the memory
        mu_cur = new double[Tot_N]; //Tot_N vector
        mu_pro = new double[Tot_N]; //Tot_N vector
        m_cur = new double[Tot_N]; //Tot_N vector
        m_pro = new double[Tot_N]; //Tot_N vector
        Beta_j_cur = new double[PP]; //p-dim vector
        Beta_j_pro = new double[PP]; //p-dim vector
        
        BeZ_j_cur = new double[Tot_N]; //Tot_N-dim vector
        BeZ_j_pro = new double[Tot_N]; //Tot_N-dim vector
        
        //cnt = 0;
        //## one OTU at a time
        for(j=0; j < JJ; j++)
        {
            ind1 = PP*j; // index for p*J dim matrix -- Be,
            ind2 = Tot_N*j;  // index for Tot_N*J dim matrix -- Y, Kth, BeZ, MMu
            
            // take all p betas for otu j
            for (p=0; p < PP; p++) {
                // Beta: (p*J) matrix
                Beta_j_cur[p] = Beta_cur[ind1 + p];     //# p-dim vector --- cov for OTU j
                Beta_j_pro[p] = Beta_j_cur[p];
            } //for (p=0; p < PP; p++) {
            
            // MMu:  Tot_N*J matrix:
            // take all Tot_N mu for otu j
            for (i_sam=0; i_sam < Tot_N; i_sam++) {
                BeZ_j_cur[i_sam] = BeZ_cur[ind2 + i_sam];            //# Tot_N-dim
                mu_cur[i_sam] = MMu_cur[ind2 + i_sam];               //# Tot_N-dim
                
                // compute the mean
                m_cur[i_sam] = RR[i_sam]*SS[j]*mu_cur[i_sam];        //# Tot_N-dim
                
            } //for (i_sam=0; i_sam < Tot_N; i_sam++) {
            
            //### one covariate at a time
            for(p=0; p < PP; p++)
            {
                //propose new beta_{pj}
                GetRNGstate();
                Beta_j_pro[p] = Beta_j_pro[p] + rnorm(0.0, 0.3);
                PutRNGstate();
                
                //#### Z: Tot_N*p covariates
                fn_compute_BeZ(ZZ, Beta_j_pro, Tot_N, PP, BeZ_j_pro);  // \sum_p beta_jp*Z_ip ==> Tot_N-dim vector
                
                // Prior
                A_cur =  - pow(Beta_j_cur[p], 2.0)/2.0/sig2[j]/phi[ind1 + p];
                A_pro =  - pow(Beta_j_pro[p], 2.0)/2.0/sig2[j]/phi[ind1 + p];
                
                for (i_sam=0; i_sam < Tot_N; i_sam++) {
                    mu_pro[i_sam] = exp(Th0[j] + KTh[ind2 + i_sam] + BeZ_j_pro[i_sam]);             //# Tot_N-dim
                    m_pro[i_sam] = RR[i_sam]*SS[j]*mu_pro[i_sam];                                   //# Tot_N-dim
                    
                    //### Y: Tot_N*J matrix: J otus and Tot_N samples
                    A_cur = A_cur + YY[ind2 + i_sam]*log(m_cur[i_sam]/(1.0+m_cur[i_sam])) - 1.0/SS[j]*log(1.0 + m_cur[i_sam]);
                    A_pro = A_pro + YY[ind2 + i_sam]*log(m_pro[i_sam]/(1.0+m_pro[i_sam])) - 1.0/SS[j]*log(1.0 + m_pro[i_sam]);
                } // for (i_sam=0; i_sam < Tot_N; i_sam++) {
                
                //Let's do M-H
                GetRNGstate();
                u = runif(0.0, 1.0);
                PutRNGstate();
                
                //### Accept?
                if(log(u) < (A_pro - A_cur))
                {
                    //cnt = cnt + 1;
                    // accept
                    Beta_j_cur[p] = Beta_j_pro[p];
                    
                    for (i_sam=0; i_sam < Tot_N; i_sam++) {
                        BeZ_j_cur[i_sam] = BeZ_j_pro[i_sam];  //# Tot_N-dim
                        mu_cur[i_sam] = mu_pro[i_sam];        //# Tot_N-dim
                        m_cur[i_sam] = m_pro[i_sam];          //# Tot_N-dim
                    }
                }else{ //if(log(u) < (prob_pro - prob_cur))
                    Beta_j_pro[p] = Beta_j_cur[p];
                }//if(log(u) < (prob_pro - prob_cur))
                
            } //## for(p in 1:PP)

            //### save the current Beta ==> so BeZ and Mu
            // take all p betas for otu j
            for (p=0; p < PP; p++) {
                // ind1 = PP*j;
                // Beta: (p*J) matrix
                Beta_cur[ind1 + p] = Beta_j_cur[p];     //# p-dim vector --- cov for OTU j
            } //for (p=0; p < PP; p++) {
            
            // MMu:  Tot_N*J matrix:
            // take all Tot_N mu for otu j
            for (i_sam=0; i_sam < Tot_N; i_sam++) {
                //ind2 = Tot_N*j;
                BeZ_cur[ind2 + i_sam] = BeZ_j_cur[i_sam];            //# Tot_N-dim
                MMu_cur[ind2 + i_sam] = mu_cur[i_sam];               //# Tot_N-dim
            } //for (i_sam=0; i_sam < Tot_N; i_sam++) {
            
        } //## for(j in 1:JJ)
        
        //printf("\n accpt single beta cnt=%f", 1.0*cnt/(1.0*JJ*PP));

        
        // Free the memory
        delete[] mu_pro; mu_pro = NULL;
        delete[] m_pro; m_pro = NULL;
        delete[] Beta_j_pro; Beta_j_pro = NULL;
        delete[] BeZ_j_pro; BeZ_j_pro = NULL;
        delete[] mu_cur; mu_cur = NULL;
        delete[] m_cur; m_cur = NULL;
        delete[] Beta_j_cur; Beta_j_cur = NULL;
        delete[] BeZ_j_cur; BeZ_j_cur = NULL;

    }

    
    /*
     ################################################################################################################################
     ## Beta : beta_{pj}: effect of covariate p on OTU j
     ## ### \Beta: p*J matrix.
     ### MMu:  Tot_N*J matrix: J otus and Tot_N samples ** Note ** some rows are repitition due to replicate.
     ################################################################################################################################
     
     
     ### Y: Tot_N*J matrix: J otus and Tot_N samples
     
     
     #### K: Tot_N * M where  M: the number of basis points
     #### Th0: J-dim vector
     #### Th: M*J matrix
     ####  ==> KTh: Tot_N*J matrix
     
     ### r: Tot_N
     #### s: J-dim vector
     
     #### Z: Tot_N*p covariates   *** Note ***: Z is the same for all K_t replicates at time t
     #### BeZ: (p*J)*(Tot_N*p)  ==> Tot_N*J
     #### Change Beta  ===> Change BeZ =====> Chage Mu
     
     ####  \beta_jp \sim N(0, \sig^2_j*phi_jp)
     ####  \phi: p*J matrix.
     */
    
    void fn_update_ALL_Beta(int JJ, int PP, int Tot_N, double *ZZ, double *Th0, double *KTh, double *YY, double *SS, double *RR, double *sig2, double *phi, double *Beta_cur, double *BeZ_cur, double *MMu_cur)
    {
        int j, p, ind1, ind2, i_sam;
        double A_cur, A_pro, u, th0_j, s_j;
        double *mu_pro, *m_pro, *m_cur, *Beta_j_pro, *BeZ_j_pro;
        
        //Allocate the memory
        mu_pro = new double[Tot_N]; //Tot_N vector
        
        m_cur = new double[Tot_N]; //Tot_N vector
        m_pro = new double[Tot_N]; //Tot_N vector
        
        Beta_j_pro = new double[PP]; //p-dim vector
        BeZ_j_pro = new double[Tot_N]; //Tot_N-dim vector
        
        //cnt = 0;
        //## one OTU at a time
        for(j=0; j < JJ; j++)
        {
            ind1 = PP*j; // index for p*J dim matrix -- Be,
            ind2 = Tot_N*j;  // index for Tot_N*J dim matrix -- Y, Kth, BeZ, MMu
            
            // th0 and s for otu j
            th0_j = Th0[j];
            s_j = SS[j];
            
            // Prior
            A_cur = 0.0;
            A_pro = 0.0;
            
            // take all p betas for otu j & propose new values of beta
            for (p=0; p < PP; p++) {
                // Beta: (p*J) matrix
                //propose new beta_{pj}
                GetRNGstate();
                Beta_j_pro[p] = Beta_cur[ind1 + p] + rnorm(0.0, 0.01);
                PutRNGstate();
                
                // Prior
                A_cur =  A_cur - pow(Beta_cur[ind1 + p], 2.0)/2.0/sig2[j]/phi[ind1 + p];
                A_pro =  A_pro - pow(Beta_j_pro[p], 2.0)/2.0/sig2[j]/phi[ind1 + p];

            } //for (p=0; p < PP; p++) {
            
            //#### Z: Tot_N*p covariates
            fn_compute_BeZ(ZZ, Beta_j_pro, Tot_N, PP, BeZ_j_pro);  // \sum_p beta_jp*Z_ip ==> Tot_N-dim vector
            
            // MMu:  Tot_N*J matrix:
            // take all Tot_N mu for otu j
            for (i_sam=0; i_sam < Tot_N; i_sam++) {
                
                // with current values
                // compute the mean
                m_cur[i_sam] = RR[i_sam]*s_j*MMu_cur[ind2 + i_sam];        //# Tot_N-dim
                
                // with proposed values
                mu_pro[i_sam] = exp(th0_j + KTh[ind2 + i_sam] + BeZ_j_pro[i_sam]);            //# Tot_N-dim
                m_pro[i_sam] = RR[i_sam]*s_j*mu_pro[i_sam];                                   //# Tot_N-dim
                
                //### Y: Tot_N*J matrix: J otus and Tot_N samples
                A_cur = A_cur + YY[ind2 + i_sam]*log(m_cur[i_sam]/(1.0+m_cur[i_sam])) - 1.0/s_j*log(1.0 + m_cur[i_sam]);
                A_pro = A_pro + YY[ind2 + i_sam]*log(m_pro[i_sam]/(1.0+m_pro[i_sam])) - 1.0/s_j*log(1.0 + m_pro[i_sam]);

            } //for (i_sam=0; i_sam < Tot_N; i_sam++) {
            
            
            //Let's do M-H
            GetRNGstate();
            u = runif(0.0, 1.0);
            PutRNGstate();
            
            //### Accept?
            if(log(u) < (A_pro - A_cur))
            {
                //cnt = cnt + 1;
                // accept******************************************
                //### save the current Beta ==> so BeZ and Mu
                // take all p betas for otu j
                for (p=0; p < PP; p++) {
                    // ind1 = PP*j;
                    // Beta: (p*J) matrix
                    Beta_cur[ind1 + p] = Beta_j_pro[p];     //# p-dim vector --- cov for OTU j
                } //for (p=0; p < PP; p++) {
                
                // MMu:  Tot_N*J matrix:
                // take all Tot_N mu for otu j
                for (i_sam=0; i_sam < Tot_N; i_sam++) {
                    //ind2 = Tot_N*j;
                    BeZ_cur[ind2 + i_sam] = BeZ_j_pro[i_sam];            //# Tot_N-dim
                    MMu_cur[ind2 + i_sam] = mu_pro[i_sam];               //# Tot_N-dim
                } //for (i_sam=0; i_sam < Tot_N; i_sam++) {

            }//if(log(u) < (prob_pro - prob_cur))
            
        } //## for(j in 1:JJ)
        
        //printf("\n accpt joint beta cnt=%f", 1.0*cnt/(1.0*JJ));

        // Free the memory
        delete[] mu_pro; mu_pro = NULL;
        delete[] m_pro; m_pro = NULL;
        delete[] Beta_j_pro; Beta_j_pro = NULL;
        delete[] BeZ_j_pro; BeZ_j_pro = NULL;
        delete[] m_cur; m_cur = NULL;
        
    }


    
    /*
     ################################################################################################################################
     ## Beta : beta_{pj}: effect of covariate p on OTU j
     ## ### \Beta: p*J matrix.
     ### MMu:  Tot_N*J matrix: J otus and Tot_N samples ** Note ** some rows are repitition due to replicate.
     ################################################################################################################################
     
     
     ### Y: Tot_N*J matrix: J otus and Tot_N samples
     
     
     #### K: Tot_N * M where  M: the number of basis points
     #### Th0: J-dim vector
     #### Th: M*J matrix
     ####  ==> KTh: Tot_N*J matrix
     
     ### r: Tot_N
     #### s: J-dim vector
     
     #### Z: Tot_N*p covariates   *** Note ***: Z is the same for all K_t replicates at time t
     #### BeZ: (p*J)*(Tot_N*p)  ==> Tot_N*J
     #### Change Beta  ===> Change BeZ =====> Chage Mu
     
     ####  \beta_jp \sim N(0, \sig^2_j*phi_jp)
     ####  \phi: p*J matrix.
     */

    
    void fn_update_ALL_Beta_Th0(int JJ, int PP, int Tot_N, double *ZZ, double *Th0_bar, double u2_th0, double *Th0_cur, double *KTh, double *YY, double *SS, double *RR, double *sig2, double *phi, double *Beta_cur, double *BeZ_cur, double *MMu_cur)
    {
        int j, p, ind1, ind2, i_sam;
        double A_cur, A_pro, u, Th0_j_pro, s_j;
        double *mu_pro, *m_pro, *m_cur, *Beta_j_pro, *BeZ_j_pro;
        
        //Allocate the memory
        mu_pro = new double[Tot_N]; //Tot_N vector
        
        m_cur = new double[Tot_N]; //Tot_N vector
        m_pro = new double[Tot_N]; //Tot_N vector
        
        Beta_j_pro = new double[PP]; //p-dim vector
        BeZ_j_pro = new double[Tot_N]; //Tot_N-dim vector
        
        //cnt = 0;
        //## one OTU at a time
        for(j=0; j < JJ; j++)
        {
            ind1 = PP*j; // index for p*J dim matrix -- Be, Phi
            ind2 = Tot_N*j;  // index for Tot_N*J dim matrix -- Y, Kth, BeZ, MMu
            
            // th0 and s for otu j
            s_j = SS[j];
            
            // Prior
            A_cur = 0.0;
            A_pro = 0.0;
            
            // take all p betas for otu j & propose new values of beta
            for (p=0; p < PP; p++) {
                // Beta: (p*J) matrix
                //propose new beta_{pj}
                GetRNGstate();
                Beta_j_pro[p] = Beta_cur[ind1 + p] + rnorm(0.0, 0.01);
                PutRNGstate();
                
                // Prior
                A_cur =  A_cur - pow(Beta_cur[ind1 + p], 2.0)/2.0/sig2[j]/phi[ind1 + p];
                A_pro =  A_pro - pow(Beta_j_pro[p], 2.0)/2.0/sig2[j]/phi[ind1 + p];
                
            } //for (p=0; p < PP; p++) {
            
            //#### Z: Tot_N*p covariates
            fn_compute_BeZ(ZZ, Beta_j_pro, Tot_N, PP, BeZ_j_pro);  // \sum_p beta_jp*Z_ip ==> Tot_N-dim vector
            
            
            //propose new th_{0j}
            GetRNGstate();
            Th0_j_pro = Th0_cur[j] + rnorm(0.0, 0.01);
            PutRNGstate();
            
            // prior
            A_cur = A_cur - pow(Th0_cur[j] - Th0_bar[j], 2.0)/2.0/u2_th0;
            A_pro = A_pro - pow(Th0_j_pro - Th0_bar[j], 2.0)/2.0/u2_th0;
            

            // MMu:  Tot_N*J matrix:
            // take all Tot_N mu for otu j
            for (i_sam=0; i_sam < Tot_N; i_sam++) {
                
                // with current values
                // compute the mean
                m_cur[i_sam] = RR[i_sam]*s_j*MMu_cur[ind2 + i_sam];        //# Tot_N-dim
                
                // with proposed values
                mu_pro[i_sam] = exp(Th0_j_pro + KTh[ind2 + i_sam] + BeZ_j_pro[i_sam]);        //# Tot_N-dim
                m_pro[i_sam] = RR[i_sam]*s_j*mu_pro[i_sam];                                   //# Tot_N-dim
                
                //### Y: Tot_N*J matrix: J otus and Tot_N samples
                A_cur = A_cur + YY[ind2 + i_sam]*log(m_cur[i_sam]/(1.0+m_cur[i_sam])) - 1.0/s_j*log(1.0 + m_cur[i_sam]);
                A_pro = A_pro + YY[ind2 + i_sam]*log(m_pro[i_sam]/(1.0+m_pro[i_sam])) - 1.0/s_j*log(1.0 + m_pro[i_sam]);
                
            } //for (i_sam=0; i_sam < Tot_N; i_sam++) {
            
            
            //Let's do M-H
            GetRNGstate();
            u = runif(0.0, 1.0);
            PutRNGstate();
            
            //### Accept?
            if(log(u) < (A_pro - A_cur))
            {
                //cnt = cnt + 1;
                // accept******************************************
                
                // accept
                Th0_cur[j] = Th0_j_pro;
                
                //### save the current Beta ==> so BeZ and Mu
                // take all p betas for otu j
                for (p=0; p < PP; p++) {
                    // ind1 = PP*j;
                    // Beta: (p*J) matrix
                    Beta_cur[ind1 + p] = Beta_j_pro[p];     //# p-dim vector --- cov for OTU j
                } //for (p=0; p < PP; p++) {
                
                // MMu:  Tot_N*J matrix:
                // take all Tot_N mu for otu j
                for (i_sam=0; i_sam < Tot_N; i_sam++) {
                    //ind2 = Tot_N*j;
                    BeZ_cur[ind2 + i_sam] = BeZ_j_pro[i_sam];            //# Tot_N-dim
                    MMu_cur[ind2 + i_sam] = mu_pro[i_sam];               //# Tot_N-dim
                } //for (i_sam=0; i_sam < Tot_N; i_sam++) {
                
            }//if(log(u) < (prob_pro - prob_cur))
            
        } //## for(j in 1:JJ)
        
        //printf("\n accpt joint th0, beta cnt=%f", 1.0*cnt/(1.0*JJ));
        
        // Free the memory
        delete[] mu_pro; mu_pro = NULL;
        delete[] m_pro; m_pro = NULL;
        delete[] Beta_j_pro; Beta_j_pro = NULL;
        delete[] BeZ_j_pro; BeZ_j_pro = NULL;
        delete[] m_cur; m_cur = NULL;
        
    }
    

    
    /*
     ################################################################################################################################
     ## Th_pj : baseline of  OTU j
     ################################################################################################################################
     
     #### K: Tot_N * M where  M: the number of basis points
     #### Th0: J-dim vector
     #### Th: M*J matrix
     ####  ==> KTh: Tot_N*J matrix
     
     #### Z: Tot_N*p covariates   *** Note ***: Z is the same for all K_t replicates at time t
     #### BeZ: (p*J)*(Tot_N*p)  ==> Tot_N*J
     
     ### MMu:  Tot_N*J matrix: J otus and Tot_N samples ** Note ** some rows are repitition due to replicate.
     ### r: Tot_N
     
     #### Change Th  ===> Change KTh =====> Chage Mu
     */
    void fn_update_ThP(int JJ, int MM, int Tot_N, double *KK, double *YY, double *SS, double *RR, double *tau2, double *BeZ, double *Th0, double *Th_cur, double *KTh_cur, double *MMu_cur)
    {
        
        int j, m, ind1, ind2, i_sam;
        double A_cur, A_pro, u;
        double *mu_cur, *mu_pro, *m_pro, *m_cur, *KTh_j_cur, *KTh_j_pro, *Th_j_cur, *Th_j_pro;
        
        //Allocate the memory
        mu_cur = new double[Tot_N]; //Tot_N vector
        mu_pro = new double[Tot_N]; //Tot_N vector
        m_cur = new double[Tot_N]; //Tot_N vector
        m_pro = new double[Tot_N]; //Tot_N vector
        Th_j_cur = new double[MM]; //M-dim vector
        Th_j_pro = new double[MM]; //M-dim vector
        
        KTh_j_cur = new double[Tot_N]; //Tot_N-dim vector
        KTh_j_pro = new double[Tot_N]; //Tot_N-dim vector

        //## One OTU at a time
        for(j=0; j < JJ; j++)
        {
            ind1 = MM*j; // index for M*J dim matrix -- Th,
            ind2 = Tot_N*j;  // index for Tot_N*J dim matrix -- Y, Kth, BeZ, MMu
            
            // take all m thetas for otu j
            for (m=0; m < MM; m++) {
                // Theta: (M*J) matrix
                Th_j_cur[m] = Th_cur[ind1 + m];     //# p-dim vector --- cov for OTU j
                Th_j_pro[m] = Th_j_cur[m];
            } //for (p=0; p < PP; p++) {
            
            
            // MMu:  Tot_N*J matrix:
            // take all Tot_N mu for otu j
            for (i_sam=0; i_sam < Tot_N; i_sam++) {
                KTh_j_cur[i_sam] = KTh_cur[ind2 + i_sam];            //# Tot_N-dim
                mu_cur[i_sam] = MMu_cur[ind2 + i_sam];               //# Tot_N-dim
                
                // compute the mean
                m_cur[i_sam] = RR[i_sam]*SS[j]*mu_cur[i_sam];        //# Tot_N-dim
            } //for (i_sam=0; i_sam < Tot_N; i_sam++) {
            
            
            //##### One basis at a time
            for(m=0; m < MM; m++)
            {
                //propose new beta_{pj}
                GetRNGstate();
                Th_j_pro[m] = Th_j_cur[m] + rnorm(0.0, 0.8);
                PutRNGstate();
                
                // printf("\n j=%d, m=%d, Th_cur=%f, Th_pro=%f", j, m, Th_j_cur[m], Th_j_pro[m]);
                
                //#### KK: Tot_N * M covariates
                fn_compute_BeZ(KK, Th_j_pro, Tot_N, MM, KTh_j_pro);  // \sum_m theta_mj*K_im ==> Tot_N-dim vector
                
                // Prior
                A_cur =  - pow(Th_j_cur[m], 2.0)/2.0/tau2[j];
                A_pro =  - pow(Th_j_pro[m], 2.0)/2.0/tau2[j];
                
                for (i_sam=0; i_sam < Tot_N; i_sam++) {
                    mu_pro[i_sam] = exp(Th0[j] + KTh_j_pro[i_sam] + BeZ[ind2 + i_sam]);               //# Tot_N-dim
                    m_pro[i_sam] = RR[i_sam]*SS[j]*mu_pro[i_sam];                                     //# Tot_N-dim
                    
                    //### Y: Tot_N*J matrix: J otus and Tot_N samples
                    A_cur = A_cur + YY[ind2 + i_sam]*log(m_cur[i_sam]/(1.0+m_cur[i_sam])) - 1.0/SS[j]*log(1.0 + m_cur[i_sam]);
                    A_pro = A_pro + YY[ind2 + i_sam]*log(m_pro[i_sam]/(1.0+m_pro[i_sam])) - 1.0/SS[j]*log(1.0 + m_pro[i_sam]);
                } // for (i_sam=0; i_sam < Tot_N; i_sam++) {
                
                //Let's do M-H
                GetRNGstate();
                u = runif(0.0,1.0);
                PutRNGstate();
                
                //### Accept?
                if(log(u) < (A_pro - A_cur))
                {
                    // accept
                    Th_j_cur[m] = Th_j_pro[m];
                    
                    for (i_sam=0; i_sam < Tot_N; i_sam++) {
                        KTh_j_cur[i_sam] = KTh_j_pro[i_sam];  //# Tot_N-dim
                        mu_cur[i_sam] = mu_pro[i_sam];        //# Tot_N-dim
                        m_cur[i_sam] = m_pro[i_sam];          //# Tot_N-dim
                    }
                }else{ //if(log(u) < (prob_pro - prob_cur))
                    
                    Th_j_pro[m] = Th_j_cur[m];
                }//if(log(u) < (prob_pro - prob_cur))
            } //## for(m=0; m < MM; m++)
            
            //### save the current Th_j ==> so KTh and Mu
            // take all m theta for otu j
            for (m=0; m < MM; m++) {
                // ind1 = PP*j;
                // Beta: (p*J) matrix
                Th_cur[ind1 + m] = Th_j_cur[m];     //# M-dim vector --- cov for OTU j
            } //for (m=0; m < MM; m++) {
            
            // MMu:  Tot_N*J matrix:
            // take all Tot_N mu for otu j
            for (i_sam=0; i_sam < Tot_N; i_sam++) {
                //ind2 = Tot_N*j;
                KTh_cur[ind2 + i_sam] = KTh_j_cur[i_sam];      //# Tot_N-dim
                MMu_cur[ind2 + i_sam] = mu_cur[i_sam];               //# Tot_N-dim
            } //for (i_sam=0; i_sam < Tot_N; i_sam++) {
            
            
        } //for(j=0; j < JJ; j++)
        
        // Free the memory
        delete[] mu_pro; mu_pro = NULL;
        delete[] m_pro; m_pro = NULL;
        delete[] Th_j_pro; Th_j_pro = NULL;
        delete[] KTh_j_pro; KTh_j_pro = NULL;
        delete[] mu_cur; mu_cur = NULL;
        delete[] m_cur; m_cur = NULL;
        delete[] Th_j_cur; Th_j_cur = NULL;
        delete[] KTh_j_cur; KTh_j_cur = NULL;
    }
    
    
/*
################################################################################################################################
## Th_0j : baseline of  OTU j
################################################################################################################################
    
#### K: Tot_N * M where  M: the number of basis points
#### Th0: J-dim vector
#### Th: M*J matrix
####  ==> KTh: Tot_N*J matrix
    
#### Z: Tot_N*p covariates   *** Note ***: Z is the same for all K_t replicates at time t
#### BeZ: (p*J)*(Tot_N*p)  ==> Tot_N*J
#### Change Beta  ===> Change BeZ =====> Chage Mu
    
### MMu:  Tot_N*J matrix: J otus and Tot_N samples ** Note ** some rows are repitition due to replicate.
### r: Tot_N
### change th0 ===> change mu
*/
    
    void fn_update_Th0(int JJ, int Tot_N, double *YY, double *SS, double *RR, double *Th0_bar, double u2_th0, double *BeZ, double *KTh, double *Th0_cur, double *MMu_cur)
    {
        int j, i_sam, ind2;
        double Th0_j_pro, u, m_cur, m_pro, A_cur, A_pro, KTh_BeZ;
        double *mu_pro;
        mu_pro = new double[Tot_N]; //Tot_N vector

        
        //## one OTU at a time
        for(j=0; j < JJ; j++)
        {
            ind2 = Tot_N*j;  // index for Tot_N*J dim matrix -- Y, Kth, BeZ, MMu

            //propose new th_{0j}
            GetRNGstate();
            Th0_j_pro = Th0_cur[j] + rnorm(0.0, 0.3); //0.4
            PutRNGstate();

            // prior
            A_cur = - pow(Th0_cur[j] - Th0_bar[j], 2.0)/2.0/u2_th0;
            A_pro = - pow(Th0_j_pro - Th0_bar[j], 2.0)/2.0/u2_th0;
            
            for (i_sam=0; i_sam < Tot_N; i_sam++) {
                
                // compute mu and m with proposed th0_j_pro
                KTh_BeZ = KTh[ind2 + i_sam] + BeZ[ind2 + i_sam];
                
                mu_pro[i_sam] = exp(Th0_j_pro + KTh_BeZ);     //# Tot_N-dim
                m_pro = RR[i_sam]*SS[j]*mu_pro[i_sam];        //# Tot_N-dim
                
                m_cur = RR[i_sam]*SS[j]*MMu_cur[ind2 + i_sam];
                
                //### Y: Tot_N*J matrix: J otus and Tot_N samples
                A_cur = A_cur + YY[ind2 + i_sam]*log(m_cur/(1.0+m_cur)) - 1.0/SS[j]*log(1.0 + m_cur);
                A_pro = A_pro + YY[ind2 + i_sam]*log(m_pro/(1.0+m_pro)) - 1.0/SS[j]*log(1.0 + m_pro);
            } // for (i_sam=0; i_sam < Tot_N; i_sam++) {
            
            //Let's do M-H
            GetRNGstate();
            u = runif(0.0,1.0);
            PutRNGstate();
            
            //### Accept?
            if(log(u) < (A_pro - A_cur))
            {
                // accept
                Th0_cur[j] = Th0_j_pro;
                
                for (i_sam=0; i_sam < Tot_N; i_sam++) {
                    MMu_cur[ind2 + i_sam] = mu_pro[i_sam];    //# Tot_N-dim
                }
            } //if(log(u) < (prob_pro - prob_cur))
        }// for(j=0; j < JJ; j++)
        
        // Free the memory
        delete[] mu_pro; mu_pro = NULL;

    }
    
/*
################################################################################################################################
### update s_j (overdispersion parameter for each OTU)
################################################################################################################################
###
 
 ### MMu:  Tot_N*J matrix: J otus and Tot_N samples ** Note ** some rows are repitition due to replicate.
 ### r: Tot_N

# s: J-dim vector
#######################
*/
    void fn_update_s(double a_s, double b_s, double *YY, double *RR, double *MMu, int JJ, int Tot_N, double *s_cur)
    {
        int j, i_sam, ind2;
        double s_j_cur, s_j_pro, inv_s_j_cur, inv_s_j_pro, A_cur, A_pro;
        double R_Mu, m_cur, m_pro, u;
        
        //### for each OTU (j)
        for(j=0; j < JJ; j++)
        {
            ind2 = Tot_N*j;  // index for Tot_N*J dim matrix -- Y, Kth, BeZ, MMu
            
            //## M-H
            s_j_cur = s_cur[j];
            
            GetRNGstate();
            s_j_pro = exp(log(s_j_cur) + rnorm(0.0, 0.8));
            PutRNGstate();
            
            inv_s_j_cur = 1.0/s_j_cur;
            inv_s_j_pro = 1.0/s_j_pro;
            
            // acceptance prob
            A_cur = a_s*log(s_j_cur) - b_s*s_j_cur; //###  Prior + Jacobian
            A_pro = a_s*log(s_j_pro) - b_s*s_j_pro; //###  Prior + Jacobian
            
            // likelihood
            for (i_sam=0; i_sam < Tot_N; i_sam++) {
                R_Mu = RR[i_sam]*MMu[ind2 + i_sam];   //## OTU j:  (K[1], K[2], ... K[n])-dim
                
                m_cur = s_j_cur*R_Mu;
                m_pro = s_j_pro*R_Mu;

                A_cur = A_cur + fn_lgamma(YY[ind2 + i_sam] + inv_s_j_cur) + YY[ind2 + i_sam]*log(m_cur/(1.0+m_cur)) - inv_s_j_cur*log(1.0 + m_cur);
                A_pro = A_pro + fn_lgamma(YY[ind2 + i_sam] + inv_s_j_pro) + YY[ind2 + i_sam]*log(m_pro/(1.0+m_pro)) - inv_s_j_pro*log(1.0 + m_pro);
            } //  for (i_sam=0; i_sam < Tot_N; i_sam++) {
            
            A_cur = A_cur - 1.0*Tot_N*lgamma(inv_s_j_cur);
            A_pro = A_pro - 1.0*Tot_N*lgamma(inv_s_j_pro);
            
            //Let's do M-H
            GetRNGstate();
            u = runif(0.0, 1.0);
            PutRNGstate();
            
            //### Accept?
            if(log(u) < (A_pro - A_cur))
            {
                s_cur[j] = s_j_pro;
            }
        } // for(j=0; j < JJ; j++)
    }
    
    
    /*
     ################################################################################################################################
     ## Sig2 : sig2_j: var for beta_{pj}
     ###   \Beta: p*J matrix.
     ####  \phi: p*J matrix.
     ################################################################################################################################
     */
    void fn_update_sig2(int JJ, int PP, double *Beta, double *phi, double a_sig, double b_sig, double *sig2_cur)
    {
        int p, j;
        
        double aa, bb;
        
        aa = a_sig + 1.0*PP/2.0;

        // for each covariate at a time
        for (j=0; j < JJ; j++) {
            
            bb = b_sig;
            for (p=0; p < PP; p++) {
                bb = bb + pow(Beta[PP*j + p], 2.0)/2.0/phi[PP*j + p];
            }
            
            GetRNGstate();
            sig2_cur[j] = 1.0/rgamma(aa, 1.0/bb);
            PutRNGstate();
            
        } //for (j=0; j < JJ; j++) {
    }
    

    /*
     ################################################################################################################################
     ## lam2 : phi_jp ~ Exp(lam2_j/2)
     ################################################################################################################################
     */
    void fn_update_lam2(int JJ, int PP, double *phi, double a_lam, double b_lam, double *lam2_cur)
    {
        int p, j;
        
        double aa, bb;
        
        aa = a_lam + 1.0*PP;
        
        // for each covariate at a time
        for (j=0; j < JJ; j++) {
            
            bb = b_lam;
            for (p=0; p < PP; p++) {
                bb = bb + phi[PP*j + p]/2.0;
            }
            
            GetRNGstate();
            lam2_cur[j] = rgamma(aa, 1.0/bb);  // mean=a/b
            PutRNGstate();
            
        } //for (j=0; j < JJ; j++) {
    }
    

    /*
     ################################################################################################################################
     ## draw a random variate from Inverse Gaussin(mu, lam)
     ################################################################################################################################
     */
    double fn_sample_InverseGau(double mu, double lam)
    {
        double y, y2, mu2, x, u;
        
       // Steps 1 & 2
        GetRNGstate();
        y = rnorm(0, 1.0);
        PutRNGstate();
        
        y = y*y;  // y \sim \chi_sq(1)
        y2 = y*y;
        mu2 = mu*mu;
        
        // Step 3
        x = mu + mu2*y/2.0/lam - mu/2.0/lam*sqrt(4.0*mu*lam*y + mu2*y2);
        
        // Step 4
        GetRNGstate();
        u = runif(0.0, 1.0);
        PutRNGstate();

        if (u < (mu/(mu+x))) {
            return(x);
        }else{
            return(mu2/x);
        }
    }

    
    /*
     ################################################################################################################################
     ## phi_jp : phi_jp ~ Exp(lam2_j/2)
     ###   \Beta: p*J matrix.
     ####  \phi: p*J matrix.
     ################################################################################################################################
     */
    void fn_update_phi(int JJ, int PP, double *Beta, double *lam2, double *sig2, double *phi_cur)
    {
        int p, j;
        
        double lam_prime, mu_prime, phi_jp;
        
        
        // for each covariate at a time
        for (j=0; j < JJ; j++) {
            
            for (p=0; p < PP; p++) {
                
                lam_prime = lam2[j];
                mu_prime = sqrt(lam2[j]*sig2[j]/pow(Beta[PP*j + p], 2.0));
                
                phi_jp = fn_sample_InverseGau(mu_prime, lam_prime);
                
                // to avoid numerical issues
                if (phi_jp > 0.0000001) {
                    phi_cur[PP*j + p] = 1.0/phi_jp;
                }else{
                    phi_cur[PP*j + p] = 1.0/0.0000001;
                    
                    printf("\n j=%d, p=%d, phi error", j, p);
                } //if (phi_jp > 0.0000001) {
            } //for (p=0; p < PP; p++) {
        } //for (j=0; j < JJ; j++) {
    }
    
    
    
    
    
    void fn_update_IG(double *a, double *b, double *sam)
        {
        int p;
        
        // for each covariate at a time
        for (p=0; p < 10000; p++) {
            GetRNGstate();
            sam[p] = 1.0/rgamma(*a, 1.0/(*b));  // IG(a, b) with mean b/(a-1)  x^{-a-1}\exp(-b/x)
            PutRNGstate();
            
        } //for (p=0; p < PP; p++) {
    }
    


    
    void fn_update_Gam(double *a, double *b, double *sam)
    {
        int p;
        
        // for each covariate at a time
        for (p=0; p < 10000; p++) {
            GetRNGstate();
            sam[p] = rgamma(*a, 1.0/(*b));  // Ga(a, b) with mean a/b
            PutRNGstate();
            
        } //for (p=0; p < PP; p++) {
    }

    
    
/*
################################################################################################################################
## tau2_0 : prior variance for th0
## tau2_j : prior variance for th_mj
################################################################################################################################
*/
//#### Th: M*J matrix
    
    void fn_update_tau2(int JJ, int MM, double aa_tau, double bb_tau, double *Th, double *tau2_cur)
    {
        int j, m;
        double aa, bb;
        
        //###  tau2_j for var of th_mj
        // for each covariate at a time
        aa = aa_tau + 1.0*MM/2.0;
        
        for (j=0; j < JJ; j++) {
            
            bb = bb_tau;
            for (m=0; m < MM; m++) {
                bb = bb + pow(Th[m*JJ + j], 2.0)/2.0;
            }
            
            GetRNGstate();
            tau2_cur[j] = 1.0/rgamma(aa, 1.0/bb);
            PutRNGstate();
            
        } //for (p=0; p < PP; p++) {
    }

    /*
     ################################################################################################################################
     ## r_ik : from LL*2 components
     ## indexes by Delta
     ## Mu: exp(th0 + K*th + beta*Z) -- Tot_N*J
     ################################################################################################################################
     ### MMu:  Tot_N*J matrix: J otus and Tot_N samples ** Note ** some rows are repitition due to replicate.
     ################################################################################################################################
     
     
     ### Y: Tot_N*J matrix: J otus and Tot_N samples
     ####  ==> KTh: Tot_N*J matrix
     #### BeZ: (p*J)*(Tot_N*p)  ==> Tot_N*J
     
     ### r: Tot_N
     #### s: J-dim vector
     */
    
    void fn_update_r(double JJ, double *rt_m, int Tot_N, double *YY, double *MMu, double *SS, double u2, double *RR_cur)
    {
        int i_sam, j;
        double r_ik_cur, r_ik_pro, rt_ik_cur, rt_ik_pro, u, A_pro, A_cur, m_cur, m_pro;
        
        for(i_sam=0; i_sam < Tot_N; i_sam++)
        {
            //##  M-H
            r_ik_cur = RR_cur[i_sam];
            rt_ik_cur = log(r_ik_cur);
            
            GetRNGstate();
            rt_ik_pro = rt_ik_cur + rnorm(0, 0.2);
            PutRNGstate();
            
            r_ik_pro = exp(rt_ik_pro);
            
            // prior
            A_cur = - pow(rt_ik_cur - rt_m[i_sam], 2.0)/2.0/u2;
            A_pro = - pow(rt_ik_pro - rt_m[i_sam], 2.0)/2.0/u2;
            
            // go over all OTUs
            for (j=0; j < JJ; j++) {
                
                m_cur = r_ik_cur*MMu[j*Tot_N + i_sam]*SS[j];
                m_pro = r_ik_pro*MMu[j*Tot_N + i_sam]*SS[j];
                
                A_cur = A_cur + YY[j*Tot_N + i_sam]*log(m_cur/(1.0+m_cur)) - 1.0/SS[j]*log(1.0 + m_cur);
                A_pro = A_pro + YY[j*Tot_N + i_sam]*log(m_pro/(1.0+m_pro)) - 1.0/SS[j]*log(1.0 + m_pro);
            }
            
            //Let's do M-H
            GetRNGstate();
            u = runif(0.0,1.0);
            PutRNGstate();
            
            //### Accept?
            if(log(u) < (A_pro - A_cur))
            {
                RR_cur[i_sam] = r_ik_pro;
            }
        } //## for(i in 1:Tot_N)
        
    }
    

/*
################################################################################################################################
## Th_0j and r_j : baseline of  OTU j
################################################################################################################################
    
### Y: Tot_N*J matrix: J otus and Tot_N samples
    
#### K: Tot_N * M where  M: the number of basis points
#### Th0: J-dim vector
#### Th: M*J matrix
####  ==> KTh: Tot_N*J matrix
    
#### Z: Tot_N*p covariates   *** Note ***: Z is the same for all K_t replicates at time t
#### BeZ: (p*J)*(Tot_N*p)  ==> Tot_N*J
    
#### s: J-dim vector
    
    
### MMu:  Tot_N*J matrix: J otus and Tot_N samples ** Note ** some rows are repitition due to replicate.
### r: Tot_N
### change th0 ===> change mu
*/
    
    void fn_update_r_and_Th0(int JJ, int Tot_N, double *YY, double *SS, double u2, double *Th0_bar, double u2_th0, double *BeZ, double *KTh, double *rt_m, double *RR_cur, double *Th0_cur, double *MMu_cur)
    {
        int i_sam, j, ind2;
        double rt_ik_cur, rt_ik_pro, r_ik_cur, r_ik_pro, m_cur, m_pro, u, mu_pro;
      
        double *RR_pro, *Th0_pro, *MMu_pro;
        
        //Allocate the memory
        RR_pro = new double[Tot_N]; //Tot_N vector
        Th0_pro = new double[JJ]; //J vector
        MMu_pro = new double[Tot_N*JJ]; //Tot_N*JJ vector
        
        double A_cur = 0.0;
        double A_pro = 0.0;
        
        for(i_sam=0; i_sam < Tot_N; i_sam++)
        {
            //##  M-H
            r_ik_cur = RR_cur[i_sam];
            rt_ik_cur = log(r_ik_cur);
            
            GetRNGstate();
            rt_ik_pro = rt_ik_cur + rnorm(0, 0.01);
            PutRNGstate();
            
            r_ik_pro = exp(rt_ik_pro);
            RR_pro[i_sam] = r_ik_pro;
            
            // prior -- rt_ik
            A_cur = A_cur - pow(rt_ik_cur - rt_m[i_sam], 2.0)/2.0/u2;
            A_pro = A_pro - pow(rt_ik_pro - rt_m[i_sam], 2.0)/2.0/u2;
        }
        
        // go over all OTUs
        for (j=0; j < JJ; j++) {
            
            ind2 = Tot_N*j;  // index for Tot_N*J dim matrix -- Y, Kth, BeZ, MMu
            
            //propose new th_{0j}
            GetRNGstate();
            Th0_pro[j] = Th0_cur[j] + rnorm(0.0, 0.01);
            PutRNGstate();
            
            // prior
            A_cur = A_cur - pow(Th0_cur[j]-Th0_bar[j], 2.0)/2.0/u2_th0;
            A_pro = A_pro - pow(Th0_pro[j]-Th0_bar[j], 2.0)/2.0/u2_th0;
        }
        
        // Evaluate the likelihood
        //## one OTU at a time
        for(j=0; j < JJ; j++)
        {
            ind2 = Tot_N*j;  // index for Tot_N*J dim matrix -- Y, Kth, BeZ, MMu
            
            for (i_sam=0; i_sam < Tot_N; i_sam++) {
                
                m_cur = RR_cur[i_sam]*SS[j]*MMu_cur[ind2 + i_sam];
                
                // compute mu and m with proposed th0_pro
                // Tot_N*J dim matrix -- Y, Kth, BeZ, MMu
                mu_pro = exp(Th0_pro[j] + KTh[ind2 + i_sam] + BeZ[ind2 + i_sam]);
                MMu_pro[ind2 + i_sam] = mu_pro;
                m_pro = RR_pro[i_sam]*SS[j]*mu_pro;        //# Tot_N-dim
                
                //### Y: Tot_N*J matrix: J otus and Tot_N samples
                A_cur = A_cur + YY[ind2 + i_sam]*log(m_cur/(1.0+m_cur)) - 1.0/SS[j]*log(1.0 + m_cur);
                A_pro = A_pro + YY[ind2 + i_sam]*log(m_pro/(1.0+m_pro)) - 1.0/SS[j]*log(1.0 + m_pro);
            } // for (i_sam=0; i_sam < Tot_N; i_sam++) {
        } //for(j=0; j < JJ; j++)
        
        //Let's do M-H
        GetRNGstate();
        u = runif(0.0, 1.0);
        PutRNGstate();
        
        //### Accept?
        if(log(u) < (A_pro - A_cur))
        {
            // accept
            for(i_sam=0; i_sam < Tot_N; i_sam++)
            {
                RR_cur[i_sam] = RR_pro[i_sam];
            }
            
            // go over all OTUs
            for (j=0; j < JJ; j++) {
                
                Th0_cur[j] = Th0_pro[j];
                
                // corrected Jan-30-18
                ind2 = Tot_N*j;  // index for Tot_N*J dim matrix -- Y, Kth, BeZ, MMu
                for (i_sam=0; i_sam < Tot_N; i_sam++) {
                    MMu_cur[ind2 + i_sam] = MMu_pro[ind2 + i_sam];
                } // for (i_sam=0; i_sam < Tot_N; i_sam++) {
            } //  for (j=0; j < JJ; j++) {
        } // if(log(u) < (A_pro - A_cur))
        
        
        // Free the memory
        delete[] RR_pro; RR_pro = NULL;
        delete[] MMu_pro; MMu_pro = NULL;
        delete[] Th0_pro; Th0_pro = NULL;
        

    }
    
    
    void fn_update_Beta_Theta(int *JJ, int *PP, int *MM, int *Tot_N, double *YY,
                              double *a_sig, double *b_sig, double *sig2,
                              double *a_tau, double *b_tau, double *tau2,
                              double *a_lam, double *b_lam, double *lam2,
                              double *Rt_m, double *u2, double *RR,
                              double *a_s, double *b_s, double *SS,
                              double *Th0_bar, double *u2_th0, double *Th0,
                              double *KK, double *Th, double *KTh,
                              double *ZZ, double *phi, double *Beta, double *BeZ,
                              double *MMu)
    {
        int i;
        
        for (i=0; i < 2; i++) {
            //#### Change Beta  ===> Change BeZ =====> Chage Mu
            fn_update_Beta(*JJ, *PP, *Tot_N, ZZ, Th0, KTh, YY, SS, RR, sig2, phi, Beta, BeZ, MMu);
            
            //#### Change ALL Beta together, one otu at a time ===> Change BeZ =====> Chage Mu
            fn_update_ALL_Beta(*JJ, *PP, *Tot_N, ZZ, Th0, KTh, YY, SS, RR, sig2, phi, Beta, BeZ, MMu);
        }
        
        //#### Change Th  ===> Change KTh =====> Chage Mu
        fn_update_ThP(*JJ, *MM, *Tot_N, KK, YY, SS, RR, tau2, BeZ, Th0, Th, KTh, MMu);

        // ### change th0 ===> change mu
        fn_update_Th0(*JJ, *Tot_N, YY, SS, RR, Th0_bar, *u2_th0, BeZ, KTh, Th0, MMu);
        
        /*
        for (i=0; i < 5; i++) {
            //#### Change ALL Beta and th0 together, one otu at a time ===> Change BeZ =====> Chage Mu
            fn_update_ALL_Beta_Th0(*JJ, *PP, *Tot_N, ZZ, Th0_bar, *u2_th0, Th0, KTh, YY, SS, RR, sig2, phi, Beta, BeZ, MMu);
        }
         */

        // ### update r - library size adjustment
        fn_update_r(*JJ, Rt_m, *Tot_N, YY, MMu, SS, *u2, RR);
        
        // ### update s_j (overdispersion parameter for each OTU)
        fn_update_s(*a_s, *b_s, YY, RR, MMu, *JJ, *Tot_N, SS);
        
        // #### update sig2_j
        fn_update_sig2(*JJ, *PP, Beta, phi, *a_sig, *b_sig, sig2);
        
        // ### update phi_pj
        fn_update_phi(*JJ, *PP, Beta, lam2, sig2, phi);
        
        // ### update lam2_j
        fn_update_lam2(*JJ, *PP, phi, *a_lam, *b_lam, lam2);

        //### update tau2_0 and tau2
        fn_update_tau2(*JJ, *MM, *a_tau, *b_tau, Th, tau2);
        
        for (i=0; i < 5; i++) {
            //### jointly change r and th0 ===> change mu
            fn_update_r_and_Th0(*JJ, *Tot_N, YY, SS, *u2, Th0_bar, *u2_th0, BeZ, KTh, Rt_m, RR, Th0, MMu);
        }
    }
    

}
