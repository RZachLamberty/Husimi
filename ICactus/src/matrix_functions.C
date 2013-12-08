#include"bethe_lapack_interface.h"
#include"matrix_functions.h"
#include"printing_functions.h"
#include"math_utilities.h"
using namespace std;

/////////////////////////////////////////////////////////////////////////////
void general_real_diagonalize(Matrix &a, 
                 std::vector< complex<double> > &eigs, 
                 Matrix &eigenvecs)
{
    int n=int(a.NRows());
    int info=0;
    std::vector<double> reigs(n),imeigs(n);
    std::vector<double> work(4*n);
    Matrix vl(n,n), b(n,n);
    
    b=a;
    eigs.clear();
    eigenvecs.resize(n,n);

    work.assign(4*n,0.);
         
    dgeev('N','V',n,&*b.begin(),n,&*reigs.begin(),&*imeigs.begin(),
          &*vl.begin(),n,&*eigenvecs.begin(),n,
           &*work.begin(),int(work.size()),info);
    
    for (int i=0;i<n;i++) eigs.push_back(complex<double>(reigs[i],imeigs[i]));

    //gram_schmidt(eigenvecs,eigs);
}
//////////////////////////////////////////////////////////////////////////////
void perform_lu(Matrix &mat, Matrix &lu, std::vector<int> &ipiv)
{
    lu.clear();
    lu.resize(mat.NRows(),mat.NCols());
    lu=mat;

    ipiv.clear();
    ipiv.resize(int(min(lu.NRows(),lu.NCols())));

    int info=0;
    dgetrf(int(lu.NRows()),int(lu.NCols()),&*lu.begin(),int(lu.NRows()),&*ipiv.begin(),info);
}


void invert_real_matrix(Matrix &mat, Matrix &inverse)
{
    int j;
    int info=0;
    int n=int(mat.NRows());

    std::vector<int> ipiv(n);

    // Perform LU decomposition
    perform_lu(mat,inverse,ipiv);

    std::vector<double> work(8*n);
    work.assign(8*n,0.);
    int lwork=work.size();

    // Now invert the LU matrix
    dgetri(n,&*inverse.begin(),n,&*ipiv.begin(),&*work.begin(),lwork,info);

}

//////////////////////////////////////////////////////////////////////////////
void symmetric_diagonalize( Matrix &a, 
                 std::vector<double> &eigs, 
                 Matrix &eigenvecs)
{
    int w,info=0;
    int n=int(a.NRows());
    std::vector<double> work(4*n);

    //cout<<"Making copy of matrix"<<endl;
    eigenvecs=a;
    work.assign(4*n,0.);

    //print_real_mat(a);

    //cout<<"Actual diagonalization"<<endl;
    dsyev('V','U',n,&*eigenvecs.begin(),n,
           &*eigs.begin(),&*work.begin(),int(work.size()),info);
    
    //cout<<"Finished diagonalizing (now returning)"<<endl;

    // Sorted the eigenvalues and corresponding eigenvectors
    // Eigenvectors are columns of the matrix named eigenvecs

}

//////////////////////////////////////////////////////////////////////////////

void gram_schmidt( Matrix &eigenvecs,std::vector< complex<double> > &eigs)
{
	std::vector< std::vector<double> > v_s;
	int size=eigenvecs.NRows();
	std::vector<double> tmp_vec(size);
	double q;

	for (int j=0;j<size;j++)
	{
		for (int m=0;m<size;m++)
		{tmp_vec[m]=eigenvecs(m,j);}

		
       		for (int k=0;k<j;k++)
       		{
            		if (abs(real(eigs[j])-real(eigs[k]))<1.0e-10)
			{
				cout<<"Orthogonalizing subspace"<<endl;
				q=ddot(size,&*tmp_vec.begin(),1,&*v_s[k].begin(),1);
				daxpy(size,-q,&*v_s[k].begin(),1,&*tmp_vec.begin(),1);
			}
		}
            	q=ddot(size,&*tmp_vec.begin(),1,&*tmp_vec.begin(),1);
		dscal(size,1.0/sqrt(q),&*tmp_vec.begin(),1);

		v_s.push_back(tmp_vec);
		for (int m=0;m<size;m++)
		{eigenvecs(m,j)=tmp_vec[m];}
	}
}

//////////////////////////////////////////////////////////////////////////////
void eta_orthogonalize( std::vector< std::vector<double> > &eigenvecs,
			std::vector<double> const &eta)

{
  	std::vector< std::vector<double> > tmp_eigenvecs;
	int 				   size=eigenvecs[0].size();
	int 				   num_vecs=eigenvecs.size();
	std::vector<double>		   overlap(num_vecs);
	std::vector<double>		   old_temp_vec;
	double				   prod_overlap=0.0;		
	double				   q=0.0;	
	int 				   tries;
        int 				   tries_max=100;
 	double 				   rnd;
	
	tries=0;
	tmp_eigenvecs=eigenvecs;
	
	do
	{
		//cout<<"Trying to orthogonalize"<<endl;
		prod_overlap=1.0;
		for (int j=0;j<num_vecs;j++)
		{
			old_temp_vec=tmp_eigenvecs[j];
			for (int k=0;k<j;k++)
			{
					//cout<<"Eta Orthogonalizing subspace"<<endl;
					q=0.0;
					for (int l=0;l<size;l++) q=q+(old_temp_vec[l]*tmp_eigenvecs[k][l]*eta[l]);
					daxpy(size,-q/overlap[k],&*tmp_eigenvecs[k].begin(),1,&*tmp_eigenvecs[j].begin(),1);
			}
			q=ddot(size,&*tmp_eigenvecs[j].begin(),1,&*tmp_eigenvecs[j].begin(),1);
			dscal(size,1.0/sqrt(q),&*tmp_eigenvecs[j].begin(),1);
			overlap[j]=0.0;
			for (int l=0;l<size;l++) overlap[j]=overlap[j]+(tmp_eigenvecs[j][l]*tmp_eigenvecs[j][l]*eta[l]);
			prod_overlap=prod_overlap*overlap[j];
			//cout<<"overlap[j]="<<overlap[j]<<endl;

			if (abs(prod_overlap)<1.0e-7) 
			{
				prod_overlap=0.0;
				j=num_vecs;
			}

		}
		
		if (pow(abs(prod_overlap),1.0/double(num_vecs))<1.0e-7)
		{
			tmp_eigenvecs=eigenvecs;
			for (int j=0;j<num_vecs;j++)
			{
				for (int k=0;k<num_vecs;k++)
				{
					rnd=uniform_rnd();
					daxpy(size,rnd,&*eigenvecs[k].begin(),1,&*tmp_eigenvecs[j].begin(),1);
				}
				q=ddot(size,&*tmp_eigenvecs[j].begin(),1,&*tmp_eigenvecs[j].begin(),1);
				dscal(size,1.0/sqrt(q),&*tmp_eigenvecs[j].begin(),1);
			}
		}

		tries=tries+1;

	} while (abs(prod_overlap)<1.0e-7 and tries<tries_max) ;

	if (tries==tries_max) {
	//cout<<"I could not Eta orthogonalize this set for SOME reason"<<endl;
	}
	
	eigenvecs=tmp_eigenvecs;
}

///////////////////////////////////////////////////////////////////////////////

void real_matrix_multiply( Matrix &a, Matrix &b, Matrix &c)
{
    double alpha=1.0;
    double beta=0.0;
    
    //Matrix a1,b1;
    //a1=a;b1=b;

    c.resize(a.NRows(),b.NCols());

    /*dgemm('N', 'N' , int(a.NRows()) ,int(b.NCols()) ,int(a.NCols()) ,alpha ,
          &*a1.begin(), int(a.NRows()), &*b1.begin(), int(b.NRows()),
          beta,&*c.begin(),int(c.NRows()));  */
    
    dgemm('N', 'N' , int(a.NRows()) ,int(b.NCols()) ,int(a.NCols()) ,alpha ,
          &*a.begin(), int(a.NRows()), &*b.begin(), int(b.NRows()),
          beta,&*c.begin(),int(c.NRows()));
}

//////////////////////////////////////////////////////////////////////////////

void real_matrix_times_vector( Matrix &a,
                              std::vector<double> &b,
                              std::vector<double> &c)
{
    double alpha=1.0;
    double beta=0.0;
    int brows=b.size();
    
    //Matrix a1;
    //std::vector<double> b1;
    //a1=a;
    //b1=b;

    c.resize(brows);

    dgemm('N', 'N' , int(a.NRows()) ,1 ,int(a.NCols()) ,alpha ,
          &*a.begin(), int(a.NRows()), &*b.begin(), brows,
          beta,&*c.begin(),brows);   
}

//////////////////////////////////////////////////////////////////////////////

void real_vector_times_vector(std::vector<double> &a,           
                              std::vector<double> &b,
                              Matrix &c)
{
    double alpha=1.0;
    double beta=0.0;
    int brows=b.size();

    c.resize(brows,brows);

    dgemm('N', 'T' , a.size(),b.size(),1,alpha ,
          &*a.begin(), a.size(), &*b.begin(), b.size(),
          beta,&*c.begin(),b.size());   
}

//////////////////////////////////////////////////////////////////////////////

void svd( Matrix &a, 
          std::vector<double> &eigs, 
          Matrix &u, Matrix &vt)
{
    int info=0;
    int m=int(a.NRows());
    int n=int(a.NRows());
    std::vector<double> work(4*m);
    Matrix a_copy;

    a_copy=a;
    work.assign(4*m,0.);

    u.resize(m,m);
    vt.resize(n,n);
    eigs.resize(min(m,n));

    cout<<"Performing SVD"<<endl;
    dgesvd('A','A',m,n,&*a_copy.begin(),m,
           &*eigs.begin(),&*u.begin(),m,
	   &*vt.begin(),n,
	   &*work.begin(),int(work.size()),info);
    
    cout<<"Finished SVD"<<endl;

}
//////////////////////////////////////////////////////////////////////////////
void matrix_lanczos(int 			      		iterations,
		    int 			      		how_many_eigenvecs, 
                    Matrix 					&a,
		    std::vector<double> 			&eigs,
		    Matrix 					&eigenvecs,
	            bool 					ipr)
{
   time_t 				start,end;
   int 					i,it,j,k,l,m;
   int 					size=a.NRows();
   bool 				orth_failed;
   std::vector<double>			tmp_eigs;

   iterations=min(iterations,size);
   how_many_eigenvecs=min(how_many_eigenvecs,iterations);
   eigs.resize(how_many_eigenvecs);
   tmp_eigs.resize(iterations);
   
   double 				dif,tmp;
   double				q,alpha,beta,norm;
   std::vector<double> 			alphas,betas;
   std::vector<double> 			h_dot_v(size),w(size);
   std::vector<double> 			v_p(size),v_o(size),v_p_old(size),tmpv;
   std::vector< std::vector<double> >   vs;
   Matrix 				t_mat(iterations,iterations);
   Matrix				t_eigenvecs(iterations,iterations);

   // Initializations
   for (i=0;i<size;i++){v_p[i]=uniform_rnd();v_o[i]=0.0;w[i]=0.0;}
       
   norm=sqrt(ddot(size,&*v_p.begin(),1,&*v_p.begin(),1));
   dscal(size,1.0/norm,&*v_p.begin(),1);
  
   beta=0.0;betas.push_back(beta);

   for (it=0;it<iterations;it++)
   {
       if (ipr) cout<<"Doing Lanczos iteration "<<it<<endl;
       vs.push_back(v_p);     
       
       time (&start);
       // Computing H*v_p - This is the bulk of the operation
       for (i=0;i<size;i++) 
       {
	    if (abs(v_p[i])>1.0e-16)                       // Preventing additions of lots of zeros //HJC June 5 2011
	    {
		    for (k=0;k<size;k++)
		    {w[k]+=(a(k,i)*v_p[i]);}
	    }
       }
       
       daxpy(size,-beta,&*v_o.begin(),1,&*w.begin(),1);
       alpha=ddot(size,&*w.begin(),1,&*v_p.begin(),1);
       alphas.push_back(alpha);
       daxpy(size,-alpha,&*v_p.begin(),1,&*w.begin(),1);
       v_o=v_p;
       beta=sqrt(ddot(size,&*w.begin(),1,&*w.begin(),1));
       v_p=w;
       dscal(size,1.0/beta,&*v_p.begin(),1);
       betas.push_back(beta);
       dscal(size,0.0,&*w.begin(),1);

       // Now reorthogonalize vectors
       if (vs.size()<iterations)
       {
	       orth_failed=false;
	       for (i=0;i<vs.size();i++)
	       {
		    q=ddot(size,&*v_p.begin(),1,&*vs[i].begin(),1);
		    if (abs(q)>1.0e-10)
		    {
			i=it;
			if (ipr)
			{
				cout<<"q (overlap) ="<<q<<endl;
				cout<<"--------------------------------------------------"<<endl;
				cout<<"Orthogonalization failed... choosing random vector"<<endl;
				cout<<"--------------------------------------------------"<<endl;
			}
			orth_failed=true;
		    }
		    else
		    {
			daxpy(size,-q,&*vs[i].begin(),1,&*v_p.begin(),1);
			dscal(size,1.0/(1.0-(q*q)),&*v_p.begin(),1);
		    }
	       }
	       
	       norm=0.0;
	       if (orth_failed)
	       {
			while (abs(norm)<1.0e-2)
			{
				for (j=0;j<size;j++)  {v_p_old[j]=(1.0+uniform_rnd());}
				norm=sqrt(ddot(size,&*v_p_old.begin(),1,&*v_p_old.begin(),1));
				dscal(size,1.0/norm,&*v_p_old.begin(),1);
				//cout<<"norm="<<norm<<endl;
				v_p=v_p_old;
				for (i=0;i<vs.size();i++)
				{
					q=ddot(size,&*v_p_old.begin(),1,&*vs[i].begin(),1);
					//cout<<"q (after orth fail )="<<q<<endl;
					daxpy(size,-q,&*vs[i].begin(),1,&*v_p.begin(),1);
				}
				norm=sqrt(ddot(size,&*v_p.begin(),1,&*v_p.begin(),1));
			}
			dscal(size,1.0/norm,&*v_p.begin(),1);
	       }
       }

       /*for (i=0;i<it;i++)
       {
            q=ddot(size,&*v_p.begin(),1,&*vs[i].begin(),1);
            daxpy(size,-q,&*vs[i].begin(),1,&*v_p.begin(),1);
            dscal(size,1.0/(1.0-(q*q)),&*v_p.begin(),1);
       }*/

       time (&end);
       dif=difftime(end,start);
       
       if (ipr) cout<<"Time to perform Lanczos iteration "<<it<<" was "<<dif<<" seconds"<<endl;
   }

   if (ipr) cout<<"Now building T matrix"<<endl;
   
   time (&start);
   for (j=0;j<iterations;j++)
   {
        t_mat(j,j)=alphas[j];
        if (j+1<iterations)
        {t_mat(j,j+1)=betas[j+1];t_mat(j+1,j)=betas[j+1];}
   }
   time (&end);
   dif=difftime(end,start);
  
   if (ipr) cout<<"Time to build T was "<<dif<<" seconds"<<endl;
 
   time (&start);
   symmetric_diagonalize(t_mat,tmp_eigs,t_eigenvecs);
   time (&end);
   dif=difftime(end,start);

   if (ipr) cout<<"Time to diagonalize T was "<<dif<<" seconds"<<endl;
	   
   eigenvecs.resize(size,how_many_eigenvecs);

   //if (ipr) cout<<"how_many_eigenvecs="<<how_many_eigenvecs<<endl;//if (ipr) cout<<"size="<<size<<endl;

   for (j=0;j<how_many_eigenvecs;j++)
   {
	      eigs[j]=tmp_eigs[j];
	      if (ipr) cout<<"Making j = "<<j<<" Ritz eigenvector"<<endl;
	      for (i=0;i<size;i++)
	      {
		eigenvecs(i,j)=0.0;
		for (k=0;k<iterations;k++){eigenvecs(i,j)+=vs[k][i]*t_eigenvecs(k,j);}
	      }
   }      
   //f (ipr) print_real_mat(eigenvecs);
   if (ipr) cout<<"Done computing Ritz eigenvectors"<<endl;
}


