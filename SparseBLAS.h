/*
* Sparse BLAS (Basic Linear Algebra Subprograms) Library
*
* This code has been modified from the NIST distribution, in that the interface is templated.
* It thus no longer conforms to the ANSI C interface specification
* of the Sparse BLAS in the BLAS Technical Forum Standard[1].
*
* Mathematical and Computational Sciences Division
* National Institute of Technology,
* Gaithersburg, MD USA
*
*
* [1] BLAS Technical Forum: www.netlib.org/blas/blast-forum/
* [2] I. S. Duff, M. A. Heroux, R. Pozo, "An Overview of the Sparse Basic
*     Linear Algebra Subprograms: The new standard of the BLAS Techincal
*     Forum,"  Vol. 28, No. 2, pp. 239-267,ACM Transactions on Mathematical 
*     Software (TOMS), 2002.
*
*
* DISCLAIMER:
*
* This software was developed at the National Institute of Standards and
* Technology (NIST) by employees of the Federal Government in the course
* of their official duties. Pursuant to title 17 Section 105 of the
* United States Code, this software is not subject to copyright protection
* and is in the public domain. NIST assumes no responsibility whatsoever for
* its use by other parties, and makes no guarantees, expressed or implied,
* about its quality, reliability, or any other characteristic.
*
*
*/

#ifndef SPARSEBLAS_H
#define SPARSEBLAS_H

/* numeric is for accumulate() below */
#include <iostream>
#include <complex>
#include <numeric>
#include <vector>
#include <utility>
/* pair defined here */

typedef int blas_sparse_matrix;

enum blas_order_type {
            blas_rowmajor = 101,
            blas_colmajor = 102 };

enum blas_trans_type {
            blas_no_trans   = 111,
            blas_trans      = 112,
            blas_conj_trans = 113 };

enum blas_uplo_type  {
            blas_upper = 121,
            blas_lower = 122 };

enum blas_diag_type {
            blas_non_unit_diag = 131,
            blas_unit_diag     = 132 };

enum blas_side_type {
            blas_left_side  = 141,
            blas_right_side = 142 };

enum blas_cmach_type {
            blas_base      = 151,
            blas_t         = 152,
            blas_rnd       = 153,
            blas_ieee      = 154,
            blas_emin      = 155,
            blas_emax      = 156,
            blas_eps       = 157,
            blas_prec      = 158,
            blas_underflow = 159,
            blas_overflow  = 160,
            blas_sfmin     = 161};

enum blas_norm_type {
            blas_one_norm       = 171,
            blas_real_one_norm  = 172,
            blas_two_norm       = 173,
            blas_frobenius_norm = 174,
            blas_inf_norm       = 175,
            blas_real_inf_norm  = 176,
            blas_max_norm       = 177,
            blas_real_max_norm  = 178 };

enum blas_sort_type {
            blas_increasing_order = 181,
            blas_decreasing_order = 182 };

enum blas_conj_type {
            blas_conj    = 191,
            blas_no_conj = 192 };

enum blas_jrot_type {
            blas_jrot_inner  = 201,
            blas_jrot_outer  = 202,
            blas_jrot_sorted = 203 };

enum blas_base_type {
            blas_zero_base = 221,
            blas_one_base  = 222 };

enum blas_symmetry_type {
            blas_general          = 231,
            blas_symmetric        = 232,
            blas_hermitian        = 233,
            blas_triangular       = 234,
            blas_lower_triangular = 235,
            blas_upper_triangular = 236,
            blas_lower_symmetric  = 237,
            blas_upper_symmetric  = 238,
            blas_lower_hermitian  = 239,
            blas_upper_hermitian  = 240  };

enum blas_size_type {
            blas_num_rows      = 251,
            blas_num_cols      = 252,
            blas_num_nonzeros  = 253  };

enum blas_handle_type{
            blas_invalid_handle = 261,
            blas_new_handle     = 262,
            blas_open_handle    = 263,
            blas_valid_handle   = 264};

enum blas_sparsity_optimization_type {
            blas_regular       = 271,
            blas_irregular     = 272,
            blas_block         = 273,
            blas_unassembled   = 274 };

#ifdef SPBLAS_ERROR_FATAL
#include <cassert>
#define ASSERT_RETURN(x, ret_val) assert(x)
#define ERROR_RETURN(ret_val)  assert(0)
#else
#define ASSERT_RETURN(x, ret_val) {if (!(x)) return ret_val;}
#define ERROR_RETURN(ret_val) return ret_val
#endif

using namespace std;

namespace SparseBLAS
{


inline const double& conj(const double &x)
{ 
  return x;
}

inline const float& conj(const float &x)
{ 
  return x;
}

/*
   Generic sparse matrix (base) class: defines only the structure 
   (size, symmetry, etc.) and maintains state during construction, 
   but does not specify the actual nonzero values, or their type. 

*/

class Sp_mat
{
  private:

    int num_rows_;
    int num_cols_;
    int num_nonzeros_;

    /* ... */

    int void_;
    int nnew_;      /* avoid using "new" since it is a C++ keyword */
    int open_;
    int valid_;

    int unit_diag_ ;
    int upper_triangular_;
    int lower_triangular_;
    int upper_symmetric_;
    int lower_symmetric_;
    int upper_hermitian_;
    int lower_hermitian_;
    int general_;

    int one_base_;


        /* optional block information */

    int Mb_;                /* matrix is partitioned into Mb x Nb blocks    */
    int Nb_;                /* otherwise 0, if regular (non-blocked) matrix */
    int k_;                 /* for constant blocks, each block is k x l     */
    int l_;                 /* otherwise 0, if variable blocks are used.   */

    int rowmajor_;          /* 1,if block storage is rowm major.  */
    int colmajor_;          /* 1,if block storage is column major. */

    /* unused optimization paramters */

    int opt_regular_;
    int opt_irregular_;
    int opt_block_;
    int opt_unassembled_;

    vector<int> K_; /* these are GLOBAL index of starting point of block     */
    vector<int> L_; /* i.e. block(i,j) starts at global location (K[i],L[i]) */
                    /* and of size (K[i+1]-K[i] x L[i+1]-L[i])               */

  public:
    Sp_mat(int M, int N);

    int& num_rows();
    int& num_cols();
    int& num_nonzeros();
    int num_rows() const;
    int num_cols() const;
    int num_nonzeros() const;
    int is_one_base() const;
    int is_zero_base() const;
    int is_void() const;
    int is_new() const;
    int is_open() const;
    int is_valid() const;
    int is_unit_diag() const;
    int is_upper_triangular() const;
    int is_lower_triangular() const;
    int is_triangular() const;
    int is_lower_symmetric() const;
    int is_upper_symmetric() const;
    int is_symmetric() const;
    int is_lower_hermitian() const;
    int is_upper_hermitian() const;
    int is_hermitian() const;
    int is_general() const;
    int is_lower_storage() const;
    int is_upper_storage() const;

    int is_opt_regular() const;
    int is_opt_irregular() const;
    int is_opt_block() const;
    int is_opt_unassembled() const;

    int K(int i) const;
    int L(int i) const;

    int is_rowmajor() const;
    int is_colmajor() const;

    void set_one_base();
    void set_zero_base();

    void set_void();
    void set_new();
    void set_open();
    void set_valid();

    void set_unit_diag();
    void set_upper_triangular();
    void set_lower_triangular();
    void set_upper_symmetric();
    void set_lower_symmetric();
    void set_upper_hermitian();
    void set_lower_hermitian();

    void set_const_block_parameters(int Mb, int Nb, int k, int l);
    void set_var_block_parameters(int Mb, int Nb, const int *k, const int *l);
    virtual int end_construction();
    virtual void print() const;
    virtual void destroy();
    virtual ~Sp_mat();
};

// The table of sparse matrices
extern vector<Sp_mat *> Table;
extern unsigned int Table_active_matrices;
int Table_insert(Sp_mat* S);
int Table_remove(unsigned int i);

template <class T>
class TSp_mat : public Sp_mat
{
  private:
    vector< vector< pair<T, int> > > S;
    /* optional diag if matrix is
       triangular. Created
       at end_construction() phase */
    vector<T> diag;   

  private:
    inline T sp_dot_product( const vector< pair<T, int> > &r, 
                             const T* x, int incx ) const
    {
      T sum(0);
      
      if (incx == 1)
        {
          for ( typename vector< pair<T,int> >::const_iterator p = r.begin(); 
                p < r.end(); p++)
            {
              //sum = sum + p->first * x[p->second];
              sum += p->first * x[p->second];
            }
        }
      else /* incx != 1 */
        {
          for ( typename vector< pair<T,int> >::const_iterator p = r.begin(); 
                p < r.end(); p++)
            {
              //sum = sum + p->first * x[p->second * incx];
              sum += p->first * x[p->second * incx];
            }
        }
    
      return sum;
    }

    inline T sp_conj_dot_product( const vector< pair<T, int> > &r, 
                                  const T* x, int incx ) const
    {
      T sum(0);
      
      if (incx == 1)
        {
          for ( typename vector< pair<T,int> >::const_iterator p = r.begin(); 
                p < r.end(); p++)
            {
              sum += conj(p->first) * x[p->second];
            }
        }
      else /* incx != 1 */
        {
          for ( typename vector< pair<T,int> >::const_iterator p = r.begin(); 
                p < r.end(); p++)
            {
              //sum = sum + p->first * x[p->second * incx];
              sum += conj(p->first) * x[p->second * incx];
            }
        }
      
      return sum;
    }


    inline void sp_axpy( const T& alpha, const vector< pair<T,int> > &r, 
                         T*  y, int incy) const
    {
      if (incy == 1)
        {
          for (typename vector< pair<T,int> >::const_iterator p = r.begin(); 
               p < r.end(); p++)
            y[p->second] += alpha * p->first;  
          
        }
      else /* incy != 1 */
        {
          for (typename vector< pair<T,int> >::const_iterator p = r.begin(); 
               p < r.end(); p++)
            y[incy * p->second] += alpha * p->first;  
        } 
    } 
    
    inline void sp_conj_axpy( const T& alpha, const vector< pair<T,int> > &r, 
                              T*  y, int incy) const
    {
      if (incy == 1)
        {
          for (typename vector< pair<T,int> >::const_iterator p = r.begin(); 
               p < r.end(); p++)
            y[p->second] += alpha * conj(p->first);  
          
        }
      else /* incy != 1 */
        {
          for (typename vector< pair<T,int> >::const_iterator p = r.begin(); 
               p < r.end(); p++)
            y[incy * p->second] += alpha * conj(p->first);  
        } 
    } 
    
    void mult_diag(const T& alpha, const T* x, int incx, T* y, int incy) 
      const
    {
      const T* X = x;
      T* Y = y;
      typename vector<T>::const_iterator d= diag.begin();
      for ( ; d < diag.end(); X+=incx, d++, Y+=incy)
        {
          *Y += alpha * *d * *X;
        }
    }
    
    void mult_conj_diag(const T& alpha, const T* x, int incx, T* y, int incy) 
      const
    {
      const T* X = x;
      T* Y = y;
      typename vector<T>::const_iterator d= diag.begin();
      for ( ; d < diag.end(); X+=incx, d++, Y+=incy)
        {
          *Y += alpha * conj(*d) * *X;
        }
    }
    
    
    void nondiag_mult_vec(const T& alpha, const T* x, int incx, 
                          T* y, int incy) const
    {
      
      int M = num_rows();
      
      if (incy == 1)
        {
          for (int i=0; i<M; i++)
            y[i] += alpha * sp_dot_product(S[i], x, incx);
        }
    else
      {
        for (int i=0; i<M; i++)
          y[i * incy] += alpha * sp_dot_product(S[i], x, incx);
      }
    }
    
    void nondiag_mult_vec_conj(const T& alpha, const T* x, int incx, 
                               T* y, int incy) const
    {
      
      int M = num_rows();
      
      if (incy == 1)
        {
          for (int i=0; i<M; i++)
            y[i] += alpha * sp_conj_dot_product(S[i], x, incx);
        }
      else
        {
          for (int i=0; i<M; i++)
            y[i * incy] += alpha * sp_conj_dot_product(S[i], x, incx);
        }
    }
    
    void nondiag_mult_vec_transpose(const T& alpha, const T* x, int incx, 
                                    T* y, int incy) const
    {
      /* saxpy: y += (alpha * x[i]) row[i]  */
      
      int M = num_rows();
      const T* X = x;
      for (int i=0; i<M; i++, X += incx)
        sp_axpy( alpha * *X, S[i], y, incy);
    }
    
    void nondiag_mult_vec_conj_transpose(const T& alpha, const T* x, int incx, 
                                         T* y, int incy) const
    {
      /* saxpy: y += (alpha * x[i]) row[i]  */
      
      int M = num_rows();
      const T* X = x;
      for (int i=0; i<M; i++, X += incx)
        sp_conj_axpy( alpha * *X, S[i], y, incy);
    }
    
    void mult_vec(const T& alpha, const T* x, int incx, T* y, int incy) 
      const
    {
      nondiag_mult_vec(alpha, x, incx, y, incy);
      
      if (is_triangular() || is_symmetric())
        mult_diag(alpha, x, incx, y, incy);
      
      if (is_symmetric())
        nondiag_mult_vec_transpose(alpha, x, incx, y, incy);
    }
    
    
    void mult_vec_transpose(const T& alpha, const T* x, int incx, T* y, 
                            int incy) const
    {
      
      nondiag_mult_vec_transpose(alpha, x, incx, y, incy);
      
      if (is_triangular() || is_symmetric())
        mult_diag(alpha, x, incx, y, incy);
      
      if (is_symmetric())
        nondiag_mult_vec(alpha, x, incx, y, incy);
    }
    
    
    void mult_vec_conj_transpose(const T& alpha, const T* x, int incx, T* y, 
                                 int incy) const
    {
      
      nondiag_mult_vec_conj_transpose(alpha, x, incx, y, incy);
      
      if (is_triangular() || is_symmetric())
        mult_conj_diag(alpha, x, incx, y, incy);
      
      if (is_symmetric())
        nondiag_mult_vec_conj(alpha, x, incx, y, incy);
    }
    
    int triangular_solve(T alpha, T* x, int incx ) const
    {
      if (alpha == (T) 0.0)
        ERROR_RETURN(1);
      
      if ( ! is_triangular() )
        ERROR_RETURN(1);
      
      int N = num_rows();
      
      if (is_lower_triangular())
        {
          for (int i=0, ii=0; i<N; i++, ii += incx)
            {
              x[ii] = (x[ii] - sp_dot_product(S[i], x, incx)) / diag[i];
            }
          if (alpha != (T) 1.0)
            {
              for (int i=0, ii=0; i<N; i++, ii += incx)
                x[ii] /= alpha; 
            }
        }
      else if (is_upper_triangular())
        {
          
          for (int i=N-1, ii=(N-1)*incx ;   0<=i ;    i--, ii-=incx)
            {
              x[ii] = (x[ii] - sp_dot_product(S[i],x, incx)) / diag[i];
            }
          if (alpha != (T) 1.0)
            {
              for (int i=N-1, ii=(N-1)*incx ;   0<=i ;    i--, ii-=incx)
                x[ii] /= alpha; 
            }
          
        }
      else
        ERROR_RETURN(1);
      
      return 0;
    }
    
    int transpose_triangular_solve(T alpha, T* x, int incx) const
    {
      if ( ! is_triangular())
        return -1;
      
      int N = num_rows();
      
      if (is_lower_triangular())
        {
          for (int j=N-1, jj=(N-1)*incx; 0<=j; j--, jj -= incx)
            {
              x[jj] /= diag[j] ;
              sp_axpy( -x[jj], S[j], x, incx);
            }
          if (alpha != (T) 1.0)
            {
              for (int jj=(N-1)*incx; 0<=jj; jj -=incx)
                x[jj] /= alpha;
            }
        }
      else if (is_upper_triangular())
        {
          for (int j=0, jj=0; j<N; j++, jj += incx)
            {
              x[jj] /= diag[j];
              sp_axpy(- x[jj], S[j], x, incx);
            }
          if (alpha != (T) 1.0)
            {
              for (int jj=(N-1)*incx; 0<=jj; jj -=incx)
                x[jj] /= alpha;
            }
        }
      else
        ERROR_RETURN(1);
      
      return 0;
    }
    
    int transpose_triangular_conj_solve(T alpha, T* x, int incx) const
    {
      if ( ! is_triangular())
        return -1;
      
      int N = num_rows();
      
      if (is_lower_triangular())
        {
          for (int j=N-1, jj=(N-1)*incx; 0<=j; j--, jj -= incx)
            {
              x[jj] /= conj(diag[j]) ;
              sp_conj_axpy( -x[jj], S[j], x, incx);
            }
          if (alpha != (T) 1.0)
            {
              for (int jj=(N-1)*incx; 0<=jj; jj -=incx)
                x[jj] /= alpha;
            }
        }
      else if (is_upper_triangular())
        {
          for (int j=0, jj=0; j<N; j++, jj += incx)
            {
              x[jj] /= conj(diag[j]);
              sp_conj_axpy(- x[jj], S[j], x, incx);
            }
          if (alpha != (T) 1.0)
            {
              for (int jj=(N-1)*incx; 0<=jj; jj -=incx)
                x[jj] /= alpha;
            }
        }
      else
        ERROR_RETURN(1);
      
      return 0;
    }

 public:

  inline T& val(pair<T, int> &VP) { return VP.first; }
  inline int& col_index(pair<T,int> &VP) { return VP.second; } 

  inline const T& val(pair<T, int> const &VP) const { return VP.first; }
  inline int col_index(pair<T,int> const &VP) const { return VP.second; } 


  TSp_mat( int M, int N) : Sp_mat(M,N), S(M), diag() {}

  void destroy()
  {
    // set vector sizes to zero
    (vector<T>(0)).swap(diag);
    (vector< vector< pair<T, int> > > (0) ).swap(S);
  }

/**
    This function is the entry point for all of the insert routines in 
    this implementation.  It fills the sparse matrix, one entry at a time.
    If matrix is declared unit_diagonal, then inserting any diagonal
    values is ignored.  If it is symmetric (upper/lower) or triangular
    (upper/lower) inconsistent values are not caught.  (That is, entries
    into the upper region of a lower triangular matrix is not reported.)

    [NOTE: the base is determined at the creation phase, and can be determined
    by testing whether  BLAS_usgp(A, blas_one_base) returns 1.  If it returns 0,
    then offsets are zero based.]

    @param val  the numeric value of entry A(i,j)
    @param i  the row index of A(i,j)  
    @param j  the column index of A(i,j)

    @return 0 if succesful, 1 otherwise
*/  
  int insert_entry(T val, int i, int j)
  {
    if (is_one_base())        
    {
      i--;
      j--;
    }

    /* make sure the indices are in range */
    ASSERT_RETURN(i >= 0, 1);
    ASSERT_RETURN(i < num_rows(), 1);
    ASSERT_RETURN(j >= 0, 1);
    ASSERT_RETURN(j < num_cols(), 1);

    /* allocate space for the diagonal, if this is the first time
     * trying to insert values.
    */
    if (is_new())
    {
      set_open();
      
      if (is_triangular() || is_symmetric())
      {
        diag.resize(num_rows());

        if (is_unit_diag())
        {
          for (unsigned int ii=0; ii< diag.size(); ii++)
              diag[ii] = T(1.0); 
        }
        else
        {
          for (unsigned int ii=0; ii< diag.size(); ii++)
              diag[ii] = (T) 0.0; 
        }
      }

    }
    if (is_open())
    {

      if (i==j && (is_triangular() || is_symmetric() || is_hermitian()) )
      {
        if (!is_unit_diag())
        {
          diag[i] += val;
        }
        else /* if unit diagonal */
        {
          if (val != (T) 1) 
            ERROR_RETURN(0);    /* tries to insert non-unit diagonal */
        }

        if (is_upper_storage() && i > j)
            ERROR_RETURN(0);    /* tries to fill lower-triangular region */
        else 
        
          if (is_lower_storage() && i < j)
            ERROR_RETURN(0);  /* tries to fill upper-triangular region */

      }
      else
      {
        S[i].push_back( make_pair(val, j) );
      }

      num_nonzeros() ++;
    }


    return 0;
  }

  int insert_entries( int nz, const T* Val, const int *I, const int *J)
  {
    for (int i=0; i<nz; i++)
    {
      insert_entry(Val[i], I[i], J[i]) ;
    }
    return 0;

  }

  int insert_row(int k, int nz, const T* Val, const int *J)
  {
    for (int i=0; i<nz; i++)
      insert_entry(Val[i], k, J[i]);  
    return 0;
  }

  int insert_col(int k, int nz, const T* Val, const int *I)
  {
    for (int i=0; i<nz; i++)
      insert_entry(Val[i], I[i], k);  
    return 0;
  }

  int insert_block(const T* Val, int row_stride, 
        int col_stride, int bi, int bj)
  {
    /* translate from block index to global indices */
    int Iend = K(bi+1);
    int Jend = L(bj+1);
    for (int i=K(bi), r=0; i<Iend; i++, r += row_stride)
      for (int j=L(bi); j<Jend; j++, r += col_stride)
        insert_entry( Val[r], i, j );

    return 0;
  }

  int end_construction()
  {
    return Sp_mat::end_construction();
  }



  int usmv(enum blas_trans_type transa, const T& alpha, const  T* x , int incx, 
    T* y, int incy) const
  {
    
  ASSERT_RETURN(is_valid(), -1);

  if (transa == blas_no_trans)
    mult_vec(alpha, x, incx, y, incy);
  else
  if (transa == blas_conj_trans)
    mult_vec_conj_transpose(alpha, x, incx, y, incy);
  else
  if ( transa == blas_trans)
    mult_vec_transpose(alpha, x, incx, y, incy);
  else
    ERROR_RETURN(1);
    
  
    return 0;

  }


  int usmm(enum blas_order_type ordera, enum blas_trans_type transa, 
    int nrhs, const T& alpha, const  T* b, int ldb, T* C, int ldC) const
  {
    if (ordera == blas_rowmajor)
    {
      /* for each column of C, perform a mat_vec */
      for (int i=0; i<nrhs; i++)
      {
        usmv( transa, alpha, &b[i], ldb, &C[i], ldC );
      }
      return 0;
    }
    else
    if (ordera == blas_colmajor)
    {
      /* for each column of C, perform a mat_vec */
      for (int i=0; i<nrhs; i++)
      {
        usmv( transa, alpha, &b[i*ldb], 1, &C[i*ldC], 1 );
      }
      return 0;
    }
    else
      ERROR_RETURN(1);
  }

  int ussv( enum blas_trans_type transa, const T& alpha,  T* x, int incx) const
  {
      if (transa == blas_trans)
        return transpose_triangular_solve(alpha, x, incx);
      else 
      if (transa == blas_conj_trans)
        return transpose_triangular_conj_solve(alpha, x, incx);
      else
      if (transa == blas_no_trans)
        return triangular_solve(alpha, x, incx);
      else
        ERROR_RETURN(1);


  }

  int ussm( enum blas_order_type ordera, enum blas_trans_type transa, int nrhs,
      const T& alpha, T* C, int ldC) const
  {
    if (ordera == blas_rowmajor)
    {
      /* for each column of C, perform a usmv */
      for (int i=0; i<nrhs; i++)
      {
        ussv( 
            transa, alpha, &C[i], ldC );
      }
      return 0;
    }
    else
    if (ordera == blas_colmajor)
    {
      /* for each column of C, perform a mat_vec */
      for (int i=0; i<nrhs; i++)
      {
        ussv( transa, alpha, &C[i*ldC], 1 );
      }
      return 0;
    }
    else
      ERROR_RETURN(1);
  } 

  void print() const
  {
    Sp_mat::print();  /* print matrix header info */

    /* if there is actual data, print out contents */
    for (int i=0; i<num_rows(); i++)
      for (unsigned int j=0; j< S[i].size(); j++)
        cout << i << "    " << col_index(S[i][j]) <<
              "        "  << val(S[i][j]) << "\n";

    /* if matrix is triangular, print out diagonals */
    if (is_upper_triangular() || is_lower_triangular())
    {
      for (unsigned int i=0; i< diag.size(); i++)
        cout << i << "    " << i << "     " << diag[i] << "\n";
    }
  }    
};

int BLAS_ussp(blas_sparse_matrix A, int pname);
int BLAS_usgp(blas_sparse_matrix A, int pname);
int BLAS_usds(int A);
void table_print();
void print(int A);

} /* namespace SparseBLAS */

using namespace std;
using namespace SparseBLAS;

/* Level 1 */

/* --------------------------- */
/*  Level 1 generic routines   */
/* --------------------------- */

/* dummy routines for real version of usdot to compile. */

template <class T>
void BLAS_xusdot( enum blas_conj_type conj_flag, int nz,
    const T *x,  const int *index,  const T *y, int incy,
    T *r, enum blas_base_type index_base)
{

  T t(0);

  if (index_base == blas_one_base)
    y -= incy;

  if (conj_flag == blas_no_conj)
  {
    for (int i=0; i<nz; i++)
      t += x[i] * y[index[i]*incy];
  }
  else
    for (int i=0; i<nz; i++)
      t += conj(x[i]) * y[index[i]*incy];


  *r = t;
}


template <class T>
void BLAS_xusaxpy(int nz, T alpha, const T *x,
    const int *index, T *y, int incy,
    enum blas_base_type index_base)
{

  if (index_base == blas_one_base)
    y -= incy;

  for (int i=0; i<nz; i++)
  {
//     y[index[i]*incy] +=  (alpha * x[i]);

  }
}

template <class T>
void BLAS_xusga( int nz, const T *y, int incy, T *x, const int *indx,
              enum blas_base_type index_base )
{
  if (index_base == blas_one_base)
    y -= incy;

  for (int i=0; i<nz; i++)
    x[i] = y[indx[i]*incy];

}

template <class T>
void BLAS_xusgz( int nz, T *y, int incy, T *x, const int *indx,
              enum blas_base_type index_base )
{
  if (index_base == blas_one_base)
    y -= incy;

  for (int i=0; i<nz; i++)
  {
    x[i] = y[indx[i]*incy];
    y[indx[i]*incy] = (T) 0.0;
  }

}


template <class T>
void BLAS_xussc(int nz, const T *x, T *y, int incy, const int *index,
    enum blas_base_type index_base)
{
  if (index_base == blas_one_base)
    y -= incy;

  for (int i=0; i<nz; i++)
    y[index[i]*incy] = x[i];

}

/* --------------------------- */
/* Level 2&3 generic precision */
/* --------------------------- */

template <class T>
int BLAS_xuscr_insert_entry(blas_sparse_matrix A,  const T& val, int i, int j)
{
  return ((TSp_mat<T> *)Table[A])->insert_entry(val, i, j);
}

template <class T> 
int BLAS_xuscr_insert_entries(blas_sparse_matrix A, int nz, const T* Val, 
      const int* I, const int *J)
{
  return ((TSp_mat<T>*) Table[A])->insert_entries(nz, Val, I, J);
}


template <class T> 
int BLAS_xuscr_insert_col(blas_sparse_matrix A, int j, int nz, const T* Val, 
      const int* indx)
{
  return ((TSp_mat<T>*) Table[A])->insert_col(j, nz, Val, indx);
}

template <class T> 
int BLAS_xuscr_insert_row(blas_sparse_matrix A, int i, int nz, const T* Val, 
      const int* indx)
{
  return ((TSp_mat<T>*) Table[A])->insert_row(i, nz, Val, indx);
}

template <class T> 
int BLAS_xuscr_insert_clique(blas_sparse_matrix A, int k, int l, const T* Val, 
      const int row_stride, const int col_stride, const int *indx, 
      const int *jndx)
{
  return ((TSp_mat<T>*) Table[A])->insert_clique(k, l, Val, row_stride, 
            col_stride, indx, jndx);
}

template <class T> 
int BLAS_xuscr_insert_block(blas_sparse_matrix A, const T* Val, 
      const int row_stride, const int col_stride, int bi, int bj )
{
  return ((TSp_mat<T>*) Table[A])->insert_block(Val, 
        row_stride, col_stride, bi, bj); 
}

inline int BLAS_xuscr_end(blas_sparse_matrix A)
{
  return (Table[A])->end_construction();
}


template <class T>
int BLAS_xusmv(enum blas_trans_type transa, const T& alpha, 
    blas_sparse_matrix A, const T *x, int incx, T *y, int incy )
{
  TSp_mat<T> *M = (TSp_mat<T> *) Table[A];

  ASSERT_RETURN(M->is_valid(), 1);

  return M->usmv(transa, alpha, x, incx, y, incy);
}


template <class T>
int BLAS_xusmm(enum blas_order_type ordera, enum blas_trans_type transa, 
    int nrhs, const T& alpha, blas_sparse_matrix A, 
    const T *B, int ldB, T* C, int ldC)
{
  TSp_mat<T> *M = (TSp_mat<T> *) Table[A];

  ASSERT_RETURN(M->is_valid(), 1);

  return M->usmm(ordera, transa, nrhs, alpha, B, ldB, C, ldC);
}


template <class T>
int BLAS_xussv(enum blas_trans_type transa, const T& alpha, 
    blas_sparse_matrix A, T *x, int incx)
{
  TSp_mat<T> *M = 
        (TSp_mat<T> *) Table[A];

  ASSERT_RETURN(M->is_valid(), 1);

  return M->ussv(transa, alpha, x, incx);
}

template <class T>
int BLAS_xussm(enum blas_order_type orderA, enum blas_trans_type transa, 
    int nrhs, const T& alpha, blas_sparse_matrix A, T *C, int ldC)
{
  TSp_mat<T> *M = 
        (TSp_mat<T> *) Table[A];

  ASSERT_RETURN(M->is_valid(), 1);

  return M->ussm(orderA, transa, nrhs, alpha, C, ldC);
}

template <class T>
int BLAS_xuscr_begin(int M, int N)
{
  TSp_mat<T> *A = new TSp_mat<T>(M, N);
  return Table_insert(A);
}

#endif
