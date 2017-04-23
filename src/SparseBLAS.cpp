#include "SparseBLAS.h"

namespace SparseBLAS
{
  vector<Sp_mat *> Table;
  unsigned int Table_active_matrices = 0;

/* 
  finds an empty slot in global sparse marix table, and fills it with
  given entry.  Returns -1 if no spot found, or table is corrupt.
*/
int Table_insert(Sp_mat* S)
{
  if (Table_active_matrices <= Table.size())
  {
    Table.push_back(S);
    Table_active_matrices++;
    return Table.size() - 1;
  }
  else
  {
    /* there is an available slot; find it. */
    for (unsigned int i=0; i<Table.size(); i++)
    {
      if (Table[i] == NULL)
      {
        Table[i] = S;
        Table_active_matrices++;
        return i;
      }
    }
  }

  return -1;
}

/* 
  removes an exisiting sparse matrix from global table.  Returns 0, if
  successfull, 1 otherwise.
*/
int Table_remove(unsigned int i)
{
  if (i < Table.size() && Table[i] != NULL)
  {
    Table[i] = NULL;
    Table_active_matrices--;
    return 0;
  }
  else
    return -1;
}

void table_print()
{
  cout << "Table has " << Table.size() << " element(s). \n";
  for (unsigned int i=0; i< Table.size(); i++)
  {
    if (Table[i] != 0)
    {
      cout << "***** Table[" << i << "]: \n";
      Table[i]->print();
      cout << "\n\n";
    }
  }
}

void print(int A)
{
  cout << "\n";
  Table[A]->print();
  cout << "\n";
}


Sp_mat :: Sp_mat(int M, int N) : 
  num_rows_(M),         /* default construction */
  num_cols_(N),
  num_nonzeros_(0),
  
  void_(0),
  nnew_(1),
  open_(0),
  valid_(0),
  
  unit_diag_(0),
  upper_triangular_(0),
  lower_triangular_(0),
  upper_symmetric_(0),
  lower_symmetric_(0),
  upper_hermitian_(0),
  lower_hermitian_(0),
  general_(0),
  one_base_(0),
  Mb_(0),
  Nb_(0),
  k_(0),
  l_(0),
  rowmajor_(0),
  colmajor_(0),
  opt_regular_(0),
  opt_irregular_(1),
  opt_block_(0),
  opt_unassembled_(0),
  K_(),
  L_() {}

int& Sp_mat :: num_rows()           { return num_rows_; }
int& Sp_mat :: num_cols()           { return num_cols_; }
int& Sp_mat :: num_nonzeros()         { return num_nonzeros_;}
int Sp_mat :: num_rows() const        { return num_rows_; }
int Sp_mat :: num_cols() const        { return num_cols_; }
int Sp_mat :: num_nonzeros() const      { return num_nonzeros_;}
int Sp_mat :: is_one_base() const     { return (one_base_ ? 1 : 0); }
int Sp_mat :: is_zero_base() const    { return (one_base_ ? 0 : 1); }
int Sp_mat :: is_void() const         { return void_; }
int Sp_mat :: is_new() const          { return nnew_; }
int Sp_mat :: is_open() const         { return open_; }
int Sp_mat :: is_valid() const        { return valid_; }
int Sp_mat :: is_unit_diag() const    { return unit_diag_; }
int Sp_mat :: is_upper_triangular() const   { return upper_triangular_;}
int Sp_mat :: is_lower_triangular() const   { return lower_triangular_;}
int Sp_mat :: is_triangular() const     { return upper_triangular_ || lower_triangular_; }
int Sp_mat :: is_lower_symmetric() const    { return lower_symmetric_; }
int Sp_mat :: is_upper_symmetric() const    { return upper_symmetric_; }
int Sp_mat :: is_symmetric() const      { return upper_symmetric_ || lower_symmetric_; }
int Sp_mat :: is_lower_hermitian() const    { return lower_hermitian_; }
int Sp_mat :: is_upper_hermitian() const    { return upper_hermitian_; }
int Sp_mat :: is_hermitian() const  { return lower_hermitian_ || upper_hermitian_; }
int Sp_mat :: is_general() const { return !( is_hermitian() || is_symmetric()) ; }
int Sp_mat :: is_lower_storage() const { return is_lower_triangular() || is_lower_symmetric() || is_lower_hermitian() ; }
int Sp_mat :: is_upper_storage() const { return is_upper_triangular() || is_upper_symmetric() || is_upper_hermitian() ; }
int Sp_mat :: is_opt_regular() const { return opt_regular_; }
int Sp_mat :: is_opt_irregular() const { return opt_irregular_; }
int Sp_mat :: is_opt_block() const { return opt_block_;} 
int Sp_mat :: is_opt_unassembled() const { return opt_unassembled_;}
int Sp_mat :: K(int i) const { return (k_ ? i*k_ : K_[i] ); }
int Sp_mat :: L(int i) const { return (l_ ? i*l_ : L_[i] ); }
int Sp_mat :: is_rowmajor() const { return rowmajor_; }
int Sp_mat :: is_colmajor() const { return colmajor_; }
void Sp_mat :: set_one_base()   { one_base_ = 1; }
void Sp_mat :: set_zero_base()  { one_base_ = 0; }
void Sp_mat :: set_void()       { void_ = 1;  nnew_ = open_ =  valid_ = 0;}
void Sp_mat :: set_new()        { nnew_ = 1;  void_ = open_ =  valid_ = 0;}
void Sp_mat :: set_open()       { open_ = 1;  void_ = nnew_  = valid_ = 0;}
void Sp_mat :: set_valid()      { valid_ = 1; void_ = nnew_ =  open_ = 0; }
void Sp_mat :: set_unit_diag()    { unit_diag_ = 1;}
void Sp_mat :: set_upper_triangular()   { upper_triangular_ = 1; }
void Sp_mat :: set_lower_triangular()   { lower_triangular_ = 1; }
void Sp_mat :: set_upper_symmetric()  { upper_symmetric_ = 1; }
void Sp_mat :: set_lower_symmetric()  { lower_symmetric_ = 1; }
void Sp_mat :: set_upper_hermitian()  { upper_hermitian_ = 1; }
void Sp_mat :: set_lower_hermitian()  { lower_hermitian_ = 1; }

void Sp_mat :: set_const_block_parameters(int Mb, int Nb, int k, int l)
{
  Mb_ = Mb;
  Nb_ = Nb;
  k_ = k;
  l_ = l;
}

void Sp_mat :: set_var_block_parameters(int Mb, int Nb, const int *k, const int *l)
{
  Mb_ = Mb;
  Nb_ = Nb;
  k_ = 0;
  l_ = 0;
  
  K_.resize(Mb+1);
  K_[0] = 0;
  for (int i=0; i<Mb; i++)
    K_[i+1] = k[i] + K_[i];
  
  L_.resize(Nb+1);
  L_[0] = 0;
  for (int j=0; j<Mb; j++)
    K_[j+1] = k[j] + K_[j];
}

int Sp_mat :: end_construction()
{
  if (is_open() || is_new())
    {
      set_valid();
      
      return 0;
    }
  else
    ERROR_RETURN(1);
}

void Sp_mat :: destroy() {}

Sp_mat :: ~Sp_mat() {}

void Sp_mat :: print() const
{
  cout << "State : " <<
    (is_void() ? "void" : 
     is_new()  ? "new" :  
     is_open() ? "open" :
     is_valid() ? "valid" : "unknown") << "\n";
  
  cout << "M = " <<  num_rows() <<  "  N = " << num_cols() <<
    "  nz = " << num_nonzeros() << "\n";
  
#define yesno(exp) ( (exp) ? "yes" : "no" )
  cout << "upper_triangular: " << yesno(is_upper_triangular()) << "\n";
  cout << "lower_triangular: " << yesno(is_lower_triangular()) << "\n";
  
  cout << "regular:    " << yesno(is_opt_regular()) << "\n";
  cout << "irregular:  " << yesno(is_opt_irregular()) << "\n";
  cout << "block:      " << yesno(is_opt_block()) << "\n";
  cout << "unassembled:" << yesno(is_opt_unassembled()) << "\n";
#undef yesno
}

/*------------------------------------*/
/* Non-precision Sparse BLAS routines */
/*------------------------------------*/

/* -------- */
/*  USSP()  */
/* -------- */
int BLAS_ussp(blas_sparse_matrix A, int pname)
{
    Sp_mat *S = Table[A];


    /* Note: these are returns, in the case */
    /* statement, so "break" is not needed.  */

    switch (pname)
    {
        case (blas_zero_base) : S->set_zero_base(); break;
        case (blas_one_base)  : S->set_one_base(); break;

        case (blas_unit_diag) : S->set_unit_diag(); break;

        case (blas_lower_triangular) : S->set_lower_triangular(); break;
        case (blas_upper_triangular) : S->set_upper_triangular(); break;

        case (blas_lower_symmetric) : S->set_lower_symmetric(); break;
        case (blas_upper_symmetric) : S->set_upper_symmetric(); break;

        case (blas_lower_hermitian) : S->set_lower_hermitian(); break;
        case (blas_upper_hermitian) : S->set_upper_hermitian(); break;


                                        /* optimizations not used */
        case (blas_regular )        :
        case (blas_irregular)       : 
        case (blas_block)           :
        case (blas_unassembled)     :    return 0;

        default:                return -1;  /* invalid property */
    }

    return 0;

}

/* -------- */
/*  USGP()  */
/* -------- */

int BLAS_usgp(blas_sparse_matrix A, int pname)
{
  Sp_mat *S = Table[A];

    switch (pname)
    {
        case (blas_num_rows)          : return S->num_rows(); 
        case (blas_num_cols)          : return S->num_cols(); 
        case (blas_num_nonzeros)      : return S->num_nonzeros(); 

        case (blas_lower_triangular) : return S->is_lower_triangular(); break;
        case (blas_upper_triangular) : return S->is_upper_triangular(); break;

        case (blas_general)           : return S->is_general();
        case (blas_symmetric)         : return S->is_symmetric();
        case (blas_hermitian)         : return S->is_hermitian();

        case (blas_zero_base) : return S->is_zero_base(); 
        case (blas_one_base) : return  S->is_one_base(); 


        case (blas_rowmajor)          : return S->is_rowmajor();
        case (blas_colmajor)          : return S->is_colmajor();
        case (blas_new_handle)        : return S->is_new();
        case (blas_valid_handle)      : return S->is_valid();
        case (blas_open_handle)       : return S->is_open();
        case (blas_invalid_handle)    : return S->is_void();

        case (blas_regular)           : return S->is_opt_regular();
        case (blas_irregular)         : return S->is_opt_irregular();
        case (blas_block)             : return S->is_opt_block();
        case (blas_unassembled)       : return S->is_opt_unassembled();

        default:                return -1;  /* invalid property */
    }
}

/* -------- */
/*  USDS()  */
/* -------- */
int BLAS_usds(int A)
{
  Sp_mat *S  = Table[A];

  S->destroy();
  Table_remove(A);

  return 0;
}

}
