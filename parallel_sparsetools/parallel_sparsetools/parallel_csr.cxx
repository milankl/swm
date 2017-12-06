#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL _scipy_sparse_sparsetools_ARRAY_API

#include "parallel_sparsetools.h"
#include "parallel_csr.h"

extern "C" {
#include "parallel_csr_impl.h"
}
