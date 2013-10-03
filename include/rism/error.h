#ifndef __RISM_ERROR_H
#define __RISM_ERROR_H

void set_std_error();
void print_error(const char *m);

#ifdef MPI
void set_mpi_error(int e);
#endif

#endif /* __RISM_ERROR_H */
