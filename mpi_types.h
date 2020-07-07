#ifdef HAVE_MPI_F08_MODULE
#define MPI_REQUEST_TYPE type(MPI_Request)
#define MPI_STATUS_TYPE type(MPI_Status)
#define MPI_COMM_TYPE type(MPI_Comm)
#else
#define MPI_REQUEST_TYPE integer
#define MPI_STATUS_TYPE integer
#define MPI_COMM_TYPE integer
#endif
