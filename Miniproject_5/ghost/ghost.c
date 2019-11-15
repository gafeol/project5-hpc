/****************************************************************
 *                                                              *
 * This file has been written as a sample solution to an        *
 * exercise in a course given at the CSCS Summer School.        *
 * It is made freely available with the understanding that      *
 * every copy of this file must include this header and that    *
 * CSCS take no responsibility for the use of the enclosed      *
 * teaching material.                                           *
 *                                                              *
 * Purpose: Exchange ghost cell in 2 directions using a topology*
 *                                                              *
 * Contents: C-Source                                           *
 *                                                              *
 ****************************************************************/

/* Use only 16 processes for this exercise
 * Send the ghost cell in two directions: left<->right and top<->bottom
 * ranks are connected in a cyclic manner, for instance, rank 0 and 12 are connected
 *
 * process decomposition on 4*4 grid
 *
 * |-----------|
 * | 0| 1| 2| 3|
 * |-----------|
 * | 4| 5| 6| 7|
 * |-----------|
 * | 8| 9|10|11|
 * |-----------|
 * |12|13|14|15|
 * |-----------|
 *
 * Each process works on a 6*6 (SUBDOMAIN) block of data
 * the D corresponds to data, g corresponds to "ghost cells"
 * xggggggggggx
 * gDDDDDDDDDDg
 * gDDDDDDDDDDg
 * gDDDDDDDDDDg
 * gDDDDDDDDDDg
 * gDDDDDDDDDDg
 * gDDDDDDDDDDg
 * gDDDDDDDDDDg
 * gDDDDDDDDDDg
 * gDDDDDDDDDDg
 * gDDDDDDDDDDg
 * xggggggggggx
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>

#define SUBDOMAIN 6
#define DOMAINSIZE (SUBDOMAIN+2)
#define T 0
#define B 1
#define L 2
#define R 3

int tag_from(int x){
    return (x^1);
}

int pos(int i, int j){
    return i*DOMAINSIZE + j;
}

int main(int argc, char *argv[])
{
    int rank, size, i, j, dims[2], periods[2], rank_top, rank_bottom, rank_left, rank_right;
    double data[DOMAINSIZE*DOMAINSIZE];
    MPI_Request request;
    MPI_Status status;
    MPI_Comm cart_comm;
    MPI_Datatype data_ghost;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size!=16) {
        printf("please run this with 16 processors\n");
        MPI_Finalize();
        exit(1);
    }

    // initialize the domain
    for (i=0; i<DOMAINSIZE*DOMAINSIZE; i++) 
        data[i]=rank;

    dims[0]=dims[1]=4;
    periods[0]=periods[1]=1;
    // TODO: Create a Cartesian communicator (4*4) with periodic boundaries (we do not allow
    // the reordering of ranks) and use it to find your neighboring
    // ranks in all dimensions in a cyclic manner.
    //int MPI_Cart_create(MPI_Comm comm_old, int ndims, int *dims, int *periods,int reorder, MPI_Comm *comm_cart)
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);

    // TODO: find your top/bottom/left/right neighbor using the new communicator, see MPI_Cart_shift()
    // rank_top, rank_bottom
    // rank_left, rank_right

    MPI_Cart_shift(cart_comm, 0, 1, &rank_top, &rank_bottom);
    MPI_Cart_shift(cart_comm, 1, 1, &rank_left, &rank_right);

    //  TODO: create derived datatype data_ghost, create a datatype for sending the column, see MPI_Type_vector() and MPI_Type_commit()
    // data_ghost
    // int MPI_Type_vector(int count, int blocklength, int stride, MPI_Datatype old_type, MPI_Datatype *newtype_p);
    MPI_Type_vector(SUBDOMAIN, 1, 0, MPI_DOUBLE, &data_ghost);
    MPI_Type_commit(&data_ghost);
    //MPI_Type_vector(1, SUBDOMAIN, 0, MPI_DOUBLE, &data_ghost); 
    //TODO: DEU MUITO ERRADO ASSIM, PQ

    //  TODO: ghost cell exchange with the neighbouring cells in all directions
    //  use MPI_Irecv(), MPI_Send(), MPI_Wait() or other viable alternatives
    //  to the top
    int cnt;

    // Irecv section

    double rcv_top[SUBDOMAIN], snd_top[SUBDOMAIN]; 
    cnt = 0; i = 1;
    for(j=1;j<=SUBDOMAIN;j++)
        snd_top[cnt++] = data[pos(i, j)];
    MPI_Status top_status;
    MPI_Request snd_top_req, rcv_top_req;
    MPI_Irecv(rcv_top, SUBDOMAIN, data_ghost, rank_top, tag_from(T), cart_comm, &rcv_top_req);
    MPI_Isend(snd_top, SUBDOMAIN, data_ghost, rank_top, T, cart_comm, &snd_top_req);

    //  to the bottom
    double rcv_bottom[SUBDOMAIN], snd_bottom[SUBDOMAIN];
    cnt = 0; i = DOMAINSIZE-2;
    for(j=1;j<=SUBDOMAIN;j++)
        snd_bottom[cnt++] = data[pos(i, j)];
    MPI_Status bottom_status;
    MPI_Request snd_bottom_req, rcv_bottom_req;
    MPI_Irecv(rcv_bottom, SUBDOMAIN, data_ghost, rank_bottom, tag_from(B), cart_comm, &rcv_bottom_req);
    MPI_Isend(snd_bottom, SUBDOMAIN, data_ghost, rank_bottom, B, cart_comm, &snd_bottom_req);

    //  to the left
    double rcv_left[SUBDOMAIN], snd_left[SUBDOMAIN];
    cnt = 0; j = 1;
    for(int i=1;i<=SUBDOMAIN;i++)
        snd_left[cnt++] = data[pos(i, j)];
    MPI_Status left_status;
    MPI_Request snd_left_req, rcv_left_req;
    MPI_Irecv(rcv_left, SUBDOMAIN, data_ghost, rank_left, tag_from(L), cart_comm, &rcv_left_req);
    MPI_Isend(snd_left, SUBDOMAIN, data_ghost, rank_left, L, cart_comm, &snd_left_req);

    //  to the right
    double rcv_right[SUBDOMAIN], snd_right[SUBDOMAIN];
    cnt = 0; j = DOMAINSIZE-2;
    for(int i=1;i<=SUBDOMAIN;i++)
        snd_right[cnt++] = data[pos(i, j)];
    MPI_Status right_status;
    MPI_Request snd_right_req, rcv_right_req;
    MPI_Irecv(rcv_right, SUBDOMAIN, data_ghost, rank_right, tag_from(R), cart_comm, &rcv_right_req);
    MPI_Isend(snd_right, SUBDOMAIN, data_ghost, rank_right, R, cart_comm, &snd_right_req);

    //Wait section

    // top section
    MPI_Wait(&rcv_top_req, &top_status);
    MPI_Wait(&snd_top_req, &top_status);
    i = 0; cnt = 0;
    for(j=1;j<=SUBDOMAIN;j++)
        data[pos(i, j)] = rcv_top[cnt++];

    // bottom section
    MPI_Wait(&rcv_bottom_req, &bottom_status);
    MPI_Wait(&snd_bottom_req, &bottom_status);
    i = DOMAINSIZE-1; cnt = 0;
    for(int j=1;j<=SUBDOMAIN;j++)
        data[pos(i, j)] = rcv_bottom[cnt++];

    // left section
    MPI_Wait(&rcv_left_req, &left_status);
    MPI_Wait(&snd_left_req, &left_status);
    j = 0; cnt = 0;
    for(int i=1;i<=SUBDOMAIN;i++)
        data[pos(i, j)] = rcv_left[cnt++];
    //right section
    MPI_Wait(&rcv_right_req, &right_status);
    MPI_Wait(&snd_right_req, &right_status);
    j = DOMAINSIZE-1; cnt = 0;
    for(int i=1;i<=SUBDOMAIN;i++)
        data[pos(i, j)] = rcv_right[cnt++];

    if (rank==9) {
        printf("data of rank 9 after communication\n");
        for (j=0; j<DOMAINSIZE; j++) {
            for (i=0; i<DOMAINSIZE; i++) {
                printf("%.1f ", data[i+j*DOMAINSIZE]);
            }
            printf("\n");
        }
    }


    MPI_Type_free(&data_ghost);
    MPI_Comm_free(&cart_comm);
    MPI_Finalize();

    return 0;
}