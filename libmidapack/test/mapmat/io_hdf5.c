/**  
 */
 
#include <hdf5.h> 

int W_indices_hdf5(int m, int nnz, int *indices, char *filename, char *id){
    hid_t       file, dataset;           // file and dataset handles //
    hid_t       datatype, dataspace;     // handles //
    hsize_t     dimsf[2];                // dataset dimensions //
    herr_t      status;                             
    // Create a new file using H5F_ACC_TRUNC access, default file creation properties, and default file access properties. 
    file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    // Describe the size of the array and create the data space for fixed size dataset.  
    dimsf[0] = nnz;
    dimsf[1] = m;
    dataspace = H5Screate_simple(2, dimsf, NULL); 
     // Define datatype for the data in the file. We will store little endian INT numbers. 
    datatype = H5Tcopy(H5T_NATIVE_INT);
    status = H5Tset_order(datatype, H5T_ORDER_LE);
    // Create a new dataset within the file using defined dataspace and datatype and default dataset creation properties. 
    dataset = H5Dcreate(file, id, datatype, dataspace,H5P_DEFAULT);
    // Write the data to the dataset using default transfer properties.*/
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, indices);
    // Close/release resources. 
    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Fclose(file);
    return 0;
}


int W_values_hdf5(int m, int nnz, double *values, char *filename, char *id){
    hid_t       file, dataset;           // file and dataset handles //
    hid_t       datatype, dataspace;     // handles //
    hsize_t     dimsf[2];                // dataset dimensions //
    herr_t      status;                             
    // Create a new file using H5F_ACC_TRUNC access, default file creation properties, and default file access properties. 
    file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    // Describe the size of the array and create the data space for fixed size dataset.  
    dimsf[0] = nnz;
    dimsf[1] = m;
    dataspace = H5Screate_simple(2, dimsf, NULL); 
     // Define datatype for the data in the file. We will store little endian INT numbers. 
    datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    status = H5Tset_order(datatype, H5T_ORDER_LE);
    // Create a new dataset within the file using defined dataspace and datatype and default dataset creation properties. 
    dataset = H5Dcreate(file, id, datatype, dataspace,H5P_DEFAULT);
    // Write the data to the dataset using default transfer properties.*/
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values);
    // Close/release resources. 
    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Fclose(file);
    return 0;
}


int R_indices_hdf5(int *m, int *nnz, int **indices, char *filename, char *id){
    hid_t       file, dataset;         // handles 
    hid_t       datatype, dataspace;   
    H5T_class_t class;                 // datatype class 
    H5T_order_t order;                 // data order 
    size_t      size;                  
    hsize_t     dims_out[2];           // dataset dimensions       
    herr_t      status;                             
    int         status_n, rank;
    // Open the file and the dataset. 
    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    dataset = H5Dopen(file, id);
    // Get datatype and dataspace handles and then query dataset class, order, size, rank and dimensions. */
    datatype  = H5Dget_type(dataset);     // datatype handle  
    class     = H5Tget_class(datatype);
    if (class == H5T_INTEGER)
      printf("Data set has INTEGER type \n");
    order     = H5Tget_order(datatype);
    if (order == H5T_ORDER_LE) 
      printf("Little endian order \n");
    size  = H5Tget_size(datatype);
    printf("Data size is %d \n", size);
    dataspace = H5Dget_space(dataset);    // dataspace handle 
    rank      = H5Sget_simple_extent_ndims(dataspace);
    status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    printf("rank %d, dimensions %lu x %lu \n", rank, (unsigned long)(dims_out[0]), (unsigned long)(dims_out[1]));
    *nnz=dims_out[0];
    *m=dims_out[1];
    *indices = (int *) malloc(dims_out[0]*dims_out[1]*sizeof(int));
    status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, *indices);
    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Fclose(file);
    return 0;
}

int R_values_hdf5(int *m, int *nnz, double **values, char *filename, char *id){
    hid_t       file, dataset;         // handles 
    hid_t       datatype, dataspace;   
    H5T_class_t class;                 // datatype class 
    H5T_order_t order;                 // data order 
    size_t      size;                  
    hsize_t     dims_out[2];           // dataset dimensions       
    herr_t      status;                             
    int         status_n, rank;
    // Open the file and the dataset. 
    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    dataset = H5Dopen(file, id);
    // Get datatype and dataspace handles and then query dataset class, order, size, rank and dimensions. */
    datatype  = H5Dget_type(dataset);     // datatype handle  
    class     = H5Tget_class(datatype);
    if (class == H5T_NATIVE_DOUBLE)
      printf("Data set has DOUBLE type \n");
    order     = H5Tget_order(datatype);
    if (order == H5T_ORDER_LE) 
      printf("Little endian order \n");
    size  = H5Tget_size(datatype);
    printf("Data size is %d \n", size);
    dataspace = H5Dget_space(dataset);    // dataspace handle 
    rank      = H5Sget_simple_extent_ndims(dataspace);
    status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    printf("rank %d, dimensions %lu x %lu \n", rank, (unsigned long)(dims_out[0]), (unsigned long)(dims_out[1]));
    *nnz=dims_out[0];
    *m=dims_out[1];
    *values = (double *) malloc(dims_out[0]*dims_out[1]*sizeof(double));
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, *values);
    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Fclose(file);
    return 0;
}


/*
int main (void){
    int		m, nnz; 
    int         i, j;
    int         *indices;
    double      *values;
    char        filename1[6]="pnt.h5";
    char        filename2[6]="sgn.h5";
    char        id1[6]="pixels";
    char        id2[6]="signal";
    m=10;
    nnz=2;
    indices = (int *) malloc(m*nnz*sizeof(int));
    values = (double *) malloc(m*nnz*sizeof(double));
    for (i = 0; i < m; i++){
      for (j = 0; j < nnz; j++){
	    indices[i*nnz+j] = i*j;
            printf(" %d ", indices[i*nnz+j]);
	    values[i*nnz+j] = 2.0*i*j;
            printf(" %lf ", values[i*nnz+j]);
      }
    }
    printf("\n");
    
    W_indices_hdf5(m, nnz, indices, filename1, id1);
    free(indices);
    R_indices_hdf5(&m, &nnz, &indices, filename1, id1);

    W_values_hdf5(m, nnz, values, filename2, id2);
    free(values);
    R_values_hdf5(&m, &nnz, &values, filename2, id2);


    for (i = 0; i < m; i++){
      for (j = 0; j < nnz; j++){
            printf(" %d ", indices[i*nnz+j]);
            printf(" %lf ", values[i*nnz+j]);
      }
    }
    printf("\n");
    free(indices);
    free(values);
    return 0;
} */    

